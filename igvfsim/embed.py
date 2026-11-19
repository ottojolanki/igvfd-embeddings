"""Embed description texts with a sentence-transformers model, with an on-disk cache.

    # sandbox (CPU): small models
    python -m igvfsim.embed --model BAAI/bge-small-en-v1.5 --variant full
    # host Mac (MPS): large models; same cache directory, so the sandbox eval picks them up
    python -m igvfsim.embed --model Qwen/Qwen3-Embedding-8B --variant full --device mps
    python -m igvfsim.embed --config configs/large.yaml       # every (model, variant) in a config

Cache: data/emb_cache/<model slug>/shard_*.npz, keyed by sha1(model | prompt | max_len | text).
Identical texts are embedded once (many Sets share a description), and reruns only encode
texts not seen before. Shards are append-only, so a host run and a sandbox run never clobber
each other.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import platform
import re
import sys
import time
import uuid

import numpy as np

from .load import DATA

CACHE = os.path.join(DATA, "emb_cache")

# Per-model prompts. "doc" is used for documents (and for doc-to-doc similarity on both sides);
# "query" for free-text search. Qwen3-Embedding is instruction-tuned on the query side only.
PROMPTS = {
    "qwen3-embedding": {"doc": "", "query": "Instruct: Given a search query, retrieve matching "
                        "genomics experiments from the IGVF data portal\nQuery: "},
    "e5": {"doc": "passage: ", "query": "query: "},
    "bge-m3": {"doc": "", "query": ""},
    "bge": {"doc": "", "query": "Represent this sentence for searching relevant passages: "},
    "nomic": {"doc": "search_document: ", "query": "search_query: "},
    "gte": {"doc": "", "query": ""},
}
# Symmetric doc-to-doc instruction for instruction-tuned models (an experiment option).
SYM_INSTRUCT = ("Instruct: Given a description of a genomics experiment, retrieve descriptions "
                "of experiments that study similar biology with similar assays\nQuery: ")
TRUST_REMOTE = ("nomic", "gte-modernbert", "jina", "stella")


def prompts_for(model: str) -> dict:
    m = model.lower()
    for k, v in PROMPTS.items():
        if k in m:
            return v
    return {"doc": "", "query": ""}


def slug(model: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "__", model)


def default_device() -> str:
    import torch
    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def _key(model, prompt, max_len, text) -> str:
    return hashlib.sha1(f"{model}\x1f{prompt}\x1f{max_len}\x1f{text}".encode()).hexdigest()


class EmbCache:
    def __init__(self, model: str):
        self.dir = os.path.join(CACHE, slug(model))
        os.makedirs(self.dir, exist_ok=True)
        self.index: dict[str, tuple[str, int]] = {}
        self._arrays: dict[str, np.ndarray] = {}
        for fn in sorted(glob.glob(os.path.join(self.dir, "shard_*.npz"))):
            with np.load(fn) as z:
                for i, k in enumerate(z["keys"]):
                    self.index[str(k)] = (fn, i)

    def encode_stats(self, keys: list[str]) -> dict:
        """Original encode time/device of the shards holding these keys (for runtime reporting)."""
        shards = {self.index[k][0] for k in keys if k in self.index}
        secs, n, devices = 0.0, 0, set()
        for fn in shards:
            with np.load(fn) as z:
                meta = json.loads(str(z["meta"]))
                secs += meta.get("seconds", 0.0)
                n += len(z["keys"])
                devices.add(meta.get("device"))
        return {"seconds": round(secs, 1), "n_encoded": n, "device": ",".join(sorted(d for d in devices if d))}

    def get(self, keys: list[str]) -> np.ndarray:
        out = []
        for k in keys:
            fn, i = self.index[k]
            if fn not in self._arrays:
                with np.load(fn) as z:
                    self._arrays[fn] = z["emb"]
            out.append(self._arrays[fn][i])
        return np.stack(out)

    def put(self, keys: list[str], emb: np.ndarray, meta: dict):
        fn = os.path.join(self.dir, f"shard_{time.strftime('%Y%m%d-%H%M%S')}_{uuid.uuid4().hex[:6]}.npz")
        np.savez(fn, keys=np.array(keys), emb=emb.astype(np.float32), meta=json.dumps(meta))
        for i, k in enumerate(keys):
            self.index[k] = (fn, i)
        self._arrays[fn] = emb.astype(np.float32)


class STModel:
    _loaded: dict = {}

    def __init__(self, model: str, device: str | None = None, max_len: int = 512, dtype: str = "auto"):
        self.model_name, self.max_len = model, max_len
        self.device = device or default_device()
        # fp16 on GPU/MPS; fp32 on CPU (bf16 checkpoints like Qwen3 otherwise load as bf16,
        # measured ~5x slower than fp32 on the sandbox's ARM CPU)
        self.dtype = dtype if dtype != "auto" else ("fp16" if self.device in ("mps", "cuda") else "fp32")
        key = (model, self.device, self.dtype)
        if key not in STModel._loaded:
            from sentence_transformers import SentenceTransformer
            kw = {}
            if any(t in model.lower() for t in TRUST_REMOTE):
                kw["trust_remote_code"] = True
            import torch
            kw["model_kwargs"] = {"dtype": {"fp16": torch.float16, "bf16": torch.bfloat16,
                                            "fp32": torch.float32}[self.dtype]}
            t0 = time.time()
            m = SentenceTransformer(model, device=self.device, **kw)
            m.max_seq_length = max_len
            print(f"loaded {model} on {self.device} ({self.dtype}) in {time.time() - t0:.1f}s", file=sys.stderr)
            STModel._loaded[key] = m
        self.st = STModel._loaded[key]

    def encode(self, texts: list[str], prompt: str = "", batch_size: int = 32) -> np.ndarray:
        return self.st.encode(texts, prompt=prompt or None, batch_size=batch_size,
                              normalize_embeddings=True, show_progress_bar=False,
                              convert_to_numpy=True).astype(np.float32)


def embed_texts(texts: list[str], model: str, *, prompt: str | None = None, device: str | None = None,
                max_len: int = 512, batch_size: int | None = None, chunk: int = 256,
                dtype: str = "auto") -> tuple[np.ndarray, dict]:
    """Embeddings for texts (row-aligned), using and filling the cache. Returns (emb, stats)."""
    prompt = prompts_for(model)["doc"] if prompt is None else prompt
    cache = EmbCache(model)
    keys = [_key(model, prompt, max_len, t) for t in texts]
    todo = {}
    for k, t in zip(keys, texts):
        if k not in cache.index and k not in todo:
            todo[k] = t
    stats = {"n": len(texts), "unique": len(set(keys)), "encoded": len(todo), "seconds": 0.0,
             "device": None}
    if todo:
        m = STModel(model, device, max_len, dtype)
        stats["device"] = m.device
        bs = batch_size or (8 if m.device == "cpu" else 32)  # measured: bs 8 fastest on CPU
        # sort by length so batches pad less; checkpoint to the cache every `chunk` texts so an
        # interrupted run resumes where it stopped
        items = sorted(todo.items(), key=lambda kv: len(kv[1]))
        t_all = time.time()
        for c0 in range(0, len(items), chunk):
            part = items[c0:c0 + chunk]
            t0 = time.time()
            emb = m.encode([t for _, t in part], prompt=prompt, batch_size=bs)
            if not np.isfinite(emb).all() or (np.linalg.norm(emb, axis=1) < 0.5).any():
                # fp16 can overflow in larger models (NaN/zero vectors). Never cache those:
                # switch this model to fp32 on the same device and redo the chunk.
                print(f"  {model}: non-finite embeddings in {m.dtype}; switching to fp32", file=sys.stderr)
                m = STModel(model, m.device, max_len, "fp32")
                t0 = time.time()
                emb = m.encode([t for _, t in part], prompt=prompt, batch_size=bs)
                if not np.isfinite(emb).all():
                    raise RuntimeError(f"{model}: non-finite embeddings even in fp32")
            secs = round(time.time() - t0, 1)
            cache.put([k for k, _ in part], emb, {"model": model, "prompt": prompt, "max_len": max_len,
                                                  "device": m.device, "dtype": m.dtype,
                                                  "host": platform.platform(),
                                                  "seconds": secs, "batch_size": bs})
            done = c0 + len(part)
            el = time.time() - t_all
            print(f"  {model}: {done}/{len(items)} encoded, {done / el:.1f} texts/s, "
                  f"ETA {(len(items) - done) * el / done / 60:.1f} min", file=sys.stderr)
        stats["seconds"] = round(time.time() - t_all, 1)
    stats["original"] = cache.encode_stats(keys)
    return cache.get(keys), stats


def embed_query(text: str, model: str, device: str | None = None, max_len: int = 512) -> np.ndarray:
    return STModel(model, device, max_len).encode([text], prompt=prompts_for(model)["query"])[0]


# ----------------------------------------------------------------------------- lexical baselines

def tfidf_matrix(texts: list[str]):
    from sklearn.feature_extraction.text import TfidfVectorizer
    v = TfidfVectorizer(ngram_range=(1, 2), sublinear_tf=True, min_df=1)
    return v.fit_transform(texts), v  # rows are L2-normalized (sparse)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", help="repeatable")
    ap.add_argument("--variant", action="append", help="description variant(s); default full")
    ap.add_argument("--config", help="YAML experiment config; embeds every (model, variant) in it")
    ap.add_argument("--device", default=None, help="cpu / mps / cuda (default: auto)")
    ap.add_argument("--max-len", type=int, default=512)
    ap.add_argument("--batch-size", type=int, default=None, help="default 8 on CPU, 32 on GPU")
    ap.add_argument("--sym-instruct", action="store_true", help="also embed with the symmetric instruction")
    ap.add_argument("--dtype", default="auto", choices=["auto", "fp16", "bf16", "fp32"],
                    help="auto = fp16 on mps/cuda, fp32 on cpu")
    ap.add_argument("--benchmark", type=int, metavar="N",
                    help="time encoding N sampled unique texts per model; no cache writes")
    a = ap.parse_args()
    from .describe import load_corpus
    jobs = []
    if a.config:
        import yaml
        cfg = yaml.safe_load(open(a.config))
        for r in cfg["runs"]:
            if r.get("embedder", "st") == "st":
                jobs.append((r["model"], r.get("variant", "full"), r.get("prompt")))
    for m in a.model or []:
        for v in a.variant or ["full"]:
            jobs.append((m, v, None))
            if a.sym_instruct:
                jobs.append((m, v, SYM_INSTRUCT))
    jobs = list(dict.fromkeys(jobs))
    corpus = load_corpus()
    if a.benchmark:
        return benchmark(corpus, jobs, a)
    for model, variant, prompt in jobs:
        if prompt == "sym_instruct":
            prompt = SYM_INSTRUCT
        texts = corpus.texts(variant)
        emb, st = embed_texts(texts, model, prompt=prompt, device=a.device, max_len=a.max_len,
                              batch_size=a.batch_size, dtype=a.dtype)
        rate = f"{st['encoded'] / st['seconds']:.1f} texts/s" if st["seconds"] else "all cached"
        print(f"{model} [{variant}{' sym' if prompt else ''}]: {st['n']} texts, {st['unique']} unique, "
              f"encoded {st['encoded']} on {st['device']} in {st['seconds']}s ({rate}); dim {emb.shape[1]}",
              flush=True)
    # Also cache the search-benchmark queries (query prompt) for every model, so search_eval can
    # score models that are too large to load in the sandbox.
    qpath = "configs/search_queries.yaml"
    if os.path.exists(qpath):
        import yaml
        queries = [q["q"] for q in yaml.safe_load(open(qpath))["queries"]]
        for model in dict.fromkeys(j[0] for j in jobs):
            embed_texts(queries, model, prompt=prompts_for(model)["query"], device=a.device,
                        max_len=a.max_len, batch_size=a.batch_size, dtype=a.dtype)
        print(f"cached {len(queries)} search queries for {len(set(j[0] for j in jobs))} model(s)", flush=True)


def benchmark(corpus, jobs, a):
    """Throughput only: encode a fixed random sample of unique texts, report load time and texts/s."""
    import random
    for model in dict.fromkeys(j[0] for j in jobs):
        variant = next(j[1] for j in jobs if j[0] == model)
        uniq = sorted(set(corpus.texts(variant)))
        sample = random.Random(0).sample(uniq, min(a.benchmark, len(uniq)))
        t0 = time.time()
        m = STModel(model, a.device, a.max_len, a.dtype)
        load_s = time.time() - t0
        bs = a.batch_size or (8 if m.device == "cpu" else 32)
        m.encode(sample[:bs], batch_size=bs)  # warm-up (kernel compilation on MPS)
        t0 = time.time()
        m.encode(sample, prompt=prompts_for(model)["doc"], batch_size=bs)
        dt = time.time() - t0
        est = len(uniq) / (len(sample) / dt)
        print(f"BENCH {model} device={m.device} dtype={m.dtype} bs={bs}: load {load_s:.1f}s, "
              f"{len(sample) / dt:.1f} texts/s -> ~{est / 60:.1f} min for {len(uniq)} unique '{variant}' texts",
              flush=True)


if __name__ == "__main__":
    main()
