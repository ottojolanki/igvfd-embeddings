# Field coverage of released Sets

Share of released objects where each field (dotted path into the embedded object) is non-empty. Fields requested from the API but never populated show 0%.

## MeasurementSet (n=7938, pulled 2026-11-18T23:04:49+00:00)

| field | % populated | example |
|---|---:|---|
| `@id` | 100.0 | ["/measurement-sets/IGVFDS0001EZIS/"] |
| `@type` | 100.0 | ["MeasurementSet", "FileSet", "Item"] |
| `accession` | 100.0 | ["IGVFDS0001EZIS"] |
| `assay_term.@id` | 100.0 | ["/assay-terms/OBI_0002764/"] |
| `assay_term.assay_slims` | 100.0 | ["chromatin accessibility"] |
| `assay_term.term_name` | 100.0 | ["single-cell ATAC-seq"] |
| `assay_titles` | 100.0 | ["single-cell ATAC-seq"] |
| `auxiliary_sets.accession` | 28.7 | ["IGVFDS7766RXOL"] |
| `auxiliary_sets.file_set_type` | 28.7 | ["lipid-conjugated oligo sequencing"] |
| `award.component` | 100.0 | ["mapping"] |
| `collections` | 27.9 | ["PD single cell multiomics"] |
| `control_file_sets.accession` | 19.5 | ["IGVFDS9173SNMO"] |
| `control_types` | 5.7 | ["non-targeting"] |
| `crispr_screen_biometric.term_name` | 27.6 | ["protein abundance"] |
| `crispr_screen_readout` | 27.6 | ["gRNA sequencing"] |
| `description` | 12.3 | ["CRISPRa Tet1-dCas9 screen with the PWS sublibrary in matSNRPN-GFP iPSCs"] |
| `donors.sex` | 99.1 | ["male"] |
| `donors.strain_background` | 6.1 | ["B6PWKF1/J", "B6NODF1/J", "C57BL/6J (B6)"] |
| `donors.taxa` | 99.1 | ["Homo sapiens"] |
| `file_set_type` | 100.0 | ["experimental data"] |
| `functional_assay_mechanisms.term_name` | 1.8 | ["protein secretion"] |
| `input_for.accession` | 97.2 | ["IGVFDS1776AENK", "IGVFDS2235HRBG", "IGVFDS5278IBSB"] |
| `lab.title` | 100.0 | ["Jason Buenrostro, Broad"] |
| `library_construction_platform.term_name` | 0.0 |  |
| `multiome_size` | 33.0 | [2] |
| `preferred_assay_titles` | 100.0 | ["SHARE-seq"] |
| `primary_crispr_screen_readout` | 0.0 |  |
| `publications.@id` | 7.4 | ["/publications/c40886f7-0c45-45e2-b00b-df7632091c29/"] |
| `related_measurement_sets.measurement_sets.accession` | 32.6 | ["IGVFDS6965XHYL"] |
| `release_timestamp` | 100.0 | ["2025-01-14T19:05:58.092766+00:00"] |
| `samples.@id` | 100.0 | ["/primary-cells/IGVFSM4247WODL/"] |
| `samples.@type` | 100.0 | ["PrimaryCell", "Biosample", "Sample"] |
| `samples.cellular_sub_pool` | 34.7 | ["SS-PKR 6"] |
| `samples.classifications` | 100.0 | ["primary cell"] |
| `samples.construct_library_sets` | 0.0 |  |
| `samples.disease_terms.term_name` | 0.0 |  |
| `samples.modifications.modality` | 30.3 | ["activation"] |
| `samples.modifications.summary` | 0.0 |  |
| `samples.sample_terms.@id` | 100.0 | ["/sample-terms/UBERON_0002371/"] |
| `samples.sample_terms.term_name` | 100.0 | ["bone marrow"] |
| `samples.sorted_from` | 0.0 |  |
| `samples.summary` | 100.0 | ["Homo sapiens (male, 37 years) bone marrow cell (cellular sub pool: SS-PKR 6)"] |
| `samples.targeted_sample_term.@id` | 13.9 | ["/sample-terms/EFO_0004905/"] |
| `samples.targeted_sample_term.term_name` | 13.9 | ["induced pluripotent stem cell"] |
| `samples.taxa` | 100.0 | ["Homo sapiens"] |
| `samples.time_post_change` | 0.0 |  |
| `samples.treatments.summary` | 10.9 | ["Treatment of 2 ng/mL interleukin-1 beta for 7 days"] |
| `samples.treatments.treatment_term_name` | 10.9 | ["interleukin-1 beta"] |
| `sequencing_library_types` | 4.4 | ["mRNA enriched"] |
| `status` | 100.0 | ["released"] |
| `summary` | 100.0 | ["single-cell ATAC-seq (SHARE-seq)"] |
| `targeted_genes.symbol` | 13.7 | ["SNRPN"] |
| `targeted_loci` | 0.0 |  |

## AnalysisSet (n=3877, pulled 2026-11-18T23:05:50+00:00)

| field | % populated | example |
|---|---:|---|
| `@id` | 100.0 | ["/analysis-sets/IGVFDS3335QAOG/"] |
| `@type` | 100.0 | ["AnalysisSet", "FileSet", "Item"] |
| `accession` | 100.0 | ["IGVFDS3335QAOG"] |
| `assay_titles` | 100.0 | ["single nucleus methylation and chromatin conformation capture sequencing"] |
| `award.component` | 99.0 | ["networks"] |
| `collections` | 52.0 | ["TF Perturb-seq", "IGVF_catalog_v1.3", "IGVF_catalog_v2.0"] |
| `description` | 44.4 | ["Identify function of all transcription factors (TFs) in cardiomyocytes: Knockd |
| `donors.sex` | 98.4 | ["female"] |
| `donors.strain_background` | 18.5 | ["C57BL/6J (B6)", "129S1/SvImJ (129)", "A/J (AJ)"] |
| `donors.taxa` | 98.4 | ["Homo sapiens"] |
| `file_set_type` | 100.0 | ["intermediate analysis"] |
| `files.content_type` | 100.0 | ["alignments", "index", "cell by bin contact matrix"] |
| `functional_assay_mechanisms.term_name` | 0.1 | ["protein secretion"] |
| `input_file_sets.@type` | 0.0 |  |
| `input_file_sets.accession` | 100.0 | ["IGVFDS9196PDRB"] |
| `input_file_sets.file_set_type` | 100.0 | ["experimental data"] |
| `input_file_sets.preferred_assay_titles` | 0.0 |  |
| `input_for.accession` | 0.0 |  |
| `lab.title` | 100.0 | ["Chongyuan Luo, UCLA"] |
| `preferred_assay_titles` | 100.0 | ["snM3C-seq"] |
| `publications.@id` | 6.3 | ["/publications/b2f69813-20ff-441a-9737-27ede20e82e8/"] |
| `release_timestamp` | 100.0 | ["2026-09-02T21:59:31.220606+00:00"] |
| `samples.@id` | 99.1 | ["/in-vitro-systems/IGVFSM3665LIKW/"] |
| `samples.@type` | 99.1 | ["InVitroSystem", "Biosample", "Sample"] |
| `samples.cellular_sub_pool` | 60.7 | ["subpool_LW322", "subpool_LW310", "subpool_LW312"] |
| `samples.classifications` | 99.1 | ["reprogrammed cell specimen"] |
| `samples.construct_library_sets` | 0.0 |  |
| `samples.disease_terms.term_name` | 0.0 |  |
| `samples.modifications.modality` | 3.9 | ["interference", "interference", "interference"] |
| `samples.modifications.summary` | 0.0 |  |
| `samples.sample_terms.@id` | 99.1 | ["/sample-terms/CL_0000057/"] |
| `samples.sample_terms.term_name` | 99.1 | ["fibroblast"] |
| `samples.sorted_from` | 0.0 |  |
| `samples.summary` | 99.1 | ["Homo sapiens (female, European) fibroblast reprogrammed cell specimen induced  |
| `samples.targeted_sample_term.@id` | 8.2 | ["/sample-terms/EFO_0004905/"] |
| `samples.targeted_sample_term.term_name` | 8.2 | ["induced pluripotent stem cell"] |
| `samples.taxa` | 99.1 | ["Homo sapiens"] |
| `samples.time_post_change` | 0.0 |  |
| `samples.treatments.summary` | 3.3 | ["Treatment of heat exposure at 35 Celsius"] |
| `samples.treatments.treatment_term_name` | 3.3 | ["heat exposure"] |
| `status` | 100.0 | ["released"] |
| `summary` | 100.0 | ["snM3C-seq"] |
| `targeted_genes.symbol` | 1.8 | ["CYP2C19"] |
| `workflows.name` | 99.9 | ["Luo Methylation and Chromatin Conformation Mapping Workflow", "Luo snmCT-seq R |

## PredictionSet (n=433, pulled 2026-11-18T23:05:05+00:00)

| field | % populated | example |
|---|---:|---|
| `@id` | 100.0 | ["/prediction-sets/IGVFDS0000USUS/"] |
| `@type` | 100.0 | ["PredictionSet", "FileSet", "Item"] |
| `accession` | 100.0 | ["IGVFDS0000USUS"] |
| `assessed_genes.symbol` | 51.5 | ["NFIA"] |
| `associated_phenotypes.term_name` | 23.3 | ["coronary artery disease"] |
| `award.component` | 100.0 | ["predictive modeling"] |
| `collections` | 97.2 | ["IGVF_catalog_v1.0", "IGVF_catalog_v1.3", "IGVF_catalog_v1.1"] |
| `description` | 100.0 | ["Predicted variant effects on TF binding affinity using the NFIA SNP effect mat |
| `donors.sex` | 55.0 | ["unspecified"] |
| `donors.strain_background` | 0.0 |  |
| `donors.taxa` | 55.0 | ["Homo sapiens"] |
| `file_set_type` | 100.0 | ["variant TF binding effects"] |
| `files.content_type` | 100.0 | ["variant binding effects"] |
| `input_file_sets.@type` | 0.0 |  |
| `input_file_sets.accession` | 0.0 |  |
| `input_for.accession` | 0.0 |  |
| `lab.title` | 100.0 | ["Alan Boyle, UMich"] |
| `large_scale_gene_list` | 0.0 |  |
| `publications.@id` | 80.4 | ["/publications/8d89db25-e479-4486-84a0-e7f270a8ef99/"] |
| `release_timestamp` | 100.0 | ["2025-06-05T17:35:19.879417+00:00"] |
| `samples.@id` | 45.0 | ["/tissues/IGVFSM3142EQHC/", "/tissues/IGVFSM4359MDYJ/", "/tissues/IGVFSM9259FDT |
| `samples.@type` | 45.0 | ["Tissue", "Biosample", "Sample"] |
| `samples.cellular_sub_pool` | 0.7 | ["igvf10_10x1"] |
| `samples.classifications` | 45.0 | ["tissue/organ", "tissue/organ", "tissue/organ"] |
| `samples.construct_library_sets` | 0.0 |  |
| `samples.disease_terms.term_name` | 0.0 |  |
| `samples.modifications.modality` | 1.6 | ["interference", "interference", "interference"] |
| `samples.modifications.summary` | 0.0 |  |
| `samples.sample_terms.@id` | 45.0 | ["/sample-terms/UBERON_0006920/", "/sample-terms/UBERON_0001114/", "/sample-term |
| `samples.sample_terms.term_name` | 45.0 | ["esophagus squamous epithelium", "right lobe of liver", "anterior lingual gland |
| `samples.sorted_from` | 0.0 |  |
| `samples.summary` | 45.0 | ["virtual Homo sapiens esophagus squamous epithelium tissue/organ", "virtual Hom |
| `samples.targeted_sample_term.@id` | 0.9 | ["/sample-terms/CL_0000115/", "/sample-terms/CL_0000115/", "/sample-terms/CL_000 |
| `samples.targeted_sample_term.term_name` | 0.9 | ["endothelial cell", "endothelial cell", "endothelial cell"] |
| `samples.taxa` | 45.0 | ["Homo sapiens", "Homo sapiens", "Homo sapiens"] |
| `samples.time_post_change` | 0.0 |  |
| `samples.treatments.summary` | 1.2 | ["Treatment of 100 ng/mL phorbol 13-acetate 12-myristate for 4 hours", "Treatmen |
| `samples.treatments.treatment_term_name` | 1.2 | ["phorbol 13-acetate 12-myristate", "T-cell surface glycoprotein CD3 gamma chain |
| `scope` | 30.0 | ["genome-wide"] |
| `small_scale_gene_list.symbol` | 0.2 | ["BAP1", "BARD1", "BRCA1"] |
| `software.title` | 100.0 | ["SEMVAR"] |
| `status` | 100.0 | ["released"] |
| `summary` | 100.0 | ["variant TF binding effects prediction for NFIA using SEMVAR in virtual Homo sa |
| `targeted_genes.symbol` | 0.0 |  |
