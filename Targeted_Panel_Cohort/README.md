# cHL_mutation_calling/Targeted_Panel_Cohort
Code used to call SNVs and CNAs in blood, LMD or flow-enriched healthy immune cells (TME samples), or LMD/flow-enriched HRS cells using samples from the cohort described here: https://pubmed.ncbi.nlm.nih.gov/41615883/

Builds on previous work from Dr. Gerben Duns (BCCancer), using his Snakemake pipeline for alignment and quality control of targeted sequencing data. All required resources have been updated to point to filepaths on the H4H compute cluster. The CNA calling also leverages his CNVkit pipeline. 

The SNV calling pipeline (part of the same Snakemake pipeline) has been updated to run in tumour-only mode for all samples (HRS, blood, TME). Performs Mutect2 varaint calling, variant filtering, and Funcotator annotation of variants.

Additionally, added in a python script to parse the Funcotated vcfs to create .tsv files with all information required to identify CHIP drivers and shared mutations between HRS and germline samples.

