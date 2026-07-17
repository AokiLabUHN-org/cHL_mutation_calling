# cHL_mutation_calling/WGS/SNV_Snakemake
Snakemake pipeline to call SNVs using Mutect2 in tumour-only mode and perform functional annotation (Funcotator) of variants. Performs the following functions:

1. Creates symlinks for read-only controlled bam files in a writeable directory (enables indexing within same dir)

2. Indexes read-only bam files and outputs next to symlinks

3. Filters out already identified read duplicates using the samtools flag -F 0x100 and pipes directly into GATK BQSR to avoid additional intermediate files

4. Runs Mutect2 calling in tumour-only mode

5. Filters out Mutect2 calls in contigs present in the original alignments that are absent in our reference files on H4H (EBV sequences)

6. Filters low-quality mutation calls and performs Funcotator annotation

Note: the first three steps pipeline should ideally be revised to remap the input BAMs to a standard GRCh38 reference and perform the standard GATK best practices pipeline (Picard MarkDuplicates -> BQSR). This requires large amounts of storage on H4H for all intermediate files or exteremly long runtimes to perform all steps in one run via piping output without writing files.


