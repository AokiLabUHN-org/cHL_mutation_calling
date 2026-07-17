# generate_phased_snps.sh
#
# Produces a population-phased BCF of common (MAF>1%) SNPs for a cohort of
# BAM/CRAM files, suitable for use as --phased-bcf in run_call_mCAs_pipeline.sh.
#
# This follows the same pattern used to generate calling inputs for UK Biobank
# (see UKB/generate_impute_ref and UKB/imputationQC): genotype the cohort at
# the common sites of a population reference panel, then statistically phase
# those genotypes against the same panel. The UKB scripts do this via an
# internal SNP-array scaffold + IMPUTE5/xcftools on the RAP; this script
# instead genotypes directly from BAM/CRAM with bcftools and phases with
# SHAPEIT5, so it only depends on tools an outside user can install.
#
# Steps:
#   1. For each chromosome, restrict the reference panel to common (MAF>1%)
#      sites using ../bin/vcf_extract_common_AC_AN (same tool used in
#      UKB/generate_impute_ref/run_generate_reference_RAP.sh), and combine
#      all chromosomes' sites into one genome-wide site list.
#   2. Genotype every cohort BAM/CRAM at those sites (bcftools mpileup +
#      call), once per sample across all requested chromosomes in a single
#      streamed pass. If --remap is given, each sample is first realigned to
#      --ref-file (see "Remapping" below) as part of that same stream.
#   3. Merge the per-sample genotype calls into one genome-wide cohort BCF.
#   4. Per chromosome: subset the cohort BCF to that chromosome and run
#      SHAPEIT5 to population-phase it against the reference panel + genetic
#      map for that chromosome.
#
# Remapping (--remap):
#   Use this when the input BAM/CRAM files were aligned to a different
#   reference than --ref-file. Each sample is streamed through
#   samtools collate -> samtools fastq -> bwa mem -> samtools sort -> bcftools
#   mpileup/call as a single pipe, so no intermediate FASTQ/SAM/BAM ever
#   touches disk as a named file. Realignment is the most expensive step in
#   this script, so it is done exactly once per sample (genotyping is done
#   genome-wide in one pass, not once per chromosome) rather than once per
#   sample per chromosome.
#   Requires a pre-built bwa index (--bwa-index PREFIX, defaults to the
#   GRCh38 index at /cluster/tools/data/genomes/human/GRCh38/iGenomes/
#   Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64) and the `bwa` (or
#   `--bwa PATH`) binary available. The index prefix does not need to match
#   --ref-file's filename -- they're independent (bcftools needs the fasta,
#   bwa needs its own amb/ann/bwt/pac/sa index).
#   Note: samtools collate/sort spill to disk internally for large inputs
#   (external-sort behavior); their spill-file prefixes are explicitly
#   pointed at $TMP_DIR (see -t/--tmp-dir / TMP_DIR_BASE) rather than left to
#   default to the current directory, since neither honors $TMPDIR itself.
#
# Limitations (kept out of scope to stay portable/simple):
#   - chrX PAR1/PAR2 splitting and ploidy-aware calling are not handled
#     automatically. Pass --sex-file to get diploid/haploid calling right for
#     chrX/chrY; otherwise everything is called as diploid.
#   - Only biallelic SNPs from the reference panel are considered (indels/SVs
#     in the panel are dropped implicitly by downstream tools expecting SNPs).

set -euo pipefail

usage() {
    cat <<-USAGE
Usage: $0 [options]

Generates a population-phased BCF of common SNPs per chromosome from a set
of BAM/CRAM files, for use as --phased-bcf in run_call_mCAs_pipeline.sh.

Any required option not supplied on the command line will be prompted for
interactively.

Required:
  -m, --manifest FILE        Manifest of BAM/CRAM files: 'ID PATH' per line
  -o, --out-dir DIR          Output directory
  -n, --name NAME            Output file prefix
  --ref-file FILE            Reference fasta (GRCh38), matching the build of
                              --ref-panel and used for bcftools mpileup.
                              Default:
                              /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
  --ref-panel PATTERN        Path/URL to phased reference panel VCF/BCF per
                              chromosome. Use {CHR} as a placeholder for the
                              chromosome. Default (local 1000G high-coverage
                              copy, only chr13/chr21 present as of writing --
                              pass --ref-panel to use a URL or another chrom):
                              /cluster/projects/aokigroup/Stueckmann/input/GRCh38/phased_1000genomes/1kGP_high_coverage_Illumina.{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  --genetic-map PATTERN      Path to SHAPEIT5-format genetic map per
                              chromosome. Use {CHR} as a placeholder. Default:
                              /cluster/projects/aokigroup/Stueckmann/input/GRCh38/shapeit5_resources/{CHR}.b38.gmap.gz
  --shapeit5 PATH            Path to the shapeit5_phase_common (or
                              phase_common_static) binary

Options:
  -c, --chrs LIST            Comma-separated chromosomes (default: chr1..chr22)
  --sex-file FILE            'ID SEX(M/F)' per line, used for ploidy-aware
                              calling on chrX/chrY (passed to bcftools call
                              --ploidy GRCh38 -S). Omit for autosomes only.
  --remap                    Realign input BAM/CRAM files to --ref-file
                              before genotyping (see "Remapping" above). Use
                              this when the inputs were aligned to a
                              different reference than --ref-file.
  --bwa PATH                 Path to the bwa binary (default: bwa on PATH).
                              Only used with --remap.
  --bwa-index PREFIX          Prefix of the bwa index to realign against
                              (files PREFIX.{amb,ann,bwt,pac,sa} must exist).
                              Only used with --remap; does not need to match
                              --ref-file's name. Default:
                              /cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64
  --old-ref-file FILE        Reference the input files were originally
                              aligned to. Only needed with --remap when an
                              input is CRAM (required to decode CRAM
                              reference-compressed sequence); ignored for BAM.
  -t, --tmp-dir DIR          Temporary directory (default: autogenerated
                              under /cluster/projects/vannergroup/Stueckmann/tmp)
  -j, --cores N               Total cores available (default: nproc)
  -J, --samples-parallel N    How many samples to genotype concurrently
                              (default: 1). Each sample's collate/fastq/bwa
                              mem/sort gets floor(--cores / this) threads, so
                              the two multiply out to roughly --cores total
                              -- raising this without lowering per-sample
                              threads oversubscribes the machine. Leave at 1
                              to put all cores into one sample at a time
                              (best when the manifest is small); raise it to
                              trade per-sample speed for cross-sample
                              throughput on large manifests.
  -k, --keep-tmp              Keep tmp dir after success
  -v, --verbose                Verbose output
  -h, --help                   Show this message
USAGE
}

# Defaults
CHR_LIST=""
SEX_FILE=""
TMP_DIR=""
TMP_DIR_BASE="/cluster/projects/vannergroup/Stueckmann/tmp"
CORES=$(nproc)
SAMPLES_PARALLEL=1
KEEP_TMP=false
VERBOSE=false
REMAP=false
BWA_BIN="bwa"
BWA_INDEX="/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64"
REF_FILE="/cluster/tools/data/genomes/human/GRCh38/iGenomes/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
GENETIC_MAP_PATTERN="/cluster/projects/aokigroup/Stueckmann/input/GRCh38/shapeit5_resources/{CHR}.b38.gmap.gz"
REF_PANEL_PATTERN="/cluster/projects/aokigroup/Stueckmann/input/GRCh38/phased_1000genomes/1kGP_high_coverage_Illumina.{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
OLD_REF_FILE=""
SHAPEIT5_BIN="/cluster/home/dstueckm/miniconda3/bin/SHAPEIT5_phase_common"

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--manifest) MANIFEST="$2"; shift 2;;
        -o|--out-dir) OUT_DIR="$2"; shift 2;;
        -n|--name) PFX="$2"; shift 2;;
        --ref-file) REF_FILE="$2"; shift 2;;
        --ref-panel) REF_PANEL_PATTERN="$2"; shift 2;;
        --genetic-map) GENETIC_MAP_PATTERN="$2"; shift 2;;
        --shapeit5) SHAPEIT5_BIN="$2"; shift 2;;
        -c|--chrs) CHR_LIST="$2"; shift 2;;
        --sex-file) SEX_FILE="$2"; shift 2;;
        --remap) REMAP=true; shift 1;;
        --bwa) BWA_BIN="$2"; shift 2;;
        --bwa-index) BWA_INDEX="$2"; shift 2;;
        --old-ref-file) OLD_REF_FILE="$2"; shift 2;;
        -t|--tmp-dir) TMP_DIR="$2"; shift 2;;
        -j|--cores) CORES="$2"; shift 2;;
        -J|--samples-parallel) SAMPLES_PARALLEL="$2"; shift 2;;
        -k|--keep-tmp) KEEP_TMP=true; shift 1;;
        -v|--verbose) VERBOSE=true; shift 1;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option $1" >&2; usage; exit 1;;
    esac
done

# Prompt interactively for anything required that wasn't passed as a flag.
prompt_if_missing() {
    local var_name="$1" prompt_text="$2"
    if [[ -z "${!var_name:-}" ]]; then
        if [[ ! -t 0 ]]; then
            echo "Error: $var_name is not set and input is not a terminal (use --help for flags)" >&2
            exit 1
        fi
        local value
        read -r -p "$prompt_text: " value
        printf -v "$var_name" '%s' "$value"
    fi
}

prompt_if_missing MANIFEST "Path to BAM/CRAM manifest (ID PATH per line)"
prompt_if_missing OUT_DIR "Output directory"
prompt_if_missing PFX "Output file prefix"
prompt_if_missing REF_FILE "Path to reference fasta (GRCh38)"
prompt_if_missing REF_PANEL_PATTERN "Reference panel VCF/BCF path/URL pattern (use {CHR} as placeholder)"
prompt_if_missing GENETIC_MAP_PATTERN "Genetic map path pattern for SHAPEIT5 (use {CHR} as placeholder)"
prompt_if_missing SHAPEIT5_BIN "Path to shapeit5_phase_common(_static) binary"
if [[ -z "$CHR_LIST" ]]; then
    if [[ -t 0 ]]; then
        read -r -p "Chromosomes to phase, comma-separated [chr1..chr22]: " CHR_LIST
    fi
    CHR_LIST=${CHR_LIST:-$(seq -f 'chr%g' 1 22 | paste -sd,)}
fi

[[ -f "$MANIFEST" ]] || { echo "Error: manifest not found: $MANIFEST" >&2; exit 1; }
[[ -f "$REF_FILE" ]] || { echo "Error: reference fasta not found: $REF_FILE" >&2; exit 1; }
[[ -x "$SHAPEIT5_BIN" ]] || command -v "$SHAPEIT5_BIN" >/dev/null 2>&1 \
    || { echo "Error: shapeit5 binary not found or not executable: $SHAPEIT5_BIN" >&2; exit 1; }

if $REMAP; then
    [[ -x "$BWA_BIN" ]] || command -v "$BWA_BIN" >/dev/null 2>&1 \
        || { echo "Error: bwa binary not found or not executable: $BWA_BIN" >&2; exit 1; }
    [[ -f "${BWA_INDEX}.bwt" ]] \
        || { echo "Error: no bwa index found at prefix $BWA_INDEX (missing ${BWA_INDEX}.bwt; use --bwa-index or run: bwa index -p $BWA_INDEX $REF_FILE)" >&2; exit 1; }
fi

mkdir -p "$OUT_DIR"

# Setup TMP_DIR
if [[ -z "$TMP_DIR" ]]; then
    mkdir -p "$TMP_DIR_BASE"
    TMP_DIR=$(mktemp -d "$TMP_DIR_BASE/generate_phased_snps.XXXXXX")
    CLEAN_TMP=true
else
    mkdir -p "$TMP_DIR"
    CLEAN_TMP=false
fi
if $VERBOSE; then echo "$(date): using TMP_DIR=$TMP_DIR"; fi

trap_to_cleanup() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "$(date): script exited with code $exit_code. Temporary dir: $TMP_DIR" >&2
    elif $KEEP_TMP || [[ "$CLEAN_TMP" == false ]]; then
        echo "$(date): preserving tmp dir: $TMP_DIR"
    else
        rm -rf "$TMP_DIR"
        if $VERBOSE; then echo "$(date): removed tmp dir"; fi
    fi
}
trap trap_to_cleanup EXIT

# Check dependencies
missing=()
for cmd in bcftools samtools parallel awk cut; do
    if ! command -v $cmd >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done
if [[ ${#missing[@]} -ne 0 ]]; then
    echo "Error: missing required commands: ${missing[*]}" >&2
    exit 1
fi

VCF_EXTRACT_COMMON=../bin/vcf_extract_common_AC_AN
[[ -x "$VCF_EXTRACT_COMMON" ]] || { echo "Error: $VCF_EXTRACT_COMMON not found/executable (run from pipeline/ or rebuild bin/)" >&2; exit 1; }

substitute_chr() {
    # replaces {CHR} in a path/URL pattern with the given chromosome
    local pattern="$1" chr="$2"
    echo "${pattern//\{CHR\}/$chr}"
}

ensure_indexed() {
    # bcftools view -r needs a .tbi/.csi index; create one for local files
    # that don't have it yet (remote http/s3 files are left alone -- we
    # can't index those ourselves)
    local file="$1"
    case "$file" in
        http://*|https://*|ftp://*|s3://*) return 0;;
    esac
    if [[ ! -f "${file}.tbi" && ! -f "${file}.csi" ]]; then
        echo "$(date): no index found for $file, indexing it now"
        bcftools index -t -f "$file" \
            || { echo "Error: failed to index $file (check write permissions on its directory, or pre-index it manually: bcftools index -t $file)" >&2; exit 1; }
    fi
}

echo "$(date): === preparing per-chromosome reference panel sites ==="
SITE_FILES=()
for CHR in ${CHR_LIST//,/ }; do
    echo "$(date): [$CHR] restricting reference panel to common (MAF>1%) sites"
    CHR_TMP="$TMP_DIR/$CHR"
    mkdir -p "$CHR_TMP"

    REF_PANEL=$(substitute_chr "$REF_PANEL_PATTERN" "$CHR")
    ensure_indexed "$REF_PANEL"

    # ensure AC/AN tags exist, then keep only common sites (mirrors
    # UKB/generate_impute_ref/add_AC_snp_scaffold_RAP_launch.sh +
    # run_generate_reference_RAP.sh)
    bcftools view "$REF_PANEL" -r "$CHR" -Ou \
        | bcftools +fill-tags -Oz -o "$CHR_TMP/refpanel.acan.vcf.gz" -- -t AC,AN
    "$VCF_EXTRACT_COMMON" "$CHR_TMP/refpanel.acan.vcf.gz" \
        | bcftools view -Ob -o "$CHR_TMP/refpanel.common.bcf"
    bcftools index -f "$CHR_TMP/refpanel.common.bcf"

    bcftools view -G "$CHR_TMP/refpanel.common.bcf" -Ob -o "$CHR_TMP/sites.bcf"
    bcftools index -f "$CHR_TMP/sites.bcf"
    SITE_FILES+=("$CHR_TMP/sites.bcf")
done

echo "$(date): combining sites across all requested chromosomes"
ALL_SITES="$TMP_DIR/sites.allchr.bcf"
bcftools concat -a "${SITE_FILES[@]}" -Ob -o "$ALL_SITES"
bcftools index -f "$ALL_SITES"

echo "$(date): genotyping cohort BAM/CRAM files at reference panel sites (one streamed pass per sample, all chromosomes)"
CALL_ARGS="-m"
if [[ -n "$SEX_FILE" ]]; then
    CALL_ARGS="$CALL_ARGS --ploidy GRCh38 -S $SEX_FILE"
fi

# THREADS_PER_SAMPLE * SAMPLES_PARALLEL ~= CORES, so raising --samples-parallel
# doesn't oversubscribe: each concurrent sample gets a proportionally smaller
# thread count for collate/fastq/bwa mem/sort.
THREADS_PER_SAMPLE=$(( CORES / SAMPLES_PARALLEL ))
if [[ "$THREADS_PER_SAMPLE" -lt 1 ]]; then THREADS_PER_SAMPLE=1; fi
# samtools sort's temp-spill-file naming has only 3 fallback attempts
# (name, name-001, name-002) when multiple threads race to flush a block
# under the same sequential number; at very high -@ counts, more than 3
# threads can collide on one block and exhaust that fallback, silently
# dropping the block ("File exists" errors, corrupt/incomplete sort). Sort
# also has much flatter thread scaling than bwa mem, so there's little
# upside to giving it as many threads as alignment anyway -- cap it.
SORT_THREADS_MAX=8
SORT_THREADS=$THREADS_PER_SAMPLE
if [[ "$SORT_THREADS" -gt "$SORT_THREADS_MAX" ]]; then SORT_THREADS=$SORT_THREADS_MAX; fi
if $VERBOSE; then
    echo "$(date): running $SAMPLES_PARALLEL sample(s) concurrently, $THREADS_PER_SAMPLE thread(s) each for collate/fastq/bwa (of $CORES total cores), $SORT_THREADS for sort"
fi

genotype_sample() {
    local id="$1" path="$2"
    local out="$TMP_DIR/$id.bcf"
    local collate_opts=()
    [[ -n "$OLD_REF_FILE" ]] && collate_opts=(--reference "$OLD_REF_FILE")

    if $REMAP; then
        # fully streamed remap: collate (name-group) -> fastq -> bwa mem ->
        # coordinate sort -> mpileup/call, no named intermediate file
        samtools collate -Ou "${collate_opts[@]}" -@ "$THREADS_PER_SAMPLE" "$path" "$TMP_DIR/$id.collate_tmp" \
            | samtools fastq -@ "$THREADS_PER_SAMPLE" - \
            | "$BWA_BIN" mem -t "$THREADS_PER_SAMPLE" -p "$BWA_INDEX" - \
            | samtools sort -@ "$SORT_THREADS" -T "$TMP_DIR/$id.sort_tmp" -O bam -u - \
            | bcftools mpileup -f "$REF_FILE" -R "$ALL_SITES" -a AD -Ou - \
            | bcftools call $CALL_ARGS -Ob -o "$out"
    else
        bcftools mpileup -f "$REF_FILE" -R "$ALL_SITES" -a AD -Ou "$path" \
            | bcftools call $CALL_ARGS -Ob -o "$out"
    fi
    bcftools index -f "$out"
}
export -f genotype_sample
export REF_FILE CALL_ARGS ALL_SITES REMAP BWA_BIN BWA_INDEX OLD_REF_FILE THREADS_PER_SAMPLE SORT_THREADS TMP_DIR

cat "$MANIFEST" | tr '\t' ' ' \
    | parallel --halt-on-error 2 --colsep ' ' --joblog /dev/stderr --retries 4 --jobs "$SAMPLES_PARALLEL" \
        genotype_sample {1} {2}

echo "$(date): merging cohort genotype calls across all chromosomes"
awk '{print $1}' "$MANIFEST" | while read -r id; do echo "$TMP_DIR/$id.bcf"; done > "$TMP_DIR/merge_list.txt"
bcftools merge -m none --file-list "$TMP_DIR/merge_list.txt" -Ob -o "$TMP_DIR/cohort.allchr.bcf"
bcftools index -f "$TMP_DIR/cohort.allchr.bcf"

echo "$(date): === phasing each chromosome ==="
for CHR in ${CHR_LIST//,/ }; do
    CHR_TMP="$TMP_DIR/$CHR"
    GENETIC_MAP=$(substitute_chr "$GENETIC_MAP_PATTERN" "$CHR")

    echo "$(date): [$CHR] subsetting cohort genotypes"
    bcftools view -r "$CHR" "$TMP_DIR/cohort.allchr.bcf" -Ob -o "$CHR_TMP/cohort.bcf"
    bcftools index -f "$CHR_TMP/cohort.bcf"

    echo "$(date): [$CHR] population-phasing cohort genotypes with SHAPEIT5"
    "$SHAPEIT5_BIN" \
        --input "$CHR_TMP/cohort.bcf" \
        --reference "$CHR_TMP/refpanel.common.bcf" \
        --map "$GENETIC_MAP" \
        --region "$CHR" \
        --thread "$CORES" \
        --output "$OUT_DIR/$PFX.$CHR.phased.bcf"
    bcftools index -f "$OUT_DIR/$PFX.$CHR.phased.bcf"

    rm -rf "$CHR_TMP"
    echo "$(date): [$CHR] done -> $OUT_DIR/$PFX.$CHR.phased.bcf"
done

echo "$(date): finished. Phased BCFs written to $OUT_DIR/$PFX.<chr>.phased.bcf"
echo "$(date): pass one of these as --phased-bcf to run_call_mCAs_pipeline.sh for the matching chromosome"
