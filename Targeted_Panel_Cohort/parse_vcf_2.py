#!/usr/bin/env python3
"""
parse_vcf.py
Parses a Funcotator-annotated Mutect2 VCF and extracts per-variant fields.

Output columns:
  gene             Hugo gene symbol
  location         CHROM:POS
  chrom / pos      Split location fields
  ref / alt        Reference and alternate alleles
  mut_type         SNV / insertion / deletion / MNV / indel
  VAF              Allele fraction (FORMAT AF)
  alt_reads        Alt-supporting reads (FORMAT AD)
  ref_reads        Ref reads (FORMAT AD)
  total_depth      Total depth (FORMAT DP)
  consequence      Simplified consequence (missense/nonsense/silent/etc.)
  raw_consequence  Full Funcotator variantClassification string
  FILTER           VCF FILTER field

Usage:
  python parse_vcf.py input.funcotated.vcf [output.tsv]
"""

import sys
import csv
import re


# ---------------------------------------------------------------------------
# Mutation type classification
# ---------------------------------------------------------------------------

def classify_mutation_type(ref, alts):
    types = []
    for alt in alts:
        if len(ref) == 1 and len(alt) == 1:
            types.append("SNV")
        elif len(ref) == 1 and len(alt) > 1:
            types.append("insertion")
        elif len(ref) > 1 and len(alt) == 1:
            types.append("deletion")
        elif len(ref) != len(alt):
            types.append("indel")
        else:
            types.append("MNV")
    return ",".join(types)


# ---------------------------------------------------------------------------
# FORMAT/SAMPLE field parser
# ---------------------------------------------------------------------------

def parse_format(fmt_str, sample_str):
    return dict(zip(fmt_str.split(":"), sample_str.split(":")))


# ---------------------------------------------------------------------------
# Funcotation field name extraction from ##INFO header
# ---------------------------------------------------------------------------

def parse_funcotation_fields(meta_lines):
    """Return ordered list of Funcotation field names from the ##INFO header."""
    for line in meta_lines:
        if "##INFO=<ID=FUNCOTATION" in line:
            m = re.search(r'Funcotation fields are: ([^"]+)"', line)
            if m:
                return [f.strip() for f in m.group(1).split("|")]
    return None


# ---------------------------------------------------------------------------
# FUNCOTATION INFO field parser
#
# The annotation is encoded as: FUNCOTATION=[field1|field2|...|fieldN]
# Complications:
#   - Multiple transcripts are comma-separated inside the brackets
#     (commas within field values are URL-encoded as %2C, so splitting
#     on literal comma is safe at the top level)
#   - ']' can appear inside URL-encoded content, so use a greedy regex
#     anchored to the semicolon or end-of-field that follows
# ---------------------------------------------------------------------------

def extract_funcotation(info_str, field_names):
    """
    Parse FUNCOTATION from INFO string.
    Returns a dict of annotation values.
    Takes the first transcript annotation if multiple are present.
    """
    empty = {
        "gene": "N/A", "consequence": "N/A", "transcript": "N/A",
        "strand": "N/A", "exon": "N/A", "transcript_pos": "N/A",
        "cDNA_change": "N/A", "codon_change": "N/A", "protein_change": "N/A",
    }

    m = re.search(r'FUNCOTATION=\[(.+)\](?:;|$)', info_str)
    if not m:
        return empty

    # Take first transcript (comma-separated at top level)
    first_tx = m.group(1).split(",")[0]
    parts = first_tx.split("|")

    if field_names and len(parts) >= 6:
        fd = dict(zip(field_names, parts))
        return {
            "gene":           fd.get("Gencode_34_hugoSymbol", "N/A").strip(),
            "consequence":    fd.get("Gencode_34_variantClassification", "N/A").strip(),
            "transcript":     fd.get("Gencode_34_annotationTranscript", "N/A").strip(),
            "strand":         fd.get("Gencode_34_transcriptStrand", "N/A").strip(),
            "exon":           fd.get("Gencode_34_transcriptExon", "N/A").strip(),
            "transcript_pos": fd.get("Gencode_34_transcriptPos", "N/A").strip(),
            "cDNA_change":    fd.get("Gencode_34_cDnaChange", "N/A").strip(),
            "codon_change":   fd.get("Gencode_34_codonChange", "N/A").strip(),
            "protein_change": fd.get("Gencode_34_proteinChange", "N/A").strip(),
        }
    else:
        # Positional fallback
        def get(i): return parts[i].strip() if i < len(parts) else "N/A"
        return {
            "gene":           get(0),
            "consequence":    get(5),
            "transcript":     get(12),
            "strand":         get(13),
            "exon":           get(14),
            "transcript_pos": get(15),
            "cDNA_change":    get(16),
            "codon_change":   get(17),
            "protein_change": get(18),
        }


# ---------------------------------------------------------------------------
# Consequence simplification map
# ---------------------------------------------------------------------------

CONSEQUENCE_MAP = {
    "SILENT":                 "silent",
    "SYNONYMOUS_VARIANT":     "silent",
    "MISSENSE":               "missense",
    "MISSENSE_VARIANT":       "missense",
    "NONSENSE":               "nonsense",
    "STOP_GAINED":            "nonsense",
    "NONSTOP_MUTATION":       "nonstop",
    "STOP_LOST":              "stop_lost",
    "FRAME_SHIFT_INS":        "frameshift",
    "FRAME_SHIFT_DEL":        "frameshift",
    "FRAMESHIFT_VARIANT":     "frameshift",
    "IN_FRAME_INS":           "inframe_indel",
    "IN_FRAME_DEL":           "inframe_indel",
    "SPLICE_SITE":            "splice_site",
    "SPLICE_REGION_VARIANT":  "splice_site",
    "RNA":                    "RNA",
    "LINCRNA":                "RNA",
    "IGR":                    "intergenic",
    "INTRON":                 "intronic",
    "FIVE_PRIME_UTR":         "5_prime_UTR",
    "THREE_PRIME_UTR":        "3_prime_UTR",
    "FIVE_PRIME_FLANK":       "5_prime_flank",
    "DE_NOVO_START_IN_FRAME": "start_codon",
    "START_CODON_SNP":        "start_codon",
    "START_CODON_DEL":        "start_codon",
    "START_CODON_INS":        "start_codon",
}

def simplify_consequence(raw):
    return CONSEQUENCE_MAP.get(raw.upper(), raw)


# ---------------------------------------------------------------------------
# Main VCF parser — expects a single Funcotator-annotated VCF
# ---------------------------------------------------------------------------

def parse_vcf(vcf_path):
    meta_lines = []
    funcotation_fields = None
    rows = []
    sample_name = "SAMPLE"

    import gzip
    opener = gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)

    with opener as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Collect metadata lines
            if line.startswith("##"):
                meta_lines.append(line)
                continue

            # Column header — parse Funcotation field names now that we have
            # all ## lines, then extract sample name
            if line.startswith("#CHROM"):
                funcotation_fields = parse_funcotation_fields(meta_lines)
                if funcotation_fields is None:
                    print("ERROR: No FUNCOTATION INFO header found. "
                          "Is this a Funcotator-annotated VCF?", file=sys.stderr)
                    sys.exit(1)
                cols = line.lstrip("#").split("\t")
                sample_name = cols[9] if len(cols) > 9 else "SAMPLE"
                continue

            # Data lines
            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom   = fields[0]
            pos     = fields[1]
            ref     = fields[3]
            alt_raw = fields[4]
            filt    = fields[6]
            info    = fields[7]
            fmt     = fields[8]
            sample  = fields[9] if len(fields) > 9 else ""

            alts     = alt_raw.split(",")
            mut_type = classify_mutation_type(ref, alts)

            fmt_dict = parse_format(fmt, sample)
            vaf      = fmt_dict.get("AF", ".")
            depth    = fmt_dict.get("DP", ".")
            ad       = fmt_dict.get("AD", ".")

            if ad != ".":
                ad_parts  = ad.split(",")
                ref_reads = ad_parts[0]
                alt_reads = ",".join(ad_parts[1:]) if len(ad_parts) > 1 else "."
            else:
                ref_reads = "."
                alt_reads = "."

            anno = extract_funcotation(info, funcotation_fields)
            raw_csq     = anno["consequence"]
            consequence = simplify_consequence(raw_csq)

            rows.append({
                "gene":            anno["gene"],
                "location":        f"{chrom}:{pos}",
                "chrom":           chrom,
                "pos":             pos,
                "ref":             ref,
                "alt":             alt_raw,
                "mut_type":        mut_type,
                "VAF":             vaf,
                "alt_reads":       alt_reads,
                "ref_reads":       ref_reads,
                "total_depth":     depth,
                "consequence":     consequence,
                "raw_consequence": raw_csq,
                "transcript":      anno["transcript"],
                "strand":          anno["strand"],
                "exon":            anno["exon"],
                "transcript_pos":  anno["transcript_pos"],
                "cDNA_change":     anno["cDNA_change"],
                "codon_change":    anno["codon_change"],
                "protein_change":  anno["protein_change"],
                "FILTER":          filt,
            })

    return rows, sample_name


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    vcf_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "parsed_variants.tsv"

    rows, sample_name = parse_vcf(vcf_path)

    fieldnames = [
        "gene",
        "location", "chrom", "pos", "ref", "alt",
        "mut_type",
        "VAF", "alt_reads", "ref_reads", "total_depth",
        "consequence", "raw_consequence",
        "transcript", "strand", "exon",
        "transcript_pos", "cDNA_change", "codon_change", "protein_change",
        "FILTER",
    ]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Sample:          {sample_name}")
    print(f"Variants parsed: {len(rows)}")
    print(f"Output:          {out_path}")

if __name__ == "__main__":
    main()
