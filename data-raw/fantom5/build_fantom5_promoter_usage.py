#!/usr/bin/env python3
"""Build the compact FANTOM5 promoter-usage table used by proMotif's promoter
usage layer (add_promoter_usage / plot_promoter_usage).

Source: FANTOM5 hg38 reprocessed CAGE peaks, RLE-normalised TPM, phase 1+2
(hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz, ~1,800 libraries) and
the sample sheet HumanSamples2.0.sdrf.xlsx (category, description per library).

Kept: CAGE peaks within 1 kb (same strand) of any annotated TSS of a human
protein-coding gene (Ensembl GRCh38, all transcripts of protein-coding genes);
libraries in the categories 'primary cells', 'tissues' and 'cell lines'. Donor
/ pool / replicate libraries of the same sample type are averaged.

Peak positions: hg38 coordinates and CAGE summit from
hg38_fair+new_CAGE_peaks_phase1and2.bed.gz (the expression table's peak ids
carry hg19 coordinates).

Usage: python3 build_fantom5_promoter_usage.py <tpm_ann.osc.txt.gz> <HumanSamples2.0.sdrf.xlsx> \
           <ensembl_protein_coding_transcripts.tsv> <out.tsv.gz> <CAGE_peaks.bed.gz>
"""
import sys, gzip, re, bisect, urllib.parse, collections
import numpy as np, openpyxl

tpm_f, sdrf_f, ens_f, out_f, bed_f = sys.argv[1:6]
peaks = {}
with gzip.open(bed_f, "rt") as fh:
    for l in fh:
        f = l.rstrip("\n").split("\t")
        peaks[f[3]] = (f[0], int(f[1]) + 1, int(f[2]), f[5], int(f[7]))   # chr, start, end, strand, summit
KEEP = {"primary cells", "tissues", "cell lines"}
FLANK = 1000

# ---- samples ------------------------------------------------------------------
rows = list(openpyxl.load_workbook(sdrf_f, read_only=True).active.iter_rows(values_only=True))[1:]
meta = {r[0]: (r[2] or "", r[4] or "") for r in rows}
def sample_type(desc):
    d = re.sub(r",\s*(donor|pool|biol_rep|tech_rep|rep|replicate|lot)\s*\S*.*$", "", desc, flags=re.I)
    return re.sub(r"\s+", " ", d).strip(" ,")

# ---- protein-coding TSSs ----------------------------------------------------------
tss = collections.defaultdict(list)
for i, l in enumerate(open(ens_f)):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    strand = "+" if f[3] == "1" else "-"
    tss[("chr" + f[2], strand)].append(int(f[4]) if strand == "+" else int(f[5]))
tss = {k: sorted(set(v)) for k, v in tss.items()}
def near(chrom, strand, pos):
    a = tss.get((chrom, strand))
    if not a: return False
    i = bisect.bisect_left(a, pos - FLANK)
    return i < len(a) and a[i] <= pos + FLANK

# ---- stream the expression table ---------------------------------------------------
with gzip.open(tpm_f, "rt") as fh:
    for line in fh:
        if line.startswith("##"): continue
        header = line.rstrip("\n").split("\t"); break
    cols = header[7:]
    lib = [c.split(".")[-3] for c in cols]                 # e.g. 13541-145H4
    types, cats, keep_idx = [], [], []
    for j, s in enumerate(lib):
        desc, cat = meta.get(s, ("", ""))
        if cat in KEEP:
            keep_idx.append(j); types.append(sample_type(desc)); cats.append(cat)
    utypes = sorted(set(zip(types, cats)))
    tindex = {t: k for k, t in enumerate(utypes)}
    grp = np.array([tindex[(t, c)] for t, c in zip(types, cats)])
    nlib = np.bincount(grp, minlength=len(utypes))
    keep_idx = np.array(keep_idx) + 7
    print(f"libraries kept: {len(keep_idx)} -> {len(utypes)} sample types", flush=True)

    out = gzip.open(out_f, "wt")
    out.write("\t".join(["chr", "summit", "start", "end", "strand", "gene_hint"] +
                        [f"{t} [{c}]" for t, c in utypes]) + "\n")
    n_in = n_out = 0
    for line in fh:
        if line.startswith("01STAT") or line.startswith("02STAT"): continue
        n_in += 1
        f = line.rstrip("\n").split("\t")
        pk = peaks.get(f[0])
        if pk is None: continue
        chrom, s1, e, strand, pos = pk                      # hg38; pos = CAGE summit
        if not near(chrom, strand, pos): continue
        v = np.array(f)[keep_idx].astype(float)
        agg = np.bincount(grp, weights=v, minlength=len(utypes)) / nlib
        if agg.max() < 0.5: continue                         # never expressed
        gh = f[1].split("@")[-1] if "@" in f[1] else ""
        out.write("\t".join([chrom, str(pos), str(s1), str(e), strand, gh] +
                            [f"{x:.2f}".rstrip("0").rstrip(".") if x else "0" for x in agg]) + "\n")
        n_out += 1
        if n_in % 50000 == 0: print(f"peaks read {n_in}, kept {n_out}", flush=True)
    out.close()
print(f"done: {n_out} peaks near protein-coding TSSs, {len(utypes)} sample types")
