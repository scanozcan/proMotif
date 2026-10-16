# data-raw/fantom5

Builds the folder `fantom5_promoter_tpm_hg38/` used by `add_promoter_usage()`
(promoter usage across cell types, normal vs cancer). The folder (~40 MB, one
`.rds` per chromosome and strand plus `samples.rds`) is **not shipped** with the
package. Build it once and point proMotif to it:

```r
options(promotif.fantom5_dir = "/path/to/fantom5_promoter_tpm_hg38")
# or copy it to tools::R_user_dir("proMotif", "data")
```

## Inputs

FANTOM5 hg38 reprocessed CAGE peaks, phase 1+2 (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/):

1. `CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz` (RLE-normalised TPM)
2. `CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz` (hg38 coordinates and summits; the
   expression table's peak ids carry hg19 coordinates)
3. `HumanSamples2.0.sdrf.xlsx` (sample sheet; https://fantom.gsc.riken.jp/5/datafiles/latest/basic/)
4. Ensembl BioMart (GRCh38), all transcripts of protein-coding genes, columns:
   gene stable ID, gene name, chromosome, strand, transcript start, transcript end
   (tab-separated, with header).

## Build

```bash
python3 build_fantom5_promoter_usage.py \
  hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz HumanSamples2.0.sdrf.xlsx \
  ensembl_protein_coding_transcripts.tsv fantom5_promoter_tpm_hg38.tsv.gz \
  hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
Rscript make_rds.R fantom5_promoter_tpm_hg38.tsv.gz fantom5_promoter_tpm_hg38
```

Kept: CAGE peaks within 1 kb (same strand) of an annotated protein-coding TSS;
libraries of the categories primary cells, tissues and cell lines (1,023
libraries -> 608 sample types; donors / replicates averaged). TPM is stored as
integer x 100.
