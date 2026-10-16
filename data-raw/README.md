# data-raw

`make_promoter_background.py` builds `inst/extdata/promoter_background_JASPAR2026_hg38.tsv`,
the promoter background used by discovery mode (`rank_by = "background"`).

Inputs:

1. `canon.tsv` - Ensembl BioMart (GRCh38), dataset `hsapiens_gene_ensembl`, filters
   `biotype = protein_coding` and `transcript_is_canonical`, attributes: gene stable ID,
   gene name, chromosome, strand, transcript start, transcript end, transcript stable ID.
2. `JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt` from
   https://jaspar.elixir.no/download/data/2026/CORE/

```bash
python3 make_promoter_background.py canon.tsv \
  JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt 2026 all
```

All 20,094 protein-coding genes on chromosomes 1-22, X, Y are used (pass a number
instead of `all` to sample, seed 1); promoter sequences
(-2500/+500 bp around the canonical TSS) come from the UCSC API (hg38). Every
matrix is converted exactly as `proMotif:::.pfm_to_logodds()` and its best-site
score (both strands) is recorded per promoter. Output: per matrix, upper-tail
quantiles `q_p` (P(best site >= q_p) = p) plus an exponential tail scale for
scores beyond the strongest background promoter. Scoring is resumable: each block of 500 promoters is saved to `chunks/`, finished
blocks are skipped on restart, and `WORKER=i NWORK=k` runs blocks in parallel.
Runtime ~1 h on 2 cores (numpy). Held-out check (2-fold): for motifs with
continuous scores, 4.9% / 0.93% / 0.07% of promoters fall at p <= 0.05 / 0.01 /
0.001; short motifs whose maximum score is common are conservative (ties).
