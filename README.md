<picture>
  <source media="(prefers-color-scheme: dark)" srcset="man/figures/logo-dark.svg">
  <img src="man/figures/logo-light.svg" alt="proMotif" width="440">
</picture>

# proMotif

[![R-CMD-check](https://github.com/scanozcan/proMotif/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scanozcan/proMotif/actions/workflows/R-CMD-check.yaml)

Scan a gene's promoter for **transcription-factor binding sites** using JASPAR
and HOCOMOCO motifs. You give it a **gene symbol** and a **TF**; proMotif looks
up the gene's transcripts from
Ensembl, fetches each promoter sequence, scans both strands for the TF's motif
from both databases, assigns each hit a statistical significance, and draws the
results per transcript and on genomic coordinates. An interactive Shiny app is
included.

It reuses the [GLproxScape](https://github.com/scanozcan/GLproxScape) package's
Ensembl / JASPAR / HOCOMOCO fetchers.

---

## Installation

From GitHub (installs GLproxScape automatically):

```r
# install.packages("remotes")
remotes::install_github("scanozcan/proMotif")
```

Two optional packages unlock extra features:

```r
install.packages("shiny")      # for the interactive app, launch_app()
install.packages("TFMPvalue")  # for the significance (p-value) plots
```

An **internet connection is required at run time** — Ensembl, JASPAR and
HOCOMOCO are queried live (HOCOMOCO caches a bundle after the first use).

---

## Quick start

```r
library(proMotif)

# scan the ATP7B promoter for MTF1 sites, write a CSV + plots to results/ATP7B_MTF1/
run_binding_site_analysis("ATP7B", "MTF1")

# or launch the point-and-click app
launch_app()
```

---

## The functions

proMotif has four functions. For most uses you only need the first or the app.

| Function | What it does |
|---|---|
| `run_binding_site_analysis(gene, tf, ...)` | One-shot: scan + plots + writes a CSV and PDF figures to disk. Start here. |
| `launch_app()` | Interactive Shiny app — type a gene + TF, view the plots, download the table. |
| `scan_binding_sites(gene, tf, ...)` | Just the scan (the slow, network part). Returns the data, no plots. |
| `plot_binding_sites(scan, ...)` | Builds the plots from a `scan_binding_sites()` result (fast, local). |

The last two let you scan once and re-plot many times (e.g. restyle without
re-querying the databases):

```r
s <- scan_binding_sites("TERT", "CTCF")       # slow: fetch + scan
p <- plot_binding_sites(s, base_size = 16)    # fast: rebuild plots
p$genomic                                     # show one plot
```

### Key arguments

| Argument | Default | Meaning |
|---|---|---|
| `gene`, `tf` | — | HGNC gene symbol and TF symbol (required) |
| `species` | `"homo_sapiens"` | Ensembl species token |
| `upstream`, `downstream` | `2500`, `500` | promoter window around each TSS (bp) |
| `threshold_frac` | `0.80` | PWM cutoff as a fraction of the matrix min–max score |
| `transcripts` | `"protein_coding"` | or `"all"` (see *Alternative promoters* below) |
| `max_transcripts` | `8` | cap on distinct TSS windows scanned (canonical always kept) |
| `hocomoco_version`, `hocomoco_species` | `"v12"`, `"human"` | HOCOMOCO bundle |
| `out_root` | `"results"` | output folder |
| `color_seed` | `NULL` | fix the two database colours reproducibly (NULL = fresh each run) |
| `base_size`, `point_size` | `12`, `2.6` | plot font and marker size |

---

## Outputs

`run_binding_site_analysis()` writes to `results/<GENE>_<TF>/`:

- **`<GENE>_<TF>_binding_sites.csv`** — one row per predicted site.
- **Score plots:** `*_canonical.pdf`, `*_by_transcript.pdf`, `*_genomic.pdf`.
- **Significance plots** (if `TFMPvalue` is installed): `*_canonical_significance.pdf`, `*_by_transcript_significance.pdf`, `*_genomic_significance.pdf`.

The three plot views:

- **canonical** — predicted sites in the canonical transcript's promoter, by position relative to the TSS.
- **by transcript** — one panel per transcript / alternative promoter (TSS).
- **genomic** — all sites on real chromosome coordinates, with each TSS marked and an arrow showing the transcription direction from the canonical TSS. A site that falls inside several transcripts' promoter windows is one genomic locus, so it is drawn once here (per database and strand) rather than once per transcript as in the per-transcript view.

The score plots put `score_frac` (fraction of the matrix's max score) on the
y-axis; the significance plots put `-log10(adjusted p)` with a dotted line at
p = 0.05.

### CSV columns

| column | meaning |
|---|---|
| `transcript_id`, `is_canonical` | transcript the window came from |
| `chr`, `gene_strand`, `tss_genomic` | locus and the TSS used |
| `database`, `matrix_id`, `matrix_name` | JASPAR or HOCOMOCO and the exact motif matrix |
| `position` | site position relative to the TSS (bp; negative = upstream) |
| `strand` | strand the motif matched (`+`/`-`) |
| `score`, `score_frac` | raw PWM log-odds score and its fraction of the matrix max |
| `pvalue`, `p_adj` | per-site p-value (GC-aware background) and Bonferroni-adjusted p (NA without `TFMPvalue`) |
| `genomic_position` | absolute chromosome coordinate of the site |
| `site` | the matched promoter sequence |

---

## The Shiny app

```r
launch_app()
```

Enter a gene and TF, adjust the promoter window / threshold, and click **Run
scan**. The scan runs once; the text/point-size slider and the colour controls
restyle the existing plots instantly (no re-scan). Tabs show the three score
plots, the three significance plots, and a table you can download as CSV.

---

## Notes

- **Alternative promoters.** Many genes have more than one TSS. The default
  `transcripts = "protein_coding"` often collapses to the main promoter because
  alternative-promoter isoforms are frequently non-coding biotypes. Use
  `transcripts = "all"` to include them — the console prints how many distinct
  TSS exist across all transcripts.
- **Resilience.** Promoter sequence comes from Ensembl; if Ensembl times out or
  errors, proMotif automatically falls back to the UCSC REST API.
- **Significance.** `p_adj` is the chance of seeing a match this strong under
  the promoter's own base composition — i.e. **match significance**, not proof
  of binding. A significant motif match is not the same as a TF actually
  binding there (that depends on chromatin, accessibility and cell type).
- **Tuning.** Lower `threshold_frac` (e.g. 0.75) to surface weaker matches;
  widen `upstream`/`downstream` for a larger window.

---

## Examples

```r
library(proMotif)

# 1. Basic run: writes the CSV + PDFs to results/ATP7B_MTF1/
run_binding_site_analysis("ATP7B", "MTF1")

# 2. Include alternative promoters and lower the threshold to catch weaker sites
run_binding_site_analysis("ATP7B", "MTF1",
                          transcripts = "all", threshold_frac = 0.75)

# 3. A wider promoter window
run_binding_site_analysis("TERT", "CTCF", upstream = 5000, downstream = 1000)

# 4. A mouse gene (gene coordinates and the HOCOMOCO bundle switch to mouse)
run_binding_site_analysis("Ripk3", "Rela",
                          species = "mus_musculus", hocomoco_species = "mouse")

# 5. Scan once, then restyle without re-querying the databases
s <- scan_binding_sites("FOXP2", "FOXP1", transcripts = "all")
p <- plot_binding_sites(s, base_size = 16, point_size = 4)
p$by_transcript            # show a plot
head(s$result)             # the hits table

# 6. Keep everything in memory (no files written) and inspect the data
res <- run_binding_site_analysis("MYC", "MAX", write_files = FALSE)
res$result[order(res$result$p_adj), ]    # most significant sites first

# 7. Lock the colour scheme so re-runs look identical
run_binding_site_analysis("SNCA", "GATA1", color_seed = 42)

# 8. Scan several TFs at one locus in a loop
for (tf in c("SP1", "MAZ", "KLF4"))
  run_binding_site_analysis("TERT", tf)

# 9. Interactive: point-and-click front end
launch_app()
```
