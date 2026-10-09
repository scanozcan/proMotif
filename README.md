<img src="man/figures/logo.svg" alt="proMotif" width="440"/>

# proMotif

Promoter transcription-factor binding-site scanner for internal lab use.

Give it a **gene** and a **TF**; proMotif enumerates the gene's transcripts from
Ensembl, fetches each promoter window (with a UCSC fallback), scans both strands
for **JASPAR** and **HOCOMOCO** motif matches, and reports predicted sites in
TSS-relative and absolute genomic coordinates with per-transcript and genomic
plots. A Shiny app is included.

It reuses the GLproxScape package's Ensembl / JASPAR / HOCOMOCO fetchers and
motif-scan conventions.

## Install

```r
# from the package folder's parent:
devtools::document("proMotif")   # generate help pages from roxygen (first time)
devtools::install("proMotif")

# or straight from a checkout:
# remotes::install_local("proMotif")
```

Dependencies (GLproxScape, ggplot2, httr; shiny for the app) install with it.
An internet connection is required at run time.

## Use

```r
library(proMotif)

# one call: scan + plot + write CSV/PDFs to results/<GENE>_<TF>/
run_binding_site_analysis("ATP7B", "MTF1", threshold_frac = 0.75)

# or the two-step form (scan once, re-style many times without re-scanning)
s <- scan_binding_sites("TERT", "CTCF")
p <- plot_binding_sites(s, base_size = 16)
p$genomic

# interactive
launch_app()
```

### Main functions

| function | role |
|---|---|
| `scan_binding_sites(gene, tf, ...)` | the slow step: network fetch + motif scan; returns data |
| `plot_binding_sites(scan, ...)` | fast: builds the three ggplots (canonical, by-transcript, genomic) |
| `run_binding_site_analysis(...)` | wrapper: scan + plot + write CSV/PDFs |
| `launch_app()` | Shiny front-end |

Key arguments: `species`, `upstream`/`downstream` (2500/500), `threshold_frac`
(0.80), `transcripts` (`"protein_coding"`/`"all"`), `max_transcripts`,
`color_seed`, `base_size`/`point_size`.

## Notes

- **Resilience.** Promoter sequence comes from Ensembl; on an Ensembl timeout or
  error it falls back to the UCSC REST API.
- **Transparency.** The exact JASPAR/HOCOMOCO matrix used is recorded in the
  output (`matrix_id`, `matrix_name`).
- **Thresholding.** `threshold_frac` is a fraction of the matrix min-max score
  (a heuristic, not a p-value). For significance-based calling consider FIMO or
  motifmatchr as the engine.
