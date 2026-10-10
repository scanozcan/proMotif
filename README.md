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

For most uses you only need the first function or the app.

| Function | What it does |
|---|---|
| `run_binding_site_analysis(gene, tf, ...)` | One-shot: scan + plots + writes a CSV and PDF figures to disk. Start here. **Omit `tf` for discovery mode** (see below). |
| `launch_app()` | Interactive Shiny app — type a gene + TF, view the plots, download the table. |
| `scan_binding_sites(gene, tf, ...)` | Just the scan (the slow, network part). Returns the data, no plots. |
| `plot_binding_sites(scan, ...)` | Builds the plots from a `scan_binding_sites()` result (fast, local). |
| `run_regulator_discovery(gene, ...)` | Discovery mode: rank a gene's candidate regulators against all of JASPAR CORE (see *Discovery mode*). |
| `discover_regulators(gene, ...)` / `plot_regulators(disc)` | The discovery scan (data) and its plots. |

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
y-axis; the significance plots put `-log10(per-site p)` with a dotted line at
p = 0.05. The per-site p-value follows the FIMO approach (see *Significance*
below), so the strongest matches rise highest.

### CSV columns

| column | meaning |
|---|---|
| `transcript_id`, `is_canonical` | transcript the window came from |
| `chr`, `gene_strand`, `tss_genomic` | locus and the TSS used |
| `database`, `matrix_id`, `matrix_name` | JASPAR or HOCOMOCO and the exact motif matrix |
| `position` | site position relative to the TSS (bp; negative = upstream) |
| `strand` | strand the motif matched (`+`/`-`) |
| `score`, `score_frac` | raw PWM log-odds score and its fraction of the matrix max |
| `pvalue`, `qvalue` | FIMO-style per-site p-value (exact, GC-aware background) and its Benjamini-Hochberg FDR q-value over all positions tested (both NA without `TFMPvalue`) |
| `genomic_position` | absolute chromosome coordinate of the site |
| `site` | the matched promoter sequence |

---

## Discovery mode — "what regulates this gene?"

Call `run_binding_site_analysis()` (or `discover_regulators()`) **with no TF**
and proMotif flips from "does TF X bind gene Y?" to "which TFs *could* regulate
gene Y?". It scans the canonical promoter against the **whole JASPAR CORE
collection** and ranks transcription factors by their single strongest match,
scored with the **same FIMO per-site p-value** as confirm mode.

```r
run_binding_site_analysis("SNCA")            # no TF -> discovery
run_regulator_discovery("SNCA", top_n = 25)  # same, explicit, more TFs
d <- discover_regulators("ATP7B")            # just the ranked table
```

**Extra packages** (discovery only; install once):

```r
install.packages(c("TFMPvalue", "RSQLite", "ggrepel"))
BiocManager::install(c("TFBSTools", "JASPAR2024", "motifmatchr", "Biostrings"))
```

`TFMPvalue` (ranking metric) and a JASPAR data package (`JASPAR2024` + `RSQLite`)
with `TFBSTools` are required; `motifmatchr` + `Biostrings` only pre-shortlist
motifs for speed; `ggrepel` tidies the needle-map labels. All are `Suggests`, so
the core (single-TF) workflow needs none of them.

**Outputs** go to `results/<GENE>_regulators/`: a CSV (`*_candidate_regulators.csv`)
and three PDFs — `*_regulators_ranking.pdf` (lollipop of top TFs by `-log10` p),
`*_regulators_map.pdf` (best-site position per TF), and `*_regulators_needle.pdf`
(a promoter needle plot: position vs significance, so site hotspots stand out).

**How ties and redundancy are handled.** A TF with several matrices is
represented by its best-scoring one. TFs whose best site falls at the **same
locus** (within `site_window`, default 20 bp) are grouped — the same element is
not listed many times. The best-p TF becomes the group representative and the
others are listed in `shares_site_with`; a `(+N)` on a plot label means N more
TFs share that site. **Nothing is dropped**: the CSV keeps every TF with a
`site_group` id and an `is_representative` flag. TF family is only an
annotation (colour), not a grouping key.

**Important — read this as prioritisation, not proof.** Discovery is a
motif-match ranking, **not** motif *enrichment* (there is no background set of
sequences, unlike AME/HOMER/oPOSSUM) and **not** evidence of binding. A motif
matches far more often than the factor actually binds; real occupancy depends on
chromatin, concentration, cofactors and cell type. The ranking also favours
longer, more informative motifs (a perfect long match is rarer by chance). Use
it to decide **which TFs are worth following up**, and confirm with the single-TF
mode plus orthogonal data (ChIP, accessibility).

Key arguments: `top_n` (how many to show), `collapse_site` / `site_window`
(shared-site grouping), `upstream`/`downstream` (promoter window),
`tax_group` (e.g. `"vertebrates"`), `species`.

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
- **Significance (FIMO approach).** proMotif scores each site with a log-odds
  PWM and assigns a p-value exactly as FIMO does: an exact dynamic-programming
  computation (via `TFMPvalue`) of the probability that a random site drawn from
  a zero-order, GC-aware background scores at least as high. Because this p-value
  is monotone in the score, the strongest matches get the smallest p — this is
  what the significance plots show (`-log10 per-site p`). `qvalue` is the
  Benjamini-Hochberg FDR over all positions tested, again following FIMO.
  This is **match significance** (how motif-like the sequence is), not proof of
  binding: a strong match is not the same as a TF actually binding there, which
  additionally depends on chromatin, accessibility, cofactors and cell type.
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
res$result[order(res$result$pvalue), ]   # most significant sites first (FIMO p-value)

# 7. Lock the colour scheme so re-runs look identical
run_binding_site_analysis("SNCA", "GATA1", color_seed = 42)

# 8. Scan several TFs at one locus in a loop
for (tf in c("SP1", "MAZ", "KLF4"))
  run_binding_site_analysis("TERT", tf)

# 9. Interactive: point-and-click front end
launch_app()

# 10. Discovery mode: rank candidate regulators of a gene (no TF given)
run_binding_site_analysis("SNCA")                 # writes results/SNCA_regulators/
d <- discover_regulators("ATP7B", top_n = 30)     # just the ranked table
head(d$top[, c("tf", "pvalue", "qvalue", "n_at_site", "shares_site_with")])
```
