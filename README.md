<picture>
  <source media="(prefers-color-scheme: dark)" srcset="man/figures/logo-dark.svg">
  <img src="man/figures/logo-light.svg" alt="proMotif" width="440">
</picture>

# proMotif

[![R-CMD-check](https://github.com/scanozcan/proMotif/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scanozcan/proMotif/actions/workflows/R-CMD-check.yaml)

**proMotif** finds and ranks **transcription factor (TF) binding sites in the
promoters of a gene**. You give it a gene symbol and a TF; proMotif retrieves the
gene's transcription start sites (TSSs), scans every promoter window on both
strands with the TF's **JASPAR** and **HOCOMOCO** matrices, and reports each
predicted site with its exact genomic coordinates and sequence. Each site gets a
**binding probability** that combines motif strength, agreement between the two
motif databases and public **ChIP-Atlas** ChIP-seq evidence, so a long list of
motif matches becomes a short, ranked list of sites to test.

proMotif also works at the level of **promoters**: it groups TSSs into
promoters (P1, P2, ...), compares TFs between the alternative promoters of a
gene, ranks a gene's candidate regulators against a genome-wide promoter
background (discovery mode), and shows which promoter is used in which cell
type with FANTOM5 CAGE data. An interactive Shiny app is included.

Typical uses:

- **Plan experiments at a promoter.** Get the coordinates and sequences of the
  most likely binding sites of a TF to design ChIP-qPCR primers around them, or
  to choose bases to mutate in a reporter construct; lower-ranked sites in the
  same window serve as negative controls.
- **Check a TF-gene hypothesis** against motif evidence and hundreds of public
  ChIP-seq experiments at once.
- **Find candidate regulators** of a gene, promoter by promoter.
- **Compare alternative promoters** of a gene for the TFs that distinguish them,
  and see where each promoter is used.

Motifs: **JASPAR 2026** (CORE, non-redundant) and **HOCOMOCO v14** by default.
proMotif uses the [GLproxScape](https://github.com/scanozcan/GLproxScape)
package's Ensembl fetchers.

---

## Contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [TF binding sites and the binding probability](#tf-binding-sites-and-the-binding-probability)
4. [ChIP-Atlas evidence and the data budget](#chip-atlas-evidence-and-the-data-budget)
5. [Alternative promoters](#alternative-promoters)
6. [Promoter usage across cell types (FANTOM5 CAGE)](#promoter-usage-across-cell-types-fantom5-cage)
7. [Discovery mode: which TFs could regulate this gene?](#discovery-mode-which-tfs-could-regulate-this-gene)
8. [Comparing the promoters of a gene](#comparing-the-promoters-of-a-gene)
9. [The Shiny app](#the-shiny-app)
10. [A complete example with all settings](#a-complete-example-with-all-settings)
11. [Function reference and default values](#function-reference-and-default-values)
12. [Output files](#output-files)
13. [Caching, options and offline use](#caching-options-and-offline-use)
14. [Notes and caveats](#notes-and-caveats)

---

## Installation

From GitHub (installs GLproxScape automatically):

```r
# install.packages("remotes")
remotes::install_github("scanozcan/proMotif")
```

Optional packages:

```r
install.packages(c("TFMPvalue",   # per-site p-values, binding probability, discovery mode
                   "ggrepel",     # labels on the top-sites and discovery plots
                   "curl",        # parallel ChIP-Atlas downloads
                   "shiny"))      # the interactive app, launch_app()
BiocManager::install(c("TFBSTools", "motifmatchr", "Biostrings"))  # faster discovery (optional)
```

`TFMPvalue` is strongly recommended: without it, sites are reported with PWM
scores only (no p-values, and the binding probability falls back to score
fractions).

proMotif needs an **internet connection** the first time a gene, motif database
or ChIP-Atlas file is used (Ensembl, UCSC, JASPAR, HOCOMOCO, EPDnew,
ChIP-Atlas). Everything it downloads is cached under
`tools::R_user_dir("proMotif", "cache")`, so repeated analyses run locally.

---

## Quick start

```r
library(proMotif)

# Predicted ESR1 sites in the TFF1 promoter, ranked by binding probability,
# with ChIP-Atlas evidence from breast cancer cells.
# Writes CSVs and PDFs to results/TFF1_ESR1/
r <- run_binding_site_analysis("TFF1", "ESR1",
                               chip = TRUE,
                               chip_cell_pattern = "MCF-7|T-47D|ZR-75|breast",
                               chip_max_mb = 1000)

r$top_sites          # likely binding sites: coordinates, sequences, probability
r$plots$top_sites    # the top sites on the genome, with the ChIP-seq track

# No TF: rank the candidate regulators of a gene (discovery mode)
run_binding_site_analysis("TFF1")

# Point-and-click front end
launch_app()
```

---

## TF binding sites and the binding probability

### What a scan does

For a gene and a TF, `scan_binding_sites()`:

1. looks up the gene's transcripts in Ensembl (or its experimentally mapped TSSs
   in EPDnew, `tss_source = "epd"`);
2. groups TSSs into promoters and takes a window around each TSS
   (`upstream = 2500`, `downstream = 500` bp);
3. fetches the window sequence from Ensembl (UCSC as fallback);
4. scans both strands with the TF's JASPAR and HOCOMOCO matrices and keeps
   positions scoring at least `threshold_frac` (0.80) of the matrix score range;
5. gives each site a FIMO-style per-site p-value (exact, with the window's base
   composition as background; `TFMPvalue`) and a Benjamini-Hochberg q-value over
   all positions tested;
6. annotates each site with its genomic coordinates and sequences, and scores
   its binding probability.

```r
s <- scan_binding_sites("KLK3", "AR", upstream = 4500)   # slow part: network + scan
s <- add_chip_peaks(s, cell_pattern = "LNCaP|VCaP|22Rv1|prostate", max_mb = 1000)
head(s$sites)                                             # ranked distinct sites
top_binding_sites(s)                                      # P >= 0.5, at least 10 sites
plot_top_sites(s, top_n = 10)                             # fast: local plotting
```

### Site annotation

Every site carries:

| Field | Meaning |
|---|---|
| `genomic_start`, `genomic_end` | motif span, 1-based, genomic plus-strand coordinates |
| `motif_strand_genomic` | strand of the motif on the genome |
| `position` / `position_rel_tss` | position relative to the TSS of its window (bp; negative = upstream) |
| `site` | sequence as read in the promoter (gene orientation) |
| `site_motif` | sequence in motif orientation (reads like the logo) |
| `site_genomic` | sequence on the genomic plus strand (paste into a genome browser) |
| `context_genomic` | `site_genomic` with `flank` (10) bases on each side, flanks in lower case |

### The binding probability

A motif match alone is weak evidence: short motifs match many positions in a
3-kb window by chance. proMotif therefore combines the evidence for each site
on the logit scale:

```
logit P = -3 + 4 x motif_strength + 0.5 x db_agreement + 2.5 x chip_support + 2 x chip_replication
```

| Term | Range | Definition |
|---|---|---|
| `motif_strength` | 0-1 | log10(p) / log10(p of a perfect match to the same matrix). Puts short and long motifs on the same scale. |
| `db_agreement` | 0/1 | 1 if the other database (JASPAR / HOCOMOCO) calls a site within 10 bp on the same strand. |
| `chip_support` | 0-1 | Summit-weighted share of the ChIP-Atlas experiments binding this locus that have a peak centred on the site (see below). |
| `chip_replication` | 0-1 | Grows with the number of experiments with a peak over the site: 1, 3, 10, 20, 60 experiments give 0.002, 0.02, 0.2, 0.5, 0.9. |

What this means in practice:

| Situation | P |
|---|---|
| Perfect motif, no ChIP-seq data | 0.73 |
| Perfect motif, both databases agree | 0.82 |
| Weak motif (strength 0.3), no ChIP-seq data | 0.14 |
| Perfect motif with full ChIP-seq support and replication | > 0.99 |

ChIP-seq evidence can raise a site, never lower it, and without ChIP-seq data
the ranking still follows the motif evidence. `motif_prob` is the same score
without the ChIP-seq terms.

**ChIP-seq support in detail.** For each experiment with a peak over the site,
the closest peak summit gives a weight `exp(-d^2 / (2 sigma^2))` (`sigma = 50`
bp; doubled when only the peak interval is known), multiplied by the site's
motif strength relative to the best motif inside that experiment's peaks, so a
weak match next to a strong one receives only part of the peak's credit. The
weights are summed and divided by the number of experiments with any peak in the
scanned windows plus one (`support = "locus"`, `prior = 1`; one experiment alone
gives at most 0.5). `support = "all"` divides by all queried experiments
instead (also reported as `chip_support_all`). The replication term is
`k^2 / (k^2 + 20^2)` for `k` supporting experiments (`rep_half = 20`), scaled
by the same relative motif strength.

**Confidence and significance.** Sites are `high` (P >= 0.8), `medium`
(0.5-0.8) or `low` confidence. A site is `significant` for its locus when its
robust z-score of logit P against all distinct sites in the scanned windows,
(logit P - median) / MAD, is at least 2 **and** P >= 0.5. In the top-sites
plots, significant sites get a black ring (inside a ChIP-Atlas peak) or a grey
ring (motif evidence only).

The weights are set by design rather than fitted to training data, so P is a
transparent evidence score, not a calibrated posterior probability. Change them
with `options(promotif.score_weights = list(intercept = -3, motif = 4,
databases = 0.5, chip = 2.5, replication = 2))` (any subset), and the
significance thresholds with `options(promotif.sig_z = 2, promotif.sig_min_prob
= 0.5)`.

### Motif logo

```r
plot_motif_logo(s)     # JASPAR and HOCOMOCO matrices, forward and reverse complement
motif_matrices(s)      # probabilities, information content and consensus per position
```

Looking at both matrices is useful: databases sometimes describe the same TF
quite differently, and the logo shows which sequence the sites are matched to.

---

## ChIP-Atlas evidence and the data budget

`chip = TRUE` (or `add_chip_peaks()`) adds public ChIP-seq peaks of the same TF
from [ChIP-Atlas](https://chip-atlas.org) (human, hg38; MACS2 q < 1e-5 by
default). For heterodimer matrices such as `MAX::MYC`, the peaks of both
partners are combined. Restrict the experiments to relevant cell types with a
regular expression (spaces, `-` and `_` are ignored):

```r
run_binding_site_analysis("AFP", "FOXA1", chip = TRUE,
                          chip_cell_pattern = "HepG2|liver")
```

**Data budget.** ChIP-Atlas merges all experiments of a TF into one file, which
for well-studied TFs is very large (ESR1 about 6.9 GB, AR about 5.9 GB in hg38).
`chip_max_mb` (default 300) is a budget, not an on/off switch:

- if the merged file is cached or within the budget, it is used as a whole
  (all experiments);
- otherwise proMotif lists the TF's experiments, keeps those matching
  `chip_cell_pattern`, orders them **newest first** (highest accession number)
  and downloads small per-experiment files (about 0.5-2 MB each, with the exact
  peak summit) until the budget is used.

```r
# up to 1 GB of ESR1 data, newest MCF-7 / T-47D experiments first
run_binding_site_analysis("TFF1", "ESR1", chip = TRUE, chip_max_mb = 1000,
                          chip_cell_pattern = "MCF-7|T-47D")

# at most 50 experiments
run_binding_site_analysis("TFF1", "ESR1", chip = TRUE, chip_max_experiments = 50)

# always per-experiment files (exact summits), or only the merged file
run_binding_site_analysis("HBB", "GATA1", chip = TRUE, chip_mode = "experiments")
run_binding_site_analysis("HBB", "GATA1", chip = TRUE, chip_mode = "assembled")

# the experiment list itself, and peaks of a region directly
list_chip_experiments("AR", cell_pattern = "LNCaP")
fetch_chip_peaks_budget("AR", "19", 50850000, 50856000, budget_mb = 300)
```

All ChIP-Atlas files are cached, so the next gene with the same TF needs no
download. Each site gets `chip_n_exp` (experiments with a peak over it),
`chip_max_score`, `chip_cell_types`, `chip_supported`, `chip_summit_dist`,
`chip_support`, `chip_support_all` and `chip_replication`, and the maps gain a
peak-coverage track.

In discovery mode, `chip = TRUE` instead reports, for each top TF, the share of
its ChIP-Atlas experiments with a peak within 1, 5 or 10 kb of the gene's TSS
(`chip_distance`), from the ChIP-Atlas Target Genes tables.

---

## Alternative promoters

Annotated transcripts often start only a few bp apart, so their windows carry
almost the same sequence. proMotif groups TSSs into **promoters**: TSSs whose
windows overlap by at least `max_overlap` (20%; TSSs within about 2.4 kb for a
-2500/+500 window) form one promoter, numbered P1, P2, ... from 5' to 3'. Every
run prints the promoter structure of the gene.

```r
run_binding_site_analysis("TP73", "TP53",
                          transcripts = "all",        # include non-coding isoforms
                          scan_unit = "promoter",     # one window per promoter
                          tss_source = "epd")         # experimentally mapped TSSs
```

- `transcripts = "protein_coding"` (default) or `"all"`: alternative-promoter
  isoforms are often non-coding biotypes, so `"all"` shows more promoters.
- `scan_unit = "tss"` (default) scans every distinct TSS; `"promoter"` one
  window per promoter.
- `tss_source = "ensembl"` (default) uses annotated transcript starts;
  `"epd"` uses EPDnew TSSs mapped by CAGE and other experiments (human).
- `max_transcripts = 8` caps the number of windows; the canonical TSS comes
  first, then one TSS per promoter, so no promoter is dropped.

`plot_promoter_map()` shows all sites across the promoters (shaded windows,
every TSS drawn, one panel per promoter when they are far apart), and
`summarise_promoters()` gives a per-promoter table.

---

## Promoter usage across cell types (FANTOM5 CAGE)

Motifs say which TFs could bind each promoter; CAGE shows which promoter is
used, and where. With `usage = TRUE` (or `add_promoter_usage()`), the FANTOM5
CAGE signal of peaks within `usage_flank` (500) bp of each promoter's TSSs is
summed in each of 608 sample types (primary cells, tissues, cell lines; donors
averaged), giving each promoter's share of the gene's expression.

```r
options(promotif.fantom5_dir = "/path/to/fantom5_promoter_tpm_hg38")
r <- run_binding_site_analysis("TP73", "TP53", transcripts = "all",
                               scan_unit = "promoter", usage = TRUE)
r$plots$usage          # sample types in which each promoter is used most
r$plots$usage_cancer   # each promoter, normal vs cancer, lineage-matched

p <- plot_promoter_track("TP73")   # CAGE normal vs cancer, promoters, transcripts
ggplot2::ggsave("TP73_track.pdf", p, width = 11, height = attr(p, "height"))
```

- `plot_promoter_usage()`: per promoter, the sample types where it carries the
  largest share, plus the highest-expressing sample types (TPM).
- `plot_promoter_usage_cancer()`: share of every promoter (and total TPM) in
  normal samples vs cancer cell lines of the same lineage, two-sided Wilcoxon.
- `plot_promoter_track()`: one gene on the genome with CAGE peaks in normal and
  cancer samples, promoter arrows labelled with expression and usage share, and
  Ensembl transcripts coloured by the promoter they start from.

The FANTOM5 folder (about 40 MB) is built with `data-raw/fantom5/` and is not
shipped with the package; point `options(promotif.fantom5_dir = ...)` to it or
copy it to `tools::R_user_dir("proMotif", "data")`.

---

## Discovery mode: which TFs could regulate this gene?

Call `run_binding_site_analysis()` without a TF, or `run_regulator_discovery()`,
and proMotif scans **every promoter of the gene** against the whole **JASPAR
CORE** collection (1,019 vertebrate matrices) and ranks TFs separately for each
promoter.

**Promoter background.** Raw p-values favour long or GC-rich motifs that match
well in almost any promoter. By default (`rank_by = "background"`), each TF is
ranked by how unusual this gene's best site is compared with the best site of
the same motif in **20,094 human protein-coding promoters** (Ensembl canonical
TSS, -2500/+500 bp, hg38): `bg_pvalue` is the share of promoters whose best site
scores at least as high (an exponential tail beyond the strongest promoter),
with BH `bg_qvalue`. The background ships with the package for JASPAR 2026;
`rank_by = "fimo"` ranks by the per-site p-value instead.

```r
run_regulator_discovery("SNCA", top_n = 25)            # writes results/SNCA_regulators/
d <- discover_regulators("ATP7B", tss_source = "epd")   # the ranked table only
d <- add_chip_targets(d, distance = 1)                  # ChIP-Atlas Target Genes support
plot_regulators(d)$ranking
plot_regulators_chip(d)

ex <- expand_regulators(d, n = 3, chip = TRUE)          # full single-TF scans of the top TFs
plot_top_sites(ex[[1]])
```

TFs whose best sites fall at the same locus (within `site_window`, 20 bp) are
grouped, so one element is not listed many times; the CSV keeps every TF with a
`site_group` and `is_representative` flag. `run_regulator_discovery()` follows
up the top `expand_top` (5) TFs with targeted scans written to
`results/<GENE>_regulators/top_TFs/<TF>/`.

Discovery is a **prioritisation**: a strong best site is not evidence of
binding, which depends on chromatin, TF levels, cofactors and cell type. Follow
up candidates with the single-TF mode and ChIP-Atlas evidence.

---

## Comparing the promoters of a gene

```r
run_promoter_comparison("TP73", transcripts = "all")
cmp <- compare_promoter_regulators("TP73", transcripts = "all", tss_source = "epd")
plot_promoter_regulators(cmp, top_n = 30)
```

Each promoter is scanned against JASPAR CORE and TFs are compared by the
**window p-value** of their best site, `1 - (1 - p)^m` (m = positions tested),
the chance that a random window of this size holds a site this strong. A TF is
called present at a promoter with window p <= `present_p` (0.01) and absent with
window p > `absent_p` (0.2), giving calls such as `P1-specific` or `shared`.

---

## The Shiny app

```r
launch_app()
```

- **Gene + TF:** top sites (binding probability, hover a point to see the site
  id, coordinates, sequence with flanks, p- and q-values, ChIP-seq evidence and
  probability; ranked table and download), motif logo, promoter map,
  promoter usage and normal vs cancer plots, canonical / per-transcript /
  genomic views (scores and significance), promoter table, site table.
- **Discover regulators:** ranking, needle and map plots, motif vs ChIP-Atlas
  plot, targeted maps of the top TFs, table.
- **Compare promoters:** promoter-specific and shared TFs, promoter table.

The analysis runs once per **Run**; text size, colours and the ring cut-off
restyle the plots instantly. ChIP-Atlas evidence, its cell types and data budget,
motif database releases and promoter usage are selectable, and a **Log** tab
shows the console messages of the run.

---

## A complete example with all settings

Every argument of the main function, with its default value:

```r
library(proMotif)

r <- run_binding_site_analysis(
  gene                 = "TFF1",       # HGNC gene symbol (required)
  tf                   = "ESR1",       # TF symbol; NULL = discovery mode
  species              = "homo_sapiens",
  # ---- promoter windows and TSSs
  upstream             = 2500,         # bp upstream of each TSS
  downstream           = 500,          # bp downstream of each TSS
  transcripts          = "protein_coding",   # or "all"
  tss_source           = "ensembl",    # or "epd" (EPDnew mapped TSSs, human)
  max_transcripts      = 8,            # max TSS windows scanned
  scan_unit            = "tss",        # or "promoter"
  max_overlap          = 0.2,          # window overlap that merges TSSs into one promoter
  merge_dist           = NULL,         # or a distance in bp (overrides max_overlap)
  # ---- motifs
  threshold_frac       = 0.80,         # PWM cutoff, fraction of the matrix score range
  jaspar_release       = "2026",       # "2024", "2022" also work
  hocomoco_version     = "v14",        # "v12", "v11" also work
  hocomoco_species     = "human",      # or "mouse"
  flank                = 10,           # bases on each side in context_genomic
  # ---- ChIP-Atlas evidence
  chip                 = FALSE,        # TRUE to add ChIP-seq evidence
  chip_threshold       = "05",         # MACS2 q class: "05", "10", "20"
  chip_cell_pattern    = NULL,         # e.g. "MCF-7|T-47D|breast"
  chip_max_mb          = 300,          # data budget (MB)
  chip_mode            = "auto",       # "auto", "assembled", "experiments"
  chip_max_experiments = Inf,          # cap on experiments
  chip_sigma           = 50,           # summit-proximity scale (bp)
  chip_support         = "locus",      # or "all"
  chip_ring_p          = 1e-3,         # promoter map: ring in-peak sites with p <= this
  chip_ring_min_exp    = 1,            # promoter map: ... and >= this many experiments
  # ---- promoter usage (FANTOM5)
  usage                = FALSE,        # TRUE to add CAGE promoter usage
  usage_flank          = 500,          # bp around a promoter's TSSs
  # ---- ranked sites
  top_sites_min_prob   = 0.5,          # top-sites table: P >= this ...
  top_sites_min_n      = 10,           # ... and at least this many sites
  label_top            = 10,           # sites labelled in the top-sites plot
  # ---- output and plotting
  out_root             = "results",
  write_files          = TRUE,
  color_seed           = NULL,         # a number fixes the database colours
  base_size            = 12,
  point_size           = 2.6,
  top_n                = 20            # discovery mode: TFs shown
)

r$top_sites                # likely binding sites
r$sites                    # all distinct sites, ranked
r$result                   # every site row (per window and database)
r$promoter_summary         # one row per promoter
r$plots$top_sites          # labelled top-sites plot (and $top_sites_nolabels)
r$plots$promoter_map       # sites across promoters (-log10 p)
r$plots$motif_logo         # JASPAR and HOCOMOCO logos
r$scan                     # the full scan object for re-plotting
```

The same analysis step by step, which lets you re-plot without scanning again:

```r
s <- scan_binding_sites("TFF1", "ESR1", threshold_frac = 0.8, transcripts = "all",
                        scan_unit = "promoter", tss_source = "epd")
s <- add_chip_peaks(s, cell_pattern = "MCF-7|T-47D", max_mb = 1000,
                    mode = "auto", sigma = 50, support = "locus", rep_half = 20)
s <- add_promoter_usage(s, flank = 500)

top_binding_sites(s, min_prob = 0.6, min_n = 5)
plot_top_sites(s, top_n = 10, base_size = 14, color_seed = 1)
plot_top_sites(s, top_n = 0)                                # no labels
plot_promoter_map(s, chip_ring_p = 1e-4)
plot_promoter_map(s, y_metric = "probability", label_top = 5, chip_track = FALSE)
plot_motif_logo(s)
plot_promoter_usage(s); plot_promoter_usage_cancer(s)
plot_binding_sites(s)$genomic
```

---

## Function reference and default values

### Main workflows

| Function | What it does |
|---|---|
| `run_binding_site_analysis(gene, tf, ...)` | Scan + ChIP-Atlas + usage + plots + files for one gene and TF; without `tf`, discovery mode. |
| `run_regulator_discovery(gene, ...)` | Discovery mode with files and targeted follow-up of the top TFs. |
| `run_promoter_comparison(gene, ...)` | Compare a gene's promoters for promoter-specific TFs, with files. |
| `launch_app()` | Interactive Shiny app. |

### Building blocks

| Function | What it does |
|---|---|
| `scan_binding_sites(gene, tf, ...)` | The scan: sites, coordinates, sequences, p-values, binding probability. |
| `add_chip_peaks(scan, ...)` | ChIP-Atlas peaks of the TF, ChIP support and replication; re-scores the sites. |
| `score_binding_sites(scan, weights)` | (Re-)compute binding probabilities, confidence and significance. |
| `top_binding_sites(scan, min_prob, min_n)` | The likely binding sites. |
| `add_promoter_usage(scan, ...)` | FANTOM5 CAGE usage of each promoter. |
| `summarise_promoters(scan)` | One row per promoter. |
| `motif_matrices(scan)` | The TF's matrices as probabilities with information content. |
| `discover_regulators(gene, ...)` | Discovery scan: ranked candidate regulators per promoter. |
| `add_chip_targets(disc, ...)` | ChIP-Atlas Target Genes support for discovery candidates. |
| `expand_regulators(disc, ...)` | Targeted single-TF scans of the top candidates. |
| `compare_promoter_regulators(gene, ...)` | TFs that distinguish a gene's promoters. |
| `list_chip_experiments(tf, ...)` | A TF's ChIP-Atlas experiments, newest first. |
| `fetch_chip_peaks(tf, chr, from, to, ...)` | ChIP-Atlas peaks of a TF in a region (merged file or budget). |
| `fetch_chip_peaks_budget(tf, chr, from, to, ...)` | Peaks from per-experiment files within a data budget. |
| `clear_sequence_cache()` | Remove cached promoter sequences and Ensembl lookups. |

### Plots

| Function | What it shows |
|---|---|
| `plot_top_sites(scan, ...)` | Binding probability on the genome, ChIP-seq track, labelled top sites, significance rings. |
| `plot_promoter_map(scan, ...)` | Sites across all promoters (-log10 p or probability), TSSs, ChIP-seq track and rings. |
| `plot_motif_logo(scan, ...)` | Sequence logos of the TF's JASPAR and HOCOMOCO matrices. |
| `plot_binding_sites(scan, ...)` | Canonical, per-transcript and genomic views (scores and significance). |
| `plot_promoter_usage(scan, ...)` | Sample types where each promoter is used. |
| `plot_promoter_usage_cancer(scan, ...)` | Promoter usage in normal vs cancer samples. |
| `plot_promoter_track(gene, ...)` | CAGE signal, promoters and transcripts of a gene. |
| `plot_regulators(disc, ...)` | Discovery ranking, needle and map plots. |
| `plot_regulators_chip(disc, ...)` | Discovery: motif vs ChIP-Atlas evidence. |
| `plot_promoter_regulators(cmp, ...)` | Promoter-specific and shared TFs. |

### Default values

**`scan_binding_sites()`**

| Argument | Default | Meaning |
|---|---|---|
| `gene`, `tf` | required | gene and TF symbols |
| `species` | `"homo_sapiens"` | Ensembl species |
| `upstream`, `downstream` | `2500`, `500` | window around each TSS (bp) |
| `threshold_frac` | `0.80` | PWM cutoff (fraction of the matrix score range) |
| `hocomoco_version`, `hocomoco_species` | `"v14"`, `"human"` | HOCOMOCO release and species |
| `jaspar_release` | `"2026"` | JASPAR release |
| `transcripts` | `"protein_coding"` | or `"all"` |
| `max_transcripts` | `8` | max TSS windows |
| `merge_dist` | `NULL` | TSS merge distance (bp); `NULL` = from `max_overlap` |
| `scan_unit` | `"tss"` | or `"promoter"` |
| `max_overlap` | `0.2` | window overlap that merges TSSs |
| `tss_source` | `"ensembl"` | or `"epd"` |
| `flank` | `10` | flanking bases in `context_genomic` |

**`add_chip_peaks()`**

| Argument | Default | Meaning |
|---|---|---|
| `threshold` | `"05"` | MACS2 q class (`"05"`, `"10"`, `"20"`) |
| `cell_pattern` | `NULL` | cell-type regular expression |
| `max_mb` | `300` | data budget (MB) |
| `mode` | `"auto"` | `"auto"`, `"assembled"`, `"experiments"` |
| `max_experiments` | `Inf` | cap on experiments |
| `sigma` | `50` | summit-proximity scale (bp) |
| `support` | `"locus"` | or `"all"` |
| `prior` | `1` | pseudo-count on locus experiments |
| `rep_half` | `20` | experiments giving half the replication term |

**`run_binding_site_analysis()`**: all arguments of `scan_binding_sites()`, plus
`chip = FALSE`, `chip_threshold = "05"`, `chip_cell_pattern = NULL`,
`chip_max_mb = 300`, `chip_mode = "auto"`, `chip_max_experiments = Inf`,
`chip_sigma = 50`, `chip_support = "locus"`, `chip_ring_p = 1e-3`,
`chip_ring_min_exp = 1`, `usage = FALSE`, `usage_flank = 500`,
`top_sites_min_prob = 0.5`, `top_sites_min_n = 10`, `label_top = 10`,
`out_root = "results"`, `write_files = TRUE`, `color_seed = NULL`,
`base_size = 12`, `point_size = 2.6`, `top_n = 20`.

**Plots**

| Function | Defaults |
|---|---|
| `plot_top_sites()` | `top_n = 10`, `base_size = 12`, `point_size = 2.6`, `color_seed = NULL`, `chip_ring_p = 1`, `chip_track = TRUE` |
| `plot_promoter_map()` | `base_size = 12`, `point_size = 2.6`, `color_seed = NULL`, `chip_ring_p = 1e-3`, `chip_ring_min_exp = 1`, `facet_regions = TRUE`, `y_metric = "p"`, `label_top = 0`, `chip_track = TRUE` |
| `plot_motif_logo()` | `base_size = 12`, `revcomp = TRUE` |
| `plot_binding_sites()` | `base_size = 12`, `point_size = 2.6`, `color_seed = NULL` |
| `plot_promoter_usage()` | `top_n = 10`, `base_size = 12` |
| `plot_promoter_usage_cancer()` | `min_n = 3`, `base_size = 12` |
| `plot_promoter_track()` | `species = "homo_sapiens"`, `tss_source = "epd"`, `transcripts = "all"`, `normal = "primary cells"`, `lineage = NULL`, `cage = TRUE`, `max_transcripts = 30`, `tx_match = 500`, `flank = 2000`, `max_overlap = 0.2`, `base_size = 11`, `dir = NULL` |
| `plot_regulators()`, `plot_regulators_chip()` | `base_size = 12`, `point_size = 3` |
| `plot_promoter_regulators()` | `top_n = 30`, `base_size = 12` |

**Other functions**

| Function | Defaults |
|---|---|
| `top_binding_sites()` | `min_prob = 0.5`, `min_n = 10` |
| `score_binding_sites()` | `weights = getOption("promotif.score_weights")` (intercept -3, motif 4, databases 0.5, chip 2.5, replication 2) |
| `add_promoter_usage()` | `flank = 500`, `dir = NULL`, `min_tpm = 1`, `dominant_frac = 0.6` |
| `fetch_chip_peaks()` | `species = "homo_sapiens"`, `threshold = "05"`, `cell_pattern = NULL`, `max_mb = 300`, `mode = "auto"`, `max_experiments = Inf` |
| `fetch_chip_peaks_budget()` | `species = "homo_sapiens"`, `threshold = "05"`, `cell_pattern = NULL`, `budget_mb = 1000`, `max_experiments = Inf`, `batch = 24` |
| `list_chip_experiments()` | `species = "homo_sapiens"`, `cell_pattern = NULL` |
| `discover_regulators()` | `species = "homo_sapiens"`, `upstream = 2500`, `downstream = 500`, `collection = "CORE"`, `tax_group = "vertebrates"`, `top_n = 20`, `collapse_site = TRUE`, `site_window = 20`, `use_motifmatchr = TRUE`, `jaspar_db = NULL`, `jaspar_release = "2026"`, `rank_by = "background"`, `promoters = "all"`, `transcripts = "protein_coding"`, `merge_dist = NULL`, `max_overlap = 0.2`, `tss_source = "ensembl"` |
| `run_regulator_discovery()` | `top_n = 20`, `out_root = "results"`, `write_files = TRUE`, `base_size = 12`, `point_size = 3`, `chip = FALSE`, `chip_distance = 1`, `chip_cell_pattern = NULL`, `chip_max_mb = 300`, `expand_top = 5`, `expand_threshold = 0.8`, `expand_tfs = NULL`, `expand_by = "auto"`, plus `discover_regulators()` arguments |
| `add_chip_targets()` | `distance = 1`, `cell_pattern = NULL`, `max_mb = 300` |
| `expand_regulators()` | `tfs = NULL`, `n = 5`, `by = "auto"`, `chip = FALSE`, `threshold_frac = 0.8`, `hocomoco_version = "v14"`, `scan_unit = "promoter"`, `chip_cell_pattern = NULL`, `chip_max_mb = 300` |
| `compare_promoter_regulators()` | `species = "homo_sapiens"`, `upstream = 2500`, `downstream = 500`, `transcripts = "protein_coding"`, `merge_dist = NULL`, `max_overlap = 0.2`, `present_p = 0.01`, `absent_p = 0.2`, `collection = "CORE"`, `tax_group = "vertebrates"`, `use_motifmatchr = TRUE`, `jaspar_db = NULL`, `jaspar_release = "2026"`, `tss_source = "ensembl"` |
| `run_promoter_comparison()` | `top_n = 30`, `out_root = "results"`, `write_files = TRUE`, `base_size = 12`, plus `compare_promoter_regulators()` arguments |
| `launch_app()` | `launch.browser = TRUE` |

---

## Output files

`run_binding_site_analysis(gene, tf)` writes to `results/<GENE>_<TF>/`:

| File | Content |
|---|---|
| `<GENE>_<TF>_top_sites.csv` | likely binding sites (P >= 0.5, at least 10), one row per distinct site |
| `<GENE>_<TF>_all_sites_ranked.csv` | all distinct sites, ranked by binding probability |
| `<GENE>_<TF>_binding_sites.csv` | every site row (per TSS window and database) |
| `<GENE>_<TF>_top_sites[_chip].pdf` | top-sites plot with labels (and ChIP-seq track) |
| `<GENE>_<TF>_top_sites[_chip]_nolabels.pdf` | the same without labels |
| `<GENE>_<TF>_promoter_map[_chip].pdf` | sites across the promoters (-log10 p) |
| `<GENE>_<TF>_genomic.pdf` | sites on genomic coordinates (score fraction) |
| `<GENE>_<TF>_motif_logo.pdf`, `_motif_matrices.csv` | the TF's matrices |
| `<GENE>_<TF>_promoters.csv` | one row per promoter |
| `<GENE>_<TF>_promoter_usage.csv`, `_promoter_usage.pdf`, `_promoter_usage_normal_vs_cancer.pdf` | with `usage = TRUE` |

All plots of one run share the same two database colours (`color_seed` fixes
them across runs).

**Columns of the top-sites and ranked-sites tables**

| Column | Meaning |
|---|---|
| `site_id`, `rank` | S1, S2, ... by binding probability |
| `gene`, `tf`, `chrom` | gene, TF, chromosome (`chr` prefix) |
| `genomic_start`, `genomic_end`, `motif_strand_genomic` | motif span and strand on the genome (1-based) |
| `site_center`, `genomic_position` | centre of the site; motif start in gene orientation |
| `promoter_id`, `in_promoter_windows` | promoter of the best window; all promoter windows containing the site |
| `transcript_id`, `position_rel_tss` | window the site was scored in; position relative to its TSS (bp) |
| `databases`, `matrix_ids`, `n_databases_nearby` | databases and matrices calling the site; databases with a site within 10 bp |
| `site_motif`, `site_genomic`, `context_genomic` | sequence in motif orientation, on the plus strand, and with flanks |
| `score_frac`, `pvalue`, `qvalue` | PWM score fraction, per-site p-value, BH q-value |
| `motif_strength`, `db_agreement`, `motif_prob` | motif terms and the motif-only probability |
| `chip_n_exp`, `chip_support`, `chip_support_all`, `chip_replication`, `chip_summit_dist`, `chip_cell_types` | ChIP-Atlas evidence |
| `binding_prob`, `confidence`, `score_z`, `significant` | binding probability, class, robust z-score, locus-level significance |

`binding_sites.csv` has one row per window and database with the window-level
columns (`transcript_id`, `is_canonical`, `transcript_label`, `tss_genomic`,
`matrix_name`, `score`, `site`) in addition to the above.

Discovery mode writes to `results/<GENE>_regulators/` (`*_candidate_regulators.csv`,
ranking, map and needle plots, and `top_TFs/<TF>/` follow-up scans); promoter
comparison writes `*_promoter_regulators.csv` and its plot.

---

## Caching, options and offline use

Everything downloaded is cached under `tools::R_user_dir("proMotif", "cache")`:

| Cache | Content |
|---|---|
| `jaspar/`, `hocomoco/` | motif collections (downloaded once) |
| `chipatlas/` | merged and per-experiment ChIP-Atlas peak files, experiment lists |
| `promoter_seqs/`, `ensembl_lookup/` | promoter sequences and Ensembl gene annotations |
| `epd/` | EPDnew TSSs |

After the first run of a gene, scans need no network; `clear_sequence_cache()`
removes cached sequences and lookups.

| Option | Default | Effect |
|---|---|---|
| `promotif.score_weights` | `NULL` | weights of the binding probability, e.g. `list(chip = 1.5)` |
| `promotif.sig_z` | `2` | robust z-score for locus-level significance |
| `promotif.sig_min_prob` | `0.5` | minimum probability for significance |
| `promotif.refresh_cache` | `FALSE` | re-fetch sequences and lookups and overwrite the cache |
| `promotif.fantom5_dir` | `NULL` | folder with the FANTOM5 promoter data |

---

## Notes and caveats

- **A motif match is not binding.** Whether a TF occupies a site depends on
  chromatin, TF levels, cofactors and cell type. The binding probability adds
  public ChIP-seq evidence, but a ChIP-seq peak (~200-500 bp) shows that the TF
  binds a region in some cell type, not necessarily at the predicted motif or in
  your cells; restrict `chip_cell_pattern` to relevant cell types.
- **The binding probability is an evidence score.** Its weights are set by
  design, not fitted; compare sites within a gene and between genes, and check
  the evidence columns behind a high score.
- **Short matrices** reach their maximal motif strength more easily than long
  ones; a high `motif_strength` from a short matrix is weaker evidence.
- **HOCOMOCO** matrices were partly trained on public ChIP-seq data, so their
  agreement with ChIP-Atlas is partly circular.
- **Mouse.** Scans work for mouse genes (`species = "mus_musculus"`,
  `hocomoco_species = "mouse"`); ChIP-Atlas peaks and FANTOM5 usage are human
  only, because ChIP-Atlas mouse coordinates (mm10) differ from Ensembl (GRCm39).
- **Significance (FIMO approach).** Per-site p-values are exact
  (dynamic programming via `TFMPvalue`) under a zero-order background with the
  base composition of the window; q-values are Benjamini-Hochberg over all
  positions tested.
- **Resilience.** Ensembl requests are retried with back-off and fall back to
  the UCSC API; cached results are used whenever available.
- **Tuning.** Lower `threshold_frac` (e.g. 0.75) to include weaker matches;
  widen `upstream` to reach distal elements (e.g. 4500 bp for the KLK3 upstream
  enhancer).

## Citation

Ozcan SC. *proMotif enables promoter-resolved transcription factor binding
analysis of alternative promoter usage.* (in preparation)

## License

MIT, see `LICENSE`.
