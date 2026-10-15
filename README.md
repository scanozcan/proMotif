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

Motifs come from **JASPAR 2026** and **HOCOMOCO v14** by default (older
releases selectable). It reuses the
[GLproxScape](https://github.com/scanozcan/GLproxScape) package's Ensembl
fetchers.

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

An **internet connection is required at run time** — Ensembl is queried live.
The JASPAR 2026 and HOCOMOCO v14 collections are downloaded once and cached
(`tools::R_user_dir("proMotif", "cache")`), as are ChIP-Atlas files.

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
| `run_promoter_comparison(gene, ...)` | TFs that distinguish a gene's alternative promoters (see *Alternative promoters*). |
| `plot_promoter_map(scan)` / `summarise_promoters(scan)` | Sites across all promoters of a gene, and a per-promoter summary. |
| `add_chip_peaks(scan)` / `add_chip_targets(disc)` | Add ChIP-Atlas evidence (also via `chip = TRUE`, see *ChIP-Atlas evidence*). |

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
| `max_transcripts` | `8` | cap on TSS windows scanned (canonical first, then one TSS per promoter) |
| `jaspar_release` | `"2026"` | JASPAR release (CORE non-redundant; `"2024"`, `"2022"` also work) |
| `hocomoco_version`, `hocomoco_species` | `"v14"`, `"human"` | HOCOMOCO release (`"v12"`/`"v11"` also work; v14 covers human and mouse in one collection) |
| `out_root` | `"results"` | output folder |
| `color_seed` | `NULL` | fix the two database colours reproducibly (NULL = fresh each run) |
| `base_size`, `point_size` | `12`, `2.6` | plot font and marker size |
| **Promoters** | | |
| `scan_unit` | `"tss"` | `"tss"` = one window per distinct TSS; `"promoter"` = one window per promoter |
| `max_overlap` | `0.2` | TSSs whose windows overlap by at least this fraction form one promoter |
| `merge_dist` | `NULL` | override: TSSs closer than this (bp) form one promoter (`NULL` = from `max_overlap`, ~2.4 kb) |
| **ChIP-Atlas** (human) | | |
| `chip` | `FALSE` | add ChIP-Atlas evidence (downloaded once, then cached) |
| `chip_threshold` | `"05"` | peak q-value class: `"05"`, `"10"`, `"20"` (q < 1e-5, 1e-10, 1e-20) |
| `chip_cell_pattern` | `NULL` | regex restricting experiments to cell types, e.g. `"HepG2\|liver"` |
| `chip_max_mb` | `300` | skip ChIP-Atlas files larger than this (MB) |
| `chip_ring_p` | `1e-3` | promoter map: ring in-peak sites with per-site p at or below this |
| `chip_ring_min_exp` | `1` | promoter map: ring only sites supported by at least this many experiments |
| `chip_distance` | `1` | discovery mode only (`run_regulator_discovery`): kb around the TSS (1, 5 or 10) |

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

## Alternative promoters

Annotated transcripts often start only a few bp apart, so their promoter windows
carry almost the same sequence. proMotif groups TSSs into **promoters**: TSSs
whose scan windows overlap by at least `max_overlap` (default 20%, i.e. TSSs
within ~2.4 kb for a -2500/+500 window) form one promoter, numbered P1, P2, ...
from 5' to 3'. Every run prints the gene's promoter structure.

```r
run_binding_site_analysis("TP73", "TP53", transcripts = "all",
                          scan_unit = "promoter")   # one window per promoter
run_promoter_comparison("TP73", transcripts = "all")  # TFs that distinguish them
```

Outputs: `*_promoters.csv` (per-promoter summary) and `*_promoter_map.pdf`
(each promoter's window labelled with the number of transcripts grouped there,
`*` = canonical; every TSS drawn; one panel per promoter when they are far
apart). `max_transcripts` keeps the canonical TSS first, then one TSS per
promoter, so no promoter is silently dropped.

`run_promoter_comparison()` scans each promoter against JASPAR CORE and compares
the **window p-value** of each TF's best site, 1 - (1 - p)^m (m = positions
tested): the chance that a random window of this size holds a site this strong.
Windows have equal length, so promoters are compared fairly. Calls
(`P1-specific`, `shared`, ...) are a prioritisation, not proof of binding.

## ChIP-Atlas evidence

`chip = TRUE` adds public ChIP-seq evidence from
[ChIP-Atlas](https://chip-atlas.org) (human, hg38). Files are downloaded once and
cached; files above `chip_max_mb` (default 300 MB) are skipped with a message.

```r
run_binding_site_analysis("ATP7B", "MTF1", chip = TRUE,
                          chip_cell_pattern = "HepG2|liver")
run_regulator_discovery("MT2A", chip = TRUE, chip_distance = 1)
```

- **Confirm mode:** each site gets `chip_n_exp`, `chip_max_score`,
  `chip_cell_types` and `chip_supported`; the promoter map gains a peak-coverage
  track, and in-peak sites with p <= `chip_ring_p` (default 1e-3) get a black ring.
- **Discovery mode:** each top TF gets the share of its ChIP-Atlas experiments
  with a peak within +/- 1/5/10 kb of the TSS, plus a motif-vs-ChIP plot.

Caveats: a peak (~200-500 bp) shows the TF binds that region in some cell type,
not necessarily at the predicted motif or in your cells; HOCOMOCO matrices were
partly trained on public ChIP-seq, so agreement with HOCOMOCO is partly circular.

---

## Discovery mode — "what regulates this gene?"

Call `run_binding_site_analysis()` (or `discover_regulators()`) **with no TF**
and proMotif flips from "does TF X bind gene Y?" to "which TFs *could* regulate
gene Y?". It scans the canonical promoter against the **whole JASPAR CORE
collection** and ranks transcription factors by their single strongest match.

**Ranking against a promoter background (default).** A raw per-site p-value
favours long or GC-rich motifs that match well in almost *any* promoter. So by
default each TF is ranked by how unusual this gene's best site is **compared with
the best site of the same motif in 2,000 random human promoters** (Ensembl
canonical TSS, same -2500/+500 window, hg38): `bg_pvalue` = P(a typical
promoter's best site scores at least this high), read from that motif's
empirical background distribution (an exponential tail is used only beyond the
strongest background promoter), with BH `bg_qvalue`. `rank_by = "fimo"` restores the FIMO per-site
p-value ranking. The background is precomputed for JASPAR 2026 and shipped with
the package (`inst/extdata/`; built by `data-raw/make_promoter_background.py`);
for other windows the background is shifted by extreme-value max-stability; for other
species or releases discovery falls back to FIMO p.

```r
run_binding_site_analysis("SNCA")            # no TF -> discovery
run_regulator_discovery("SNCA", top_n = 25)  # same, explicit, more TFs
d <- discover_regulators("ATP7B")            # just the ranked table
```

**Extra packages** (discovery only; install once):

```r
install.packages(c("TFMPvalue", "ggrepel"))
BiocManager::install(c("TFBSTools", "motifmatchr", "Biostrings"))  # optional, speed
```

`TFMPvalue` (ranking metric) is required. The motif collection is JASPAR 2026
CORE vertebrates (1,019 non-redundant matrices), downloaded once from
jaspar.elixir.no and cached, so no JASPAR data package is needed.
`motifmatchr` + `Biostrings` + `TFBSTools` only pre-shortlist motifs for speed;
`ggrepel` tidies the needle-map labels. To use a Bioconductor JASPAR package
instead (e.g. `JASPAR2024` + `RSQLite`), pass `jaspar_release = NULL` or a
`jaspar_db`. All are `Suggests`, so the core (single-TF) workflow needs none of
them.

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

**Important — read this as prioritisation, not proof.** The promoter background
answers "is this gene's best site unusually strong for this motif?", not "is the
motif enriched across a gene set" (use AME/HOMER/Pscan for sets), and a motif
match is **not** evidence of binding: real occupancy depends on chromatin,
concentration, cofactors and cell type. Use it to decide **which TFs are worth
following up**, and confirm with the single-TF mode plus orthogonal data
(`chip = TRUE` for ChIP-Atlas, accessibility).

Key arguments: `top_n` (how many to show), `rank_by` (`"background"` or
`"fimo"`), `collapse_site` / `site_window` (shared-site grouping), `upstream`/`downstream` (promoter window),
`tax_group` (e.g. `"vertebrates"`), `species`.

---

## The Shiny app

```r
launch_app()
```

Three modes:

- **Gene + TF** — binding sites: promoter map, canonical / per-transcript /
  genomic plots (score and significance), per-promoter summary, table + CSV.
  Optional ChIP-Atlas evidence (peak track and black rings; the ring p cut-off
  re-styles instantly).
- **Discover regulators** — gene only: ranking, needle and map plots (promoter
  background or FIMO ranking), optional ChIP-Atlas Target Genes plot, table.
- **Compare promoters** — gene only: TFs specific to / shared between the gene's
  promoters, promoter table.

The analysis runs once per **Run**; the text/point-size slider and colour
controls restyle the plots instantly. JASPAR / HOCOMOCO releases are selectable,
and a **Log** tab shows the run's console messages (promoter structure, matrices
used, downloads).

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

# 4. A mouse gene (mouse gene coordinates; HOCOMOCO v14 covers mouse TFs)
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
