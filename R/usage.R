## ===========================================================================
## 10. PROMOTER USAGE ACROSS CELL TYPES (FANTOM5 CAGE)
## ---------------------------------------------------------------------------
## Which of a gene's promoters is actually used, and where? FANTOM5 CAGE maps
## transcription start sites genome-wide in ~1,000 human primary cells, tissues
## and cell lines (608 sample types after averaging donors / replicates). For
## each promoter (P1, P2, ...) the CAGE signal (TPM) of peaks whose summit lies
## within `flank` bp of the promoter's TSSs (same strand) is summed per sample
## type; usage = that promoter's share of the gene's total CAGE signal.
## A promoter is "dominant" in a sample type when it carries >= 60% of the
## signal and the gene is expressed (>= 1 TPM). Different dominant promoters in
## different cell types = promoter switching.
## Data: folder fantom5_promoter_tpm_hg38/ (one .rds per chromosome and strand,
## ~40 MB; not shipped with the package - built once by data-raw/fantom5/).
## Looked up in: dir argument, getOption("promotif.fantom5_dir"),
## tools::R_user_dir("proMotif", "data"), working directory. Human only.
## CAGE misses weak promoters: "not used" means "below CAGE detection".
## ===========================================================================

.f5_env <- new.env(parent = emptyenv())

## folder fantom5_promoter_tpm_hg38/ (one .rds per chromosome and strand)
.fantom5_dir <- function(dir = NULL) {
  fn   <- "fantom5_promoter_tpm_hg38"
  cand <- c(dir, getOption("promotif.fantom5_dir"),
            file.path(tools::R_user_dir("proMotif", "data"), fn), file.path(".", fn))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  d <- cand[dir.exists(cand)][1]
  if (is.na(d)) message("  promoter usage: folder ", fn, " not found (set ",
                        "options(promotif.fantom5_dir = ...) or pass dir; see ",
                        "data-raw/fantom5/README.md) - skipped.")
  if (is.na(d)) NULL else d
}

.fantom5_chunk <- function(chrom, strand, dir = NULL) {
  dd <- .fantom5_dir(dir)
  if (is.null(dd)) return(NULL)
  key <- paste0(chrom, "_", strand)
  if (!is.null(.f5_env[[key]])) return(.f5_env[[key]])
  f <- file.path(dd, paste0(key, ".rds"))
  if (!file.exists(f)) { message("  promoter usage: no FANTOM5 data for ", key, "."); return(NULL) }
  if (is.null(.f5_env$samples)) .f5_env$samples <- readRDS(file.path(dd, "samples.rds"))
  ch <- readRDS(f)
  ch$samples <- .f5_env$samples
  .f5_env[[key]] <- ch
  ch
}

#' Promoter usage across cell types (FANTOM5 CAGE)
#'
#' For each promoter (P1, P2, ...) of a confirm-mode scan, sums the FANTOM5 CAGE
#' signal (TPM) of peaks whose summit lies within \code{flank} bp of the
#' promoter's TSSs (same strand) in each of 608 sample types (primary cells,
#' tissues, cell lines; donors averaged), and reports each promoter's share of
#' the gene's signal. Human (hg38) only.
#'
#' @param scan Result of \code{\link{scan_binding_sites}}.
#' @param flank Window (bp) around each promoter's TSSs (default 500).
#' @param dir Folder \code{fantom5_promoter_tpm_hg38} (one .rds per chromosome
#'   and strand). Not shipped with the package: build it with the scripts in
#'   \code{data-raw/fantom5/}. If \code{NULL}, \code{getOption("promotif.fantom5_dir")},
#'   \code{tools::R_user_dir("proMotif", "data")} and the working directory are searched.
#' @param min_tpm Gene counted as expressed at >= this total TPM (default 1).
#' @param dominant_frac Share (default 0.6) at which a promoter is "dominant".
#' @return \code{scan} with an added \code{usage} element (\code{tpm} matrix,
#'   per-sample-type \code{table}, settings). CAGE misses weak promoters:
#'   "not used" means "below CAGE detection".
#' @export
add_promoter_usage <- function(scan, flank = 500, dir = NULL, min_tpm = 1,
                               dominant_frac = 0.6) {
  m <- scan$meta; prom <- scan$promoters
  if (!identical(m$species, "homo_sapiens")) {
    message("  promoter usage: FANTOM5 hg38 table is human only - skipped.")
    return(scan)
  }
  if (is.null(prom)) return(scan)
  ch <- .fantom5_chunk(paste0("chr", sub("^chr", "", m$chr)),
                       if (m$strand == 1) "+" else "-", dir)
  if (is.null(ch)) return(scan)
  pk <- ch$peaks; ns <- ncol(ch$tpm100)
  hits <- lapply(seq_len(nrow(prom)), function(i)
    which(pk$summit >= prom$tss_min[i] - flank & pk$summit <= prom$tss_max[i] + flank))
  tpm <- matrix(vapply(hits, function(h)
    if (length(h)) colSums(ch$tpm100[h, , drop = FALSE]) / 100 else numeric(ns),
    numeric(ns)), nrow = length(hits), byrow = TRUE)
  rownames(tpm) <- prom$promoter_id; colnames(tpm) <- ch$samples$label
  d <- ch
  tot  <- colSums(tpm)
  frac <- sweep(tpm, 2, pmax(tot, 1e-9), "/")
  top  <- rownames(tpm)[apply(tpm, 2, which.max)]
  topf <- apply(frac, 2, max)
  dominant <- ifelse(tot < min_tpm, "not expressed",
                     ifelse(topf >= dominant_frac, top, "mixed"))
  tab <- data.frame(sample_type = d$samples$sample_type, category = d$samples$category,
                    .classify_samples(d$samples$sample_type, d$samples$category),
                    total_tpm = round(tot, 2), dominant = dominant,
                    stringsAsFactors = FALSE)
  for (p in rownames(tpm)) {
    tab[[paste0(p, "_tpm")]]  <- round(tpm[p, ], 2)
    tab[[paste0(p, "_frac")]] <- round(frac[p, ], 3)
  }
  tab <- tab[order(-tab$total_tpm), ]; rownames(tab) <- NULL
  cnt <- table(factor(dominant, levels = c(rownames(tpm), "mixed", "not expressed")))
  message(sprintf("  promoter usage (FANTOM5 CAGE, %d sample types): %s",
                  ncol(tpm), paste(sprintf(ifelse(names(cnt) %in% rownames(tpm), "%s dominant in %d", "%s %d"),
                                           names(cnt), as.integer(cnt)), collapse = " | ")))
  if (sum(cnt[rownames(tpm)] > 0) >= 2)
    message("    -> promoter switching: different promoters dominate in different cell types.")
  scan$usage <- list(tpm = tpm, frac = frac, table = tab, flank = flank,
                     n_peaks = vapply(hits, length, integer(1)), min_tpm = min_tpm,
                     dominant_frac = dominant_frac)
  scan
}

## NEW plot: where each promoter is used. For every promoter, the sample types
## in which it carries the largest share of the gene's CAGE signal (gene
## expressed >= min_tpm), plus the sample types with the highest total
## expression. Bars = share of each promoter (fill colours as in the promoter
## map); numbers = total TPM of the gene in that sample type.
#' Plot promoter usage across cell types
#'
#' For each promoter, the sample types in which it carries the largest share of
#' the gene's CAGE signal, plus the sample types with the highest expression.
#' @param scan A scan with usage (\code{\link{add_promoter_usage}}).
#' @param top_n Sample types per group (default 10).
#' @param base_size Base font size.
#' @return A ggplot, or \code{NULL}.
#' @export
plot_promoter_usage <- function(scan, top_n = 10, base_size = 12) {
  u <- scan$usage; m <- scan$meta
  if (is.null(u)) return(NULL)
  tab <- u$table; pids <- rownames(u$tpm)
  ex  <- tab[tab$total_tpm >= u$min_tpm, , drop = FALSE]
  if (!nrow(ex)) return(NULL)
  pick <- list()
  for (p in pids) {
    e <- ex[ex[[paste0(p, "_tpm")]] >= u$min_tpm, , drop = FALSE]
    e <- e[order(-e[[paste0(p, "_frac")]], -e$total_tpm), , drop = FALSE]
    if (nrow(e)) pick[[p]] <- data.frame(sample_type = utils::head(e$sample_type, top_n),
                                         category = utils::head(e$category, top_n),
                                         group = sprintf("%s used most", p),
                                         stringsAsFactors = FALSE)
  }
  hi <- utils::head(ex, top_n)
  pick$top <- data.frame(sample_type = hi$sample_type, category = hi$category,
                         group = "highest expression", stringsAsFactors = FALSE)
  sel <- do.call(rbind, pick)
  sel <- sel[!duplicated(paste(sel$sample_type, sel$category)), ]
  rows <- merge(sel, tab, by = c("sample_type", "category"), sort = FALSE)
  cat_short <- c("primary cells" = "primary cell", "tissues" = "tissue", "cell lines" = "cell line")
  rows$label <- sprintf("%s [%s]", rows$sample_type, cat_short[rows$category] %||% rows$category)
  long <- do.call(rbind, lapply(pids, function(p) data.frame(
    label = rows$label, group = rows$group, promoter = p,
    frac = rows[[paste0(p, "_frac")]], total = rows$total_tpm, stringsAsFactors = FALSE)))
  first <- pids[1]
  ord <- rows$label[order(rows$group, rows[[paste0(first, "_frac")]])]
  long$label    <- factor(long$label, levels = unique(ord))
  long$promoter <- factor(long$promoter, levels = pids)
  long$group    <- factor(long$group, levels = c(sprintf("%s used most", pids), "highest expression"))
  pr <- scan$promoters
  plab <- setNames(sprintf("%s%s (%d tx)", pr$promoter_id, ifelse(pr$has_canonical, "*", ""),
                           pr$n_transcripts), pr$promoter_id)
  totals <- unique(long[, c("label", "group", "total")])
  ggplot(long, aes(frac, label)) +
    geom_col(aes(fill = promoter), width = 0.75) +
    geom_text(data = totals, aes(x = 1.02, y = label, label = sprintf("%.0f", total)),
              hjust = 0, size = base_size / 4.2, color = "grey30", inherit.aes = FALSE) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_brewer(palette = "Set2", labels = plab, name = "Promoter") +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                       breaks = c(0, .25, .5, .75, 1), limits = c(0, 1.12),
                       expand = c(0, 0)) +
    labs(x = "Share of the gene's CAGE signal", y = NULL,
         title = sprintf("Where the %s promoters are used (FANTOM5 CAGE)", m$gene),
         subtitle = sprintf("CAGE peaks within %d bp of each promoter's TSSs | numbers = total TPM | %d of %d sample types express %s (>= %g TPM) | * = canonical",
                            u$flank, nrow(ex), nrow(tab), m$gene, u$min_tpm)) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.72)),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey92", color = NA))
}

## ---------------------------------------------------------------------------
## Normal vs cancer promoter usage (FANTOM5). "Cancer" = cancer cell lines
## (carcinoma, leukemia, glioma, ... in the sample name); "normal" = primary
## cells and tissues; non-cancer cell lines (HEK293, iPS, ...) are "other".
## Samples are given a lineage by keywords so that cancer lines are compared
## with normal samples of the same lineage (colon lines vs colon / intestine,
## leukemias vs blood cells, gliomas vs brain, ...), not with all tissues.
## Caveat: FANTOM5 cancer samples are cell lines, not tumours - culture
## adaptation is part of any difference.
## ---------------------------------------------------------------------------
.cancer_words <- paste0("carcinoma|cancer|leuka?emia|lymphoma|myeloma|sarcoma|melanoma|",
                        "glioma|glioblastoma|astrocytoma|blastoma|mesothelioma|tumou?r|",
                        "seminoma|carcinoid|hodgkin|malignant|neoplasm|myelodysplastic")
.lineages <- list(   # first match wins - specific before broad
  blood      = "leuka?emia|lymphoma|myeloma|hodgkin|myelodysplastic|blood|monocyte|macrophage|t cell|b cell|nk cell|natural killer|neutrophil|lymphocyte|erythro|dendritic|basophil|eosinophil|mast cell|bone marrow|thymus|spleen|lymph node|tonsil|myeloid|megakaryo|cd[0-9]+",
  brain      = "glioma|glioblastoma|astrocytoma|neuroblastoma|medulloblastoma|meningioma|brain|astrocyte|neuron|neural|cortex|gyrus|cerebell|hippocamp|thalamus|amygdala|nucleus|lobe|medulla|spinal|substantia|putamen|pons|caudate|globus|insula|schwann|locus coeruleus|optic|pineal|mening|olfactory|corpus callosum|pituitary",
  liver      = "liver|hepat|cholangio|bile duct|gall",
  colon      = "colon|colorectal|rectal|rectum|intestin|duoden|caecum|cecum|appendix|ileum|jejun",
  stomach    = "stomach|gastric",
  esophagus  = "esophag|oesophag",
  lung       = "lung|bronch|alveol|trachea|pleura|mesothelioma",
  breast     = "breast|mammary",
  kidney     = "kidney|renal|wilms",
  prostate   = "prostat",
  pancreas   = "pancrea|islet",
  skin       = "skin|melanoma|melanocyte|keratinocyte|epiderm|dermal",
  female_reproductive = "ovar|uter|endometri|cervi|vagin|placenta|chorio|trophoblast|amnion",
  testis     = "testi|seminoma|germ cell|teratocarcinoma|sertoli",
  bone_soft_tissue = "osteosarcoma|ewing|rhabdo|leiomyo|liposarcoma|fibrosarcoma|chondro|osteo|bone|muscle|myoblast|myotube|skeletal|smooth muscle|synovi|tendon|fibroblast|adipo|preadipo|mesenchymal",
  thyroid    = "thyroid",
  bladder    = "bladder|urothel",
  head_neck  = "oral|tongue|pharyn|laryn|saliva|gingiv|tonsil|nasal",
  eye        = "retin|eye|lens|cornea|uveal",
  endothelium = "endotheli|huvec|vascular|vein|artery|aort|lymphatic")

.classify_samples <- function(sample_type, category) {
  st <- tolower(sample_type)
  condition <- ifelse(grepl(.cancer_words, st), "cancer",
               ifelse(category %in% c("primary cells", "tissues"), "normal", "other"))
  lineage <- rep(NA_character_, length(st))
  for (nm in names(.lineages)) {
    hit <- is.na(lineage) & grepl(.lineages[[nm]], st)
    lineage[hit] <- nm
  }
  data.frame(condition = condition, lineage = lineage, stringsAsFactors = FALSE)
}

## NEW plot: normal vs cancer use of every promoter, lineage-matched. One row
## per promoter: share of the gene's CAGE signal from that promoter;
## bottom: the gene's total expression. Points = sample types (gene >= min_tpm);
## p = two-sided Wilcoxon test (shown when both groups have >= 3 samples).
#' Normal vs cancer promoter usage (lineage-matched)
#'
#' Share of every promoter (and total expression) in normal samples (primary
#' cells, tissues) vs cancer cell lines of the same lineage, with two-sided
#' Wilcoxon tests. FANTOM5 cancer samples are cell lines, not tumours.
#' @param scan A scan with usage (\code{\link{add_promoter_usage}}).
#' @param min_n Minimum sample types per group and lineage (default 3).
#' @param base_size Base font size.
#' @return A ggplot, or \code{NULL} if no lineage has enough samples.
#' @export
plot_promoter_usage_cancer <- function(scan, min_n = 3, base_size = 12) {
  u <- scan$usage; m <- scan$meta
  if (is.null(u) || is.null(u$table$condition)) return(NULL)
  pr   <- scan$promoters
  alts <- pr$promoter_id                                  # all promoters
  plab <- setNames(sprintf("share of %s%s (%d tx)", pr$promoter_id,
                           ifelse(pr$has_canonical, "*", ""), pr$n_transcripts),
                   pr$promoter_id)
  d <- u$table[u$table$total_tpm >= u$min_tpm & u$table$condition %in% c("normal", "cancer") &
               !is.na(u$table$lineage), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  ## lineages with >= min_n samples in both groups
  cnt  <- table(d$lineage, d$condition)
  keep <- rownames(cnt)[cnt[, "normal"] >= min_n & cnt[, "cancer"] >= min_n]
  if (!length(keep)) return(NULL)
  d <- d[d$lineage %in% keep, , drop = FALSE]
  d_all <- d; d_all$lineage <- "ALL (matched)"
  d <- rbind(d, d_all)
  long <- do.call(rbind, lapply(alts, function(p) data.frame(
    lineage = d$lineage, condition = d$condition, sample_type = d$sample_type,
    panel = plab[[p]], value = d[[paste0(p, "_frac")]], stringsAsFactors = FALSE)))
  long <- rbind(long, data.frame(lineage = d$lineage, condition = d$condition,
                                 sample_type = d$sample_type, panel = "total TPM (log10)",
                                 value = log10(d$total_tpm), stringsAsFactors = FALSE))
  lv <- c(sort(keep), "ALL (matched)")
  long$lineage   <- factor(long$lineage, levels = lv)
  long$condition <- factor(long$condition, levels = c("normal", "cancer"))
  long$panel     <- factor(long$panel, levels = c(unname(plab[alts]), "total TPM (log10)"))
  ## Wilcoxon per lineage x panel
  st <- do.call(rbind, lapply(split(long, list(long$lineage, long$panel), drop = TRUE), function(g) {
    a <- g$value[g$condition == "normal"]; b <- g$value[g$condition == "cancer"]
    p <- if (length(a) >= min_n && length(b) >= min_n)
      tryCatch(suppressWarnings(stats::wilcox.test(a, b)$p.value), error = function(e) NA) else NA
    data.frame(lineage = g$lineage[1], panel = g$panel[1], p = p,
               y = max(g$value, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
  st <- st[!is.na(st$p), ]
  st$label <- ifelse(st$p < 0.001, "p<0.001", sprintf("p=%.2g", st$p))
  ggplot(long, aes(lineage, value, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, position = position_dodge(0.8), width = 0.7) +
    geom_point(aes(color = condition), position = position_jitterdodge(jitter.width = 0.15,
               dodge.width = 0.8, seed = 1), size = 1.2, alpha = 0.8) +
    geom_text(data = st, aes(lineage, y, label = label), inherit.aes = FALSE,
              vjust = -0.6, size = base_size / 4.2, color = "grey25") +
    facet_grid(panel ~ ., scales = "free_y") +
    scale_fill_manual(values = c(normal = "#4C9F70", cancer = "#C8553D")) +
    scale_color_manual(values = c(normal = "#2E6B47", cancer = "#8E3324")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = NULL, y = NULL, fill = NULL, color = NULL,
         title = sprintf("%s promoter use: normal vs cancer, lineage-matched (FANTOM5 CAGE)", m$gene),
         subtitle = sprintf("normal = primary cells + tissues; cancer = cancer cell lines | points = sample types expressing %s (>= %g TPM) | lineages with >= %d per group | Wilcoxon | * = canonical",
                            m$gene, u$min_tpm, min_n)) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.7)),
          strip.text.y = element_text(angle = 0),
          strip.background = element_rect(fill = "grey92", color = NA),
          legend.position = "top")
}
