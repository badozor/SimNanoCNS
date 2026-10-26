#!/usr/bin/env Rscript
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# Exports, for each target class (GBM_RTK_II, A_IDH, O_IDH, A_IDH_HG), N_PER_CLASS binarised
# methylation profiles taken from REAL, DISTINCT samples of the Capper et al. cohort
# (not always the same representative as in 01_export_probe_targets.R): this is what gives
# a simulated cohort real inter-individual diversity rather than identical copies.
#
# Usage: Rscript 02_export_cohort_probe_targets.R <data_dir> <out_dir> <n_per_class> [seed]
suppressPackageStartupMessages({library(GenomicRanges)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript 02_export_cohort_probe_targets.R <data_dir> <out_dir> <n_per_class> [seed]", call. = FALSE)
data_dir <- args[1]; out_dir <- args[2]
n_per_class <- as.integer(args[3])
if (is.na(n_per_class) || n_per_class < 1) stop("<n_per_class> must be a positive integer", call. = FALSE)
seed <- if (length(args) >= 4) as.integer(args[4]) else 2026
set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(data_dir, "top_probes_hm450.Rdata"))                 # top_probes
load(file.path(data_dir, "capper_top_100k_betas_binarised.Rdata"))  # betas
load(file.path(data_dir, "HM450.hg38.manifest.gencode.v22.Rdata"))  # hm450
# NB: updateObject(hm450) fails here (see 01_export_probe_targets.R); direct subsetting, with
# no impact since only the coordinates/strands/names are needed, not the Seqinfo slot.
hm450 <- hm450[!as.character(seqnames(hm450)) %in% c("chrX", "chrY")]

probes <- intersect(top_probes, names(hm450))
probes <- intersect(probes, colnames(betas))
cat("Common probes (top_probes ∩ autosomal hg38 manifest ∩ betas):", length(probes), "\n")

gr <- hm450[probes]
pos_tbl <- data.frame(probe = names(gr), chrom = as.character(seqnames(gr)),
                       pos0 = start(gr) - 1L, pos1 = end(gr), strand = as.character(strand(gr)))

targets <- list(gbm = "GBM_RTK_II", astro = "A_IDH", oligo = "O_IDH",
                 astro_hg = "A_IDH_HG")  # grade 4 subgroup (--cdkn2a-del), see build_cohort.py
manifest <- data.frame()
for (prof in names(targets)) {
  dx <- targets[[prof]]
  sub <- as.matrix(betas[betas$Dx == dx, probes, drop = FALSE])
  storage.mode(sub) <- "numeric"
  n_avail <- nrow(sub)
  if (n_avail == 0) {
    cat(sprintf("WARNING: no sample of class %s in the training set: %s skipped\n", dx, prof))
    next
  }
  if (n_per_class > n_avail) {
    cat(sprintf("WARNING: %d individuals requested for %s but only %d real samples available -> reuse (sampling with replacement) beyond %d.\n",
                n_per_class, dx, n_avail, n_avail))
    idx <- c(sample.int(n_avail), sample.int(n_avail, n_per_class - n_avail, replace = TRUE))
  } else {
    idx <- sample.int(n_avail, n_per_class)
  }
  for (k in seq_len(n_per_class)) {
    real_id <- rownames(sub)[idx[k]]
    vals <- sub[idx[k], ]
    out <- cbind(pos_tbl, target = as.integer(vals[match(pos_tbl$probe, probes)]))
    out <- out[!is.na(out$target), ]
    fname <- sprintf("probe_targets_%s_%02d.tsv", prof, k)
    write.table(out, file.path(out_dir, fname), sep = "\t", quote = FALSE, row.names = FALSE)
    manifest <- rbind(manifest, data.frame(profile = prof, index = k, dx_class = dx,
                                            real_sample_id = real_id, file = fname,
                                            n_probes = nrow(out), pct_methyl = round(100 * mean(out$target))))
  }
  cat(sprintf("%-8s: %d files written (class %s, %d real samples available)\n",
              prof, n_per_class, dx, n_avail))
}
write.table(manifest, file.path(out_dir, "cohort_probe_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Manifest:", file.path(out_dir, "cohort_probe_manifest.tsv"), "\n")
