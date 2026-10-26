#!/usr/bin/env Rscript
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# Exports, for each target class (GBM_RTK_II, A_IDH, O_IDH), the binarised methylation
# profile of ONE real representative sample of the Capper et al. cohort (the one whose
# profile is closest to the centroid of its class), at the hg38 positions of the probes
# used by methylation_classification.R (top_probes ∩ HM450 manifest, autosomes
# only: the classification script excludes chrX/chrY).
#
# Usage: Rscript 01_export_probe_targets.R <data_dir> <out_dir>
suppressPackageStartupMessages({library(GenomicRanges)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 01_export_probe_targets.R <data_dir> <out_dir>", call. = FALSE)
data_dir <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(data_dir, "top_probes_hm450.Rdata"))                 # top_probes
load(file.path(data_dir, "capper_top_100k_betas_binarised.Rdata"))  # betas
load(file.path(data_dir, "HM450.hg38.manifest.gencode.v22.Rdata"))  # hm450
# NB: updateObject(hm450) fails here (package 'Seqinfo' not found, probably a Bioconductor
# version mismatch between the serialised object and the installed apt packages).
# The GRanges accessors work directly without this call.
hm450 <- hm450[!as.character(seqnames(hm450)) %in% c("chrX", "chrY")]

probes <- intersect(top_probes, names(hm450))
probes <- intersect(probes, colnames(betas))
cat("Probes with an autosomal hg38 coordinate AND present in top_probes/betas:", length(probes), "/", length(top_probes), "\n")

gr <- hm450[probes]
pos_tbl <- data.frame(
  probe  = names(gr),
  chrom  = as.character(seqnames(gr)),
  pos0   = start(gr) - 1L,   # BED 0-based
  pos1   = end(gr),
  strand = as.character(strand(gr))
)

targets <- list(gbm = "GBM_RTK_II", astro = "A_IDH", oligo = "O_IDH")
for (prof in names(targets)) {
  dx <- targets[[prof]]
  sub <- betas[betas$Dx == dx, probes, drop = FALSE]
  if (nrow(sub) == 0) {
    cat(sprintf("WARNING: no sample of class %s in the training set: %s skipped\n", dx, prof))
    next
  }
  sub <- as.matrix(sub)
  storage.mode(sub) <- "numeric"
  centroid <- colMeans(sub, na.rm = TRUE)
  # real sample closest to the centroid of its class (L1 distance; on binary profiles,
  # equivalent to the Hamming distance to the majority profile)
  d <- rowSums(abs(sweep(sub, 2, centroid, "-")), na.rm = TRUE)
  best <- which.min(d)
  cat(sprintf("%-6s class=%-12s n_samples=%-4d representative=%s (distance to centroid=%.0f/%d probes)\n",
              prof, dx, nrow(sub), rownames(sub)[best], d[best], length(probes)))
  vals <- sub[best, ]
  out <- cbind(pos_tbl, target = as.integer(vals[pos_tbl$probe]))
  out <- out[!is.na(out$target), ]
  write.table(out, file.path(out_dir, paste0("probe_targets_", prof, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  -> %d probes written (%.0f%% methylated)\n", nrow(out), 100 * mean(out$target)))
}
