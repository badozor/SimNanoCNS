# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
##### Rscript to perform methylation classification using a random forest ad-hoc model developed under Rapid-CNS2

# list all packages necessary
pkgs <- c(
  "optparse",
  "GenomicRanges",
  "ranger",
  "matrixStats",
  "data.table",
  # required for the serialised model to load, even though no glmnet function
  # is called directly here - do not remove
  "glmnet"
)

## load each one, hiding startup banners
invisible(
  lapply(
    pkgs,
    \(p) suppressPackageStartupMessages(
           library(p, character.only = TRUE, warn.conflicts = FALSE)
         )
  )
)

#Parse arguments
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample name, used as the prefix of the output files (required)", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="output directory (required)", metavar="character"),
  make_option(c("-i", "--in_file"), type="character", default=NULL,
              help="bedMethyl file (modkit format) of the sample (required)", metavar="character"),
  make_option(c("-p", "--probes_file"), type="character", default="top_probes_hm450.Rdata",
              help="top probes Rdata [default %default]", metavar="character"),
  make_option(c("-a", "--array_file"), type="character", default="HM450.hg38.manifest.gencode.v22.Rdata",
              help="HM450 hg38 manifest Rdata [default %default]", metavar="character"),
  make_option(c("-d", "--training_data"), type="character", default="capper_top_100k_betas_binarised.Rdata",
              help="Capper et al. training set (data/ directory of Rapid-CNS2_nf) [default %default]", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=16,
              help="number of threads [default %default]", metavar="integer"),
  make_option(c("-n", "--num_trees"), type="integer", default=20000,
              help="number of random forest trees (20000 in production; fewer for quick tests) [default %default]", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

for (req in c("sample", "out_dir", "in_file")) {
  if (is.null(opt[[req]])) {
    print_help(opt_parser)
    stop(sprintf("option --%s is required", req), call. = FALSE)
  }
}

#Create output directory
dir.create(file.path(opt$out_dir), showWarnings = FALSE, recursive = TRUE)

#Read methylation values, keeping only the relevant columns
# modkit bedMethyl (all-tab delimited) column layout:
#   1-3 chrom/start/end, 6 strand, 10 Nvalid_cov, 11 percent modified
# fread(select=) reads only these 6 of the 18 columns: a whole-genome bedMethyl has tens of
# millions of rows, and reading all columns with read.delim needed several GB of memory.
meth_filter <- as.data.frame(fread(opt$in_file, header = FALSE, sep = "\t",
                                   select = c(1, 2, 3, 10, 6, 11)))
# explicit order (chr, start, end, cov, strand, %), whatever the data.table version
meth_filter <- meth_filter[, c("V1", "V2", "V3", "V10", "V6", "V11")]

#Set colnames
colnames(meth_filter) <- c("chr","start","end","cov","strand","methylation_percent")

#Remove X and Y chromosomes
meth_filter <- subset(meth_filter, !(chr %in% c("chrX","chrY")))
meth_filter <- meth_filter[complete.cases(meth_filter), ]

## convert str to numeric
meth_filter$cov <- as.numeric(meth_filter$cov)
meth_filter$methylation_percent <- as.numeric(meth_filter$methylation_percent)

#Convert to decimal for comparison to array beta values
meth_filter$methylation <- meth_filter$methylation_percent/100
print("GRanges")
#Make GRanges object
cur38 <- makeGRangesFromDataFrame(meth_filter,keep.extra.columns=TRUE)
print("Load array")
#Load hg38 450K array sites
genome(cur38) <- "hg38"
load(opt$array_file)
# NB: neither the object as loaded nor updateObject(hm450) work in this
# environment for operations touching its Seqinfo slot (package 'Seqinfo'
# not found: Bioconductor version mismatch between the serialised object and the installed
# apt packages). A fresh GRanges is rebuilt with the current classes, from the only
# useful data (coordinates, strand, probe identifier) - checked separately.
hm450_names <- names(hm450)
hm450 <- GRanges(seqnames = as.character(seqnames(hm450)),
                  ranges = IRanges::IRanges(start(hm450), end(hm450)),
                  strand = as.character(strand(hm450)))
names(hm450) <- hm450_names
load(opt$probes_file)
#Identify overlaps with 450K array CpG probes
print("Identify overlaps")
index <- findOverlaps(hm450,cur38)
hm450.matched <- hm450[queryHits(index)]
mcols(hm450.matched) <- cbind.data.frame(mcols(hm450.matched),mcols(cur38[subjectHits(index)]))
methylation_sample <- hm450.matched[,c("cov","methylation")]
# the genome-wide CpG tables are no longer needed: free them before loading the training set
rm(meth_filter, cur38, index, hm450.matched)
invisible(gc())
print("Write file")
#Write GRanges object as RDS file
fileout <- paste0(opt$out_dir,"/",opt$sample,"_methylation_hg38_HM450.RDS")

saveRDS(methylation_sample,fileout)
print("Start methylation classification")
#####################################################################
# Methylation classification#
case <- methylation_sample
case$cpg <- names(case)
probes <- intersect(unique(names(case)), top_probes)
message(paste(length(probes),"  overlapping CpG sites between sample and reference set.",sep=""))
write.table(probes,file = paste0(opt$out_dir,"/",opt$sample,"_probes_for_training.csv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Read training data
load(opt$training_data)

probes <- intersect(unique(names(case)), colnames(betas))
probes <- c(probes,"Dx")
betas_train <- betas[,match(probes, colnames(betas))]
rm(betas)

#Select most variable probes
max_CpG=10000
sds <- colSds(as.matrix(betas_train[,-ncol(betas_train)]), na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n = min(ncol(betas_train)-1, max_CpG))
ts <- betas_train[,c(maxSDs,ncol(betas_train))]
rm(betas_train)
cols <- colnames(ts[,-ncol(ts)])

#Get fractions
ts$Dx <- as.factor(ts$Dx)
Dx_fractions <- min(summary(ts$Dx,maxsum=10000)) / summary(ts$Dx,maxsum=10000)

cat(sprintf("[%s] Starting ranger: %d trees, %d probes, %d samples, %d threads\n", Sys.time(), opt$num_trees, length(cols), nrow(ts), opt$threads))
rf <- ranger(dependent.variable.name = "Dx", data = ts[,c(cols,"Dx")], num.trees=opt$num_trees, probability = T, sample.fraction = Dx_fractions,verbose=TRUE,num.threads=opt$threads)
cat(sprintf("[%s] ranger finished\n", Sys.time()))

probs <- predict(rf, ts, predict.all = F)$predictions
# predicted class of each training sample (first class in case of a tie, so that
# pred always has one value per sample)
pred <- colnames(probs)[max.col(probs, ties.method = "first")]
classes <- levels(ts$Dx)

# A GLM is trained for each class individually (1 vs all)
glm_models <- lapply(classes, function(type){
  scores <- probs[, type]
  scale_this <- data.frame(class = ifelse(pred == type & pred == ts$Dx, 1, 0), score = scores)
  scaled_scores <- glm(class ~ score, scale_this, family = binomial)
  return(scaled_scores)
})

###### For paper
details <- paste(paste0("450k_overlap: ", length(unique(names(case)))),
                           paste0("Top_100k_overlap: ", length(probes)),
                           paste0("Brier_score: ", rf$prediction.error),
                           paste0("Training_probes: ", length(cols)),
                           sep = "\n")
write.table(details, file = paste0(opt$out_dir,"/",opt$sample,"_rf_details.tsv"), row.names = F, col.names = F, quote = F)


### predict case

# Prepare the case data for prediction
case_df <- as.data.frame(unique(case))
case_df <- case_df[match(cols, case_df$cpg), ]
case_df <- case_df[, c("methylation", "cpg")]
case_df$methylation <- ((case_df$methylation >= 0.5) + 0)
methylation_matrix <- matrix(case_df$methylation, nrow = 1, dimnames = list(NULL, cols))

# Predict probabilities for each class
predicted_probabilities <- predict(rf, rbind(methylation_matrix[, cols], methylation_matrix[, cols]), predict.all = FALSE)$predictions[1, ]
predicted_class <- names(predicted_probabilities)[which.max(predicted_probabilities)]
max_probability <- predicted_probabilities[predicted_class]

# Prepare votes data frame
votes_df <- data.frame(Frequency = predicted_probabilities)
votes_df$Frequency <- votes_df$Frequency / sum(votes_df$Frequency) * 100

# Apply recalibration to the predicted probabilities
recalibrated_scores_list <- lapply(levels(ts$Dx), function(class_label) {
  class_score <- predicted_probabilities[class_label]
  glm_model <- glm_models[[which(classes == class_label)]]
  recalibrated_score <- predict(glm_model, newdata = data.frame(score = class_score), type = "response")
  return(recalibrated_score)
})
recalibrated_scores <- unlist(recalibrated_scores_list)
names(recalibrated_scores) <- levels(ts$Dx)
normalized_recalibrated_scores <- recalibrated_scores / sum(recalibrated_scores)
final_recalibrated_score <- normalized_recalibrated_scores[predicted_class]

votes_df$Recalibrated_Frequency <- normalized_recalibrated_scores[rownames(votes_df)]
votes_df$Recalibrated_Frequency <- votes_df$Recalibrated_Frequency / sum(votes_df$Recalibrated_Frequency) * 100

votes_df <- votes_df[order(votes_df$Frequency), , drop = FALSE]

### Save calibration report

report <- paste(paste0("Number of features: ", rf$num.independent.variables),
                paste0("Predicted Class: ", predicted_class),
                paste0("Initial Score: ", max_probability),
                paste0("Calibrated Score: ", final_recalibrated_score),
                sep = "\n")
write.table(report, file = paste0(opt$out_dir,"/",opt$sample,"_calibrated_classification.tsv"), row.names = F, col.names = F, quote = F)

write.table(votes_df, file = paste0(opt$out_dir,"/",opt$sample,"_votes.tsv"), quote = F)
