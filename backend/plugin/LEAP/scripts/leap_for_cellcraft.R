suppressPackageStartupMessages(library(LEAP, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))

args <- commandArgs(trailingOnly = T)

fpath_exp_data <- args[1]  # expression file path
fpath_trj_data <- args[2]  # pseudotime file path
fpath_branch_data <- args[3]  # cell select file path

maxLag <- as.numeric(args[4])  # default 0.33

fdr <- as.numeric(args[5])  # OPTIONAL, default: 0, min: 0
links <- as.numeric(args[6])  # default 1000, min: 0, max: 1

fname_grn <- args[7]  # output file name of grn
fname_odgs <- args[8]  # output file name of outdegrees

exprMatrix <- read.delim(fpath_exp_data, row.names = 1, sep = ',', check.names=FALSE)  # cell x gene
pseudotimeData <- read.delim(fpath_trj_data, header=FALSE)
cellSelectData <- read.delim(fpath_branch_data, header=FALSE)

geneNames <- colnames(exprMatrix)
cellNames <- rownames(exprMatrix)

exprMatrix <- data.frame(lapply(exprMatrix, as.numeric))

rownames(exprMatrix) <- cellNames
nonzeroIndices <- which(cellSelectData[, 1] != 0)
exprMatrix <- exprMatrix[nonzeroIndices, ]
pseudotimeData <- pseudotimeData[nonzeroIndices, ]
exprMatrix <- exprMatrix[order(pseudotimeData), ]  # cell x gene

exprMatrix <- t(exprMatrix)
exprMatrix <- as.data.frame(exprMatrix)  # gene x cell

# # input expression data
geneNames <- rownames(exprMatrix)
rownames(exprMatrix) <- c()

# Run LEAP's compute Max. Absolute Correlation
# MAC_cutoff is set to zero to get a score for all TFs
# max_lag_prop is set to the max. recommended value from the paper's supplementary file
# Link to paper: https://academic.oup.com/bioinformatics/article/33/5/764/2557687

MAC_results = MAC_counter(data = exprMatrix, max_lag_prop=maxLag, MAC_cutoff = 0,
                          file_name = "temp", lag_matrix = FALSE, symmetric = FALSE)

# Write output to a file
source <- geneNames[MAC_results[,'Row gene index']]
target <- geneNames[MAC_results[,'Column gene index']]
weight <- MAC_results[,'Correlation']
weight <- abs(weight)


edges <- data.frame(source, weight, target)

res <- edges %>%
  group_by(source, target) %>%
  filter(weight == max(weight)) %>%
  ungroup()

finalDF <- res %>%
  arrange(desc(weight))

if (fdr != 0) {
  finalDF <- as.data.frame(finalDF)
  zscore <- (finalDF$weight - mean(finalDF$weight)) / sd(finalDF$weight)
  pval <- 1 - pnorm(zscore)
  computed_fdr <- p.adjust(pval, method = "fdr")

  inds_cutoff <- computed_fdr < fdr
  finalDF <- finalDF[inds_cutoff, ]


} else {
  finalDF <- finalDF[1:links, ]
}

outdegree_count <- as.data.frame(table(finalDF$source))

write.table(finalDF, fname_grn, sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

colnames(outdegree_count) <- c("Node", "Outdegree")

outdegree_count <- outdegree_count %>%
  arrange(desc(Outdegree))

write.table(outdegree_count, file = fname_odgs,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
