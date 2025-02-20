suppressPackageStartupMessages(library(monocle, warn.conflicts = FALSE , quietly = TRUE))
suppressPackageStartupMessages(library(Scribe, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (igraph, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (dplyr, warn.conflicts = FALSE, quietly = TRUE))

cal_ncenter <- function(ncells, ncells_limit = 100){
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}

args <- commandArgs(trailingOnly = TRUE)

fpath_exp_data <- args[1]  # expression file path
fpath_trj_data <- args[2]  # pseudotime file path
fpath_branch_data <- args[3]  # cell select file path

fpath_newcelldataset <- FALSE  # OPTIONAL, .RDS file path, default: FALSE

lower_limit <- as.numeric(args[4])  # float, defualt 0.0
expression_family <- args[5]  # default: 'uninormal' if not use negbinomial.size() automatically
method <- args[6]  # default: ucRDI

delay <- args[7]  # default: 1, min: 1
log <- args[8]  # default: FALSE, type: chracter

# ignore_pt <- args[10]  # default: FALSE, type: chracter

links <- as.numeric(args[9])  # default: 1000, min: 0
fdr <- as.numeric(args[10])  # OPTIONAL, default 0, min: 0, max: 1

fname_edges <- args[11]  # output file name of ranked edges, default: ranked_edges.txt
fname_grn <- args[12]  # output file name of grn
fname_odgs <- args[13]  # output file name of outdegrees

# if the new cell dataset file is not provided create one
if (fpath_newcelldataset == FALSE){

  if (length(fpath_exp_data) == 0){
      stop("Please enter the path to expression file")
    }

  if (length(fpath_trj_data) == 0){
      stop("Please enter the path to pseudotime file.")
    }

  if (length(fpath_branch_data) == 0){
    stop("Please enter the path to cell select file")
  }

  # Read data
  # Set check.names to False to avoid R adding an 'X' to the beginning of columns that start with an integer
  exprMatrix <- read.delim(fpath_exp_data, row.names = 1, sep = ',', check.names=FALSE)  # cell x gene
  pseudotimeData <- read.delim(fpath_trj_data, header=FALSE)
  cellSelectData <- read.delim(fpath_branch_data, header=FALSE)

  geneNames <- colnames(exprMatrix)
  cellNames <- rownames(exprMatrix)

  exprMatrix <- data.frame(lapply(exprMatrix, as.numeric))

  rownames(exprMatrix) <- cellNames
  colnames(exprMatrix) <- geneNames

  geneData <- data.frame(gene_name=geneNames, stringsAsFactors = FALSE)
  colnames(geneData)[colnames(geneData) == "gene_name"] <- "gene_short_name"
  rownames(geneData) <- geneNames

  nonzeroIndices <- which(cellSelectData[, 1] != 0)
  exprMatrix <- exprMatrix[nonzeroIndices, ]
  pseudotimeData <- pseudotimeData[nonzeroIndices, ]
  cellData <- as.data.frame(pseudotimeData, header=FALSE)

  cellNames <- rownames(exprMatrix)
  rownames(cellData) <- cellNames
  colnames(cellData)[colnames(cellData) == "pseudotimeData"] <- "Time"

  exprMatrix <- t(exprMatrix)
  exprMatrix <- as.data.frame(exprMatrix)  # gene x cell

  # # Min-Max 정규화 함수 적용
  # exprMatrix <- t(apply(exprMatrix, 1, function(x) {
  #   (x - min(x)) / (max(x) - min(x))
  # }))

  exprMatrix <- as.data.frame(exprMatrix)

  cd <- new("AnnotatedDataFrame", data = cellData)
  gd <- new("AnnotatedDataFrame", data = geneData)


  # Use uninormal if it is simulated data
  if (expression_family == 'uninormal'){
    cat("Using uninormal() as expression family.\n")
    CDS <- newCellDataSet(as(as.matrix(exprMatrix), "sparseMatrix"),
                        phenoData = cd,
                        featureData = gd,
                        lowerDetectionLimit = lower_limit,
                        expressionFamily = uninormal())

  sizeFactors(CDS) <- 1 # Same as the one in neuronal_sim_cCDS
  } else{
    # For scRNA-Seq data (counts, RPKM/FPKM)
    cat("Using negbinomial.size() as expression family.\n")
    CDS <- newCellDataSet(as(as.matrix(exprMatrix), "sparseMatrix"),
                          phenoData = cd,
                          featureData = gd,
                          lowerDetectionLimit = lower_limit,
                          expressionFamily = negbinomial.size())
    CDS <- estimateSizeFactors(CDS)
    CDS <- estimateDispersions(CDS)
    disp_table <- dispersionTable(CDS)
    ordering_genes <- disp_table$gene_id
    CDS <- setOrderingFilter(CDS, ordering_genes)
  }

  CDS$Pseudotime <- CDS$Time
  CDS@phenoData@data$Pseudotime <- CDS@phenoData@data$Time
  CDS@phenoData@data$State <- 1
  CDS$State <- 1

  # if (ignore_pt == TRUE){
  #   cat("Using experimental time instead of PseudoTime computed using
  #       monocle.\n")
  #   CDS$Pseudotime <- CDS$Time
  #   CDS@phenoData@data$Pseudotime <- CDS@phenoData@data$Time
  #   CDS@phenoData@data$State <- 1
  #   CDS$State <- 1
  #
  # }else{
  #
  # cat("Computing pseudotime.\n")
  # CDS <- reduceDimension(CDS, norm_method ="none")
  # CDS <- orderCells(CDS)
  #
  # saveRDS(CDS, file= 'dataset.RDS')
  # write.csv(CDS@phenoData@data, file= 'PseudoTime.csv', quote = FALSE)
  # }
}

# If newcelldataset RDS is available already
if (fpath_newcelldataset != FALSE){
  CDS <- readRDS(fpath_newcelldataset)
}

### Run scribe

cat("Computing",method,"\n")


delay <- as.numeric(strsplit(delay, ",")[[1]])

if (method == 'uRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = TRUE, log = log)
  netOut <- net$max_rdi_value
  # computes CLR if we use uRDI
  # TODO: Make this an optional
  netOut <- clr(netOut)
} else if (method == 'ucRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = TRUE, log = log)
  netOut <- calculate_conditioned_rdi(CDS, rdi_list = net, uniformalize = TRUE, log = log)
} else if (method == 'RDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = FALSE, log = log)
  netOut <- net$max_rdi_value
  # computes CLR if we use RDI
  # TODO: Make this optional
  netOut <- clr(netOut)
} else if (method == 'cRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = FALSE, log = log)
  netOut <- calculate_conditioned_rdi(CDS, rdi_list = net, uniformalize = FALSE, log = log)
} else{
  stop("Method must be one of RDI, cRDI, uRDI, or ucRDI.
       Run Rscript runScribe.R -h for more details.")
}

outGraph <- graph_from_adjacency_matrix(netOut, mode = 'directed', weighted = TRUE)

edges <- as.data.frame(as_edgelist(outGraph, names = TRUE))
weights <- E(outGraph)$weight
edges <- cbind(edges, weights)
colnames(edges) <- c("source", "target", "weight")
edges <- edges[, c("source", "weight", "target")]

write.table(edges, file = fname_edges,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

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

  outdegree_count <- as.data.frame(table(finalDF$source))
} else {
  finalDF <- finalDF[1:links, ]
  outdegree_count <- as.data.frame(table(finalDF$source))
}

write.table(finalDF, file = fname_grn,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

colnames(outdegree_count) <- c("Node", "Outdegree")

outdegree_count <- outdegree_count %>%
  arrange(desc(Outdegree))

write.table(outdegree_count, file = fname_odgs,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Done.\n")