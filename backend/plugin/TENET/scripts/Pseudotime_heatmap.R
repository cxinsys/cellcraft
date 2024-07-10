#!/usr/bin/env Rscript

library(ggplot2)
library(pheatmap)
library(plotly)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
input_expression_file <- args[1]
input_trajectory_file <- args[2]
output_file_path <- args[3]

pseudotime_heatmap <- function(matrix, selected_gene, gene_list,
                               pseudotime, span=0.5, use_pseudotime_origin = T, use_z_score = T,
                               order_pseudo_score=F, p_min=-0.7, p_max=0.7, p_legend=T,
                               max_min = T, out_result = F, target_average = F, target, pseudo_decrease=T) {
  #====================subset & order by pseudotime =================
  
  sub_matrix <- matrix[, selected_gene]
  sub_matrix <- cbind(pseudotime, sub_matrix)
  sub_matrix <- subset(sub_matrix, sub_matrix[, 1] != "Inf" & sub_matrix[, 1] != "NA")
  sub_matrix <- as.data.frame(sub_matrix)
  
  sub_matrix_order <- order(sub_matrix[, 1])
  
  sub_matrix_sorted <- sub_matrix[sub_matrix_order, ]
  
  sub_matrix_sorted <- as.matrix(sub_matrix_sorted)
  
  #====================use origin pseudotime or use order pseudotime =================
  if(use_pseudotime_origin == T) {
  } else {
    expanded_pseudotime <- data.frame(1:dim(sub_matrix_sorted)[1])
    sub_matrix_sorted[, 1] <- expanded_pseudotime[, 1]
  }
  
  #=================== target average heatmap =================
  if (target_average == T) {
    average_matrix <- sub_matrix_sorted[, 1]
    for (j in 1:length(selected_gene)) {
      target_list <- subset(gene_list, gene_list[, 1] == gsub("[.]", "-", selected_gene[j]))[, 3]
      target_list <- gsub("-", ".", target_list)
      temp_matrix <- as.data.frame(matrix[, target_list])
      if (length(target_list) == 1) {
        colnames(temp_matrix) <- target_list
      }
      temp_matrix <- cbind(pseudotime, temp_matrix)
      temp_matrix <- subset(temp_matrix, temp_matrix[, 1] != "Inf" & temp_matrix[, 1] != "NA")
      temp_matrix <- as.data.frame(temp_matrix)
      temp_matrix_order <- order(temp_matrix[, 1])
      temp_matrix_sorted <- temp_matrix[temp_matrix_order, ]
      temp_matrix_sorted <- as.matrix(temp_matrix_sorted)
      if (length(target_list) != 1) {
        temp_average <- apply(temp_matrix_sorted[, -1], 1, mean)
      } else {
        temp_average <- temp_matrix_sorted[, 2]
      }
      temp_average <- as.data.frame(temp_average)
      colnames(temp_average) <- paste0(selected_gene[j], " (", dim(temp_matrix_sorted)[2] - 1, ")")
      average_matrix <- cbind(average_matrix, temp_average)
    }
    sub_matrix_sorted <- average_matrix
  }
  
  #====================z-score =================
  if(use_z_score == T) {
    for (i in 2:dim(sub_matrix_sorted)[2]) {
      sub_matrix_sorted[, i] <- (sub_matrix_sorted[, i] - mean(sub_matrix_sorted[, i])) / sd(sub_matrix_sorted[, i])
    }
  }
  
  #====================normalization =================
  temp <- cbind(sub_matrix_sorted[, 1], apply(sub_matrix_sorted[, -1], 2, function(j) {
    k <- sub_matrix_sorted[, 1]
    lo <- loess(j ~ k, span = span)
    xl <- seq(min(k), max(k), (max(k) - min(k)) / (length(k) - 1))
    predict(lo, xl)
  }))
  colnames(temp)[1] <- "pseudotime"
  temp <- t(temp)
  
  normalized_sub_matrix_sorted <- temp
  
  #====================pseudo_score =================
  if(order_pseudo_score == T) {
    pseudo_score <- data.frame()
    
    for (i in 2:dim(normalized_sub_matrix_sorted)[1]) {
      pseudo_score[i, 1] <- sum(normalized_sub_matrix_sorted[1, ] * normalized_sub_matrix_sorted[i, ])
    }
    
    normalized_sub_matrix_sorted_pseudo_ordered <- cbind(pseudo_score, normalized_sub_matrix_sorted)
    rownames(normalized_sub_matrix_sorted_pseudo_ordered) <- rownames(normalized_sub_matrix_sorted)
    temp <- normalized_sub_matrix_sorted_pseudo_ordered[-1, ]
    temp <- temp[order(temp[, 1], decreasing = pseudo_decrease), ]
    normalized_sub_matrix_sorted_pseudo_ordered <- rbind(normalized_sub_matrix_sorted_pseudo_ordered[1, ], temp)
    final_matrix <- normalized_sub_matrix_sorted_pseudo_ordered
    pheat_start <- 2
  } else {
    final_matrix <- normalized_sub_matrix_sorted
    pheat_start <- 1
  }
  
  #====================make pheatmap =================
  if(max_min == T) {
    s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]],
                   breaks = seq(p_min, p_max, length.out=100), show_colnames=F, cluster_rows = F, cluster_cols = F, legend = p_legend)
  } else {
    s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]],
                   show_colnames=F, cluster_rows = F, cluster_cols = F, legend = p_legend)
  }
  
  #====================Plotly 변환====================
  # Ensure 'z' is a matrix
  z_matrix <- as.matrix(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]])
  
  fig <- plot_ly(
    z = z_matrix, 
    x = colnames(final_matrix)[pheat_start:dim(final_matrix)[2]], 
    y = rownames(final_matrix)[2:dim(final_matrix)[1]],
    type = "heatmap",
    colorscale = "Viridis"
  )
  
  # 올바른 title 형식 지정
  layout <- list(
    title = "Heatmap",
    xaxis = list(title = "X Axis Title"),
    yaxis = list(title = "Y Axis Title"),
    scene = list(zaxis = list(title = "Z Axis Title"))
  )
  
  fig <- fig %>% layout(layout)
  
  # JSON 형식으로 변환
  fig_json <- plotly_json(fig, jsonedit = FALSE)
  
  # frame 키를 name으로 변경하고, title 값을 문자열로 변환
  fig_json <- gsub('"frame":', '"name":', fig_json)
  fig_json <- gsub('"title":\\[\\]', '"title":""', fig_json)
  
  write(fig_json, output_file_path)
}

# 입력 파일 경로
expression_data <- read.csv(input_expression_file)
trajectory <- read.table(input_trajectory_file, quote="\"", comment.char="")

# 고정값으로 사용되는 예시 데이터
selected_gene <- c("Gene1", "Gene2", "Gene3")  # 예시 값, 실제 사용 시 수정 필요
gene_list <- data.frame(Gene = c("Gene1", "Gene2", "Gene3"), TF = c("TF1", "TF2", "TF3"))  # 예시 값, 실제 사용 시 수정 필요

# 함수 호출
pseudotime_heatmap(matrix = expression_data, gene_list = gene_list,
                   selected_gene = selected_gene, pseudotime = trajectory,
                   span=0.7, use_pseudotime_origin = T, use_z_score = T,
                   order_pseudo_score=T, max_min = T, p_min=-0.7, p_max=0.7, p_legend=T, 
                   out_result =F, target_average = T, pseudo_decrease=T)
