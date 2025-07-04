#!/usr/bin/env Rscript

library(ggplot2)
library(pheatmap)
library(plotly)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

# 구분자 감지 함수 정의
detect_separator <- function(file_path) {
    first_line <- readLines(file_path, n=1)  # 첫 번째 줄만 읽기
    if (grepl(",", first_line)) {
        return(",")
    } else if (grepl("\t", first_line)) {
        return("\t")
    } else if (grepl(" ", first_line)) {
        return(" ")
    } else {
        return(NULL)  # 감지 실패 시 NULL 반환
    }
}

# 파일에서 적절한 구분자 감지
sep_detected <- detect_separator(args[3])

input_expression_file <- args[1]
input_trajectory_file <- args[2]
# 감지된 구분자를 기반으로 데이터 로드
if (!is.null(sep_detected)) {
    tenet_sif <- read.table(args[3], header=FALSE, sep=sep_detected, stringsAsFactors=FALSE)
} else {
    stop("구분자를 감지할 수 없습니다. 파일 형식을 확인하세요.")
}
output_file_path <- args[4]
top_num <- as.numeric(args[5])
sampled_num <- as.numeric(args[6])

# 파일 생성 후 JSON에서 빈 리스트([])를 빈 문자열("")로 변경하는 함수
replace_empty_titles_in_file <- function(file_path) {
  # 파일 읽기
  json_content <- readLines(file_path, warn = FALSE)
  
  # 파일 내용을 하나의 문자열로 합침
  json_content <- paste(json_content, collapse = "\n")
  
  # 빈 리스트를 빈 문자열로 대체 (title 값에 대해서만 처리)
  json_content <- gsub('"title":\\[\\]', '"title":""', json_content)
  
  # zaxis.title도 빈 리스트가 있으면 대체
  json_content <- gsub('"zaxis":\\{"title":\\[\\]\\}', '"zaxis":{"title":""}', json_content)
  json_content <- gsub('"xaxis":\\{"title":\\[\\]\\}', '"xaxis":{"title":""}', json_content)
  json_content <- gsub('"yaxis":\\{"title":\\[\\]\\}', '"yaxis":{"title":""}', json_content)

  # 변경된 내용을 다시 파일에 저장
  writeLines(json_content, file_path)
}

pseudotime_heatmap <- function(matrix, selected_gene, gene_list,
                               pseudotime, span=0.5, use_pseudotime_origin = T, use_z_score = T,
                               order_pseudo_score=F, p_min=-0.7, p_max=0.7, p_legend=T,
                               max_min = T, out_result = F, target_average = F, target, pseudo_decrease=T) {
  #====================subset & order by pseudotime =================
  
  # 데이터 길이 확인 및 조정
  n_samples <- length(pseudotime)
  
  sub_matrix <- matrix[, selected_gene]
  sub_matrix <- cbind(data.frame(pseudotime=pseudotime), sub_matrix)
  
  # NA와 Inf 값 제거
  valid_rows <- which(!is.infinite(sub_matrix$pseudotime) & !is.na(sub_matrix$pseudotime))
  sub_matrix <- sub_matrix[valid_rows, ]
  
  # pseudotime으로 정렬
  sub_matrix <- sub_matrix[order(sub_matrix$pseudotime), ]
  
  # 모든 열의 길이가 동일하도록 보장
  n_rows <- nrow(sub_matrix)
  if (!use_pseudotime_origin) {
    sub_matrix$pseudotime <- 1:n_rows
  }
  
  #=================== target average heatmap =================
  if (target_average == T) {
    average_matrix <- data.frame(pseudotime = sub_matrix$pseudotime)
    
    for (j in 1:length(selected_gene)) {
      current_gene <- selected_gene[j]
      
      # target_list 구성
      if (ncol(gene_list) >= 3) {
        target_list <- subset(gene_list, gene_list[,1] == gsub("[.]", "-", current_gene))
        if (nrow(target_list) > 0) {
          target_list <- gsub("-", ".", target_list[,3])
        } else {
          target_list <- current_gene
        }
      } else {
        target_list <- current_gene
      }
      
      # 평균 계산
      if (length(target_list) == 1) {
        temp_values <- sub_matrix[[target_list]]
      } else {
        temp_values <- rowMeans(sub_matrix[, target_list, drop=FALSE])
      }
      
      average_matrix[[paste0(current_gene, " (", length(target_list), ")")]] <- temp_values
    }
    
    sub_matrix <- average_matrix
  }
  
  #====================z-score =================
  if(use_z_score == T) {
    for (i in 2:dim(sub_matrix)[2]) {
      sub_matrix[, i] <- (sub_matrix[, i] - mean(sub_matrix[, i])) / sd(sub_matrix[, i])
    }
  }
  
  #====================normalization =================
  normalized_data <- list()
  normalized_data[[1]] <- sub_matrix[, 1]  # pseudotime 값 저장
  
  # 각 유전자에 대해 정규화 수행
  for (i in 2:ncol(sub_matrix)) {
    x <- sub_matrix[, 1]  # pseudotime
    y <- sub_matrix[, i]  # expression values
    
    # NA, NaN, Inf 제거
    valid_idx <- which(is.finite(x) & is.finite(y))
    if (length(valid_idx) > 2) {  # loess를 위해 최소 3개 이상의 포인트 필요
      tryCatch({
        lo <- loess(y[valid_idx] ~ x[valid_idx], span = span)
        normalized_data[[i]] <- predict(lo, newdata = x)
      }, error = function(e) {
        # loess 실패 시 원본 데이터 사용
        normalized_data[[i]] <- y
      })
    } else {
      # 유효한 데이터가 부족한 경우 원본 데이터 사용
      normalized_data[[i]] <- y
    }
  }
  
  # 리스트를 행렬로 변환
  normalized_matrix <- do.call(cbind, normalized_data)
  colnames(normalized_matrix) <- colnames(sub_matrix)
  normalized_sub_matrix_sorted <- t(normalized_matrix)
  
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
    xaxis = list(
      title = "Pseudotime",
      titlefont = list(size = 14)
    ),
    yaxis = list(
      title = "Genes",
      titlefont = list(size = 14)
    ),
    colorbar = list(
      title = list(
        text = "Expression Level",
        font = list(size = 12)
      )
    )
  )
  
  fig <- fig %>% layout(layout)
  
  # JSON으로 직접 변환
  fig_json <- plotly_json(fig, jsonedit = FALSE, pretty = TRUE)
  
  # JSON 파일 저장
  write(fig_json, output_file_path)
}

count_degree <- function(data, decreasing=T, degree="out") {
    # 데이터 형식 확인 (열 개수)
    num_cols <- ncol(data)
    
    if (num_cols == 3) {  # 3열 형식 (source, weight, target)
        if (degree == "out") {
            att <- as.data.frame(table(data[,1]))[
                order(as.data.frame(table(data[,1]),)[,2], 
                decreasing = decreasing),]
        }
        if (degree == "in") {
            att <- as.data.frame(table(data[,3]))[
                order(as.data.frame(table(data[,3]),)[,2], 
                decreasing = decreasing),]
        }
    } else if (num_cols == 2) {  # 2열 형식 (gene, weight)
        # 2열 형식에서는 degree 파라미터 무시
        att <- data.frame(
            Var1 = data[,1],
            Freq = data[,2]
        )
        # 가중치(Freq)로 정렬
        att <- att[order(att$Freq, decreasing=decreasing),]
    } else {
        stop("입력 데이터는 2열 또는 3열 형식이어야 합니다.")
    }
    
    return(att)
}

# 입력 파일 경로
expression_data <- read.csv(input_expression_file)
trajectory <- read.table(input_trajectory_file, quote="\"", comment.char="")

# 고정값으로 사용되는 예시 데이터
count_TENET <- count_degree(tenet_sif, degree = "out")
selected_gene <- count_TENET$Var1

# 수치형 열만 선택
numeric_columns <- sapply(expression_data, is.numeric)
numeric_expression_data <- expression_data[, numeric_columns]

# 상위 {top_num}개 유전자 계산
top_genes <- head(order(rowSums(numeric_expression_data), decreasing = TRUE), top_num)

# 상위 {top_num}개의 유전자만 샘플링
sampled_expression_data <- numeric_expression_data[top_genes, ]

# 상위 {sampled_num}개의 샘플만 샘플링
sampled_trajectory <- trajectory[1:sampled_num, ]

# selected_gene 중 expression_data의 열 이름에 있는 유전자만 선택
valid_genes <- selected_gene[selected_gene %in% colnames(sampled_expression_data)]

# 함수 호출
pseudotime_heatmap(matrix = sampled_expression_data, gene_list = tenet_sif,
                   selected_gene = as.character(valid_genes), pseudotime = sampled_trajectory,
                   span=0.7, use_pseudotime_origin = T, use_z_score = T,
                   order_pseudo_score=T, max_min = T, p_min=-0.7, p_max=0.7, p_legend=T, 
                   out_result =F, target_average = T, pseudo_decrease=T)
