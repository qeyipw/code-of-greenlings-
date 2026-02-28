#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# 设置参数
autosomal_file <- ifelse(length(args) > 0, args[1], "autosomal_fpkm.txt")
affected_file <- ifelse(length(args) > 1, args[2], "affected_fpkm.txt")

# 加载必要的包
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos="https://cloud.r-project.org")
  library(dplyr)
}

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos="https://cloud.r-project.org")
  library(ggplot2)
}

if (!require("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr", repos="https://cloud.r-project.org")
  library(ggpubr)
}

# 读取2种类型的FPKM矩阵
cat("Reading FPKM matrices for 2 gene types...\n")

read_and_label <- function(file_path, label) {
  if (!file.exists(file_path)) {
    cat("Warning: File", file_path, "does not exist. Skipping.\n")
    return(NULL)
  }
  data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  
  data$GeneType <- label
  data$GeneID <- rownames(data)
  return(data)
}

autosomal_data <- read_and_label(autosomal_file, "auto")
affected_data <- read_and_label(affected_file, "sex")

# 检查数据是否成功读取
if (is.null(autosomal_data) || is.null(affected_data)) {
  stop("One or more input files could not be read. Please check file paths and names.")
}

# 找出两个文件共有的样本
common_samples <- intersect(colnames(autosomal_data), colnames(affected_data))
common_samples <- common_samples[!common_samples %in% c("GeneType", "GeneID")]

cat("Common samples found:", length(common_samples), "\n")
cat("Common sample names:", paste(common_samples, collapse = ", "), "\n")

# 只保留共有的样本进行分析
autosomal_data_common <- autosomal_data[, c("GeneType", "GeneID", common_samples)]
affected_data_common <- affected_data[, c("GeneType", "GeneID", common_samples)]

# 合并数据
all_data <- rbind(autosomal_data_common, affected_data_common)
cat("Combined data dimensions:", dim(all_data), "\n")

# 解析样本信息 - 适配X_Liver_1, Y_Muscle_2格式
parse_sample_info <- function(sample_names) {
  sex <- sapply(strsplit(sample_names, "_"), function(x) x[1])
  tissue <- sapply(strsplit(sample_names, "_"), function(x) x[2])
  individual <- sapply(strsplit(sample_names, "_"), function(x) x[3])
  
  sex_full <- ifelse(sex == "X", "Female", ifelse(sex == "Y", "Male", sex))
  
  return(data.frame(
    Sample = sample_names,
    Sex = sex_full,
    Tissue = tissue,
    Individual = individual,
    stringsAsFactors = FALSE
  ))
}

sample_info <- parse_sample_info(common_samples)
cat("Sample information:\n")
print(table(sample_info$Sex, sample_info$Tissue))

# 定义组织类型
tissues <- unique(sample_info$Tissue)

# 存储所有组织的分析结果
all_results <- data.frame()

# 对每个组织进行2种基因类型的雌雄差异分析
for (tissue in tissues) {
  cat("\nAnalyzing sex differences in", tissue, "for 2 gene types...\n")
  
  # 提取该组织的样本
  tissue_samples <- sample_info[sample_info$Tissue == tissue, ]
  female_samples <- tissue_samples[tissue_samples$Sex == "Female", "Sample"]
  male_samples <- tissue_samples[tissue_samples$Sex == "Male", "Sample"]
  
  # 检查是否有足够的样本，如果没有也继续分析
  if (length(female_samples) < 1 || length(male_samples) < 1) {
    cat("Warning: Not enough samples for", tissue, "analysis. Female:", length(female_samples), "Male:", length(male_samples), "\n")
    # 即使样本不足也继续，使用可用样本
    if (length(female_samples) < 1 && length(male_samples) < 1) {
      cat("No samples available for", tissue, ". Skipping.\n")
      next
    }
  }
  
  cat("  Female samples:", length(female_samples), "-", paste(female_samples, collapse = ", "), "\n")
  cat("  Male samples:", length(male_samples), "-", paste(male_samples, collapse = ", "), "\n")
  
  # 提取该组织的表达数据
  tissue_data <- all_data[, c("GeneType", "GeneID", female_samples, male_samples)]
  
  # 对每种基因类型进行分析
  for (gene_type in c("auto", "sex")) {
    cat("  Analyzing", gene_type, "genes...\n")
    
    # 提取该基因类型的数据
    type_data <- tissue_data[tissue_data$GeneType == gene_type, ]
    
    if (nrow(type_data) == 0) {
      cat("    No data for", gene_type, "genes in", tissue, ". Creating empty record.\n")
      # 创建空记录而不是跳过
      empty_data <- data.frame(
        GeneID = NA,
        Female_Mean = NA,
        Male_Mean = NA,
        Ratio = NA,
        Tissue = tissue,
        GeneType = gene_type,
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, empty_data)
      next
    }
    
    # 提取表达矩阵
    expr_matrix <- type_data[, c(female_samples, male_samples)]
    rownames(expr_matrix) <- type_data$GeneID
    
    # 计算雌雄平均表达量
    female_mean <- apply(expr_matrix[, female_samples, drop = FALSE], 1, mean)
    male_mean <- apply(expr_matrix[, male_samples, drop = FALSE], 1, mean)
    
    # 过滤条件1：雌雄平均FPKM都小于1的基因
    keep_genes <- (female_mean >= 1) | (male_mean >= 1)
    
    # 应用过滤1
    female_mean_filtered <- female_mean[keep_genes]
    male_mean_filtered <- male_mean[keep_genes]
    expr_matrix_filtered <- expr_matrix[keep_genes, , drop = FALSE]
    
    cat("    Genes retained after FPKM filtering:", sum(keep_genes), "/", length(keep_genes), "\n")
    
    # 如果没有基因保留，创建空记录
    if (sum(keep_genes) == 0) {
      cat("    All genes filtered out for", gene_type, "in", tissue, ". Creating empty record.\n")
      empty_data <- data.frame(
        GeneID = NA,
        Female_Mean = NA,
        Male_Mean = NA,
        Ratio = NA,
        Tissue = tissue,
        GeneType = gene_type,
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, empty_data)
      next
    }
    
    # 计算雌雄表达量比值 (Female/Male)，不进行log转换
    ratio <- female_mean_filtered / (male_mean_filtered + 1e-10)
    
    # 过滤条件2：雌雄比值大于2.7的基因（避免极端值）
    keep_ratio <- ratio >= 0.3 & ratio <= 2.5
    ratio_filtered <- ratio[keep_ratio]
    female_mean_final <- female_mean_filtered[keep_ratio]
    male_mean_final <- male_mean_filtered[keep_ratio]
    gene_ids_final <- rownames(expr_matrix_filtered)[keep_ratio]
    
    cat("    Genes retained after ratio filtering:", sum(keep_ratio), "/", length(keep_ratio), "\n")
    
    # 如果没有基因保留，创建空记录
    if (sum(keep_ratio) == 0) {
      cat("    All genes filtered out due to ratio > 2.5 for", gene_type, "in", tissue, ". Creating empty record.\n")
      empty_data <- data.frame(
        GeneID = NA,
        Female_Mean = NA,
        Male_Mean = NA,
        Ratio = NA,
        Tissue = tissue,
        GeneType = gene_type,
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, empty_data)
      next
    }
    
    # 创建该类型的详细比值数据
    type_ratio_data <- data.frame(
      GeneID = gene_ids_final,
      Female_Mean = female_mean_final,
      Male_Mean = male_mean_final,
      Ratio = ratio_filtered,
      Tissue = tissue,
      GeneType = gene_type,
      stringsAsFactors = FALSE
    )
    
    # 添加到总结果
    all_results <- rbind(all_results, type_ratio_data)
  }
}

# 输出总结报告和可视化
if (nrow(all_results) > 0) {
  # 确保因子顺序正确
  all_results$Tissue <- factor(all_results$Tissue, levels = tissues)
  all_results$GeneType <- factor(all_results$GeneType, 
                                 levels = c("auto", "sex"))
  
  # 创建箱型图 - 使用原始比值，不进行log转换
  # 移除NA值用于绘图
  plot_data <- all_results[!is.na(all_results$Ratio), ]
  
  if (nrow(plot_data) > 0) {
    boxplot <- ggplot(plot_data, aes(x = GeneType, y = Ratio, fill = GeneType)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      facet_wrap(~ Tissue, nrow = 1) +
      # 限制y轴范围，避免极端值
      coord_cartesian(ylim = c(0, 3)) +
      # 修改x轴标签
      scale_x_discrete(labels = c("auto" = "ref/alt\n(autosome)", 
                                 "sex" = "sex_linked\ngenes")) +
      labs(
        x = "",
        y = "X/Y expression"
      ) +
      theme_classic() +
      theme(
        plot.title = element_blank(),
        text = element_text(family = "Arial", size = 18, color = "black"),
        axis.text = element_text(family = "Arial", size = 18, color = "black"),
        axis.text.x = element_text(family = "Arial", size = 14),
        axis.title = element_text(family = "Arial", size = 18, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),  # Y轴标题字体变小
        legend.position = "none",
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Arial", size = 20, color = "black")
      ) +
      scale_fill_manual(values = c("auto" = "lightblue", 
                                  "sex" = "lightcoral"))
    
    # 保存箱型图
    boxplot_output <- "two_gene_types_boxplot.png"
    ggsave(boxplot_output, boxplot, width = 6, height = 3, dpi = 900, bg = "white")
    cat("Boxplot saved to:", boxplot_output, "\n")
  } else {
    cat("No valid data for plotting. Skipping boxplot generation.\n")
  }
  
  # 计算并保存每个组织的雌雄表达量总结
  expression_summary <- all_results %>%
    group_by(Tissue, GeneType) %>%
    summarise(
      Female_Mean_Expression = mean(Female_Mean, na.rm = TRUE),
      Male_Mean_Expression = mean(Male_Mean, na.rm = TRUE),
      Median_Ratio = median(Ratio, na.rm = TRUE),
      Mean_Ratio = mean(Ratio, na.rm = TRUE),
      Genes_Count = sum(!is.na(Ratio)),  # 计算非NA的基因数量
      Total_Genes_Attempted = n(),
      .groups = 'drop'
    )
  
  # 保存表达量总结
  expression_output <- "two_gene_types_expression_summary.txt"
  write.table(expression_summary, expression_output, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Expression summary saved to:", expression_output, "\n")
  
  # 保存所有组织的详细比值数据
  ratio_output <- "two_gene_types_ratio_data.txt"
  write.table(all_results, ratio_output, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("All ratio data saved to:", ratio_output, "\n")
  
  # 打印比值统计
  cat("\nRatio statistics by tissue and gene type:\n")
  ratio_stats <- all_results %>%
    group_by(Tissue, GeneType) %>%
    summarise(
      Median_Ratio = median(Ratio, na.rm = TRUE),
      Mean_Ratio = mean(Ratio, na.rm = TRUE),
      Genes_Count = sum(!is.na(Ratio)),
      Total_Genes_Attempted = n(),
      .groups = 'drop'
    )
  print(ratio_stats)
} else {
  cat("No data available after filtering. Please check your input files and filtering criteria.\n")
}

cat("\nThree gene types sex difference analysis completed successfully!\n")