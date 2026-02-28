#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# 设置参数
autosomal_file <- ifelse(length(args) > 0, args[1], "autosomal_fpkm.txt")
affected_file <- ifelse(length(args) > 1, args[2], "affected_fpkm.txt")
unaffected_file <- ifelse(length(args) > 2, args[3], "unaffected_fpkm.txt")

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

# 读取三种类型的FPKM矩阵
cat("Reading FPKM matrices for three gene types...\n")

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

autosomal_data <- read_and_label(autosomal_file, "aut")
affected_data <- read_and_label(affected_file, "aff")
unaffected_data <- read_and_label(unaffected_file, "rem")

# 检查数据是否成功读取
if (is.null(autosomal_data) || is.null(affected_data) || is.null(unaffected_data)) {
  stop("One or more input files could not be read. Please check file paths and names.")
}

# 合并数据
all_data <- rbind(autosomal_data, affected_data, unaffected_data)
cat("Combined data dimensions:", dim(all_data), "\n")

# 提取样本信息
samples <- setdiff(colnames(all_data), c("GeneType", "GeneID"))

# 解析样本信息 - 适配F_Liver_1, M_Muscle_2格式
parse_sample_info <- function(sample_names) {
  sex <- sapply(strsplit(sample_names, "_"), function(x) x[1])
  tissue <- sapply(strsplit(sample_names, "_"), function(x) x[2])
  individual <- sapply(strsplit(sample_names, "_"), function(x) x[3])
  
  sex_full <- ifelse(sex == "F", "Female", ifelse(sex == "M", "Male", sex))
  
  return(data.frame(
    Sample = sample_names,
    Sex = sex_full,
    Tissue = tissue,
    Individual = individual,
    stringsAsFactors = FALSE
  ))
}

sample_info <- parse_sample_info(samples)
cat("Sample information:\n")
print(table(sample_info$Sex, sample_info$Tissue))

# 定义组织类型
tissues <- unique(sample_info$Tissue)

# 存储所有组织的分析结果
all_results <- data.frame()

# 对每个组织进行三种基因类型的雌雄差异分析
for (tissue in tissues) {
  cat("\nAnalyzing sex differences in", tissue, "for three gene types...\n")
  
  # 提取该组织的样本
  tissue_samples <- sample_info[sample_info$Tissue == tissue, ]
  female_samples <- tissue_samples[tissue_samples$Sex == "Female", "Sample"]
  male_samples <- tissue_samples[tissue_samples$Sex == "Male", "Sample"]
  
  # 检查是否有足够的样本
  if (length(female_samples) < 1 || length(male_samples) < 1) {
    cat("Not enough samples for", tissue, "analysis. Skipping.\n")
    next
  }
  
  # 提取该组织的表达数据
  tissue_data <- all_data[, c("GeneType", "GeneID", female_samples, male_samples)]
  
  # 对每种基因类型进行分析
  for (gene_type in c("aut", "aff", "rem")) {
    cat("  Analyzing", gene_type, "genes...\n")
    
    # 提取该基因类型的数据
    type_data <- tissue_data[tissue_data$GeneType == gene_type, ]
    
    if (nrow(type_data) == 0) {
      cat("    No data for", gene_type, "genes in", tissue, ". Skipping.\n")
      next
    }
    
    # 提取表达矩阵
    expr_matrix <- type_data[, c(female_samples, male_samples)]
    rownames(expr_matrix) <- type_data$GeneID
    
    # 计算雌雄平均表达量
    female_mean <- apply(expr_matrix[, female_samples], 1, mean)
    male_mean <- apply(expr_matrix[, male_samples], 1, mean)
    
    # 过滤条件1：雌雄平均FPKM都小于0.1的基因
    keep_genes <- (female_mean >= 0.1) | (male_mean >= 0.1)
    
    # 应用过滤1
    female_mean <- female_mean[keep_genes]
    male_mean <- male_mean[keep_genes]
    expr_matrix <- expr_matrix[keep_genes, ]
    
    cat("    Genes retained after FPKM filtering:", sum(keep_genes), "/", length(keep_genes), "\n")
    
    # 如果没有基因保留，跳过
    if (sum(keep_genes) == 0) {
      cat("    All genes filtered out for", gene_type, "in", tissue, ". Skipping.\n")
      next
    }
    
    # 计算雌雄表达量比值 (Female/Male)，不进行log转换
    ratio <- female_mean / (male_mean + 1e-10)
    
    # 过滤条件2：雌雄比值大于0.2小于2.5的基因
    keep_ratio <- ratio >= 0.2 & ratio <= 2.3
    ratio <- ratio[keep_ratio]
    female_mean <- female_mean[keep_ratio]
    male_mean <- male_mean[keep_ratio]
    
    cat("    Genes retained after ratio filtering:", sum(keep_ratio), "/", length(keep_ratio), "\n")
    
    # 如果没有基因保留，跳过
    if (sum(keep_ratio) == 0) {
      cat("    All genes filtered out due to ratio > 3 for", gene_type, "in", tissue, ". Skipping.\n")
      next
    }
    
    # 创建该类型的详细比值数据
    type_ratio_data <- data.frame(
      GeneID = rownames(expr_matrix)[keep_ratio],
      Female_Mean = female_mean,
      Male_Mean = male_mean,
      Ratio = ratio,
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
                                 levels = c("aut", "aff", "rem"))
  
  # 创建箱型图 - 使用原始比值，不进行log转换
  boxplot <- ggplot(all_results, aes(x = GeneType, y = Ratio, fill = GeneType)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~ Tissue, nrow = 1) +
    # 限制y轴范围，避免极端值
    coord_cartesian(ylim = c(0, 3)) +
    # 简化主题和标签
    labs(
      x = "",
      y = "Female/Male expression"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial"),  # 设置所有文本为Arial
      plot.title = element_blank(),
      axis.text = element_text(size = 18, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),  # 进一步增大x轴文字
      axis.title = element_text(size = 18, color = "black"),  # 进一步增大坐标轴标题
      axis.title.y = element_text(size = 14, color = "black"),  # Y轴标题字体变小
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 20, color = "black")  # 进一步增大分面标签
    ) +
    scale_fill_manual(values = c("aut" = "lightblue", 
                                "aff" = "lightcoral", 
                                "rem" = "lightgreen"))
  
  # 使用stat_compare_means添加p值
  comparisons_list <- list(c("aut", "aff"), c("aut", "rem"), c("aff", "rem"))
  
  boxplot <- boxplot +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.format",
      size = 4,  # 进一步增大p值标签
      vjust = 0.05,
      step.increase = 0.12,
      tip.length = 0,
      family = "Arial"  # 设置p值标签字体为Arial
    )
  
  # 保存箱型图 - 减小图片尺寸但保持高分辨率
  boxplot_output <- "three_gene_types_boxplot.png"
  ggsave(boxplot_output, boxplot, width = 6, height = 3, dpi = 900, bg = "white")  # 减小图片尺寸
  cat("Boxplot with significance markers saved to:", boxplot_output, "\n")
  
  # 计算并保存每个组织的雌雄表达量总结
  expression_summary <- all_results %>%
    group_by(Tissue, GeneType) %>%
    summarise(
      Female_Mean_Expression = mean(Female_Mean, na.rm = TRUE),
      Male_Mean_Expression = mean(Male_Mean, na.rm = TRUE),
      Median_Ratio = median(Ratio, na.rm = TRUE),
      Mean_Ratio = mean(Ratio, na.rm = TRUE),
      Genes_Count = n(),
      .groups = 'drop'
    )
  
  # 保存表达量总结
  expression_output <- "three_gene_types_expression_summary.txt"
  write.table(expression_summary, expression_output, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Expression summary saved to:", expression_output, "\n")
  
  # 保存每个组织的详细雌雄表达量
  for (tissue in tissues) {
    tissue_data <- all_results[all_results$Tissue == tissue, ]
    
    # 保存该组织的详细表达量数据
    tissue_output <- paste0("expression_", tissue, ".txt")
    write.table(tissue_data, tissue_output, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Detailed expression data for", tissue, "saved to:", tissue_output, "\n")
  }
  
  # 保存所有组织的详细比值数据
  ratio_output <- "three_gene_types_ratio_data.txt"
  write.table(all_results, ratio_output, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("All ratio data saved to:", ratio_output, "\n")
  
  # 打印比值统计
  cat("\nRatio statistics by tissue and gene type:\n")
  ratio_stats <- all_results %>%
    group_by(Tissue, GeneType) %>%
    summarise(
      Median_Ratio = median(Ratio, na.rm = TRUE),
      Mean_Ratio = mean(Ratio, na.rm = TRUE),
      Genes_Count = n(),
      .groups = 'drop'
    )
  print(ratio_stats)
} else {
  cat("No data available after filtering. Please check your input files and filtering criteria.\n")
}

cat("\nThree gene types sex difference analysis completed successfully!\n")