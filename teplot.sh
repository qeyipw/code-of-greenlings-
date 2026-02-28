#!/bin/bash
# te_plot_batch.sh
# 批量运行转座子密度绘图脚本

echo "开始批量绘制转座子密度图..."
echo "========================================"

# 创建日志目录
log_dir="te_plot_logs"
mkdir -p "$log_dir"

# 记录开始时间
start_time=$(date)
echo "开始时间: $start_time"
echo ""

# 定义任务参数数组 - 特别注意Y任务的高亮区域字符串
declare -A tasks=(
    ["btc"]="/mnt/DATA/SDB/home/SJW/btc/btc_chr14_plot.gff:chr14:4200000-6800000;chr14:6800000-15150000;chr14:15150000-19900000"
    ["btx"]="/mnt/DATA/SDB/home/SJW/BTM/chr05_PLOT.gff:chr05:3914066-6500000;chr05:6500000-14500000;chr05:14500000-19339089"
    ["dlc"]="/mnt/DATA/SDB/home/SJW/dlc_new/plot.gff:chr14:4800000-6960000;chr14:6960000-13750000;chr14:13750000-14600000"
    ["dlx"]="/mnt/DATA/SDB/home/SJW/dlm/dlm_plot.gff:chr14:5133725-7430000;chr14:7430000-14490000;chr14:14490000-15303667"
    ["Y"]="/mnt/DATA/SDB/home/SJW/bt_dan/h1/Y_plot.gff:chr14:4500000-7100000;chr14:7100000-14700000;chr14:14700000-19600000"  # 修正引号问题
)

# 运行所有任务
for output in "${!tasks[@]}"; do
    # 解析参数
    IFS=':' read -r input_file highlight_regions <<< "${tasks[$output]}"
    
    echo "========================================"
    echo "运行任务: $output"
    echo "输入文件: $input_file"
    echo "高亮区域: $highlight_regions"
    echo "----------------------------------------"
    
    # 创建单独的日志文件
    log_file="$log_dir/${output}_plot.log"
    
    # 运行命令并记录输出到日志文件
    {
        echo "========================================"
        echo "任务: $output"
        echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "输入文件: $input_file"
        echo "高亮区域: $highlight_regions"
        echo "----------------------------------------"
        
        # 使用正确的Python脚本路径
        if [ "$output" = "Y" ]; then
            python_script="/mnt/DATA/SDB/home/SJW/te_density_wux_plot.py"
        else
            python_script="/mnt/DATA/SDB/home/SJW/te_density_plot.py"
        fi
        
        echo "使用Python脚本: $python_script"
        
        # 直接传递高亮区域参数，不使用引号包裹
        python "$python_script" \
            -i "$input_file" \
            -w 50000 \
            -o "$output" \
            --highlight "$highlight_regions"
        
        echo "----------------------------------------"
        echo "任务: $output"
        echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "========================================"
    } 2>&1 | tee "$log_file"
    
    # 检查命令执行状态
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "✅ 任务 $output 完成"
    else
        echo "❌ 任务 $output 失败，请查看日志: $log_file"
    fi
    
    echo ""
    sleep 1
done

# 记录结束时间
end_time=$(date)
echo "========================================"
echo "所有任务完成!"
echo "开始时间: $start_time"
echo "结束时间: $end_time"
echo ""
echo "生成的文件:"
echo "----------------------------------------"

# 列出生成的文件
for output in "${!tasks[@]}"; do
    echo "对于 $output:"
    ls -la ${output}_combined.png 2>/dev/null || echo "  未找到 combined.png 文件"
    ls -la ${output}_legend.png 2>/dev/null || echo "  未找到 legend.png 文件"
    ls -la ${output}_data.tsv 2>/dev/null || echo "  未找到 data.tsv 文件"
    ls -la ${output}_total_lengths.tsv 2>/dev/null || echo "  未找到 total_lengths.tsv 文件"
    echo ""
done

echo "日志文件保存在: $log_dir/"
ls -la "$log_dir/"
echo ""
echo "批量任务完成!"