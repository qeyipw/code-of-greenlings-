import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse
import os
import re
import matplotlib as mpl
from matplotlib.patches import Rectangle

# 设置系统可用的字体
# 首选 Arial，如果不存在则使用其他无衬线字体
available_fonts = ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans', 
                   'sans-serif', 'STIXGeneral']

# 检查系统字体
from matplotlib import font_manager
import matplotlib

# 获取所有可用字体
font_list = [f.name for f in font_manager.fontManager.ttflist]

# 选择可用的字体
selected_font = 'Arial'  # 默认首选
for font in available_fonts:
    if font in font_list:
        selected_font = font
        break

print(f"使用字体: {selected_font}")

# 设置matplotlib字体
plt.rcParams['font.family'] = selected_font
plt.rcParams['mathtext.fontset'] = 'stix'  # 数学字体
plt.rcParams['axes.unicode_minus'] = False

# 如果系统中没有 Arial，打印提示
if 'Arial' not in font_list:
    print("警告: 系统中未找到 Arial 字体，使用替代字体")
    print("可用的类似字体:", [f for f in font_list if 'Arial' in f or 'Sans' in f])

class TransposonDensityPlotter:
    def __init__(self, gff_file, window_size=100000, output_prefix="transposon_density", 
                 chromosomes=None, regions=None, te_types=None, style="default",
                 highlight_regions=None, max_x_position=35):
        self.gff_file = gff_file
        self.window_size = window_size
        self.output_prefix = output_prefix
        self.chromosome_data = defaultdict(lambda: defaultdict(list))
        self.target_chromosomes = chromosomes
        self.target_regions = regions
        self.te_types = te_types
        self.style = style
        self.highlight_regions = highlight_regions  # 新增：高亮区域
        self.max_x_position = max_x_position  # 新增：最大X轴位置
        
        # 新增：存储每个转座子类型的总长度
        self.te_total_lengths = defaultdict(int)
        
        # 设置绘图风格
        self.set_plot_style(style)
        
        # 存储选中的字体
        self.selected_font = selected_font
        
        # 默认高亮区域颜色和透明度 - 三个不同颜色的阴影
        self.highlight_colors = ['#FFE4E1', '#E0FFFF', '#FFDAB9']  # 浅粉色、浅青色、浅橙色
        self.highlight_alpha = 0.7
    
    def parse_highlight_regions(self, highlight_str):
        """解析高亮区域字符串，支持多个不同颜色的区域"""
        if not highlight_str:
            return []
        
        regions = []
        try:
            # 格式: "chr:start-end;chr:start-end;chr:start-end"
            parts = highlight_str.split(';')
            for i, part in enumerate(parts):
                if not part.strip():
                    continue
                    
                chrom, pos_range = part.strip().split(':')
                start, end = pos_range.split('-')
                
                # 转换为Mb单位
                start_mb = int(start) / 1e6
                end_mb = int(end) / 1e6
                
                # 分配颜色（循环使用颜色列表）
                color_idx = i % len(self.highlight_colors)
                
                regions.append({
                    'chromosome': chrom,
                    'start_mb': start_mb,
                    'end_mb': end_mb,
                    'start': int(start),
                    'end': int(end),
                    'color': self.highlight_colors[color_idx]  # 添加颜色属性
                })
                
                print(f"高亮区域 {i+1}: {chrom}:{start}-{end} ({start_mb:.2f}Mb - {end_mb:.2f}Mb), 颜色: {self.highlight_colors[color_idx]}")
        except Exception as e:
            print(f"解析高亮区域时出错: {e}")
            print("高亮区域格式应为: chr:start-end;chr2:start-end;chr3:start-end")
        
        return regions
    
    def set_plot_style(self, style):
        """设置绘图风格"""
        if style == "dark":
            plt.style.use('dark_background')
            self.bg_color = '#1a1a1a'
            self.grid_color = '#333333'
            self.text_color = 'white'
        elif style == "seaborn":
            plt.style.use('seaborn-v0_8-whitegrid')
            self.bg_color = 'white'
            self.grid_color = '#dddddd'
            self.text_color = 'black'
        elif style == "classic":
            plt.style.use('classic')
            self.bg_color = 'white'
            self.grid_color = '#cccccc'
            self.text_color = 'black'
        else:  # default
            plt.style.use('default')
            self.bg_color = 'white'
            self.grid_color = '#e0e0e0'
            self.text_color = 'black'
    
    def parse_region_string(self, region_str):
        """解析区域字符串，格式: chr:start-end 或 chr"""
        if ':' in region_str:
            chrom, pos = region_str.split(':')
            if '-' in pos:
                start, end = map(int, pos.split('-'))
                return chrom, start, end
            else:
                return chrom, None, None
        else:
            return region_str, None, None
    
    def standardize_te_type(self, te_type):
        """标准化转座子类型名称"""
        te_type = str(te_type).strip()
        
        # 将常见的转座子类型名称标准化
        type_mapping = {
            'hAT': 'hAT',
            'cacta': 'CACTA',
            'mutator': 'Mutator',
            'pif_harbinger': 'PIF_Harbinger',
            'pif': 'PIF_Harbinger',
            'harbinger': 'PIF_Harbinger',
            'tc1_mariner': 'Tc1_Mariner',
            'tc1': 'Tc1_Mariner',
            'mariner': 'Tc1_Mariner',
            'helitron': 'Helitron',
            'line': 'LINE',
            'ltr': 'LTR',
            'ltr_retrotransposon': 'LTR',
            'unknown': 'Unknown',
            'simple_repeat': 'Simple_repeat',
            'simple': 'Simple_repeat'
        }
        
        # 尝试匹配已知类型
        for key, value in type_mapping.items():
            if key in te_type.lower():
                return value
        
        # 如果无法匹配，返回原始类型（首字母大写）
        return te_type.capitalize()
    
    def merge_overlapping_intervals(self, intervals):
        """合并重叠的区间"""
        if not intervals:
            return []
        
        # 按起始位置排序
        intervals.sort(key=lambda x: x[0])
        
        merged = []
        current_start, current_end = intervals[0]
        
        for start, end in intervals[1:]:
            if start <= current_end:  # 有重叠
                current_end = max(current_end, end)  # 合并区间
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        
        merged.append((current_start, current_end))
        return merged
    
    def parse_gff(self):
        """解析GFF文件，直接使用第三列的转座子类型，并合并重叠的转座子"""
        print(f"正在解析GFF文件: {self.gff_file}")
        
        # 解析目标区域
        parsed_regions = []
        if self.target_regions:
            for region in self.target_regions:
                chrom, start, end = self.parse_region_string(region)
                parsed_regions.append((chrom, start, end))
                print(f"目标区域: {chrom}:{start if start else 'start'}-{end if end else 'end'}")
        
        # 解析高亮区域
        highlight_regions = []
        if self.highlight_regions:
            highlight_regions = self.parse_highlight_regions(self.highlight_regions)
        
        # 定义我们感兴趣的TE类型
        target_te_types = ['hAT', 'CACTA', 'Mutator', 'PIF_Harbinger', 
                          'Tc1_Mariner', 'Helitron', 'LINE', 'LTR', 
                          'Unknown', 'Simple_repeat']
        
        # 如果用户没有指定TE类型，使用我们定义的列表
        if self.te_types is None:
            self.te_types = target_te_types
            print(f"使用预定义的TE类型: {', '.join(self.te_types)}")
        
        try:
            # 读取GFF文件
            df = pd.read_csv(self.gff_file, sep='\t', comment='#', header=None,
                           names=['seqid', 'source', 'type', 'start', 'end', 
                                  'score', 'strand', 'phase', 'attributes'])
            
            # 不再从属性中提取类型，直接使用第三列的类型
            print(f"找到 {len(df)} 个特征")
            
            # 按染色体和类型分组存储位置信息和长度，应用过滤
            type_counts = defaultdict(int)
            total_length = 0
            
            # 临时存储所有转座子，用于合并重叠
            all_intervals = defaultdict(lambda: defaultdict(list))
            
            for _, row in df.iterrows():
                chrom = row['seqid']
                start = int(row['start'])
                end = int(row['end'])
                te_type = self.standardize_te_type(row['type'])
                
                # 应用TE类型过滤 - 只保留我们感兴趣的TE类型
                if te_type not in self.te_types:
                    continue
                
                # 应用染色体过滤
                if self.target_chromosomes and chrom not in self.target_chromosomes:
                    continue
                
                # 应用区域过滤
                keep = True
                if parsed_regions:
                    keep = False
                    for target_chrom, target_start, target_end in parsed_regions:
                        if chrom == target_chrom:
                            if target_start is None and target_end is None:
                                keep = True  # 整个染色体
                                break
                            elif target_start <= start <= target_end or target_start <= end <= target_end:
                                keep = True  # 与目标区域有重叠
                                break
                
                if keep:
                    # 存储转座子的起始和结束位置
                    all_intervals[chrom][te_type].append((start, end))
                    type_counts[te_type] += 1
            
            # 打印找到的所有TE类型
            print("\n找到的TE类型:")
            for te_type, count in type_counts.items():
                print(f"  {te_type}: {count}")
            
            # 合并重叠的转座子
            print("\n正在合并重叠的转座子...")
            for chrom in all_intervals:
                for te_type in all_intervals[chrom]:
                    intervals = all_intervals[chrom][te_type]
                    merged_intervals = self.merge_overlapping_intervals(intervals)
                    self.chromosome_data[chrom][te_type] = merged_intervals
                    
                    # 计算合并后的总长度
                    for start, end in merged_intervals:
                        length = end - start + 1
                        self.te_total_lengths[te_type] += length
                        total_length += length
            
            # 打印合并前后的统计信息
            print("\n转座子类型统计 (合并重叠前):")
            for te_type in self.te_types:
                if te_type in type_counts:
                    print(f"  {te_type}: {type_counts[te_type]}")
            
            print("\n转座子类型统计 (合并重叠后):")
            for te_type in self.te_types:
                count = sum(len(self.chromosome_data[chrom][te_type]) for chrom in self.chromosome_data)
                if count > 0:
                    print(f"  {te_type}: {count} (总长度: {self.te_total_lengths[te_type]} bp)")
            
            print(f"\n合并后总转座子长度: {total_length} bp")
            
            # 检查哪些染色体有数据
            print("\n有转座子数据的染色体:")
            for chrom in self.chromosome_data:
                te_types_in_chrom = [te_type for te_type in self.chromosome_data[chrom] if self.chromosome_data[chrom][te_type]]
                if te_types_in_chrom:
                    print(f"  {chrom}: {', '.join(te_types_in_chrom)}")
                
        except Exception as e:
            print(f"解析GFF文件时出错: {e}")
            import traceback
            traceback.print_exc()
            return False
            
        return True
    
    def calculate_density(self):
        """计算每个窗口的转座子长度密度，按类型分别计算"""
        print("正在计算转座子长度密度...")
        
        density_data = {}
        parsed_regions = []
        
        # 解析目标区域
        if self.target_regions:
            for region in self.target_regions:
                chrom, start, end = self.parse_region_string(region)
                parsed_regions.append((chrom, start, end))
        
        for chrom, type_positions in self.chromosome_data.items():
            if not any(type_positions.values()):
                print(f"染色体 {chrom} 没有转座子数据，跳过")
                continue
                
            # 查找该染色体的目标区域
            chrom_regions = [(start, end) for c, start, end in parsed_regions 
                           if c == chrom] if parsed_regions else []
            
            if chrom_regions:
                # 有特定区域，为每个区域单独计算
                for region_idx, (region_start, region_end) in enumerate(chrom_regions):
                    if region_start is None or region_end is None:
                        # 整个染色体
                        region_start = 1
                        region_end = max([max(end for start, end in positions) for positions in type_positions.values() if positions])
                        region_name = chrom
                    else:
                        # 特定区域
                        region_name = f"{chrom}_{region_start//1000000}Mb_{region_end//1000000}Mb"
                    
                    # 为每个TE类型计算长度密度
                    te_densities = {}
                    for te_type, positions in type_positions.items():
                        if not positions:
                            continue
                            
                        num_windows = (region_end - region_start) // self.window_size + 1
                        density = np.zeros(num_windows)
                        
                        for start_pos, end_pos in positions:
                            # 计算转座子与目标区域的交集
                            overlap_start = max(start_pos, region_start)
                            overlap_end = min(end_pos, region_end)
                            
                            if overlap_start > overlap_end:
                                continue  # 没有重叠
                            
                            # 计算转座子覆盖的窗口范围
                            first_window = (overlap_start - region_start) // self.window_size
                            last_window = (overlap_end - region_start) // self.window_size
                            
                            # 对于每个被覆盖的窗口，计算转座子在该窗口中的长度
                            for window_idx in range(first_window, last_window + 1):
                                if window_idx < 0 or window_idx >= num_windows:
                                    continue
                                
                                window_start = region_start + window_idx * self.window_size
                                window_end = window_start + self.window_size - 1
                                
                                # 计算转座子与当前窗口的重叠长度
                                overlap_start_window = max(overlap_start, window_start)
                                overlap_end_window = min(overlap_end, window_end)
                                
                                if overlap_start_window <= overlap_end_window:
                                    overlap_length = overlap_end_window - overlap_start_window + 1
                                    density[window_idx] += overlap_length
                        
                        # 转换为每100kb的转座子长度 (bp/100kb)
                        density_per_100kb = density / (self.window_size / 100000)
                        
                        # 计算位置坐标
                        window_starts = np.arange(num_windows) * self.window_size + region_start
                        
                        te_densities[te_type] = {
                            'positions': window_starts / 1e6,  # 转换为Mb
                            'density': density_per_100kb
                        }
                    
                    if te_densities:
                        density_data[region_name] = te_densities
                        total_te_length = sum(end - start + 1 for positions in type_positions.values() for start, end in positions)
                        print(f"区域 {region_name}: 有 {len(te_densities)} 种TE类型")
                        for te_type in te_densities:
                            print(f"    {te_type}: {len(type_positions[te_type])} 个区域")
            else:
                # 整个染色体
                # 找出染色体的最大位置
                all_positions = []
                for positions in type_positions.values():
                    if positions:
                        all_positions.extend(positions)
                
                if not all_positions:
                    print(f"染色体 {chrom} 没有转座子位置数据")
                    continue
                
                chrom_start = 1
                chrom_end = max([end for start, end in all_positions])
                
                # 为每个TE类型计算长度密度
                te_densities = {}
                for te_type, positions in type_positions.items():
                    if not positions:
                        continue
                        
                    num_windows = (chrom_end - chrom_start) // self.window_size + 1
                    density = np.zeros(num_windows)
                    
                    for start_pos, end_pos in positions:
                        # 计算转座子覆盖的窗口范围
                        first_window = (start_pos - chrom_start) // self.window_size
                        last_window = (end_pos - chrom_start) // self.window_size
                        
                        # 对于每个被覆盖的窗口，计算转座子在该窗口中的长度
                        for window_idx in range(first_window, last_window + 1):
                            if window_idx < 0 or window_idx >= num_windows:
                                continue
                            
                            window_start = chrom_start + window_idx * self.window_size
                            window_end = window_start + self.window_size - 1
                            
                            # 计算转座子与当前窗口的重叠长度
                            overlap_start = max(start_pos, window_start)
                            overlap_end = min(end_pos, window_end)
                            
                            if overlap_start <= overlap_end:
                                overlap_length = overlap_end - overlap_start + 1
                                density[window_idx] += overlap_length
                    
                    # 转换为每100kb的转座子长度 (bp/100kb)
                    density_per_100kb = density / (self.window_size / 100000)
                    
                    te_densities[te_type] = {
                        'positions': np.arange(num_windows) * self.window_size / 1e6,
                        'density': density_per_100kb
                    }
                
                if te_densities:
                    density_data[chrom] = te_densities
                    total_te_length = sum(end - start + 1 for positions in type_positions.values() for start, end in positions)
                    print(f"染色体 {chrom}: 有 {len(te_densities)} 种TE类型")
                    for te_type in te_densities:
                        print(f"    {te_type}: {len(type_positions[te_type])} 个区域，总长度 {self.te_total_lengths[te_type]} bp")
        
        print(f"\n最终密度数据包含 {len(density_data)} 个区域/染色体")
        for region, te_densities in density_data.items():
            print(f"  {region}: {len(te_densities)} 种TE类型 - {', '.join(te_densities.keys())}")
            
        return density_data
    
    def add_highlight_regions(self, ax, region_name):
        """在指定区域添加不同颜色的阴影高亮"""
        if not self.highlight_regions:
            return
        
        highlight_regions = self.parse_highlight_regions(self.highlight_regions)
        if not highlight_regions:
            return
        
        # 获取染色体名称（如果region_name包含其他信息）
        chrom_name = region_name
        if '_' in region_name:
            chrom_name = region_name.split('_')[0]
        
        # 添加高亮区域
        for i, region in enumerate(highlight_regions):
            if region['chromosome'] == chrom_name:
                # 创建矩形补丁
                rect = Rectangle(
                    (region['start_mb'], 0),  # 左下角坐标
                    region['end_mb'] - region['start_mb'],  # 宽度
                    100000,  # 高度（覆盖整个Y轴）
                    facecolor=region['color'],  # 使用区域特定的颜色
                    alpha=self.highlight_alpha,
                    edgecolor='none',
                    zorder=0  # 确保在背景层
                )
                ax.add_patch(rect)
                print(f"在 {region_name} 添加高亮区域 {i+1}: {region['start_mb']:.2f}Mb - {region['end_mb']:.2f}Mb, 颜色: {region['color']}")
    
    def plot_combined_density(self, density_data, max_chromosomes=5):
        """绘制所有类型转座子的组合图，没有图例"""
        print("正在生成组合密度图...")
        
        if not density_data:
            print("没有找到转座子数据可绘制")
            return
        
        # 定义TE类型的颜色 - 使用高对比度、明显区分的颜色方案
        # LINE改为黑色，Unknown改为绿色
        te_colors = {
            'hAT': '#FF0000',        # 纯红色 - 最明显
            'LTR': '#0000FF',        # 纯蓝色 - 与红色对比最强
            'LINE': '#000000',       # 纯黑色 - 改为黑色
            'Helitron': '#FFA500',   # 橙色 - 鲜艳
            'CACTA': '#800080',      # 紫色 - 与红蓝绿都不同
            'Mutator': '#008080',    # 青色 - 与橙色对比
            'PIF_Harbinger': '#FF1493', # 深粉色 - 鲜艳
            'Tc1_Mariner': '#8B4513',  # 棕色 - 与其他颜色都不同
            'Unknown': '#00AA00',    # 深绿色 - 改为绿色
            'Simple_repeat': '#4B0082' # 靛蓝色 - 深色
        }
        
        # 创建子图
        n_plots = min(len(density_data), max_chromosomes)
        
        # 调整图形大小
        fig_height = 2.1 * n_plots
        fig_width = 7.1
        
        fig, axes = plt.subplots(n_plots, 1, figsize=(fig_width, fig_height))
        if n_plots == 1:
            axes = [axes]
        
        for idx, (region_name, te_densities) in enumerate(list(density_data.items())[:n_plots]):
            ax = axes[idx]
            
            # 添加高亮区域（在绘制数据之前）
            self.add_highlight_regions(ax, region_name)
            
            # 调试信息：打印当前区域有多少TE类型
            print(f"绘制区域 {region_name}: 有 {len(te_densities)} 种TE类型")
            
            # 为每个TE类型绘制密度曲线 - 使用更细的线条但颜色更鲜艳
            for te_type, data in te_densities.items():
                color = te_colors.get(te_type, '#000000')
                print(f"  绘制 {te_type}: 颜色 {color}, 数据点 {len(data['positions'])} 个")
                ax.plot(data['positions'], data['density'], 
                       color=color, linewidth=1.0, alpha=1.0, label=te_type)  # 线宽1.0，完全不透明
            
            # 设置Y轴标题
            ax.set_ylabel('TE density', fontsize=14)
            ax.set_xlabel('', fontsize=14)

            # 设置Y轴范围：最大值100000，间隔20000
            ax.set_ylim(0, 100000)

            # 定义Y轴刻度
            y_ticks = np.arange(0, 100001, 20000)  # 从0到100000，间隔20000
            ax.set_yticks(y_ticks)
            
            # 格式化Y轴标签：将20000显示为20k，40000显示为40k，60000显示为60k
            y_tick_labels = []
            for tick in y_ticks:
                if tick == 0:
                    y_tick_labels.append('0')
                elif tick == 100000:
                    y_tick_labels.append('100k')
                else:
                    y_tick_labels.append(f'{int(tick/1000)}k')
            ax.set_yticklabels(y_tick_labels)
            
            # 设置X轴范围从0开始，使用固定的最大位置
            ax.set_xlim(left=0, right=self.max_x_position)
            
            # 设置固定的X轴刻度：0, 5, 10, 15, 20, 25, 30, 35
            x_ticks = np.arange(0, self.max_x_position + 1, 5)
            ax.set_xticks(x_ticks)
            
            # 设置刻度标签字体大小
            ax.tick_params(axis='both', which='major', labelsize=14, labelrotation=0)
            
            # 设置刻度线宽度和长度
            ax.tick_params(axis='both', which='major', width=1, length=3)
            ax.tick_params(axis='both', which='minor', width=1, length=3)
            
            # 移除网格线
            ax.grid(False)
            
            # 设置背景色
            ax.set_facecolor(self.bg_color)
            
            # 设置坐标轴边框线宽和颜色
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # 设置左框线和下框线
            ax.spines['left'].set_visible(True)
            ax.spines['left'].set_linewidth(1)
            ax.spines['bottom'].set_visible(True)
            ax.spines['bottom'].set_linewidth(1)
            
            # 设置框线颜色
            ax.spines['left'].set_color('black')
            ax.spines['bottom'].set_color('black')
            
            # 移除子图之间的空隙
            if idx < n_plots - 1:
                ax.set_xlabel('')  # 只保留最下面的x轴标签
        
        # 调整布局，移除子图之间的空隙
        plt.subplots_adjust(hspace=0.05)
        
        # 只保存PNG格式（移除SVG输出）
        output_file = f"{self.output_prefix}_combined.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                   facecolor=self.bg_color, edgecolor='none')
        print(f"组合密度图已保存为: {output_file}")
        
        # 显示图形以供调试
        if os.environ.get('DISPLAY') or os.name == 'nt':
            plt.show()
        else:
            print("警告: 没有图形显示环境，无法显示图片预览")
        
        # 关闭图形，避免在非交互式环境中显示
        plt.close()
    
    def create_legend_plot(self, density_data):
        """创建单独的图例图片"""
        print("正在生成图例图片...")
        
        # 定义TE类型的颜色 - 与主图保持一致
        te_colors = {
            'hAT': '#FF0000',        # 纯红色
            'LTR': '#0000FF',        # 纯蓝色
            'LINE': '#000000',       # 纯黑色
            'Helitron': '#FFA500',   # 橙色
            'CACTA': '#800080',      # 紫色
            'Mutator': '#008080',    # 青色
            'PIF_Harbinger': '#FF1493', # 深粉色
            'Tc1_Mariner': '#8B4513',  # 棕色
            'Unknown': '#00AA00',    # 深绿色
            'Simple_repeat': '#4B0082' # 靛蓝色
        }
        
        # 收集所有出现的TE类型
        all_te_types = set()
        for te_densities in density_data.values():
            all_te_types.update(te_densities.keys())
        
        print(f"找到的TE类型用于图例: {all_te_types}")
        
        if not all_te_types:
            print("没有找到TE类型，不生成图例")
            return
        
        # 创建图例图片
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis('off')
        
        # 创建图例元素 - 使用鲜艳的颜色
        legend_elements = []
        for te_type in sorted(all_te_types):
            color = te_colors.get(te_type, '#000000')
            legend_elements.append(plt.Line2D([0], [0], color=color, lw=3.0, label=te_type))  # 线宽增加到3.0
        
        # 添加图例
        ax.legend(handles=legend_elements, 
                 loc='center', 
                 ncol=min(5, len(legend_elements)),  # 改为最多5列
                 fontsize=14,  # 减小字体大小
                 frameon=True,
                 fancybox=False,  # 改为简单框
                 shadow=False,   # 去掉阴影
                 framealpha=1.0,
                 edgecolor='black')  # 黑色边框
        
        # 只保存PNG格式（移除SVG输出）
        output_file = f"{self.output_prefix}_legend.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight',
                   facecolor=self.bg_color, edgecolor='none')
        print(f"图例图片已保存为: {output_file}")
        
        plt.close()
    
    def save_density_data(self, density_data):
        """保存密度数据到文件"""
        output_file = f"{self.output_prefix}_data.tsv"
        
        with open(output_file, 'w') as f:
            # 修改列名以反映长度密度
            f.write("Chromosome\tTE_Type\tWindow_Start\tWindow_End\tLength_Density_per_100kb\n")
            
            for chrom, te_densities in density_data.items():
                for te_type, data in te_densities.items():
                    positions = data['positions']
                    density = data['density']
                    
                    for i in range(len(positions)):
                        window_start = int(positions[i] * 1e6)
                        window_end = window_start + self.window_size
                        f.write(f"{chrom}\t{te_type}\t{window_start}\t{window_end}\t{density[i]:.2f}\n")
        
        print(f"密度数据已保存为: {output_file}")
    
    def save_total_lengths(self):
        """保存每个转座子类型的总长度到文件"""
        output_file = f"{self.output_prefix}_total_lengths.tsv"
        
        with open(output_file, 'w') as f:
            f.write("TE_Type\tTotal_Length(bp)\n")
            for te_type, length in sorted(self.te_total_lengths.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{te_type}\t{length}\n")
        
        print(f"转座子类型总长度已保存为: {output_file}")
        
        # 同时在控制台输出
        print("\n各转座子类型总长度统计:")
        print("=" * 40)
        total_all_types = sum(self.te_total_lengths.values())
        for te_type, length in sorted(self.te_total_lengths.items(), key=lambda x: x[1], reverse=True):
            percentage = (length / total_all_types) * 100 if total_all_types > 0 else 0
            print(f"{te_type:<20} {length:>10} bp ({percentage:.2f}%)")
        print("=" * 40)
        print(f"{'总计':<20} {total_all_types:>10} bp (100.00%)")
    
    def run_analysis(self):
        """运行完整的分析流程"""
        if not self.parse_gff():
            return False
            
        density_data = self.calculate_density()
        
        if not density_data:
            print("没有计算出密度数据")
            return False
            
        # 生成组合密度图和单独的图例图片
        self.plot_combined_density(density_data)
        self.create_legend_plot(density_data)
        self.save_density_data(density_data)
        
        # 新增：保存转座子类型总长度
        self.save_total_lengths()
        
        return True

def main():
    parser = argparse.ArgumentParser(description='绘制转座子密度图', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='''
使用示例:
  # 分析整个基因组
  python te_density_plot.py -i repeats.gff
  
  # 分析特定染色体
  python te_density_plot.py -i repeats.gff -c chr1,chr2,chr3
  
  # 分析特定区域
  python te_density_plot.py -i repeats.gff -r chr1:1000000-5000000
  
  # 分析多个区域
  python te_density_plot.py -i repeats.gff -r chr1:1000000-5000000 chr2:2000000-8000000
  
  # 只分析特定类型的转座子
  python te_density_plot.py -i repeats.gff -t LTR,DNA,LINE,SINE
  
  # 使用不同的绘图风格
  python te_density_plot.py -i repeats.gff --style seaborn
  
  # 在指定区域添加不同颜色的阴影（支持3个紧挨的区域）
  python te_density_plot.py -i repeats.gff --highlight "chr1:1000000-2000000;chr1:2000000-3000000;chr1:3000000-4000000"
  
  # 设置X轴最大位置
  python te_density_plot.py -i repeats.gff --max-x 40
  
  # 组合使用参数
  python te_density_plot.py -i repeats.gff -c chr1 -t LTR,DNA -w 50000 -o chr1_LTR_DNA_density --highlight "chr1:1500000-2500000;chr1:2500000-3500000;chr1:3500000-4500000" --max-x 30
''')
    
    parser.add_argument('-i', '--input', required=True, help='输入GFF文件路径')
    parser.add_argument('-w', '--window', type=int, default=100000, 
                       help='窗口大小 (bp)，默认: 100000')
    parser.add_argument('-o', '--output', default='transposon_density',
                       help='输出文件前缀，默认: transposon_density')
    parser.add_argument('-c', '--chromosomes', 
                       help='指定分析的染色体，用逗号分隔，例如: chr1,chr2,chr3')
    parser.add_argument('-r', '--regions', nargs='+',
                       help='''指定分析的基因组区域，格式: 染色体:起始-结束
                       例如: chr1:1000000-5000000 或 chr2:2000000-8000000
                       可以指定多个区域''')
    parser.add_argument('-t', '--te_types', 
                       help='''指定分析的转座子类型，用逗号分隔
                       例如: LTR,DNA,LINE,SINE,Helitron
                       如果不指定，将使用预定义的10种类型''')
    parser.add_argument('--style', choices=['default', 'seaborn', 'dark', 'classic'],
                       default='default', help='绘图风格，默认: default')
    parser.add_argument('--highlight', 
                       help='''在指定区域添加不同颜色的阴影高亮，格式: 染色体:起始-结束;染色体:起始-结束
                       例如: chr1:1000000-2000000;chr1:2000000-3000000;chr1:3000000-4000000
                       支持3个紧挨的区域，会自动分配不同颜色''')
    parser.add_argument('--max-x', type=int, default=35,
                       help='X轴最大位置 (Mb)，默认: 35')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 {args.input} 不存在")
        return
    
    # 解析染色体列表
    chromosomes = None
    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(',')]
        print(f"目标染色体: {chromosomes}")
    
    # 解析区域列表
    regions = None
    if args.regions:
        regions = args.regions
        print(f"目标区域: {regions}")
    
    # 解析TE类型列表
    te_types = None
    if args.te_types:
        te_types = [t.strip() for t in args.te_types.split(',')]
        print(f"目标TE类型: {te_types}")
    else:
        print("使用预定义的10种TE类型: hAT, CACTA, Mutator, PIF_Harbinger, Tc1_Mariner, Helitron, LINE, LTR, Unknown, Simple_repeat")
    
    # 创建绘图器并运行分析
    plotter = TransposonDensityPlotter(
        gff_file=args.input,
        window_size=args.window,
        output_prefix=args.output,
        chromosomes=chromosomes,
        regions=regions,
        te_types=te_types,
        style=args.style,
        highlight_regions=args.highlight,  # 新增高亮区域参数
        max_x_position=args.max_x  # 新增最大X轴位置参数
    )
    
    plotter.run_analysis()

if __name__ == "__main__":
    main()