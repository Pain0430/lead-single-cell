#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
铅神经毒性轨迹分析模块
Lead Neurotoxicity Trajectory Analysis

功能：
1. 细胞类型组成分析
2. 铅暴露相关的伪时间轨迹分析
3. 差异表达分析
4. 轨迹可视化

作者: Pain's AI Assistant
日期: 2026-02-23
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from scipy import stats

# 设置
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

# ============================================================================
# 配置
# ============================================================================

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(PROJECT_DIR, "data")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 脑细胞类型标记基因
CELL_TYPE_MARKERS = {
    "Neuron": ["SNAP25", "SYN1", "MAP2", "NEUN", "RBFOX3", "SLC17A6", "GAD1"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B", "CD44"],
    "Microglia": ["CX3CR1", "ITGAM", "IBA1", "TMEM119", "P2RY12", "CD68"],
    "Oligodendrocyte": ["OLIG2", "MBP", "PLP1", "SOX10", "MOG", "CNP"],
    "Endothelial": ["CLDN5", "FLT1", "PECAM1", "CDH5", "VWF", "CAV1"],
    "Pericyte": ["PDGFRB", "COL1A1", "COL3A1", "ACTA2", "RGS5"],
}

# 铅毒性相关基因
LEAD_TOXICITY_GENES = {
    "oxidative_stress": ["SOD1", "SOD2", "CAT", "GPX1", "GPX4", "NQO1", "HMOX1", "GSTA1"],
    "inflammation": ["IL1B", "IL6", "TNF", "NFKB1", "PTGS2", "CXCL8"],
    "neurotoxicity": ["APP", "MAPT", "BDNF", "MAPK1", "MAPK3", "CASP3", "MAPK8"],
    "apoptosis": ["BCL2", "BAX", "CASP3", "CASP9", "TP53", "PUMA", "NOXA"],
    "metal_transport": ["MT1A", "MT2A", "SLC11A2", "SLC39A8", "ATP13A2"],
    "synapse": ["SNAP25", "SYN1", "PSD95", "NLGN1", "NRXN1", "DLG4"],
}


# ============================================================================
# 数据加载
# ============================================================================

def load_single_cell_data():
    """加载单细胞数据"""
    print("加载单细胞数据...")
    
    # 尝试加载真实数据
    h5_files = [f for f in os.listdir(DATA_DIR) if f.endswith('.h5')]
    
    if h5_files:
        # 有h5文件，尝试用scanpy加载（如果可用）
        try:
            import scanpy as sc
            adata = sc.read(os.path.join(DATA_DIR, h5_files[0]))
            print(f"  已加载真实数据: {h5_files[0]}")
            return adata
        except ImportError:
            pass
    
    # 加载CSV数据
    expr_file = os.path.join(DATA_DIR, "sc_expression.npy")
    obs_file = os.path.join(DATA_DIR, "sc_obs.csv")
    var_file = os.path.join(DATA_DIR, "sc_var.csv")
    
    if os.path.exists(expr_file) and os.path.exists(obs_file):
        expression = np.load(expr_file)
        obs = pd.read_csv(obs_file)  # 不设置index_col，保留原始列
        var = pd.read_csv(var_file)
        
        n_cells = expression.shape[0]
        n_genes = expression.shape[1]
        
        print(f"  已加载CSV数据: {n_cells} cells x {n_genes} genes")
        
        # 处理obs与expression行数不匹配的情况
        if len(obs) > n_cells:
            print(f"  警告: obs行数({len(obs)}) > expression行数({n_cells})，取前{n_cells}行")
            obs = obs.head(n_cells)
        elif len(obs) < n_cells:
            # 用现有obs扩展
            print(f"  警告: obs行数({len(obs)}) < expression行数({n_cells})，补充缺失")
            extra_rows = n_cells - len(obs)
            extra_obs = pd.DataFrame({
                'cell_type': ['unknown'] * extra_rows,
                'cell_id': [f'cell_{i}' for i in range(len(obs), n_cells)]
            })
            obs = pd.concat([obs, extra_obs], ignore_index=True)
        
        # 确保cell_type列存在
        if 'cell_type' not in obs.columns:
            obs['cell_type'] = 'unknown'
        
        # 获取基因名
        gene_names = var['gene_name'].tolist() if 'gene_name' in var.columns else None
        if gene_names is None or len(gene_names) < n_genes:
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
        gene_names = gene_names[:n_genes]
        
        # 创建简单的数据框
        data = {
            'expression': expression,
            'obs': obs,
            'var': var,
            'cell_ids': obs['cell_id'].tolist() if 'cell_id' in obs.columns else [f"Cell_{i}" for i in range(n_cells)],
            'gene_ids': gene_names
        }
        return data
    
    # 创建模拟数据
    print("  未找到数据文件，创建模拟数据...")
    return create_simulated_data()


def create_simulated_data(n_cells=1000, n_genes=200):
    """创建模拟单细胞数据"""
    np.random.seed(42)
    
    # 定义细胞类型比例
    cell_types = ['Neuron', 'Astrocyte', 'Microglia', 'Oligodendrocyte', 'Endothelial']
    cell_type_probs = [0.4, 0.25, 0.15, 0.12, 0.08]
    
    # 分配细胞类型
    cell_type_labels = np.random.choice(cell_types, size=n_cells, p=cell_type_probs)
    
    # 为每种细胞类型生成特异性表达
    expression = np.zeros((n_cells, n_genes))
    
    # 获取标记基因
    all_markers = []
    for markers in CELL_TYPE_MARKERS.values():
        all_markers.extend(markers)
    
    # 基因名
    gene_ids = all_markers + [f"Gene_{i}" for i in range(n_genes - len(all_markers))]
    gene_ids = gene_ids[:n_genes]
    
    # 为每种细胞类型设置标记基因高表达
    for ct_idx, cell_type in enumerate(cell_types):
        ct_mask = cell_type_labels == cell_type
        ct_indices = np.where(ct_mask)[0]
        
        # 获取该细胞类型的标记基因
        markers = CELL_TYPE_MARKERS.get(cell_type, [])
        
        for i, gene in enumerate(gene_ids[:len(markers)]):
            # 标记基因高表达
            expression[ct_indices, i] = np.random.lognormal(2, 0.5, len(ct_indices))
        
        # 其他基因背景表达
        for i in range(len(markers), n_genes):
            expression[ct_indices, i] = np.random.lognormal(0, 0.5, len(ct_indices))
    
    # 添加铅暴露效应（模拟）
    lead_exposure = np.random.lognormal(1, 0.5, n_cells)
    
    # 对铅毒性基因添加暴露效应
    toxicity_genes = []
    for genes in LEAD_TOXICITY_GENES.values():
        toxicity_genes.extend(genes)
    
    for gene in toxicity_genes:
        if gene in gene_ids:
            gene_idx = gene_ids.index(gene)
            expression[:, gene_idx] += lead_exposure * np.random.uniform(0.1, 0.3)
    
    # 元数据
    obs = pd.DataFrame({
        'cell_type': cell_type_labels,
        'lead_exposure': lead_exposure,
        'cell_id': [f"Cell_{i}" for i in range(n_cells)]
    }, index=[f"Cell_{i}" for i in range(n_cells)])
    
    var = pd.DataFrame({'gene_name': gene_ids}, index=gene_ids)
    
    data = {
        'expression': expression,
        'obs': obs,
        'var': var,
        'cell_ids': obs.index.tolist(),
        'gene_ids': gene_ids
    }
    
    print(f"  模拟数据创建完成: {n_cells} cells x {n_genes} genes")
    return data


# ============================================================================
# 细胞类型分析
# ============================================================================

def analyze_cell_type_composition(data):
    """分析细胞类型组成"""
    print("\n分析细胞类型组成...")
    
    obs = data['obs']
    
    # 细胞类型计数
    type_counts = obs['cell_type'].value_counts()
    type_proportions = obs['cell_type'].value_counts(normalize=True)
    
    results = {
        'counts': type_counts.to_dict(),
        'proportions': type_proportions.to_dict(),
        'total_cells': len(obs)
    }
    
    # 可视化
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 饼图
    colors = plt.cm.Set3(np.linspace(0, 1, len(type_counts)))
    axes[0].pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', colors=colors)
    axes[0].set_title('Cell Type Composition', fontsize=14)
    
    # 条形图
    bars = axes[1].bar(type_counts.index, type_counts.values, color=colors)
    axes[1].set_xlabel('Cell Type', fontsize=12)
    axes[1].set_ylabel('Cell Count', fontsize=12)
    axes[1].set_title('Cell Type Distribution', fontsize=14)
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cell_type_composition.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  细胞类型组成图已保存")
    
    return results


# ============================================================================
# 降维分析
# ============================================================================

def perform_dimensionality_analysis(data, n_pcs=20):
    """执行降维分析 (PCA + t-SNE)"""
    print("\n执行降维分析...")
    
    expression = data['expression']
    
    # 标准化
    expression_std = (expression - expression.mean(axis=0)) / (expression.std(axis=0) + 1e-8)
    
    # PCA
    print("  计算PCA...")
    pca = PCA(n_components=min(n_pcs, min(expression_std.shape) - 1))
    pca_result = pca.fit_transform(expression_std)
    
    # t-SNE (基于PCA结果)
    print("  计算t-SNE...")
    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    tsne_result = tsne.fit_transform(pca_result[:, :10])  # 使用前10个PC
    
    # 保存结果到obs
    data['obs']['PC1'] = pca_result[:, 0]
    data['obs']['PC2'] = pca_result[:, 1]
    data['obs']['tSNE1'] = tsne_result[:, 0]
    data['obs']['tSNE2'] = tsne_result[:, 1]
    
    # 可视化
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # PCA图
    for cell_type in data['obs']['cell_type'].unique():
        mask = data['obs']['cell_type'] == cell_type
        axes[0].scatter(
            data['obs'].loc[mask, 'PC1'],
            data['obs'].loc[mask, 'PC2'],
            label=cell_type, alpha=0.6, s=20
        )
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    axes[0].set_title('PCA of Single Cells', fontsize=14)
    axes[0].legend()
    
    # t-SNE图
    for cell_type in data['obs']['cell_type'].unique():
        mask = data['obs']['cell_type'] == cell_type
        axes[1].scatter(
            data['obs'].loc[mask, 'tSNE1'],
            data['obs'].loc[mask, 'tSNE2'],
            label=cell_type, alpha=0.6, s=20
        )
    axes[1].set_xlabel('t-SNE 1')
    axes[1].set_ylabel('t-SNE 2')
    axes[1].set_title('t-SNE of Single Cells', fontsize=14)
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'dimensionality_reduction.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # PCA方差解释
    print(f"  PCA前10个PC解释方差: {sum(pca.explained_variance_ratio_[:10])*100:.1f}%")
    
    return {
        'pca_variance_explained': pca.explained_variance_ratio_.tolist(),
        'pca_result': pca_result,
        'tsne_result': tsne_result
    }


# ============================================================================
# 铅暴露分析
# ============================================================================

def analyze_lead_exposure_effects(data):
    """分析铅暴露对不同细胞类型的影响"""
    print("\n分析铅暴露效应...")
    
    obs = data['obs']
    expression = data['expression']
    gene_ids = data['gene_ids']
    
    # 检查是否有铅暴露数据
    if 'lead_exposure' not in obs.columns:
        print("  未找到铅暴露数据，跳过分析")
        return None
    
    results = []
    
    # 按细胞类型分析
    for cell_type in obs['cell_type'].unique():
        ct_mask = obs['cell_type'] == cell_type
        ct_exposure = obs.loc[ct_mask, 'lead_exposure']
        
        # 分析铅毒性基因表达与暴露的关系
        for category, genes in LEAD_TOXICITY_GENES.items():
            category_genes = [g for g in genes if g in gene_ids]
            
            if not category_genes:
                continue
            
            # 获取基因索引
            gene_indices = [gene_ids.index(g) for g in category_genes]
            
            # 计算该类别的平均表达
            gene_expr = expression[ct_mask.values, :][:, gene_indices].mean(axis=1)
            
            # 相关性分析
            if len(ct_exposure) > 10 and len(gene_expr) > 10:
                r, p = stats.spearmanr(ct_exposure, gene_expr)
                
                results.append({
                    'cell_type': cell_type,
                    'gene_category': category,
                    'n_genes': len(category_genes),
                    'correlation': r,
                    'p_value': p,
                    'significant': p < 0.05
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # 可视化
        fig, ax = plt.subplots(figsize=(12, 6))
        
        pivot = results_df.pivot_table(
            index='cell_type', 
            columns='gene_category', 
            values='correlation',
            aggfunc='mean'
        )
        
        sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                    vmin=-0.5, vmax=0.5, ax=ax)
        ax.set_title('Lead Exposure vs Gene Expression Correlation\nby Cell Type', fontsize=14)
        ax.set_xlabel('Gene Category')
        ax.set_ylabel('Cell Type')
        
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, 'lead_exposure_effects.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  显著相关结果:")
        sig_results = results_df[results_df['significant']]
        for _, row in sig_results.iterrows():
            print(f"    {row['cell_type']} - {row['gene_category']}: r={row['correlation']:.3f}")
    
    return results_df


# ============================================================================
# 伪时间轨迹分析
# ============================================================================

def pseudo_time_trajectory(data):
    """简化的伪时间轨迹分析"""
    print("\n执行伪时间轨迹分析...")
    
    obs = data['obs']
    expression = data['expression']
    gene_ids = data['gene_ids']
    
    # 选择关键基因进行轨迹分析
    trajectory_genes = []
    for genes in LEAD_TOXICITY_GENES.values():
        trajectory_genes.extend(genes)
    trajectory_genes = [g for g in trajectory_genes if g in gene_ids]
    
    if len(trajectory_genes) < 5:
        print("  铅毒性基因数据不足，使用PCA进行轨迹分析")
        # 使用PCA作为替代
        pca_result = obs[['PC1', 'PC2']].values if 'PC1' in obs.columns else None
        
        if pca_result is None:
            expression_std = (expression - expression.mean(axis=0)) / (expression.std(axis=0) + 1e-8)
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(expression_std)
    else:
        gene_indices = [gene_ids.index(g) for g in trajectory_genes]
        traj_expr = expression[:, gene_indices]
        
        # PCA用于轨迹
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(traj_expr)
    
    # 添加伪时间（基于第一个PC）
    pseudo_time = pca_result[:, 0]
    obs['pseudo_time'] = pseudo_time
    
    # 可视化
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # 伪时间分布
    axes[0].hist(pseudo_time, bins=50, color='steelblue', alpha=0.7)
    axes[0].set_xlabel('Pseudotime', fontsize=12)
    axes[0].set_ylabel('Cell Count', fontsize=12)
    axes[0].set_title('Pseudotime Distribution', fontsize=14)
    
    # 伪时间 vs 细胞类型
    cell_types = obs['cell_type'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
    
    for i, ct in enumerate(cell_types):
        mask = obs['cell_type'] == ct
        axes[1].scatter(
            obs.loc[mask, 'pseudo_time'],
            np.random.randn(mask.sum()) * 0.1,  # 添加抖动
            label=ct, alpha=0.5, s=10, c=[colors[i]]
        )
    axes[1].set_xlabel('Pseudotime', fontsize=12)
    axes[1].set_ylabel('Cell Type (jittered)', fontsize=12)
    axes[1].set_title('Pseudotime by Cell Type', fontsize=14)
    axes[1].legend()
    
    # 伪时间 vs 铅暴露
    if 'lead_exposure' in obs.columns:
        axes[2].scatter(obs['pseudo_time'], obs['lead_exposure'], 
                        alpha=0.3, s=10, c='steelblue')
        
        # 添加趋势线
        from scipy.ndimage import uniform_filter1d
        sorted_idx = np.argsort(obs['pseudo_time'])
        trend = uniform_filter1d(obs['lead_exposure'].values[sorted_idx], size=50)
        axes[2].plot(obs['pseudo_time'].values[sorted_idx], trend, 
                     color='red', linewidth=2, label='Trend')
        axes[2].set_xlabel('Pseudotime', fontsize=12)
        axes[2].set_ylabel('Lead Exposure', fontsize=12)
        axes[2].set_title('Pseudotime vs Lead Exposure', fontsize=14)
        axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'pseudo_time_trajectory.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # 早期vs晚期细胞差异分析
    if 'lead_exposure' in obs.columns:
        early_mask = pseudo_time < np.percentile(pseudo_time, 33)
        late_mask = pseudo_time > np.percentile(pseudo_time, 67)
        
        early_lead = obs.loc[early_mask, 'lead_exposure'].mean()
        late_lead = obs.loc[late_mask, 'lead_exposure'].mean()
        
        print(f"  早期伪时间铅暴露: {early_lead:.3f}")
        print(f"  晚期伪时间铅暴露: {late_lead:.3f}")
        
        if early_lead != late_lead:
            direction = "增加" if late_lead > early_lead else "减少"
            print(f"  铅暴露随伪时间{direction}")
    
    return {
        'pseudo_time': pseudo_time,
        'pca_result': pca_result
    }


# ============================================================================
# 差异表达分析
# ============================================================================

def differential_expression_analysis(data, group_col='cell_type'):
    """差异表达分析（简化版）"""
    print("\n执行差异表达分析...")
    
    obs = data['obs']
    expression = data['expression']
    gene_ids = data['gene_ids']
    
    results = []  # 比较每种细胞类型与其他类型
    
    cell_types = obs[group_col].unique()
    
    for cell_type in cell_types:
        # 该类型 vs 其他
        group1_mask = obs[group_col] == cell_type
        group2_mask = ~group1_mask
        
        # 对每个基因进行t检验
        for gene_idx, gene in enumerate(gene_ids):
            expr1 = expression[group1_mask.values, gene_idx]
            expr2 = expression[group2_mask.values, gene_idx]
            
            if len(expr1) < 5 or len(expr2) < 5:
                continue
            
            # t检验
            try:
                stat, p = stats.ttest_ind(expr1, expr2, equal_var=False)
                
                # 效应量 (Cohen's d)
                mean1, mean2 = expr1.mean(), expr2.mean()
                std1, std2 = expr1.std(), expr2.std()
                pooled_std = np.sqrt((std1**2 + std2**2) / 2)
                cohens_d = (mean1 - mean2) / pooled_std if pooled_std > 0 else 0
                
                results.append({
                    'gene': gene,
                    'cell_type': cell_type,
                    'mean_expr': mean1,
                    'mean_other': mean2,
                    'log2_fold_change': np.log2(mean1 / (mean2 + 1e-8) + 1e-8),
                    'p_value': p,
                    'cohens_d': cohens_d,
                    'significant': p < 0.05
                })
            except:
                continue
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # 保存结果
        results_df.to_csv(os.path.join(OUTPUT_DIR, 'differential_expression.csv'), index=False)
        
        # 绘制热图（Top基因）
        top_genes = results_df[results_df['significant']].groupby('gene')['p_value'].min()
        top_genes = top_genes.nsmallest(20).index.tolist()
        
        if len(top_genes) > 0:
            top_df = results_df[results_df['gene'].isin(top_genes)]
            pivot = top_df.pivot_table(index='gene', columns='cell_type', values='log2_fold_change')
            
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax, 
                       annot=True, fmt='.1f', annot_kws={'size': 8})
            ax.set_title('Differential Expression (Log2FC)', fontsize=14)
            
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, 'differential_expression_heatmap.png'), 
                       dpi=150, bbox_inches='tight')
            plt.close()
        
        print(f"  找到 {len(results_df[results_df['significant']])} 个显著差异表达基因")
    
    return results_df


# ============================================================================
# 主函数
# ============================================================================

def main():
    """主函数"""
    print("=" * 70)
    print("铅神经毒性单细胞轨迹分析")
    print("Lead Neurotoxicity Single-Cell Trajectory Analysis")
    print("=" * 70)
    
    # 加载数据
    data = load_single_cell_data()
    
    # 1. 细胞类型组成分析
    comp_results = analyze_cell_type_composition(data)
    print(f"\n细胞类型组成: {comp_results['counts']}")
    
    # 2. 降维分析
    dim_results = perform_dimensionality_analysis(data)
    
    # 3. 铅暴露效应分析
    lead_results = analyze_lead_exposure_effects(data)
    
    # 4. 伪时间轨迹分析
    traj_results = pseudo_time_trajectory(data)
    
    # 5. 差异表达分析
    diff_results = differential_expression_analysis(data)
    
    print("\n" + "=" * 70)
    print("分析完成! 结果已保存至 output/")
    print("=" * 70)


if __name__ == "__main__":
    main()
