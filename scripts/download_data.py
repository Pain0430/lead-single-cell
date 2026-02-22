#!/usr/bin/env python3
"""
单细胞数据下载脚本
Single-Cell Data Download Script

下载公开的脑细胞单细胞数据

作者: Pain's AI Assistant
日期: 2026-02-22
"""

import os
import subprocess
import requests
from urllib.parse import urljoin

# 配置
PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(PROJECT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

# 可用数据集
DATASETS = {
    "tasic_mouse_cortex": {
        "name": "Adult Mouse Cortical Cell Taxonomy",
        "source": "GEO",
        "geo_id": "GSE71585",
        "description": "Tasic et al. 2016 - 1800+ cells from mouse cortex",
        "cell_types": ["neuron", "astrocyte", "microglia", "oligodendrocyte"],
    },
    "zeisel_brain": {
        "name": "Brain Structure",
        "source": "GEO", 
        "geo_id": "GSE60361",
        "description": "Zeisel et al. 2015 - Mouse brain cells",
        "cell_types": ["neuron", "astrocyte", "microglia", "oligodendrocyte", "endothelial"],
    },
    "allenMouseBrain": {
        "name": "Allen Mouse Brain Atlas",
        "source": "Allen Institute",
        "url": "https://celltypes.brain-map.org/",
        "description": "Comprehensive mouse brain cell atlas",
        "cell_types": ["neuron", "astrocyte", "microglia", "oligodendrocyte", "endothelial", "pericyte"],
    },
}


def check_dependencies():
    """检查必要的工具"""
    print("检查依赖...")
    
    # 检查python包
    try:
        import scanpy
        print("  ✓ scanpy 已安装")
    except ImportError:
        print("  ✗ scanpy 未安装")
        print("  安装: pip install scanpy")
    
    # 检查wget/curl
    result = subprocess.run(["which", "wget"], capture_output=True)
    if result.returncode == 0:
        print("  ✓ wget 可用")
    else:
        result = subprocess.run(["which", "curl"], capture_output=True)
        if result.returncode == 0:
            print("  ✓ curl 可用")


def download_from_url(url, filename):
    """下载文件"""
    print(f"下载: {filename}...")
    
    try:
        response = requests.get(url, timeout=300)
        response.raise_for_status()
        
        filepath = os.path.join(DATA_DIR, filename)
        with open(filepath, 'wb') as f:
            f.write(response.content)
        
        size = os.path.getsize(filepath) / 1024 / 1024  # MB
        print(f"  ✓ 完成 ({size:.1f} MB)")
        return filepath
        
    except Exception as e:
        print(f"  ✗ 失败: {e}")
        return None


def create_reference_data():
    """
    创建参考数据用于分析
    
    由于真实单细胞数据很大(GB级别)，创建一个参考数据集
    包含已知的细胞类型标记基因和铅毒性相关基因
    """
    print("\n创建参考数据...")
    
    # 脑细胞类型标记基因
    cell_type_markers = {
        "Excitatory_Neuron": ["Slc17a7", "Cux2", "Rorb", "Fezf2", "Tle4"],
        "Inhibitory_Neuron": ["Gad1", "Gad2", "Pvalb", "Sst", "Htr3a", "Vip"],
        "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1", "S100b"],
        "Microglia": ["Cx3cr1", "Itgam", "Aif1", "Tmem119", "P2ry12"],
        "Oligodendrocyte": ["Olig2", "Mbp", "Plp1", "Sox10", "Mog"],
        "OPC": ["Pdgfra", "Cspg4", "Olig1", "Nkx2-2"],
        "Endothelial": ["Cldn5", "Flt1", "Pecam1", "Cdh5", "Vwf"],
        "Pericyte": ["Pdgfrb", "Col1a1", "Acta2", "Rgs5"],
    }
    
    # 铅毒性相关基因
    lead_genes = {
        "oxidative_stress": ["Sod1", "Cat", "Gpx1", "Nqo1", "Hmox1", "Gclc"],
        "inflammation": ["Il6", "Tnf", "Cxcl8", "Ccl2", "Nfkb1", "Ptgs2"],
        "apoptosis": ["Bcl2", "Bax", "Casp3", "Trp53", "Mdm2"],
        "synaptic": ["Snap25", "Syn1", "Dlg4", "Grin1", "Gria1"],
        "metabolism": ["Mtor", "Ampk", "Sirt1", "Ppargc1a", "Hif1a"],
    }
    
    # 肿瘤机制基因 (从肿瘤研究借鉴)
    cancer_mechanism_genes = {
        "metabolic_reprogramming": ["Hk2", "Pkm", "Ldha", "Pdha1", "Idh1"],
        "immune_evasion": ["Cd274", "Pdcd1", "Ctla4", "Tgfb1", "Il10"],
        "dna_damage": ["Atm", "ATR", "Brca1", "Brca2", "Xrcc1"],
        "cellular_senescence": ["Cdkn1a", "Cdkn2a", "Glb1", "Il6", "Il8"],
        "angiogenesis": ["Vegfa", "Kdr", "Angpt1", "Angpt2", "Tie1"],
    }
    
    # 保存为JSON
    import json
    
    reference_data = {
        "cell_type_markers": cell_type_markers,
        "lead_genes": lead_genes,
        "cancer_mechanism_genes": cancer_mechanism_genes,
    }
    
    filepath = os.path.join(DATA_DIR, "reference_genes.json")
    with open(filepath, 'w') as f:
        json.dump(reference_data, f, indent=2)
    
    print(f"  ✓ 参考基因数据已保存: {filepath}")
    
    return reference_data


def create_analysis_pipeline():
    """创建分析流程说明"""
    
    pipeline = """
# 铅神经毒性单细胞分析流程

## 1. 数据获取

### 方案A: 使用已有数据集
```python
import scanpy as sc

# 从Hemberg Lab数据集
# 需要安装: Bioconductor scRNAseq包

# 或直接下载:
# GSE71585: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585
# GSE60361: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361
```

### 方案B: Allen Brain Atlas
```python
# 使用abc_atlas_access包
pip install abc-atlas-access
```

## 2. 基本分析流程

```python
import scanpy as sc

# 读取数据
adata = sc.read_h5ad('data.h5ad')

# 质控
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 归一化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(而降维)

# 降维
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 细胞类型注释
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
```

## 3. 铅毒性分析

```python
# 标记细胞类型
cell_type_genes = {
    'Astrocyte': ['Gfap', 'Aqp4', 'Aldh1l1'],
    'Microglia': ['Cx3cr1', 'Itgam', 'Tmem119'],
    ...
}

# 计算铅毒性基因集得分
sc.tl.score_genes(adata, gene_list=lead_genes['oxidative_stress'], 
                   score_name='ox_stress')
```

## 4. 细胞间通讯分析

```python
# 使用CellChat或CellPhoneDB
import cellchat as cc
cc.plot_net_embedding()
```
"""
    
    filepath = os.path.join(PROJECT_DIR, "ANALYSIS_PIPELINE.md")
    with open(filepath, 'w') as f:
        f.write(pipeline)
    
    print(f"\n✓ 分析流程文档已保存: {filepath}")


def main():
    print("=" * 60)
    print("📥 单细胞数据下载")
    print("=" * 60)
    
    # 检查依赖
    check_dependencies()
    
    # 创建参考数据
    reference_data = create_reference_data()
    
    # 创建分析流程
    create_analysis_pipeline()
    
    # 打印可用数据集信息
    print("\n" + "=" * 60)
    print("📊 可用数据集")
    print("=" * 60)
    
    for key, info in DATASETS.items():
        print(f"\n{key}:")
        print(f"  名称: {info['name']}")
        print(f"  来源: {info['source']}")
        print(f"  描述: {info['description']}")
        print(f"  细胞类型: {', '.join(info['cell_types'])}")
    
    print("""
    
⚠️ 注意:
1. 真实单细胞数据很大(通常几GB)
2. 建议使用服务器环境下载
3. 推荐先从以下来源获取:
   - GEO: GSE71585 (小鼠皮层)
   - Allen Brain Atlas (全脑)
   - 10x Genomics官方案例数据

📝 下一步:
1. 手动下载数据文件
2. 使用scanpy进行分析
3. 将铅毒性基因映射到不同细胞类型

""")
    
    print("=" * 60)
    print("✅ 准备完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
