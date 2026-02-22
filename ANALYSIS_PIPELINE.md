
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
