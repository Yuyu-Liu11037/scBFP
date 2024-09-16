# Single-cell RNA-seq data imputation using Bi-level Feature Propagation

<p align="center">
    <a href="https://pytorch.org/" alt="PyTorch">
    <img src="https://img.shields.io/badge/PyTorch-%23EE4C2C.svg?e&logo=PyTorch&logoColor=white" /></a>
<img src="https://img.shields.io/badge/-Briefings%20in%20Bioinformatics-blue" />

The official source code for [**Single-cell RNA-seq data imputation using Bi-level Feature Propagation**](https://academic.oup.com/bib/article/25/3/bbae209/7665119), accepted at Briefings in Bioinformatics (Volume 25, May 2024).

## Environment
```
conda env create -f environment.yml
```

## Download data
```
mkdir dataset
pip install gdown
gdown https://drive.google.com/uc?id=1raqlykXvm5wHjam1Up0SHYT-7gq7coz4
```

## Hyperparameters

`--name:`
Name of the dataset.  
usage example :`--dataset baron_mouse`

`--gene_k:`
Number of neighbors in gene-gene graph  
usage example :`--k 10`

`--cell_k:`
Number of neighbors in cell-cell graph  
usage example :`--k 10`

`--gene_iter:`
Number of iterations in feature propagation using gene-gene graph  
usage example :`--iter 10`

`--cell_iter:`
Number of iterations in feature propagation using cell-cell graph  
usage example :`--iter 40`

Using above hyper-parmeters, you can run our model with following codes  

```
python main.py --name baron_mouse --gene_k 10 --cell_k 10 --gene_iter 10 --cell_iter 40
```

