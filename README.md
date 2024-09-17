# Single-cell RNA-seq data imputation using Bi-level Feature Propagation

<p align="center">
    <a href="https://pytorch.org/" alt="PyTorch">
    <img src="https://img.shields.io/badge/PyTorch-%23EE4C2C.svg?e&logo=PyTorch&logoColor=white" /></a>
<img src="https://img.shields.io/badge/-Briefings%20in%20Bioinformatics-blue" />

The official source code for [**Single-cell RNA-seq data imputation using Bi-level Feature Propagation**](https://academic.oup.com/bib/article/25/3/bbae209/7665119), accepted at Briefings in Bioinformatics (Volume 25, May 2024).

## Environment
```
conda create -n scBFP python=3.10
source activate scBFP
conda install pytorch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 pytorch-cuda=12.1 -c pytorch -c nvidia
conda install pytorch-sparse -c pyg
pip install munkres anndata scanpy leidenalg
```

## Download data
```
mkdir dataset
cd dataset
pip install gdown
gdown https://drive.google.com/uc?id=1raqlykXvm5wHjam1Up0SHYT-7gq7coz4
gdown https://drive.google.com/uc?id=1pilLsl2N1HX_US_Y6X6eAwmXaPso_1Mu
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
python main.py --name citeseq --gene_k 20 --cell_k 20 --gene_iter 100 --cell_iter 100
```
```
python main.py --name multiome --gene_k 20 --cell_k 20 --gene_iter 100 --cell_iter 100
```
