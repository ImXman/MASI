# MASI

MASI: marker assisted single-cell expression data integration

Manuscript: Fast model-free inetgration and transfer learning via MASI for single-cell expression data

### 1. Brief description
MASI utilizes robust marker idenfication to identify marker genes from reference data and transfers cell-type labels to target data through MACA.

![alt text](https://github.com/ImXman/MASI/blob/main/MASI/Figure%201.jpg?raw=true)

### 2. Install requirement packages
    pip install scanpy cosg rpy2 sccaf
    pip install fa2##install if doing integrative lineage analysis
    
    ##install Seurat and RobustRankAggreg separately in R
    install.packages('Seurat')
    install.packages('RobustRankAggreg')
    
### 3. Usage
    import MASI as masi
    ##step 1 transform gene expression matrix to cell-type score matrix
    ##scores can further be used for visualization and other downstream analyses
    scores, labels = masi.gene2cell(ad=ad,cell_markers=cell_markers,use_weight=True)
    ##step 2 clustering and parallel annotation
    annotation= masi.parallel(scores=scores,labels=labels,batch_size=50000)

### 4. Reproduce results in manuscript
Please see tutorials at https://github.com/ImXman/MASI/tree/main/tutorial
    

### 5. References
    1. Xu, Y., et al. "MACA: Marker-based automatic cell-type annotation for single cell expression data." (2021).
    2. Kolde, Raivo, et al. "Robust rank aggregation for gene list integration and meta-analysis." Bioinformatics 28.4 (2012): 573-580.
    3. Wolf, F. Alexander, Philipp Angerer, and Fabian J. Theis. "SCANPY: large-scale single-cell gene expression data analysis." Genome biology 19.1 (2018): 1-5.
    4. Butler, Andrew, et al. "Integrating single-cell transcriptomic data across different conditions, technologies, and species." Nature biotechnology 36.5 (2018): 411-420.
