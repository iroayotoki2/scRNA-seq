# Analysis of single-cell RNA-seq (scRNA-seq) data to investigate the effects of influenza infection on the mouse respiratory tract

## Biological Background

Single cell Transcriptomics(scRNA-seq)is the process of measuring the gene expression of many individual cells within a tissue population to study their functions, unique identities and interactions within the tissue (What Is Single Cell Transcriptomics? - Single Cell Discoveries, n.d.). The purpose of this project is to research the best methods, packages, and parameters for single cell gene expression analysis, cluster and annotate and perform differential expression of specific cell types to study the effects of the Influenza virus on these cells.

The Influenza virus refers to a family of contagious viruses collectively referred to as Orthomyxoviridae, they affect a wide range of animals from birds to other mammals and humans often leading to respiratory illnesses characterized by fever, cough, fatigue, and muscle aches, often leading to pneumonia or severe system illness (Bouvier & Palese, 2008; Short et al., 2015).

Influenza infection affects multiple cell types within the respiratory tract, leading to cellular damage, tissue disruption, and systemic consequences for the host. Understanding how different cell populations respond to infection is essential for characterizing disease pathology and identifying targets for more effective therapeutic strategies.

One of the cell types that could be interesting to study would be the secretory epithelial cells on the surface of the respiratory tract since they primary site of viral attachment, replication, and host defense, they determine the transmission severity and immune response to the virus (Denney & Ho, 2018). 

To investigate these effects at high resolution, single-cell transcriptomic data provides a powerful approach to capture cell-type-specific responses. However, extracting meaningful insights from such data requires the use of appropriate analytical tools and computational frameworks to ensure accurate identification of cell populations, differential responses, and intercellular communication patterns.

## Method comparison
Seurat is widely regarded as a leading framework for single-cell RNA sequencing (scRNA-seq) analysis due to its comprehensive end-to-end workflow. It enables preprocessing steps such as quality control, normalization, scaling, and dimensionality reduction, followed by clustering and visualization of distinct cell populations. In this analysis, Seurat is preferred over alternatives such as Scanpy because the dataset is of moderate size, making computational efficiency less of a limiting factor. Additionally, Seurat’s integration within the R ecosystem facilitates compatibility with downstream tools required for further analysis (Zhang et al., 2024).

A benchmarking study conducted in 2021 comparing single-cell annotation tools found that Seurat, SingleR, CP, SingleCellNet, and RPC performed well overall in terms of annotation accuracy. Notably, SingleR demonstrated increased robustness when distinguishing between highly similar cell types (Huang et al., 2021). This makes SingleR particularly suitable for this study, where subtle transcriptional differences between closely related cell populations are expected.

For differential gene expression analysis, multiple methods have been evaluated across metrics such as false discovery rate (FDR), sensitivity, specificity, and area under the curve (AUC). Recent comparative studies indicate that DESeq2 performs consistently well across these metrics, particularly in maintaining a balance between sensitivity and false discovery control (Li et al., 2026). As such, DESeq2 is selected as a reliable method for identifying differentially expressed genes in this analysis.

Finally, in order to compare Cell-Cell communications between our cells of interest and other cells, CellChat is preferred due to its comprehensive visualization and its ability to identify condition specific interactions across cells (Jin et al., 2024). 

# Methods
## Data Acquisition

The dataset used in this study was obtained from a previously published study available on the NCBI database (PMC11324402). It comprises single-cell RNA sequencing data generated from mouse respiratory tissues, including the respiratory mucosa, olfactory mucosa, and lateral nasal gland. Samples were collected from three mice at multiple time points during influenza infection, enabling the analysis of temporal changes in cellular responses.
Cell calling and initial preprocessing were performed prior to data distribution for this assignment by the course instructor. The processed data was provided as a Seurat object and used directly for downstream analysis within the R environment.

## Quality control
Quality control filtering was performed to remove low-quality and potentially damaged cells prior to downstream analysis. Cells with fewer than 200 detected genes (nFeature_RNA < 200) were excluded to eliminate likely empty droplets or low-complexity libraries. Additionally, cells with a high proportion of mitochondrial gene expression (percent.mt > 20%) were removed, as this is indicative of cellular stress or apoptosis.

This filtering step ensured that only high-quality cells with sufficient transcriptomic information were retained for subsequent normalization, clustering, and downstream analyses.

## Normalization and Scaling  
The filtered dataset was normalized and scaled using the `NormalizeData()` and `ScaleData()` functions from Seurat (v5.4.0), respectively. Log-normalization was applied to account for differences in sequencing depth across cells, ensuring comparability of gene expression values (Hao et al., 2023).

## Clustering  
Cell clustering was performed using the `FindClusters()` function from Seurat (v5.4.0). A resolution of 0.5 was selected to achieve a balance between cluster granularity and interpretability, allowing for the identification of distinct cell populations, including less abundant cell types (Hao et al., 2023).

## Automatic Annotation  
Automated cell type annotation was carried out using SingleR (v2.10.0). Cluster expression profiles were compared against curated mouse reference datasets, including ImmGen and MouseRNAseq, to assign cell identities and improve classification accuracy.

## Manual Annotation  
To validate and refine automated annotations, marker genes for each cluster were identified using the `FindAllMarkers()` function in Seurat. The top differentially expressed genes were cross-referenced with published literature and publicly available databases to infer cell type identities. This step was used to resolve discrepancies and improve annotation reliability (Hao et al., 2023).

## Differential Expression  
Pseudobulk aggregation was performed using the `AggregateExpression()` function, grouping cells by disease status, sample identity, and cell type (cluster). Secretory epithelial cells were selected as the cell type of interest. Differential gene expression analysis between influenza-infected and control samples was conducted using DESeq2 (v1.48.2), which models count data to identify significantly altered genes (Love et al., 2014).

## Overrepresentation Analysis (ORA)  
Gene Ontology (GO) overrepresentation analysis was conducted to identify enriched biological processes and cellular components associated with the differentially expressed genes. This analysis was performed using clusterProfiler (v4.16.0) (Xu et al., 2024). The following parameters were applied:

- **OrgDb** = org.Mm.eg.db (mouse annotation database)  
- **pAdjustMethod** = "BH" (Benjamini–Hochberg false discovery rate correction)  
- **pvalueCutoff** = 0.05  
- **qvalueCutoff** = 0.2  

## Cell–Cell Communication  
Cell–cell communication (CCC) analysis was performed using CellChat (v1.6.1) to investigate signaling interactions between secretory epithelial cells and neighboring cell populations. This approach enabled the identification of potential ligand–receptor interactions and signaling pathways involved in the cellular response to influenza infection (Jin, 2026).

The full R script used for this analysis can be found in [ScRNA-seq script](Seurat.R).

# Results

# Discussion

# References

Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, Chak S, Naikawadi RP,
Wolters PJ, Abate AR, Butte AJ, Bhattacharya M (2019). “Reference-based
analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.” _Nat. Immunol._, *20*, 163-172. doi:10.1038/s41590-018-0276-y
  
Bouvier, N. M., & Palese, P. (2008). THE BIOLOGY OF INFLUENZA VIRUSES. Vaccine, 26(Suppl 4), D49. https://doi.org/10.1016/J.VACCINE.2008.07.039

Denney, L., & Ho, L. P. (2018). The role of respiratory epithelium in host defence against influenza virus infection. Biomedical Journal, 41(4), 218. https://doi.org/10.1016/J.BJ.2018.08.004

Hao et al. Dictionary learning for integrative, multimodal and scalable
  single-cell analysis. Nature Biotechnology (2023) [Seurat V5]

Huang, Q., Liu, Y., Du, Y., & Garmire, L. X. (2021). Evaluation of Cell Type Annotation R Packages on Single-cell RNA-seq Data. Genomics, Proteomics & Bioinformatics, 19(2), 267–281. https://doi.org/10.1016/J.GPB.2020.07.004

Jin, S., Plikus, M. V., & Nie, Q. (2024). CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nature Protocols 2024 20:1, 20(1), 180–219. https://doi.org/10.1038/s41596-024-01045-4

Jin S (2026). _CellChat: Inference and analysis of cell-cell communication from single-cell and spatial transcriptomics data_. R package version 1.6.1


Li, D., Liu, P., Rahman, I., Zand, M., Pryhuber, G., Dye, T., Goniewicz, M., Gurkar, A. U., Königshoff, M., Eickelberg, O., Mora, A., Rojas, M., Ma, Q., Lugo-Martinez, J., Bar-Joseph, Z., Lanna, S., Finkel, T., & Xie, Z. (2026). Evaluation of statistical differential analysis methods for identification of senescent cells using single-cell transcriptomics. Cell Reports Methods, 6(2), 101264. https://doi.org/10.1016/J.CRMETH.2025.101264

Love, M.I., Huber, W., Anders, S.(2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550. 


Short, K. R., Richard, M., Verhagen, J. H., van Riel, D., Schrauwen, E. J. A., van den Brand, J. M. A., Mänz, B., Bodewes, R., & Herfst, S. (2015). One health, multiple challenges: The inter-species transmission of influenza A virus. One Health, 1, 1. https://doi.org/10.1016/J.ONEHLT.2015.03.001

What is Single Cell Transcriptomics? - Single Cell Discoveries. (n.d.). Retrieved April 12, 2026, from https://www.scdiscoveries.com/blog/knowledge/single-cell-transcriptomics/

S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320


Zhang, H., Zhang, W., Zhao, S., Xu, G., Shen, Y., Jiang, F., Qin, A., & Cui, L. (2024). easySCF: a tool for enhancing interoperability between R and Python for efficient single-cell data analysis. Bioinformatics, 40(12). https://doi.org/10.1093/BIOINFORMATICS/BTAE710
