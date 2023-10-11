# scGSVA: GSVA for single cell RNA seq analysis    
[Github](https://github.com/guokai8/scGSVA)   
   
scGSVA provides wrap functions to do GSVA analysis for single cell data. And scGSVA includes functions to build annotation for almost all species. scGSVA also provides function to generate figures based on the GSVA results.  
scGSVA provides functions to generate annotation data which can be used in the analysis.    
  
# GSVA: gene set variation analysis for microarray and RNA-seq data   
[Github](https://github.com/rcastelo/GSVA)   
The GSVA package allows one to perform a change in coordinate systems of molecular measurements, transforming the data from a gene by sample matrix to a gene-set by sample matrix, thereby allowing the evaluation of pathway enrichment for each sample. This new matrix of GSVA enrichment scores facilitates applying standard analytical methods such as functional enrichment, survival analysis, clustering, CNV-pathway analysis or cross-tissue pathway analysis, in a pathway-centric manner. For citing GSVA as a software package, please use the following reference:    
  
# VISION
[Github](https://github.com/YosefLab/VISION)
VISION aids in the interpretation of single-cell RNA-seq (scRNA-seq) data by selecting for gene signatures which describe coordinated variation between cells. While the software only requires an expression matrix and a signature library (available in online databases), it is also designed to integrate into existing scRNA-seq analysis pipelines by taking advantage of precomputed dimensionality reductions, trajectory inferences or clustering results. The results of this analysis are made available through a dynamic web-app which can be shared with collaborators without requiring them to install any additional software.  
  
# PAGODA2  
[Github](https://github.com/kharchenkolab/pagoda2)  
Pagoda2 is an R package for analyzing and interactively exploring large-scale single-cell RNA-seq datasets. The methods were optimized to rapidly process modern scRNAseq datasets, which are both large (approximately 1e6 cells or greater) and sparse. The package provides methods for quality control, filtering, clustering, visualization, differential expression, cross-cutting aspects/states, and geneset/pathway overdispersion analysis. The companion frontend application allows users to figure out which gene expression patterns give rise to different subpopulations within the data. The application allows users to inspect the gene expression patterns of subpopulations through annotated gene sets and pathways, including Gene Ontology (GO) categories. Users may also highlight certain clusters and perform differential expression from their browsers via the frontend application.  

# pathview  
[Github](https://github.com/datapplab/pathview)  
Pathview is a leading tool for pathway based data integration and visualization. It maps, integrates and renders a wide variety of biological data on relevant pathway graphs. Pathview has 3 important features:
 + Interpretable graphs with publication quality: KEGG view for easy interpretation and Graphviz view for better graphical control.
 + Strong data integration capacity. It works with: 1) all data mappable to pathways, 2) 30 of molecular ID types (genes/protein, compound/metabolite), 3) 5000+ species, 4) various data attributes and formats.
 + Simple and powerful: fully automated and error-resistant, seamlessly integrates with a wide range of pathway and gene set (enrichment) analysis tools.

# SBGNview  
[Github](https://github.com/datapplab/SBGNview)  
Here we introduce the SBGNview package, which adopts Systems Biology Graphical Notation (SBGN) and greatly extends the Pathview project by supporting multiple major pathway databases beyond KEGG.  
Key features:  
 + Pathway diagram and definition by the widely adopted SBGN standard formats;  
 + Supports multiple major pathway databases beyond KEGG (Reactome, MetaCyc, SMPDB, PANTHER, METACROP etc) and user defined pathways;  
 + Covers 5,200 reference pathways and over 3,000 species by default;  
 + Extensive graphics controls, including glyph and edge attributes, graph layout and sub-pathway highlight;  
 + SBGN pathway data manipulation, processing, extraction and analysis.

Supplementary Table 3. Comparison between SBGNview and other publicly available pathway based data visualization tools.
|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
|Tool|Description|Pathway Supported|Graph Style|*Full Graphic Output|Data Download|#Data Mapping|Data Integration|Automated Analysis|Interface|OS|
|SBGNview|Data mapping, integration, visualization and analysis|Major pathway databases (5200 reference pathways, 3000 species) and custom SBGN pathways|SBGN|Yes|Yes|Any data mappable|Strong|Yes|R/Bioconductor|Multiple|
|Pathview|Data mapping, integration and visualization|KEGG (all pathways, species, KO)|KEGG, Graphviz|Yes|Yes|Any data mappable|Strong|Yes|R/Bioconductor|Multiple|
|ChiBE|Pathway editing, visualization|BioPAX,  basic SIF format|Compound graphs, SBGN|No|Yes|Gene and compound data|Limited|No|Java GUI|Multiple|
|Cytoscape and plugins|Network integration, analysis, and visualization|Pathway and molecular interaction|Various layout algorithms|Yes|Yes|Molecular data|Strong|Limited|Java platform and GUI|Multiple|
|Escher|visualize data on biological pathways|Escher pathways (12 pathways for E coli, Human, Budding yeast)|Escher JSON|Yes|No|No|No|No|Web GUI|Multiple|
|g-language|Data mapping|KEGG (predefined subset)|KEGG|No color key|No|Gene and compound data|Limited|No|Web GUI|Multiple|
|KEGG Mapper|Search and color mapped pathway components|KEGG (only human and reference pathways)|KEGG|No|Yes|No|No|No|Web GUI|Multiple|
|MetaboMAPS|multi-omics data visualization in metabolic pathways|Metabolic pathways in SVG format (25 pathways, 16 species)|MetaboMAPS SVG|Yes|No|No|No|No|Web GUI|Multiple|
|PaintOmics|integrative visualization of multiple omic datasets|KEGG (~70 species)|KEGG|Yes|Yes|Gene and metabolite data|Limited|No|Web GUI|Multiple|
|Pathway Tools|Genome data management, systems biology, and omics data analysis|BioCyc|BioCyc (BioPAX, SBML)|Yes|Yes|Gene and metabolite data|Limited|No|Standalone or Web GUI|Multiple|
|Reactome Analysis Tools|Map and visualize data|Reactome|Reactome|Yes|Yes|Gene and compound data|Limited|No|Web GUI|Multiple|
|VANTED|Extendable network visualisation and analysis tool|KEGG, MetaCrop, RIMAS |Common network exchange formats including SBML, BioPAX, KGML, SBGN etc (Add-on required).|Yes|Yes|Experimental data|Limited|No|Java GUI|Multiple|
|VisANT|Visualization, editing, prediction and construction|KEGG|Compound graphs|No|Yes|Gene data|No|Yes|Java GUI and command-line|Multiple|
