

# Pseudotime and trajectory inference  

- [CellRank](https://github.com/theislab/cellrank/tree/main)  
  [Python]  
  Estimate differentiation direction based on a varied number of biological priors, including RNA velocity (La Manno et al. (2018), Bergen et al. (2020)), any pseudotime or developmental potential, experimental time points, metabolic labels, and more.  
  Compute initial, terminal and intermediate macrostates.  
  Infer fate probabilities and driver genes.  
  Visualize and cluster gene expression trends.  
  Documentation: https://cellrank.readthedocs.io/  
  https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/kernels/200_rna_velocity.html
  
- [wot](https://github.com/broadinstitute/wot)  
  [Python]  
  wot (Waddington-OT) is a tool to analyze Optimal-Transport of single-cell expression in RNA-seq data sets. The wot algorithm applies time-course data to compute the probability distributions of cells.
  Documentation: https://broadinstitute.github.io/wot/
  Home page: https://github.com/broadinstitute/wot  
  Data: https://singlecell.broadinstitute.org/single_cell
  Publications: Schiebinger G, Shu J, Tabaka M, Cleary B, Subramanian V, Solomon A, Gould J, Liu S, Lin S, Berube P, Lee L, Chen J, Brumbaugh J, Rigollet P, Hochedlinger K, Jaenisch R, Regev A, Lander ES "Optimal-Transport Analysis of Single-Cell Gene Expression Identifies Developmental Trajectories in Reprogramming." Cell. **2019** Feb 7;176(4):928-943.e22. Jan 31. https://doi.org/10.1016/j.cell.2019.01.006

- [CytoTRACE](https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz)
  [R]  
  Home: https://cytotrace.stanford.edu/  

- [Palantir](https://github.com/dpeerlab/Palantir)  
  [Python](https://www.nature.com/articles/s41587-019-0068-4)
  A tool for computation of cellular differentiation trajectories in single-cell RNA-seq and Mass cytometry data. The Palantir algorithm models differentiation as a stochastic process to estimate cell fates.  
  Documentation: https://nbviewer.jupyter.org/github/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
  Publications: Setty M, Kiseliovas V, Levine J, Gayoso A, Mazutis L, Pe'er D "Characterization of cell fate probabilities in single-cell data with Palantir." Nat Biotechnol. **2019** Apr;37(4):451-460. Mar 21. https://doi.org/10.1038/s41587-019-0068-4
  PMID: 30899105

- [PhenoPath](https://github.com/kieranrcampbell/phenopath)  
  [R] v1.26  
  Single-cell pseudotime with heterogeneous genetic and environmental backgrounds, including Bayesian significance testing of iteractions.
  phenopath is an R tool to estimate pseud-time trajectories in single-cell RNA-seq data. The phenopath algorithm uses a new statistical approach, a mixture regression-latent variable model, which can detect known and unique covariate-pseudo-time interaction influences.
  Publications: Campbell KR, Yau C "Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data." Nat Commun. **2018** Jun 22;9(1):2442. https://doi.org/10.1038/s41467-018-04696-6
  PMID: 29934517
  PMCID: PMC6015076
  Documentation: https://bioconductor.org/packages/release/bioc/vignettes/phenopath/inst/doc/introduction_to_phenopath.html


- [scVelo](https://github.com/theislab/scvelo)   
  [Python]  
  scVelo is a scalable toolkit for RNA velocity analysis in single cells. It generalizes the concept of RNA velocity by relaxing previously made assumptions with a dynamical model. It allows to identify putative driver genes, infer a latent time, estimate reaction rates of transcription, splicing and degradation, and detect competing kinetics.

- [slingshot](https://github.com/kstreet13/slingshot)  
  [R]  
  Functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data.

- [VELOCYTO](http://velocyto.org/)  
  [Python, R]  
  Estimating RNA velocity in single cell RNA sequencing datasets.

- [CytoRuter](https://github.com/edroaldo/cellrouter)
  [R]
  CellRouter is a tool for the identification of cell-state transition trajectories in single-cell RNA-seq data. The CellRouter algorithm uses flow networks to identify subpopulations.
  Publications:  Lummertz da Rocha E, Rowe RG, Lundin V, Malleshaiah M, Jha DK, Rambo CR, Li H, North TE, Collins JJ, Daley GQ "Reconstruction of complex single-cell trajectories using CellRouter." Nat Commun. 2018 Mar 1;9(1):892. https://doi.org/10.1038/s41467-018-03214-y
  PMID: 29497036
  PMCID: PMC5832860

- [CALISTA](https://github.com/CABSEL/CALISTA) - [R] - CALISTA provides a user-friendly toolbox for the analysis of single cell expression data. CALISTA accomplishes three major tasks: 1)	Identification of cell clusters in a cell population based on single-cell gene expression data, 2)	Reconstruction of lineage progression and produce transition genes, and 3)	Pseudotemporal ordering of cells along any given developmental paths in the lineage progression.
- [CoSpar](https://cospar.readthedocs.io/) - [python] - CoSpar is a toolkit for dynamic inference by integrating state and lineage information. It gains superior robustness and accuracy by exploiting both the local coherence and sparsity of differentiation transitions, i.e., neighboring initial states share similar yet sparse fate outcomes. When only state information is available, CoSpar also improves upon existing dynamic inference methods by imposing sparsity and coherence.
- [DensityPath](https://doi.org/10.1101/276311) - [.] - DensityPath: a level-set algorithm to visualize and reconstruct cell developmental trajectories for large-scale single-cell RNAseq data
- [dynverse](https://github.com/dynverse/dynverse/) - [R] - [A comparison of single-cell trajectory inference methods: towards more accurate and robust tools](https://www.nature.com/articles/s41587-019-0071-9)
- [ECLAIR](https://github.com/GGiecold/ECLAIR) - [python] - ECLAIR stands for Ensemble Clustering for Lineage Analysis, Inference and Robustness. Robust and scalable inference of cell lineages from gene expression data.
- [K-Branches](https://github.com/theislab/kbranches) - [R] - The main idea behind the K-Branches method is to identify regions of interest (branching regions and tips) in differentiation trajectories of single cells. So far, K-Branches is intended to be used on the diffusion map representation of the data, so the user should either provide the data in diffusion map space or use the destiny package perform diffusion map dimensionality reduction.
- [MERLoT](https://github.com/soedinglab/merlot) - [R/python] - Reconstructing complex lineage trees from scRNA-seq data using MERLoT.
- [ouija](https://github.com/kieranrcampbell/ouija) - [R] - [A descriptive marker gene approach to single-cell pseudotime inference](https://doi.org/10.1101/060442)
- [ouijaflow](http://www.github.com/kieranrcampbell/ouijaflow) - [python] - [A descriptive marker gene approach to single-cell pseudotime inference](https://doi.org/10.1101/060442)
- [pseudodynamics](https://github.com/theislab/pseudodynamics) - [MATLAB] - [Inferring population dynamics from single-cell RNA-sequencing time series data](https://www.nature.com/articles/s41587-019-0088-0)
- [psupertime](https://github.com/wmacnair/psupertime) - [R] - psupertime is an R package which uses single cell RNAseq data, where the cells have labels following a known sequence (e.g. a time series), to identify a small number of genes which place cells in that known order. It can be used for discovery of relevant genes, for exploration of unlabelled data, and assessment of one dataset with respect to the labels known for another dataset. - [preprint](https://www.biorxiv.org/content/10.1101/622001v1)
- [SCDIFF](https://github.com/phoenixding/scdiff) - [Python, JavaScript] - SCDIFF is a single-cell trajectory inference method with interactive visualizations powered by D3.js. SCDIFF utilized the TF regulatory information to mitigate the impact of enormous single-cell RNA-seq noise (such as drop-out). With the TF regulatory information, SCDIFF is also able to predict the TFs (and their activation time), which drive the cells to different cell fates. Such predictive power has been [experimentally validated](https://genome.cshlp.org/content/28/3/383).
- [SCIMITAR](https://github.com/dimenwarper/scimitar) - [Python] - Single Cell Inference of Morphing Trajectories and their Associated Regulation module (SCIMITAR) is a method for inferring biological properties from a pseudotemporal ordering. It can also be used to obtain progression-associated genes that vary along the trajectory, and genes that change their correlation structure over the trajectory; progression co-associated genes.
- [SCORPIUS](https://cran.r-project.org/package=SCORPIUS) - [R] - An accurate and easy tool for performing linear trajectory inference on single cells using single-cell RNA sequencing data. In addition, SCORPIUS provides functions for discovering the most important genes with respect to the reconstructed trajectory, as well as nice visualisation tools. Cannoodt et al. (2016) [doi:10.1101/079509](https://doi.org/10.1101/079509).
- [SCUBA](https://github.com/gcyuan/SCUBA) - [matlab/R] - SCUBA stands for "Single-cell Clustering Using Bifurcation Analysis." SCUBA is a novel computational method for extracting lineage relationships from single-cell gene expression data, and modeling the dynamic changes associated with cell differentiation.
- [SLICER](https://github.com/jw156605/SLICER) - [R] - Selective Locally linear Inference of Cellular Expression Relationships (SLICER) algorithm for inferring cell trajectories.
- [SPADE](http://www.nature.com/nprot/journal/v11/n7/full/nprot.2016.066.html) - [R] - Visualization and cellular hierarchy inference of single-cell data using SPADE.
- [TASIC](https://www.andrew.cmu.edu/user/sabrinar/TASIC) - [matlab] - TASIC is a new method for determining temporal trajectories, branching and cell assignments in single cell time series experiments. Unlike prior approaches TASIC uses on a probabilistic graphical model to integrate expression and time information making it more robust to noise and stochastic variations.
- [TopSLAM](https://github.com/mzwiessele/topslam) - [python] - Extracting and using probabilistic Waddington's landscape recreation from single cell gene expression measurements
- [TSCAN](https://github.com/zji90/TSCAN) - [R] - Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
