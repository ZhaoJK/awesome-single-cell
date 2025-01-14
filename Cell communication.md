

## CellphoneDB {:cellphonedb}
 
Avaliable: <python, R, terminal>  
update v4.1.0 2023-03-09    
[Github source](https://github.com/Teichlab/cellphonedb)    
[Web interface](https://www.cellphonedb.org/)   
[Pipy](https://pypi.org/project/CellphoneDB/)   
 
CellphoneDB is a suite to study cell-cell communication from single-cell transcriptomics data. Identifying ligand–receptor interactions from scRNA-seq requires both the annotation of complex ligand–receptor relationships from the literature (i.e. database) and a statistical method that integrates the resource with scRNA-seq data and selects relevant interactions from the dataset (ie.e tool). CellphoneDB is composed of two units, a database and a tool.  

CellphoneDB database is a publicly available repository of curated receptors, ligands and their interactions. Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies. It integrates new manually reviewed interactions with evidenced roles in cell-cell communication together with existing datasets that pertain to cellular communication (such as *Shilts et al. 2022* and *Kanemaru et al. 2023*). Recently, the database expanded to include non-protein molecules acting as ligands.  

The database can be used to search for a particular ligand/receptor or in combination with the tool to interrogate your own single-cell transcriptomics data.  

## Liana  {:liana}  

Avaliable: <python, R>  
[Github source](https://github.com/saezlab/liana-py)   
   
LIANA is a Ligand-Receptor inference framework that enables the use of any LR method with any resource. This is its faster and memory efficient Python implementation, an R version is also available [here](https://github.com/saezlab/liana).  
The methods implemented in this repository are: CellPhoneDBv2, NATMI, Connectome, SingleCellSignalR, CellChat (+), 1-vs-rest expression LogFC score, Geometric Mean - ligand-receptor geometric mean with pvalues obtained via the permutation approach implemented by CellPhoneDBv2, rank_aggregate of the predictions calculated with the RobustRankAggregate method.  
  
  
   
