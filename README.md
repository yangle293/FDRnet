# FDRnet 1.0.0 [![DOI](https://zenodo.org/badge/196633557.svg)](https://zenodo.org/badge/latestdoi/196633557)
FDRnet is a method for identifying significantly mutated (or differentially expressed) subnetworks in human diseases. 

## Setup

The setup process for FDRnet includes the following steps:

### Download

Download FDRnet. The following command clones the current FDRnet repository from GitHub:

`git clone https://github.com/yangle293/FDRnet.git`

### Installation

The following software is required for installing FDRnet:

- Linux/Unix
- [Python (2.7 or 3.5)](www.python.org)
- [NumPy (1.14)](https://www.numpy.org)
- [NetworkX (2.0)](https://networkx.github.io/)
- [Matplotlib (3.1)](https://matplotlib.org/)
- CPLEX (12.7)

#### Installing CPLEX
Academic users can obtain a free, complete version of CPLEX via [IBM Academic Initiative](https://my15.digitalexperience.ibm.com/b73a5759-c6a6-4033-ab6b-d9d4f9a6d65b/dxsites/151914d1-03d2-48fe-97d9-d21166848e65/home).

After the installation of CPLEX, you also need to install the CPLEX-Python modules. Run

    cd yourCPLEXhome/python/VERSION/PLATFORM

    python setup.py install

## Use

### Input
There are three input files, which together define a network used in FDRnet. For example, the following files define a network containing two connected nodes A and B, which have a score of 0.5 and 0.2, respectively. 
#### Index-to-gene file
This file associates each node with an index:

    1  A
    2  B
    
##### Edge list file
This file defines a network by using the indices in the index-to-gene file:

    1    2
    
##### Gene-to-score file
This file associates each gene with a local FDR score:

    A 0.5
    B 0.2
    
If you have a list of p-values for individual genes, you can calculate local FDR scores using either the original R implementation (http://cran.r-project.org/web/packages/locfdr/), or the python-implementation of the locfdr method (https://github.com/leekgroup/locfdr-python). For ease of use, we provide a copy of the python implementation (`locfdr-python`) and a wrapper function `locfdr_compute.py` in `example` directory. 

### Running
1. Compute local FDRs by running `example/locfdr_compute.py` script.

2. Run FDRnet by running `src/FDRnet_main.py` script. 

See the `Examples` section for a full minimal working example of FDRnet.

### Output
The output file is a .csv file organized as follows:

    Seed Gene, Running Time, Optimization Status, Subnetwork
### Usage

    usage: FDRnet_main.py [-h] -igi INPUT_GENE_INDEX -iel INPUT_EDGE_LIST -igl
                      INPUT_GENE_LFDR -ofn OUTPUT_FILE_NAME -se SEED
                      [-bd BOUND] [-al ALPHA] [-sz SIZE]
                      [-tl TIME_LIMIT] [-rg RELATIVE_GAP]
                      
### Program argument

    -h             Show help message and exit
    -igi           File name of the index-to-gene file
    -iel           File name of the input edge list file
    -igl           File name of the gene-to-score file
    -ofn           File name of output, a .csv file
    -se            Seed gene names, either specify a seed gene name (e.g., TP53) OR set 'all' to use all the genes with local FDRs less than FDR bound 
                      as seeds
    -bd            FDR bound, default 0.1
    -al            Random-walk parameter, default 0.85
    -sz            Local graph size, default 400
    -tl            Time limit for each seed for solving MILP problem, default 100
    -rg            Relative gap in a MILP problem, default 0.01
    
### Examples
We provide three files in the `example` directory: an index-to-gene file `irefindex9_index_gene`, an edge list file `irefindex9_edge_list` for the iRefIndex9.0 PPI network, and a gene-to-score file (p-value) `TCGA_BRCA_SNV.txt` generated by MutSig2CV applied to the TCGA breast cancer mutation data. 

We can calculate the local FDRs from p-values by running the following code:

    python example/locfdr_compute.py example/TCGA_BRCA_SNV.txt

There are two ways to use our algorithm. First, users can try to identify a subnetwork around any gene regardless of its local FDR score by using `-se gene_name` option (possibly no solution) to show a local picture of perturbation. Second, users can obtain a landscape of all perturbed regions in a PPI network by using `-se all`. Running with this option will return a set of significantly perturbed subnetworks around all the seeds (i.e., genes with local FDR score less than the given bound B). See the supplemental information of our paper for details. Here, we use the famous cancer gene TP53 as an example, that is, to identify a subnetwork around TP53 gene, by running the following code:

    python src/FDRnet_main.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl example/TCGA_BRCA_SNV_lfdr.txt -ofn example/test.csv -se TP53

The output file look like this one:

| Seed Gene | Running Time  | Optimization Status | Subnetwork|
|:-------:|:-------:|:-----:|:------:|
| TP53	|3.54608988762	| MIP_optimal	|SP3 PTEN IRF9 GPS2 TP53 NCOR1|

We also provide a simple script `plot_subnetworks.py` to visualize the identified subnetwork. The inputs are the three files in the `example` directory and the output file from above. For example, we can visualize the identified subnetwork around TP53 by running the code:

    python src/plot_subnetworks.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl     example/TCGA_BRCA_SNV_lfdr.txt -ofn example/test.csv

The result should be:
![alt text](https://github.com/yangle293/FDRnet/blob/master/example/seed_TP53.png)
## Additional information
### Reproducibility
The `simulation` directory contains the files and codes used to generate synthetic data and plot figures in the paper. The `preprocessing` directory contains the preprocessing scripts used in the breast cancer study and lymphoma study. The code used to plot the results for breast cancer and lymphoma studies is same as that used in the simulation.
### Support
Please send us an email if you have any question when using FDRnet.
### License
See `LICENSE.txt` for license information.
### Citation
If you use FDRnet in your work, please cite the following manuscript:
L. Yang, R. Chen, S. Goodison, Y. Sun. An efficient and effective method to identify significantly perturbed subnetworks in cancer, <em>Nature Computational Science</em>, 1, pages 79–88 (2021).
