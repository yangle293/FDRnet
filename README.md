# FDRnet
FDRnet is a method for identifying significantly mutated (or differentially expressed) subnetworks in human diseases. 

## Setup

The setup process for FDRnet requires the following steps:

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
Academic users can obtain a free but complete version of CPLEX via [IBM Academic Initiative](https://my15.digitalexperience.ibm.com/b73a5759-c6a6-4033-ab6b-d9d4f9a6d65b/dxsites/151914d1-03d2-48fe-97d9-d21166848e65/home).

After the installation of CPLEX, you also need to install the CPLEX-Python modules. Run

    cd yourCPLEXhome/python/VERSION/PLATFORM

    python setup.py install

## Use

### Input
There are three input files for FDRnet that together define a network with scores on the nodes of the network. For example, the following example defines a network with an edge between the nodes A and B, which have scores 0.5 and 0.2, respectively. 
#### Index-to-gene file
This file associates each gene with an index, which we use for an edge list as well as a similarity matrix:

    1  A
    2  B
    
##### Edge list file
This file defines a network using the indices in an index-to-gene file:

    1    2
    
##### Gene-to-score file
This file associates each gene with a local FDR score:

    A 0.5
    B 0.2
    
If you have a p-value for each gene, we provide a wrapping function to calculate the local FDRs in `locfdr` directory using the [`locfdr`](https://github.com/leekgroup/locfdr-python) package.
### Running
1. Compute local FDRs by running `locfdr/locfdr_compute.py` script.

2. Run FDRnet by running `src/FDRnet_main.py` script. 

See the `Examples` section for a full minimal working example of FDRnet.
### Output
The output file is a .csv file organized as follows:

    Seed Gene, Running Time, Optimization Status, Subnetwork
### Usage

    usage: FDRnet_main.py [-h] -igi INPUT_GENE_INDEX -iel INPUT_EDGE_LIST -igl
                      INPUT_GENE_LFDR -ofn OUTPUT_FILE_NAME [-se SEED]
                      [-bd BOUND] [-al ALPHA] [-sz SIZE]
                      [-tl TIME_LIMIT] [-rg RELATIVE_GAP]
                      
### Optional argument

    -h             Show help message and exit
    -igi           File name of the input gene index file
    -iel           File name of the input edge list file
    -igl           File name of the input local FDRs file
    -ofn           File name of output, a .csv file
    -se            Seed gene names, 'all' for setting all the genes with local FDRs less than FDR bound 
                      as seeds OR a specified seed gene name (e.g. TP53)
    -bd            FDR bound, default 0.1
    -al            Random-walk parameter, default 0.85
    -sz            Local graph size, default 400
    -tl            Time limit for each seed in MILP problem, default 100
    -rg            Relative gap in MILP problem, default 0.01
    
### Examples
We provide three files in the `example` directory: an index-to-gene file `irefindex9_index_gene`, an edge list file `irefindex9_edge_list` of the iRefIndex9.0 PPI network, and a gene-to-score file `BRCA_fdr.txt` calculated from TCGA breast cancer data. 

We can identify a subnetwork around PSMB3 gene by running the following code:

    python src/FDRnet_main.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl example/BRCA_fdr.txt -ofn example/test.csv -se PSMB3

The output file should be:

| Seed Gene | Running Time  | Optimization Status | Subnetwork|
|:-------:|:-------:|:-----:|:------:|
| PSMB3	|4.573695183	| MIP_optimal	|PSMD11 PSMD12 GFPT2 CDC6 PSMD3 PSMC4 PSMC5 VDAC3 PSMA7 PSMB4 PSMB3|

We also provide a simple script `plot_subnetworks.py` to visualize the identified subnetwork. The input is the three files in the `example` directory and the output file from above. For example, we can visualize the subnetwork around PSMB3 by running the code:

    python src/plot_subnetworks.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl     example/BRCA_fdr.txt -ofn example/test.csv

The result should be:
![alt text](https://github.com/yangle293/FDRnet/blob/master/example/seed_PSMB3.png)
## Additional information
### Support
Please send us an email if you meet any problem in using FDRnet.
### License
See `LICENSE.txt` for license information.
### Citation
If you use FDRnet in your work, please cite the following manuscript:
L. Yang, R. Chen, S. Goodison, Y. Sun. FDRnet:  A novel method to identify significantly mutated subnetworks in cancer.
