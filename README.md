# FDRnet


## Setup

The setup process for FDRnet requires the following steps:

### Download

Download FDRnet. The following command clones the current FDRnet repository from GitHub:

`git clone https://github.com/yangle293/FDRnet.git`

### Installation

The following software is required for FDRnet.

- Linux/Unix
- [Python (2.7 or 3.5)](www.python.org)
- [NumPy (1.14)](https://www.numpy.org)
- [NetworkX (2.0)](https://networkx.github.io/)
- [Matplotlib (3.1)](https://matplotlib.org/)
- CPLEX (12.7)

#### Installing CPLEX
Academic users can obtain a free but complete version of CPLEX via [IBM Academic Initiative](https://my15.digitalexperience.ibm.com/b73a5759-c6a6-4033-ab6b-d9d4f9a6d65b/dxsites/151914d1-03d2-48fe-97d9-d21166848e65/home).

After the installation of CPLEX, you also need to install the CPLEX-Python modules. Run

`cd yourCPLEXhome/python/VERSION/PLATFORM`

`python setup.py install`

### Testing


## Use

### Input
There are three input files for Hierarchical HotNet that together define a network with scores on the nodes of the network. For example, the following example defines a network with an edge between the nodes ABC and DEF, which have scores 0.5 and 0.2, respectively. For convenience, these files use the same format as the input files for HotNet2.
#### Index-to-gene file
This file associates each gene with an index, which we use for the edge list as well as a similarity matrix:

    1  ABC
    2  DEF
    
##### Edge list file
This file defines a network using the indices in the index-to-gene file:

    1    2
    
##### Gene-to-score file
This file associates each gene with a score:

    ABC 0.5
    DEF 0.2
    
### Running

### Output

## Additional information

### Examples
We provide three files in the `example` directory as an example. 

`python src/FDRnet_main.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl example/BRCA_fdr.txt -ofn example/test.csv -se PSMB3`

We also provide a simple program `plot_subnetworks.py` to visualize the identified subnetwork. The input is the three files in the `example` directory and the output file from above. For example, run the code below:

`python src/plot_subnetworks.py -igi example/irefindex9_index_gene -iel example/irefindex9_edge_list -igl example/BRCA_fdr.txt -ofn example/test.csv`

We can see the subnetwork around PSMB3:
![alt text](https://github.com/yangle293/FDRnet/blob/master/example/seed_PSMB3.png)
### Support
For support with Hierarchical HotNet, please
### License
See `LICENSE.txt` for license information.
### Citation
If you use FDRnet in your work, then please cite the following manuscript:
L. Yang, R. Chen, S. Goodison, Y. Sun. FDRnet:  A Novel Method to Identify SignificantlyMutated Subnetworks in Cancer.
