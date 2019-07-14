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
- NumPy (1.14)
- Scipy (0.19)
- NetworkX (2.0)
- CPLEX (12.7)

### Testing


## Use

### Input

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

### License

### Citation
