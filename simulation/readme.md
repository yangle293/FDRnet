# Simulation
The files in this directory are used in simulation study in the paper.

## Files
- `[network name]_index_gene` and `[network name]_edge_list` are the gene files and edge files for different networks. 
- `coreComplexes.txt` is the CORUM database download from CORUM website.
- `Cancer_complex_from_CORUM.xlsx` is the file including CORUM protein complexes that are breast cancer-related.
- `script_generate_simulation_data.m` is the code used to generate the synthetic data. It loads the CORUM database and a PPI network, outputs 10 sets of simulated p-values for each of signal levels (from 0.01 to 0.1).
- `script_plot_simulation_results.m` is the code for generating the figures shown in the paper. It plots F-score, F-sub score, FDR bubble graph and running time. 
