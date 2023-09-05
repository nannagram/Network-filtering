# Network-filtering
Repository for the filtering of complex networks. It includes codes to generate network from Erdós-Renyi and Barabasi-Albert models, the algorithms of four filter and wide range of tool to compare and analyze their performance. 

The repository is divided in 4 branches: 

# BA_networks
This branch containes the code that create BA networks, it is possible to create them undirected, directed, binary and weighted.

# ER_network
This repository is analogous to the previous one but for the ER networks

# Filtering_algorithms
This repository containes the functions that implement 4 different filters: 
- The Hypergeometric filter ( Mantegna et al.)
- The Disparity filter ( Serrano, Bogu ̃n ́a, & Vespignani)
- The Pólya filter ( Marcaccioli & Livan)
- The Max Entropy filter ( Squartini, MAstrandrea & Garlaschelli )

# Filter_analysis
This branch contains some files to analyse the performance of filtering algorithm: 
- comparison_heterogeneity.m and comparison_plot.m are files with metrics implemented to compare the backbones extracted by the filters
- gtest_reshuffling.m gtest3.m and gtest4.m are ground-truth-based tests implemented to test the abilities of filters to detect anomalies
- pvaluesapprox.m and test_hypergeom.m are dedicated to the study of the Hyeprgeometric filter
- plot_net.m and write_map.m are dedicated to the visualisation of the network of the US air traffic 

