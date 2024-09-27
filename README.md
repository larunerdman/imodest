# Purpose
This repo has been created to share the code and details of the modes of gene expression regulation project for the iModEst paper and tool development. 

The modes of regulation project seeks to identify the relationship between different gene expression regulators and gene expression itself. This has been done by fitting elastic net models associating each of a set of regulators with the genes they are proposed to target. Details on how the regulators were chosen and an in-depth description of the methodology can be found in the [iModEst manuscript](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing). 

# Organization of this repo 
This repo contains the scripts used to create the figures and conduct the analyses detailed in the iModEst manuscript \
   (1) **data_preprocessing**: Scripts to preprocess data used to run models and generate PRESS R&#x00B2.  \
   (2) **process_model_results**: Scripts to generate the added the PRESS R&#x00B2; for a given regulator group and gene as well as the average and weighted average coefficients for a given individual regulator and gene.\
   (3) **graphs**: Scripts to reproduce the graphs in the paper with the same color. \
   (4) **manuscript_figure_replication**: Scripts to replicate each figure in the manuscript. \
