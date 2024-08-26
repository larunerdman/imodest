# Purpose
This repo has been created to share the code and details of the modes of gene expression regulation project for the iModEst paper and tool development. 

The modes of regulation project seeks to identify the relationship between different gene expression regulators and gene expression itself. This has been done by fitting elastic net models associating each of a set of regulators with the genes they are proposed to target. Details on how the regulators were chosen and an in-depth description of the methodology can be found in the [iModEst manuscript draft](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing). 

# Organization of this repo 
This repo contains the scripts used to create the figures and conduct the analyses detailed in the [iModEst manuscript draft](https://docs.google.com/document/d/1Cr4qanLEvsm4gJR2_z3cMIneZFDzSYjcFGg4n0TVO1s/edit?usp=sharing). \
   (1) **generate_model_results**: Scripts to generate the PRESS R&#x00B2; and model coefficients. These scripts reproduce the fundamental analyses of our project and can be used for reference but don't need to be run (they take a long time since they implementent a leave-one-out procedure on all of our data). \
   (2) **process_model_results**: Scripts to generate the added the PRESS R&#x00B2; for a given regulator group and gene as well as the average and weighted average coefficients for a given individual regulator and gene.\
   (3) **graphs**: Scripts to reproduce the graphs in the paper with the same color. \
   (4) **annotate_regs**: Scripts to annotate regulators. \
   (5) To potentially be added: preprocessing scripts for reference
