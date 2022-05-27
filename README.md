# Scaling-issues-of-neutral-theory-reveal-violations-of-ecological-equivalence-for-dominant-Amazonian-
R_scripts accompanying our Ecology Letters paper titled: "Scaling issues of neutral theory reveal violations of ecological equivalence for dominant Amazonian tree species". 

Repository includes all the necessary files and scripts to recreate figures and data from the main manuscript. Scripts are based on an example dataset from Guyana/Suriname with details as specified in the main manuscript. A list of packages used can be found in the main R script. 

When using this script ---> cite original paper including the DOI

Directory includes the following files:

- Script   ----------------------> Script to run simulation in parallel, including example dataset and a list of packages used.
- Script to create all figures --> Script to create all figures (both main manuscript and supplementary material)
- Functions_3Dmodel.R -----------> Script with necessary functions for both above mentioned scripts

- species_data.csv --> an example dataset based on Guyana/Suriname as detailed in the manuscript
- plotdata.csv ------> the plot metadata belonging to species_data.csv

- Forest_OSB_8000ha_m0_1185.rds --> resulting .rds file from simulation using migration parameter settings as explained in the main                                           manuscript to recreate figures
- Forest_OSB_m_uni.rds -----------> resulting .rds file from simulation using migration parameter near unity settings
- Forest_OSB_mnull.rds -----------> resulting .rds file from simulation using migration parameter near null settings 
