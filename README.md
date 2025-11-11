# DAA-CRC-Chemotoxicity-analysis

Comparative evaluation of six commonly used differential abundance analysis methods applied to microbiome data associated with chemotherapy-induced toxicity in colorectal cancer patients.



```
DAA-CRC-Chemotoxicity-analysis/
│
├── data/                     # Folder containing the phyloseq RDS objects
│   ├── derep_ps_16SGlobalpaper.RDS   # Original subset phyloseq object from Conde-Pérez et al. (2024) 
│   ├── tox_ps_final.RDS   # Unfiltered processed phyloseq object
│   └── physeq_filtered.RDS   # Filtered processed phyloseq object
│   └── non-ffpe_metadata_tox_2024.csv  # metadata
│
├── scripts/                  # Folder containing the scripts
│   ├── install_R_libraries.R   # Script for installing required R packages
│   └── DAA_functions.R      # Auxiliary functions used in main analysis pipeline (R file)
│   └── DAA_Chemotoxicity_workflow.Rmd      # Main analysis pipeline (R Markdown file)
│
├── results/                   
│   └── DAA_Benchmarking_CRC_ChemoToxicity_Report.html   # Full report with explanations, code and visualizations
│
└── README.md
```

### Software Requirements

- R version 4.5.1


- R packages: phyloseq, ALDEx2, ANCOMBC, DESeq2, LEfSe (microbiomeMarker), LinDA, ZicoSeq (GUniFrac), dplyr and tidyverse, among others. The RMD file contains session information at the end, specifying all the libraries and packages used.

  
- A reproducible installation script (install_libraries.R) is provided to assist with setup. An additional script with useful functions is also included (scripts/DAA_functions.R).

  

# Reproducibility and Workflow

The analysis pipeline was developed and executed using an R Markdown (.Rmd) workflow. An additional script is included to facilitate the installation of libraries and avoid potential issues. 

To ensure reproducibility, the intermediate `phyloseq` objects used in the analyses will be provided in this repository, including:
- `derep_ps_16SGlobalpaper.RDS`- The initial phyloseq object.
- `tox_ps_final.RDS`- The unfiltered phyloseq object.
- `physeq_filtered.RDS` - The filtered phyloseq object by a prevalence of 5%.

These objects allow users to directly inspect the processed microbiome data and reproduce the main DAA comparisons without re-running the complete preprocessing workflow. For this reason, an HTML report is also provided, containing both the analysis code and the corresponding visualizations.


## Dataset Information

The dataset evaluated in this study corresponds to a subset of a previously analysed cohort from Conde-Pérez et al (2024) [1, 2]. The citations for both studies are provided below. The specific subset used here is described in the Methods section of the manuscript and was selected to focus on chemotherapy-induced toxicity in colorectal cancer patients.

# Manuscript Under Review

This repository contains supplementary code and materials for a manuscript currently under peer review. 
Additional details and citation will be added after the review process is complete.


# References

[1] Conde-Pérez et al. (2024). *The multispecies microbial cluster of Fusobacterium, Parvimonas, Bacteroides, and Faecalibacterium as a precision biomarker for colorectal cancer diagnosis.* Molecular Oncology.  
DOI: [https://doi.org/10.1002/1878-0261.13604](https://doi.org/10.1002/1878-0261.13604)

[2] Conde-Pérez, K., Buetas, E., Aja-Macaya, P., Martin-De Arribas, E., Iglesias-Corrás, I., Trigo-Tasende, N., Nasser-Ali, M., Estévez, L. S., Rumbo-Feal, S., Otero-Alén, B., Noguera, J. F., Concha, Á., Pardiñas-López, S., Carda-Diéguez, M., Gómez-Randulfe, I., Martínez-Lago, N., Ladra, S., Aparicio, L. A., Bou, G., . . .  Poza, M. (2024). *Parvimonas micra can translocate from the subgingival sulcus of the human oral cavity to colorectal adenocarcinoma.* Molecular Oncology, 18(5), 1143-1173. 
DOI: [https://doi.org/10.1002/1878-0261.13506](https://doi.org/10.1002/1878-0261.13506)
