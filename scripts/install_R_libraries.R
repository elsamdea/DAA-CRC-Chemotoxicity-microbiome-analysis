# ===== Automated Installation of R Packages for Microbiome Analysis =====
# Run with:
#   Rscript install_R_libraries.R
# or in interactive mode:
#   source("install_R_libraries.R")

# ==================
# ---- Helpers ----
# ==================
options(repos = c(CRAN = "https://cloud.r-project.org"))
NCPUS <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

quietly <- function(expr) try(suppressPackageStartupMessages(suppressWarnings(expr)), silent = TRUE)

# check if installed
is_installed <- function(pkg) {
  is.element(pkg, installed.packages()[, "Package"])
}

# Helper: install from CRAN
install_cran_pkgs <- function(pkgs) {
  pkgs <- unique(pkgs[!vapply(pkgs, is_installed, logical(1))])
  if (length(pkgs)) {
    message("Installing from CRAN: ", paste(pkgs, collapse = ", "))
    install.packages(pkgs, dependencies = TRUE, Ncpus = NCPUS)
  }
}

# Helper: install from Bioconductor
install_bioc_pkgs <- function(pkgs) {
  pkgs <- unique(pkgs[!vapply(pkgs, is_installed, logical(1))])
  if (!length(pkgs)) return(invisible())
  if (!is_installed("BiocManager")) install.packages("BiocManager", Ncpus = NCPUS)
  message("Installing from Bioconductor: ", paste(pkgs, collapse = ", "))
  BiocManager::install(pkgs, update = FALSE, ask = FALSE)
}

# Helper: install from GitHub
install_github_pkgs <- function(repo_vec) {
  if (!length(repo_vec)) return(invisible())
  if (!is_installed("remotes") && !is_installed("devtools")) install.packages("remotes", Ncpus = NCPUS)
  gh_fun <- if (is_installed("remotes")) remotes::install_github else devtools::install_github
  for (repo in repo_vec) {
    pkg <- sub(".*/", "", repo)
    if (!is_installed(pkg)) {
      message("Installing from GitHub: ", repo)
      gh_fun(repo, dependencies = TRUE, upgrade = "never")
    }
  }
}


# ====================
# ---- Libraries ----
# ====================


# CRAN
cran_pkgs <- c(
  "kableExtra","rmarkdown","markdown","readxl", "knitr","devtools","xtable","data.table",
  "reshape2","formatR","tidyr","purrr","plyr","dplyr","ggrepel","janitor","stringr",
  "gtools","ggh4x","patchwork","readr","VennDiagram","RColorBrewer","egg",
  "ggalluvial","reprtree","camcorder","plotly","randomcoloR",
  "colorblindcheck","ggvenn","glue","tidytree","vegan","microbiome","GUniFrac",
  "pheatmap","tibble","viridis","ggraph","igraph","scales", "modeest", "lmerTest", 
  "foreach", "parallel", "ggplot2", "ggrepel", "GUniFrac", "tibble"
)

bio_pkgs <- c("phyloseq", "ALDEx2", "ANCOMBC", "DESeq2", "microbiomeMarker")


github_pkgs <- c("jbisanz/qiime2R","zhouhj1994/LinDA")

# ======================
# Package Installation
# ======================
install_cran(cran_pkgs)
install_bioc(bioc_pkgs)
install_github_safe(github_repos)

# ====================================
# Optional: Missing packages Check
# ====================================

all_pkgs <- c(cran_pkgs, bioc_pkgs, sub(".*/", "", github_repos))
missing <- all_pkgs[!vapply(all_pkgs, is_installed, logical(1))]

if (length(missing)) {
  message("\nSome packages could not be installed (check above for errors): ",
          paste(missing, collapse = ", "))
} else {
  message("\nAll packages installed successfully.")
}

# ====================================
# Optional: Load all packages quietly
# ====================================
quietly(lapply(all_pkgs, require, character.only = TRUE))
message("\nPackages loaded successfully (silent mode).")