# Install bioconductor
BiocManager::valid()
if (!requireNamespace('BiocManager', quietly = TRUE))
  BiocManager::install()

# Install bioconductor package
BiocManager::install("SeuratDisk")

# Install R package with pak
pak::pkg_install("devtools")
