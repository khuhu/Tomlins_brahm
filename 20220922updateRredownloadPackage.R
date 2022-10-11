lib_loc <- "/home/kevhu/R/x86_64-pc-linux-gnu-library/4.1/"
to_install <- unname(installed.packages(lib.loc = lib_loc)[, "Package"])
to_install
install.packages(pkgs = to_install)

### needed for devtools
###  sudo apt-get install -y libfribidi-dev 
install.packages("devtools")

BiocManager::install(c("GenomicRanges"))
devtools::install_github("dami82/mutSignatures", force = TRUE, build_vignettes = TRUE)
BiocManager::install("preprocessCore")
devtools::install_github("andrewparkermorgan/argyle", dependencies = TRUE)

BiocManager::install("BSgenome")
