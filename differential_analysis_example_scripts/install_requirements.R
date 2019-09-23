if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Check that you have write permissions to these directories:
installed.packages()[, c("Package", "LibPath")]
# e.g., I had to:
# $ sudo chown -R myuser:myuser /usr/lib64/R/
# $ sudo chown -R myuser:myuser /usr/share/doc/R

BiocManager::install("limma")
# (updated all packages in list of requirements when prompted)
BiocManager::install("affy")
# "annotate" requires Rcurl, which requires libcurl at the OS level:
# $ dnf install libcurl-devel
# "annotate" requires XML, which requires libxml at the OS level:
# $ dnf install libxml-devel
BiocManager::install("annotate")
BiocManager::install("hgu133plus2.db")
BiocManager::install("rat2302.db")
BiocManager::install("stringr")
