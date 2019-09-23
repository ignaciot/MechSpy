library(hgu133plus2.db)

x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
write.table(xx, file="./hgu133plus2_all_gene_symbols.txt", row.names = FALSE, col.names = FALSE, sep="\n")
