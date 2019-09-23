# To automate this:
# - Pass targetsdir as variable
# - Read the chip type directly from the Attributes file (and install the library as needed)
# - Come up with the control/treatment variable names from the Attributes file

library(limma)
library(affy)
library(annotate)
library(hgu133plus2.db)
library(stringr)

significance_threshold = 0.05

args = commandArgs(trailingOnly=TRUE)

#targetsdir = "time/series/dir/containing/tsv/descriptions/and/celfiles/subdir"
targetsdir = args[1]
#outputdir = "./microarray"
outputdir = args[2]
#metadata = "cyclosporin-A-20uM.tsv"
metadata = args[3]
attr_data = readLines(paste(targetsdir, metadata, sep="/"))
Sys.setlocale(locale="C")

# can parse from any of the rows
chem_short_name = unlist(strsplit(attr_data[15], "\t"))[3]
chem_long_name = unlist(strsplit(attr_data[15], "\t"))[3]
chem_short_name = gsub("[^0-9a-zA-Z]", "-", chem_short_name)
chem_long_name = gsub("[^0-9a-zA-Z]", "-", chem_long_name)

# parse the concentrations used, always in microMolar units
chem_concentr = unlist(strsplit(attr_data[15], "\t"))[6]
chem_concentr_unit = unlist(strsplit(attr_data[15], "\t"))[7]
array_design = "hgu133plus2"

datadir = paste(targetsdir, "celfiles", sep="/")
targets <- readTargets(metadata, path=targetsdir, sep="\t", row.names="Source Name")
# massage the barcode strings before passing them to ReadAffy to add the leading 0s and the .CEL extension
ab <- ReadAffy(filenames=paste(gsub(" ", "_", targets$Subject.ID), "CEL", sep="."), celfile.path=datadir)
eset <- rma(ab)
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, paste(array_design, "db", sep="."))
fData(eset) <- data.frame(Symbol=Symbol)
treatments <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8),
                     labels=c("ctrl_12h", "ctrl_24h", "ctrl_48h", "ctrl_72h", 
                              "treatment_12h", "treatment_24h", "treatment_48h", "treatment_72h"))
contrasts(treatments) <- cbind(Time=c(0,1,2,3,0,1,2,3),
                               treatment_12h=c(0,0,0,0,1,0,0,0),
                               treatment_24h=c(0,0,0,0,0,1,0,0),
                               treatment_48h=c(0,0,0,0,0,0,1,0),
                               treatment_72h=c(0,0,0,0,0,0,0,1))
design <- model.matrix(~treatments)
colnames(design) <- c("Intercept","Time",
                      "treatment_12h","treatment_24h","treatment_48h","treatment_72h", 
                      "treatments","treatments")
fit <- lmFit(eset,design)

# contrast control 12hr with treatment 12hr
cont.matrix <- cbind(ctrl2=c(1,0,0,0,0,0,0,0),treat2=c(0,0,0,0,1,0,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl2=results[,1],treat2=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat2",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
# the output file format will be CHEMNAME_CONCENTRATION-uM_TIMEPOINT_10k_genes.txt" and will contain
# up to the top 10k genes with a significance less than significance_threshold
write.table(diff_exp_brief, 
            file=paste(outputdir, 
                       paste(chem_long_name, 
                             paste(chem_concentr, 
                                   paste(chem_concentr_unit, "12h", "liver-HepG2_top_10k_genes.txt", sep="_"), 
                                   sep="-"), 
                             sep="_"), 
                       sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 24hr with treatment 24hr
cont.matrix <- cbind(ctrl2=c(0,1,0,0,0,0,0,0),treat2=c(0,0,0,0,0,1,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl2=results[,1],treat2=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat2",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
# the output file format will be CHEMNAME_CONCENTRATION-uM_TIMEPOINT_10k_genes.txt" and will contain
# up to the top 10k genes with a significance less than significance_threshold
write.table(diff_exp_brief, 
            file=paste(outputdir, 
                       paste(chem_long_name, 
                             paste(chem_concentr, 
                                   paste(chem_concentr_unit, "24h", "liver-HepG2_top_10k_genes.txt", sep="_"), 
                                   sep="-"), 
                             sep="_"), 
                       sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 48hr with treatment 48hr
cont.matrix <- cbind(ctrl2=c(0,0,1,0,0,0,0,0),treat2=c(0,0,0,0,0,0,1,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl2=results[,1],treat2=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat2",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
# the output file format will be CHEMNAME_CONCENTRATION-uM_TIMEPOINT_10k_genes.txt" and will contain
# up to the top 10k genes with a significance less than significance_threshold
write.table(diff_exp_brief, 
            file=paste(outputdir, 
                       paste(chem_long_name, 
                             paste(chem_concentr, 
                                   paste(chem_concentr_unit, "48h", "liver-HepG2_top_10k_genes.txt", sep="_"), 
                                   sep="-"), 
                             sep="_"), 
                       sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 72hr with treatment 72hr
cont.matrix <- cbind(ctrl2=c(0,0,0,1,0,0,0,0),treat2=c(0,0,0,0,0,0,0,1))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl2=results[,1],treat2=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat2",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
# the output file format will be CHEMNAME_CONCENTRATION-uM_TIMEPOINT_10k_genes.txt" and will contain
# up to the top 10k genes with a significance less than significance_threshold
write.table(diff_exp_brief, 
            file=paste(outputdir, 
                       paste(chem_long_name, 
                             paste(chem_concentr, 
                                   paste(chem_concentr_unit, "72h", "liver-HepG2_top_10k_genes.txt", sep="_"), 
                                   sep="-"), 
                             sep="_"), 
                       sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

sessionInfo()

