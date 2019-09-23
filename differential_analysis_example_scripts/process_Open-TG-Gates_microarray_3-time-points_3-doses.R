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
attr_data = readLines(paste(targetsdir, "Attribute.tsv", sep="/"))
Sys.setlocale(locale="C")

# can parse from any of the rows
chem_short_name = unlist(strsplit(attr_data[2], "\t"))[9]
chem_long_name = unlist(strsplit(attr_data[2], "\t"))[8]
chem_long_name = gsub(" ", "-", chem_long_name)

# parse the concentrations used, always in microMolar units
chem_concentr_low = unlist(strsplit(attr_data[8], "\t"))[19]
chem_concentr_mid = unlist(strsplit(attr_data[14], "\t"))[19]
chem_concentr_high = unlist(strsplit(attr_data[20], "\t"))[19]
array_design = gsub("[^[:alnum:] ]", "", tolower(unlist(strsplit(attr_data[20], "\t"))[2]))

datadir = paste(targetsdir, "celfiles", sep="/")
targets <- readTargets("Attribute.tsv", path=targetsdir, sep="\t", row.names="BARCODE")
# massage the barcode strings before passing them to ReadAffy to add the leading 0s and the .CEL extension
ab <- ReadAffy(filenames=paste(str_pad(targets$BARCODE, 12, pad="0"), "CEL", sep="."), celfile.path=datadir)
eset <- rma(ab)
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, paste(array_design, "db", sep="."))
fData(eset) <- data.frame(Symbol=Symbol)
treatments <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12),
                     labels=c("ctrl_2h", "ctrl_8h", "ctrl_24h", 
                              "low_dose_2h", "low_dose_8h", "low_dose_24h", 
                              "mid_dose_2h", "mid_dose_8h", "mid_dose_24h", 
                              "high_dose_2h", "high_dose_8h", "high_dose_24h"))
contrasts(treatments) <- cbind(Time=c(0,1,2,0,1,2,0,1,2,0,1,2),
                               low_dose_2h=c(0,0,0,1,0,0,0,0,0,0,0,0),
                               low_dose_8h=c(0,0,0,0,1,0,0,0,0,0,0,0),
                               low_dose_24h=c(0,0,0,0,0,1,0,0,0,0,0,0),
                               mid_dose_2h=c(0,0,0,0,0,0,1,0,0,0,0,0),
                               mid_dose_8h=c(0,0,0,0,0,0,0,1,0,0,0,0),
                               mid_dose_24h=c(0,0,0,0,0,0,0,0,1,0,0,0),
                               high_dose_2h=c(0,0,0,0,0,0,0,0,0,1,0,0),
                               high_dose_8h=c(0,0,0,0,0,0,0,0,0,0,1,0),
                               high_dose_24h=c(0,0,0,0,0,0,0,0,0,0,0,1))
design <- model.matrix(~treatments)
colnames(design) <- c("Intercept","Time","low_dose_2h","low_dose_8h","low_dose_24h",
                      "mid_dose_2h","mid_dose_8h","mid_dose_24h",
                      "high_dose_2h","high_dose_8h","high_dose_24h", 
                      "treatments")
fit <- lmFit(eset,design)

# contrast control 2hr with treatment 2hr (low concentration)
cont.matrix <- cbind(ctrl2=c(1,0,0,0,0,0,0,0,0,0,0,0),treat2=c(0,0,0,1,0,0,0,0,0,0,0,0))
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
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_low, "uM_2h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 8hr with treatment 8hr (low concentration)
cont.matrix <- cbind(ctrl8=c(0,1,0,0,0,0,0,0,0,0,0,0),treat8=c(0,0,0,0,1,0,0,0,0,0,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl8=results[,1],treat8=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat8",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_low, "uM_8h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 24hr with treatment 24hr (low concentration)
cont.matrix <- cbind(ctrl24=c(0,0,1,0,0,0,0,0,0,0,0,0),treat24=c(0,0,0,0,0,1,0,0,0,0,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl24=results[,1],treat24=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat24",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_low, "uM_24h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 2hr with treatment 2hr (mid concentration)
cont.matrix <- cbind(ctrl2=c(1,0,0,0,0,0,0,0,0,0,0,0),treat2=c(0,0,0,0,0,0,1,0,0,0,0,0))
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
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_mid, "uM_2h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 8hr with treatment 8hr (mid concentration)
cont.matrix <- cbind(ctrl8=c(0,1,0,0,0,0,0,0,0,0,0,0),treat8=c(0,0,0,0,0,0,0,1,0,0,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl8=results[,1],treat8=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat8",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_mid, "uM_8h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 24hr with treatment 24hr (mid concentration)
cont.matrix <- cbind(ctrl24=c(0,0,1,0,0,0,0,0,0,0,0,0),treat24=c(0,0,0,0,0,0,0,0,1,0,0,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl24=results[,1],treat24=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat24",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_mid, "uM_24h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 2hr with treatment 2hr (high concentration)
cont.matrix <- cbind(ctrl2=c(1,0,0,0,0,0,0,0,0,0,0,0),treat2=c(0,0,0,0,0,0,0,0,0,1,0,0))
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
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_high, "uM_2h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 8hr with treatment 8hr (high concentration)
cont.matrix <- cbind(ctrl8=c(0,1,0,0,0,0,0,0,0,0,0,0),treat8=c(0,0,0,0,0,0,0,0,0,0,1,0))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl8=results[,1],treat8=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat8",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_high, "uM_8h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

# contrast control 24hr with treatment 24hr (high concentration)
cont.matrix <- cbind(ctrl24=c(0,0,1,0,0,0,0,0,0,0,0,0),treat24=c(0,0,0,0,0,0,0,0,0,0,0,1))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
i <- grep("AFFX",featureNames(eset))
options(digits=10)
summary_output <- summary(fit2$F.p.value[i])
p_value_cutoff <- signif(as.numeric(summary_output)[1], digits=3) *0.99
results <- classifyTestsF(fit2, p.value=p_value_cutoff)
summary(results)
table(ctrl24=results[,1],treat24=results[,2])
vennDiagram(results,include="up")
vennDiagram(results,include="down")
options(digits=3)
diff_exp <- topTable(fit2,coef="treat24",n=10000, p=significance_threshold)
diff_exp_brief <- data.frame(diff_exp$Symbol, diff_exp$logFC, diff_exp$adj.P.Val)
write.table(diff_exp_brief, 
            file=paste(outputdir, paste(chem_long_name, paste(chem_concentr_high, "uM_24h_top_10k_genes.txt", sep="-"), sep="_"), sep="/"), 
            quote=FALSE, sep='\t', col.names = NA)

sessionInfo()

