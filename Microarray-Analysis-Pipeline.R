# #####################
# General configuration
# #####################

baseDir <- getwd()
print(baseDir)
targetsFile <- 'targets.tsv'

options(scipen = 99) # prevent scientific notation

library(BiocManager)
require(limma)
require(statmod)

require(stringr)

library(kableExtra)
library(gplots)
library(tidyr)



# ################
# Make valid names
# ################
sh_sy5y <- make.names("SH-SY5Y")
sk_n_be_2 <- make.names("SK-N-BE(2)")
rencell_vm <- make.names("ReNcell VM")
kelly <- make.names("Kelly")
adult_brain <- make.names("Adult Brain")
hippocampus <- make.names("Hippocampus")
insula <- make.names("Insula")



# ##################################
# Read or create the annotation file
# ##################################

# Create a folder for annotation files
dir_annotation <- paste0(baseDir, "/annotation-files")

if (!dir.exists(dir_annotation)) {
  dir.create(dir_annotation)
}

setwd(dir_annotation)

tryCatch(
  expr = {
    currentMonth <- format(Sys.Date(), format="%Y-%m")

    annotationFile <- paste0(
      'Human_agilent_sureprint_g3_ge_8x60k_v2_', 
      gsub("-", "_", as.character(currentMonth)), 
      '.tsv'
    )
    
    if (!file.exists(annotationFile)) {
      source(paste0(baseDir, '/', 'Annotation-File-Generator.R'), local = TRUE)
    }
  },
  finally = {
    setwd(baseDir)
  }
)



# ################
# Read in the data
# ################

# Targets
targetinfo <- readTargets(targetsFile, row.names = 'Name')
targetinfo[order(targetinfo$Name),]

# Converts the raw data to an EListRaw object
wtAgilent.GFilter <- function(qta) { qta[,"gIsPosAndSignif"] }
eset <- read.maimages(
  targetinfo,
  source = 'agilent.median',
  green.only = TRUE,
  path = "data",
  other.columns = 'gIsWellAboveBG',
  wt.fun = wtAgilent.GFilter
)

colnames(eset) <- row.names(targetinfo)
eset <- eset[, order(colnames(eset))]

# Add the spot type
spotTypes <- readSpotTypes(file = 'SpotTypes.tsv')
eset$genes$Status <- controlStatus(spotTypes, eset)



# ###################
# Annotate the probes
# ###################

annotLookup <- read.csv(
  paste0(dir_annotation, '/', annotationFile),
  header = TRUE,
  sep = '\t',
  stringsAsFactors = FALSE
)

colnames(annotLookup)[1] <- 'AgilentID'

annotLookup <- annotLookup[which(annotLookup$AgilentID %in% eset$genes$ProbeName),]
annotLookup <- annotLookup[match(eset$genes$ProbeName, annotLookup$AgilentID),]
table(eset$genes$ProbeName == annotLookup$AgilentID) # check that annotations are aligned

eset$genes$AgilentID <- annotLookup$AgilentID
eset$genes$wikigene_description <- annotLookup$wikigene_description
eset$genes$ensembl_gene_id <- annotLookup$ensembl_gene_id
eset$genes$entrezgene <- annotLookup$entrezgene
eset$genes$gene_biotype <- annotLookup$gene_biotype
eset$genes$external_gene_name <- annotLookup$external_gene_name



# #############################
# Perform background correction
# #############################

eset <- backgroundCorrect(eset, method = 'normexp')
plotDensities(eset, legend = FALSE, main = 'Density plot with background-corrected data')



# ##############
# Normalize data
# ##############

eset <- normalizeBetweenArrays(eset, method = 'quantile')
plotDensities(eset, legend = FALSE, main = 'Density plot with normalized, background-corrected data')



# #######
# MA-Plot
# #######

plotMDS(eset)



# ##########################
# Filter out control probes, 
# those with no symbol, 
# and those that fail
# ##########################

Control <- eset$genes$ControlType != 0
NoSymbol <- is.na(eset$genes$external_gene_name)
IsExpr <- rowSums(eset$other$gIsWellAboveBG > 0) >= 3

eset <- eset[!Control & !NoSymbol & IsExpr, ]



# ##########################################
# Remove annotation columns no longer needed
# ##########################################

eset$genes <- eset$genes[,c(
  'ProbeName','wikigene_description','ensembl_gene_id','entrezgene','gene_biotype','external_gene_name'
)]



# ############################################
# Replace replicate probes with the mean value
# ############################################

eset <- avereps(eset, ID = eset$genes$ProbeName)



# #####################
# Array quality weights
# #####################

array.weights <- arrayWeights(eset)
barplot(array.weights, xlab = "Array", ylab = "Weight", main = "Array weights", col = "white", las = 2)
abline(h = 1, lwd = 1, lty = 2)



# ###############################
# Create a linear fit of the data
# ###############################

colnames <- c(sh_sy5y, sk_n_be_2, rencell_vm, kelly, insula, hippocampus, adult_brain)
f <- factor(targetinfo$CellType, levels = colnames)

design <- model.matrix(~0+f)
colnames(design) <- colnames

fit <- lmFit(eset, design = design, weights = array.weights)
fit <- eBayes(fit)



# ###########################
# Create a Student's Q-Q Plot
# ###########################

qqt(fit$t, df = fit$df.prior + fit$df.residual, pch = 16, cex = 0.2)
abline(0,1)



# ################
# Create a boxplot
# ################

cols <- eset$targets$CellType
cols[cols == sh_sy5y] <- "steelblue1"
cols[cols == sk_n_be_2] <- "steelblue2"
cols[cols == rencell_vm] <- "steelblue3"
cols[cols == kelly] <- "steelblue4"
cols[cols == adult_brain] <- "springgreen1"
cols[cols == hippocampus] <- "springgreen2"
cols[cols == insula] <- "springgreen3"

boxplot(eset$E~col(eset$E), names = rownames(eset$targets), col = cols, xlab = "Cell Type", ylab = "E-values")



# ##############################
# Pairwise comparison of FC data
# ##############################

# Adult brain vs cell lines
difference.adult_brain.sh_sy5y <- paste(adult_brain,"-",sh_sy5y, sep="")
difference.adult_brain.sk_n_be_2 <- paste(adult_brain,"-",sk_n_be_2, sep="")
difference.adult_brain.rencell_vm <- paste(adult_brain,"-",rencell_vm, sep="")
difference.adult_brain.kelly <- paste(adult_brain,"-",kelly, sep="")

# Hippocampus vs cell lines
difference.hippocampus.sh_sy5y <- paste(hippocampus,"-",sh_sy5y, sep="")
difference.hippocampus.sk_n_be_2 <- paste(hippocampus,"-",sk_n_be_2, sep="")
difference.hippocampus.rencell_vm <- paste(hippocampus,"-",rencell_vm, sep="")
difference.hippocampus.kelly <- paste(hippocampus,"-",kelly, sep="")

# Insula vs cell lines
difference.insula.sh_sy5y <- paste(insula,"-",sh_sy5y, sep="")
difference.insula.sk_n_be_2 <- paste(insula,"-",sk_n_be_2, sep="")
difference.insula.rencell_vm <- paste(insula,"-",rencell_vm, sep="")
difference.insula.kelly <- paste(insula,"-",kelly, sep="")

# Tests for the all the cell lines
contrast.matrix.cell_lines <- makeContrasts(
  difference.adult_brain.sh_sy5y,
  difference.hippocampus.sh_sy5y,
  difference.insula.sh_sy5y,
  difference.adult_brain.sk_n_be_2,
  difference.hippocampus.sk_n_be_2,
  difference.insula.sk_n_be_2,
  difference.adult_brain.kelly,
  difference.hippocampus.kelly,
  difference.insula.kelly,
  difference.adult_brain.rencell_vm,
  difference.hippocampus.rencell_vm,
  difference.insula.rencell_vm,
  levels = design
)

fit.cell_lines <- contrasts.fit(fit, contrast.matrix.cell_lines)
fit.cell_lines <- eBayes(fit.cell_lines, trend = TRUE, robust = TRUE)
results.cell_lines <- decideTests(fit.cell_lines, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary(results.cell_lines)

# Tests for the SH-SY5Y cell line
contrast.matrix.sh_sy5y <- makeContrasts(
  difference.adult_brain.sh_sy5y,
  difference.hippocampus.sh_sy5y,
  difference.insula.sh_sy5y,
  levels = design
)

fit.sh_sy5y <- contrasts.fit(fit, contrast.matrix.sh_sy5y)
fit.sh_sy5y <- eBayes(fit.sh_sy5y, trend = TRUE, robust = TRUE)
results.sh_sy5y <- decideTests(fit.sh_sy5y, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary(results.sh_sy5y)

# Tests for the SK-N-BE(2) cell line
contrast.matrix.sk_n_be_2 <- makeContrasts(
  difference.adult_brain.sk_n_be_2, 
  difference.hippocampus.sk_n_be_2, 
  difference.insula.sk_n_be_2,
  levels = design
)

fit.sk_n_be_2 <- contrasts.fit(fit, contrast.matrix.sk_n_be_2)
fit.sk_n_be_2 <- eBayes(fit.sk_n_be_2, trend = TRUE, robust = TRUE)
results.sk_n_be_2 <- decideTests(fit.sk_n_be_2, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary(results.sk_n_be_2)

# Tests for the ReNcell VM cell line
contrast.matrix.rencell_vm <- makeContrasts(
  difference.adult_brain.rencell_vm, 
  difference.hippocampus.rencell_vm, 
  difference.insula.rencell_vm,
  levels = design
)

fit.rencell_vm <- contrasts.fit(fit, contrast.matrix.rencell_vm)
fit.rencell_vm <- eBayes(fit.rencell_vm, trend = TRUE, robust = TRUE)
results.rencell_vm <- decideTests(fit.rencell_vm, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary(results.rencell_vm)

# Tests for the Kelly cell line
contrast.matrix.kelly <- makeContrasts(
  difference.adult_brain.kelly, 
  difference.hippocampus.kelly,
  difference.insula.kelly,
  levels = design
)

fit.kelly <- contrasts.fit(fit, contrast.matrix.kelly)
fit.kelly <- eBayes(fit.kelly, trend = TRUE, robust = TRUE)
results.kelly <- decideTests(fit.kelly, method = "separate", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary(results.kelly)



# ###############################################
# Create graphs to visualize the analysis results
# ###############################################

# Create a folder for graphs
dir_graphics <- paste0(baseDir, "/graphics")

if (!dir.exists(dir_graphics)) {
  dir.create(dir_graphics)
}

setwd(dir_graphics)

tryCatch(
  expr = {
    # Create graphs for the SH-SY5Y cell line
    pdf(file = 'plotSA.sh_sy5y.pdf')
    plotSA(fit.sh_sy5y)
    dev.off()
    
    pdf(file = 'vennDiagram.sh_sy5y.up.pdf')
    vennDiagram(results.sh_sy5y, include = "up", main = "Up-regulated")
    dev.off()
    
    pdf(file = 'vennDiagram.sh_sy5y.down.pdf')
    vennDiagram(results.sh_sy5y, include = "down", main = "Down-regulated")
    dev.off()
    
    # Create graphs for the SK-N-BE(2) cell line
    pdf(file = 'plotSA.sk_n_be_2.pdf')
    plotSA(fit.sk_n_be_2)
    dev.off()
    
    pdf(file = 'vennDiagram.sk_n_be_2.up.pdf')
    vennDiagram(results.sk_n_be_2, include = "up", main = "Up-regulated")
    dev.off()
    
    pdf(file = 'vennDiagram.sk_n_be_2.down.pdf')
    vennDiagram(results.sk_n_be_2, include = "down", main = "Down-regulated")
    dev.off()
    
    # Create graphs for the ReNcell VM cell line
    pdf(file = 'plotSA.rencell_vm.pdf')
    plotSA(fit.rencell_vm)
    dev.off()
    
    pdf(file = 'vennDiagram.rencell_vm.up.pdf')
    vennDiagram(results.rencell_vm, include = "up", main = "Up-regulated")
    dev.off()
    
    pdf(file = 'vennDiagram.rencell_vm.down.pdf')
    vennDiagram(results.rencell_vm, include = "down", main = "Down-regulated")
    dev.off()
    
    # Create graphs for the Kelly cell line
    pdf(file = 'plotSA.kelly.pdf')
    plotSA(fit.kelly)
    dev.off()
    
    pdf(file = 'vennDiagram.kelly.up.pdf')
    vennDiagram(results.kelly, include = "up", main = "Up-regulated")
    dev.off()
    
    pdf(file = 'vennDiagram.kelly.down.pdf')
    vennDiagram(results.kelly, include = "down", main = "Down-regulated")
    dev.off()
  },
  finally = {
    setwd(baseDir)
  }
)


# ############################################
# Replace replicate arrays with the mean value
# ############################################

eset_ave <- avearrays(eset, ID = eset$targets$CellType, weights = eset$weights)
eset_ave <- eset_ave[rowSums(is.na(eset_ave$E))==0,]
colnames(eset_ave$E) <- c('Adult Brain', 'Hippocampus', 'Insula', 'Kelly', 'ReNcell VM', 'SH-SY5Y', 'SK-N-BE(2)')



# ################
# Create a heatmap
# ################

coolmap(eset_ave, cluster.by="de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol=45, labRow='')

setEPS()
postscript(file = paste0(dir_graphics, '/heatmap_z-value.eps'), width = 5, height = 4)
coolmap(eset_ave, cluster.by="de pattern", show.dendrogram = "column", margins = c(7, 1), srtCol=45, labRow='')
dev.off()



coolmap(eset_ave, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='')

setEPS()
postscript(file = paste0(dir_graphics, '/heatmap_expression-level.eps'), width = 5, height = 4)
coolmap(eset_ave, cluster.by="expression level", show.dendrogram = "column", col = "redblue", margins = c(7, 1), srtCol=45, labRow='')
dev.off()


plotMDS(eset_ave, labels = NULL, pch = c(rep(15, 3), rep(16, 4)), col = 1:7)
legend(-3.5, 4, legend = colnames(eset_ave$E), col = 1:7, pch = c(rep(15, 3), rep(16, 4)), bty = 'n', cex = 0.9)

setEPS()
postscript(file = paste0(dir_graphics, '/multidimensional_scaling.eps'), width = 5, height = 4)
plotMDS(eset_ave, labels = NULL, pch = c(rep(15, 3), rep(16, 4)), col = 1:7)
legend(-3.5, 4, legend = colnames(eset_ave$E), col = 1:7, pch = c(rep(15, 3), rep(16, 4)), bty = 'n', cex = 0.9)
dev.off()
