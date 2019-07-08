########################################################################################
## Processing Illumina EPIC array data, identifying DMPs between islands and villages ##
########################################################################################

# Read in libraries 
library(DMRcate)
library(ggplot2)
library(IlluminaHumanMethylationEPICmanifest)
library(lumi)
library(limma)
library(methylationArrayAnalysis)
library(minfi)
library(missMethyl)
library(reshape)

########################
## Process EPIC array ##
########################

# Setting the working directory
setwd("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/Differential_Indonesia/")

# Set the directory containing the raw data and the sample sheet
baseDir_batch1 = ("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/171110-171011AM-01_Guy Jacobs_80/(1)Raw.data/idat/")
baseDir_batch2 = ("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/180727-180710AM-01_Guy_Jacobs_40/(1)Raw.data/idat/")

# Reads in .csv file found in the baseDir
targets_batch1 <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/171110-171011AM-01_Guy Jacobs_80/(1)Raw.data/Sample.Table.txt", 
                             stringsAsFactors = FALSE, sep = "\t", header = TRUE)
head(targets_batch1)
targets_batch2 <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/180727-180710AM-01_Guy_Jacobs_40/(1)Raw.data/Sample.Table.txt", 
                             stringsAsFactors = FALSE, sep = "\t", header = TRUE)
head(targets_batch2)

# Add basename
targets_batch1$Basename <- file.path(baseDir_batch1, paste0(targets_batch1$Sentrix.Barcode, "_",
                                                            targets_batch1$Sample.Section))
targets_batch2$Basename <- file.path(baseDir_batch2, paste0(targets_batch2$Sentrix.Barcode, "_",
                                                            targets_batch2$Sample.Section))

targets_batch1$id2 <- paste0(targets_batch1$Sample.ID, "_methyl1")
targets_batch2$id2 <- paste0(targets_batch2$Sample.ID, "_methyl2")

# Identifies .idat files and creats an RGset object 
RGset_batch1 <- read.metharray.exp(targets = targets_batch1)
RGset_batch2 <- read.metharray.exp(targets = targets_batch2)

# Combining arrays
RGset_batch1_batch2 <- combineArrays(RGset_batch1, RGset_batch2, outType = "IlluminaHumanMethylationEPIC", verbose = TRUE)

# Gets phenotype info from data
pd_batch1_batch2 <- pData(RGset_batch1_batch2)

# Makes a pdf of various QC steps
qcReport(RGset_batch1_batch2, sampNames = pd_batch1_batch2$id2, sampGroups = pd_batch1_batch2$Sample.Group, pdf = "/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/Plots/batch1_batch2_qcReport.pdf")

# Corrects background from the array
datas.bg_batch1_batch2 <- bgcorrect.illumina(RGset_batch1_batch2)

# Performs other preprocessing on the array
datas.bg_batch1_batch2.raw <- preprocessRaw(datas.bg_batch1_batch2)

# Keep only probes that do not overlap SNPs segregating in these populations
# Dataframe with probe IDs and sites
probes <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/epic_probes.tsv", header=TRUE)
head(probes)
# Dataframe with segregating sites
sites_to_remove <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/unique_segregaring_sites.tsv", header=FALSE)
head(sites_to_remove)
sites_to_remove <- sites_to_remove$V1
length(sites_to_remove)
probes_to_remove <- subset(probes, (position %in% sites_to_remove))
dim(probes_to_remove)
probes_to_keep <- subset(probes, !(position %in% sites_to_remove))
probes.to.keep = probes_to_keep$probeid
length(probes.to.keep)

# Makes a vector of the probes to keep
unique <- featureNames(datas.bg_batch1_batch2.raw) %in% as.matrix(probes.to.keep)

# Identifies and removes probes with low signal across arrays
signal_batch1_batch2 <- detectionP(datas.bg_batch1_batch2)
signal.sufficient <- signal_batch1_batch2 < 0.01
signal.keep <- rowMeans(signal.sufficient) > 0.75

# Number of probes to keep
length(signal.keep)

# Filtering out bad quality probes
datas.bg_batch1_batch2.raw.filtered <- datas.bg_batch1_batch2.raw[unique & signal.keep, ]
dim(datas.bg_batch1_batch2.raw.filtered)

# QC on filtered probes
indoQC <- getQC(datas.bg_batch1_batch2.raw.filtered)
head(indoQC)
pdf(file="/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/Plots/indo_raw_qc.pdf")
plotQC(indoQC)
dev.off()

# Performs SWAN normalization (Subset-quantile Within Array Normalization) for an array
set.seed(456)
datas.bg_batch1_batch2.filtered.swan <- preprocessSWAN(RGset_batch1_batch2, mSet = datas.bg_batch1_batch2.raw.filtered)

# Get red and green channel (Meth and Unmeth)
M <- minfi::getMeth(datas.bg_batch1_batch2.filtered.swan)
U <- getUnmeth(datas.bg_batch1_batch2.filtered.swan)

# Quantile normalize red and green channels
M.normalized <- lumiN(M, method = 'quantile')
U.normalized <- lumiN(U, method = 'quantile')

# Creates a MethylSet object using normalized data
meth.normalized <- MethylSet(Meth = as.matrix(M.normalized), Unmeth = as.matrix(U.normalized))

# Get annotations for final set of probes
annotations <- read.csv("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/MethylationEPIC_v-1-0_B4.csv", skip = 7, header = TRUE)
annotations_to_keep <- subset(annotations, IlmnID %in% probes.to.keep)

# Get a data frame of betas for all of the probes
meth.final <- as.data.frame(getBeta(meth.normalized))

# Renames columns with individual ID
colnames(meth.final) <- c(targets_batch1[,"id2"], targets_batch2[,"id2"])
colnames(meth.final) <- gsub("_", "-", colnames(meth.final))
colnames(meth.final) <- gsub("-new", "", colnames(meth.final))

# Calculating M-values
mval <- log2(meth.final/(1-meth.final))

# Writing beta and m values to a file
beta_path <- "/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/beta_vals.tsv"
m_path <- "/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/m_vals.tsv"
write.table(meth.final, beta_path, sep = "\t")
write.table(mval, m_path, sep = "\t")

##############################
## Differential methylation ##
##############################

meth.final <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/beta_vals.tsv", sep = "\t", header = TRUE)
head(meth.final)
mval <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/m_vals.tsv")
head(mval)

# Read covariate matrix
covariates <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/metadata.tsv", header=TRUE, sep="\t")

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
covariates$Age[which(is.na(covariates$Age) == T)]=45

# Sorting to match the order of samples in the methylation dataframe
samplenames <- colnames(mval)
samplenames <- gsub(".methyl1", "", samplenames)
samplenames <- gsub(".methyl2", "", samplenames)
samplenames <- gsub("\\.", "-", samplenames)
length(samplenames)
covariates <- covariates[match(samplenames, covariates$Sample.ID),]
dim(covariates)

covariates$Island <- gsub("West Papua", "Mappi", covariates$Island)
covariates$Island <- as.factor(covariates$Island)
covariates$Sequencing.Batch <- as.factor(covariates$Sequencing.Batch)
covariates$methyl_id <- colnames(mval)
covariates$methyl_exp_id <- with(covariates, paste0(methyl_id, "-rna", Sequencing.Batch))

# Adding blood deconvolution estimates
blood <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/predictedCellCounts_DeconCell_edit.txt", sep="\t", as.is=T, header=T)
head(blood)
covariates[,16:22] <- blood[match(covariates$Sample.ID, blood$ID),1:7]
head(covariates)

# Adding methylation batch
covariates$methyl_batch <- as.numeric(grepl('methyl1', covariates$methyl_id, ignore.case=T))

# Removing MPI-296 batch 2 sample 
# covariates <- covariates[!(covariates$methyl_id=="MPI-296-methyl2"),]
# colnames(mval)
# mval <- mval[,!(colnames(mval)=="MPI-296-methyl2")]

# Replacing column names
colnames(mval) <- covariates$Sample.ID

# Creating a design matrix for islands
design <- model.matrix(~0 + covariates$Island + covariates$Age + covariates$methyl_batch + covariates$CD8T + covariates$CD4T + covariates$NK + covariates$Bcell + covariates$Mono + covariates$Gran)
colnames(design) <- c("Mappi", "Mentawai", "Sumba", "Age", "Methyl_batch", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# Defining pairwise comparisons
contrasts <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
all_comparisons <-colnames(data.frame(contrasts))

# Fitting the linear model
fit <- lmFit(mval, design)
vfit <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(vfit)

# Getting significant probes between islands
coef=1

for (i in all_comparisons){
  toptable <- topTable(efit, adjust="BH", coef=coef, num=Inf)
  sig_probes <- toptable[which(toptable$adj.P.Val<=1 & abs(toptable$logFC)>=0),]
  path <- paste("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMP_", i, "_adjp1_logFC0_DeconCell_new_allsamples.txt", sep = "")
  write.table(sig_probes, file=path, sep="\t", quote=TRUE)
  coef=coef+1
}

# Finding differentially methylated regions
mval_matrix <- as.matrix(mval)
coef=1

for (i in all_comparisons){
  myAnnotation <- cpg.annotate(object = mval_matrix, datatype = "array", what = "M",
                               annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19"),
                               analysis.type = "differential", design = design,
                               contrasts = TRUE, cont.matrix = contrasts,
                               coef = i, arraytype = "EPIC", fdr = 0.01)
  
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, betacutoff = 0, min.cpgs=2)
  
  # Converting the regions to annotated genomic ranges
  DMRs_results_ranges <- extractRanges(DMRs, genome = "hg19")
  
  # Saving annotated ranges to a file
  DMRs_results_ranges_df = as(DMRs_results_ranges, "data.frame")
  path <- paste("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMRs_", i, "_fdr001_beta0_pcutoffdefault_mincpgs2_DeconCell_new_allsamples.txt", sep = "")
  write.table(DMRs_results_ranges_df, path, sep = "\t")
}



# Creating a design matrix for villages
design <- model.matrix(~0 + covariates$Sampling.Site + covariates$Age + covariates$methyl_batch + covariates$CD8T + covariates$CD4T + covariates$NK + covariates$Bcell + covariates$Mono + covariates$Gran)
colnames(design) <- c("Anakalung", "Bilarenge", "HupuMada", "Madobag", "Mappi", "PadiraTana", "PatialaBawa", "Rindi", "Taileleu", "Wunga", "WuraHomba", "Age", "Methyl_batch", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# Defining pairwise comparisons
# contrasts <- makeContrasts(ANKvsWNG=Anakalung - Wunga, ANKvsRIN=Anakalung - Rindi, ANKvsHPM=Anakalung - HupuMada, ANKvsPDT=Anakalung - PatialaBawa,
#                            WNGvsRIN=Wunga - Rindi, WNGvsHPM=Wunga - HupuMada, WNGvsPDT=Wunga - PadiraTana, RINvsHPM=Rindi - HupuMada,
#                            RINvsPDT=Rindi -PadiraTana, HPMvsPDT=HupuMada - PadiraTana, ANKvsMDB=Anakalung - Madobag, ANKvsTLL=Anakalung - Taileleu, 
#                            WNGvsMDB=Wunga - Madobag, WNGvsTLL=Wunga - Taileleu, RINvsMDB=Rindi - Madobag, RINvsTLL=Rindi - Taileleu,
#                            HPMvsMDB=HupuMada - Madobag, HPMvsTLL=HupuMada - Taileleu, PDTvsMDB=PadiraTana - Madobag, PDTvsTLL=PadiraTana - Taileleu,
#                            ANKvsMPI=Anakalung - Mappi, WNGvsMPI=Wunga - Mappi, RINvsMPI=Rindi - Mappi, HPMvsMPI=HupuMada - Mappi,
#                            PDTvsMPI=PadiraTana - Mappi, MDBvsMPI=Madobag - Mappi, TLLvsMPI=Taileleu - Mappi, MDBvsTLL=Madobag - Taileleu,
#                            levels=colnames(design))

contrasts <- makeContrasts(ANKvsWNG=Anakalung - Wunga, ANKvsMDB=Anakalung - Madobag, ANKvsTLL=Anakalung - Taileleu, 
                           WNGvsMDB=Wunga - Madobag, WNGvsTLL=Wunga - Taileleu, 
                           ANKvsMPI=Anakalung - Mappi, WNGvsMPI=Wunga - Mappi,
                           MDBvsMPI=Madobag - Mappi, TLLvsMPI=Taileleu - Mappi, MDBvsTLL=Madobag - Taileleu,
                           levels=colnames(design))

all_comparisons <- colnames(data.frame(contrasts))

# Fitting the linear model
fit <- lmFit(mval, design)
vfit <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(vfit)

# Getting significant probes between villages
coef=1

for (i in all_comparisons){
  toptable <- topTable(efit, adjust="BH", coef=coef, num=Inf)
  sig_probes <- toptable[which(toptable$adj.P.Val<=0.01 & abs(toptable$logFC)>=1),]
  path <- paste("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMP_", i, "_adjp001_logFC1_DeconCell_new_allsamples.txt", sep = "")
  write.table(sig_probes, file=path, sep="\t", quote=TRUE)
  coef=coef+1
}

# Finding differentially methylated regions
mval_matrix <- as.matrix(mval)
coef=1

for (i in all_comparisons){
  myAnnotation <- cpg.annotate(object = mval_matrix, datatype = "array", what = "M",
                               annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19"),
                               analysis.type = "differential", design = design,
                               contrasts = TRUE, cont.matrix = contrasts,
                               coef = i, arraytype = "EPIC", fdr = 0.01)
  
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, betacutoff = 0, min.cpgs=2)
  
  # Converting the regions to annotated genomic ranges
  DMRs_results_ranges <- extractRanges(DMRs, genome = "hg19")
  
  # Saving annotated ranges to a file
  DMRs_results_ranges_df = as(DMRs_results_ranges, "data.frame")
  path <- paste("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMRs_", i, "_fdr001_beta0_pcutoffdefault_mincpgs2_DeconCell_new_allsamples.txt", sep = "")
  write.table(DMRs_results_ranges_df, path, sep = "\t")
}
