# TGCA Adeno(LUAD) Data Filtering and Analysis- revised groups (<3cm)
# Kelly Cagin kelly_M_Cagin@rush.edu     10/10/18      v1
# References: 


# Read in Data from local files, based on minimum time elapsed post-treatment
cagin.clin.tcga.adeno <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/Cagin_1yr_clin-tcga-adeno_selection on 09212018.csv")
cagin.clin.tcga.adeno <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/Cagin_2yr_clin-tcga-adeno_selection on 09212018.csv")


# Sample checks for group size
length(which(cagin.clin.tcga.adeno$Group == "T=1,N=0,Disease Free"))
length(which(cagin.clin.tcga.adeno$Group == "T=1,N=0,Progressed"))
length(which(cagin.clin.tcga.adeno$Group == "T=1,N>0"))
length(which(cagin.clin.tcga.adeno$Group == "T=1,N>0,Disease Free")) #If retained
length(which(cagin.clin.tcga.adeno$Group == "T=1,N>0,Progressed")) #If retained


# Using Fresh download data and TCGA subject IDs (replace dashes, etc)
manifest.key <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/manifest.key.csv", row.names = 1)
manifest.key <- plyr::arrange(manifest.key,proj.id) #  MUST USE SORTED BY ProjID
manifest<-manifest.key

# Matching manifest and clinical Data formats for TCGA Barcode (ProjIDs) using dots and omitting sample letter
manifest$proj.id <- gsub("-",".",manifest$proj.id) 
manifest$proj.id <- substr(manifest$proj.id, start = 1, nchar(manifest$proj.id)-3)
cagin.clin.tcga.adeno$ProjID <- as.character(cagin.clin.tcga.adeno$ProjID)
clinical.index <- substr(cagin.clin.tcga.adeno$ProjID, start = 1, nchar(cagin.clin.tcga.adeno$ProjID)-2)

# Matching Project IDs for both
set <- match(clinical.index, manifest$proj.id)
set <- na.omit(set)
manifest <- manifest.key[set,]

# Preparing Data to grab and import expression set
cagin.clin.tcga.adeno <- cagin.clin.tcga.adeno[-which(cagin.clin.tcga.adeno$ProjID == "TCGA.55.1595.01"), ] #  Not present in data set from TCGA
get.files <- cbind(cagin.clin.tcga.adeno[,1:2],manifest[,6:7]) #  Make a good working file, two conventions of IDs, group and filename
get.files <- `colnames <-`(get.files, c("group", "dec.proj.id", "dash.proj.id", "files")) 

# pick save files based on clinical threshold you are using
write.csv(get.files, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.1yr.csv")#  Can save if you like, good working file to keep
write.csv(get.files, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.2yr.csv")#  Can save if you like, good working file to keep

# SORT IF NEEDED
# if starting from saved, pick desired
get.files <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.1yr.csv", row.names = 1)
get.files <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.2yr.csv", row.names = 1)

# If desired, filter as simple 2 group tests- change files for 1yr or 2yr threshold
disease.prog <- get.files[which(get.files$group == "T=1,N=0,Disease Free"|get.files$group == "T=1,N=0,Progressed"), ]
write.csv(disease.prog, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/n0df.v.progressed.3cm.groups.filenames.2yr.csv")#  Can save if you like, good working file to keep
get.files <- disease.prog #  Set for input below

df.ns <- get.files[which(get.files$group == "T=1,N=0,Disease Free"|get.files$group == "T=1,N>0"), ]
write.csv(df.ns, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/n0df.v.nodes.3cmGroups.filenames.2yr.csv")
get.files <- df.ns

progressed.ns <- get.files[which(get.files$group == "T=1,N=0,Progressed"|get.files$group=="T=1,N>0"), ]
write.csv(progressed.ns, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/progressed.v.nodes.3cm.groups.filenames.2yr.csv")
get.files <- progressed.ns

# Quick and dirty removal of unwanted levels which obscure subsetting analysis in edgeR
get.files <- as.data.frame(as.matrix(get.files)) 

# Change saving filenames as appropriate below

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)


setwd("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Read Counts/Unzipped Read Counts")
luad.counts <- readDGE(files = get.files$files, path = "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Read Counts/Unzipped Read Counts", group = get.files$group)

raw.export <- edgeR::as.matrix.DGEList(luad.counts)
write.table(raw.export, file = "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/disease.progressed.luad.raw.counts.3cm.2yr.txt")

dims <- edgeR::dimnames.DGEList(luad.counts) #  checking dim names
dims$Samples <- get.files$dash.proj.id #  changing sample names from files to IDs

# Map to HGNC From Ensembl IDs
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

# Get HGNC gene symbols 
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ens.ids <- strsplit(x = dims$Tags, fixed = TRUE, split = ".") #  remove versions from Ensembl IDs
ens.ids <- unlist(ens.ids)[c(seq(from = 1, to = length(unlist(ens.ids)), by = 2))]
valid.tags <- length(dims$Tags)-5 #  quick and dirty adjustment for preserving meta tags present in DGElist, artifact of splitting
meta.tags <- dims$Tags[-c(1:valid.tags)]
ens.ids <- c(ens.ids[1:valid.tags], meta.tags)

gene.ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"), filters = "ensembl_gene_id", values = ens.ids, mart = ensembl)

# Get gene lengths for needed for rpkm later 
gene.length <- gene.ids
gene.length$bp.length <- gene.length$end_position - gene.length$start_position
with.rpkm <- match(gene.length$ensembl_gene_id, ens.ids)
write.csv(gene.length, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/gene.length.3cm.2yr.csv")

# Replace Ensmbl with gene symbol where applicable
gene.ids <- gene.ids[-which(gene.ids$hgnc_symbol == ""), ] #  Remove valid ensemble IDS lacking HGNC gene symbols
gene.index <- match(gene.ids$ensembl_gene_id, ens.ids)
gene.index <- na.omit(gene.index)
ens.ids[gene.index] <- gene.ids$hgnc_symbol
dims$Tags <- ens.ids

# fix dimensions to IDs and HGNC Symbols as desired
luad.counts <- edgeR::`dimnames<-.DGEList`(luad.counts,dims) 

luad.cpm <- cpm.DGEList(luad.counts)
write.csv(luad.cpm,"C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/all.groups.luad.cpm.3cm.2yr.csv")
 
keep <- rowSums(luad.cpm > 1) >= 30
luad.counts <- luad.counts[keep, , keep.lib.sizes=FALSE]
  
luad.counts <- calcNormFactors.DGEList(luad.counts, method = "TMM") #normalize by Trimmed Means of M
luad.exp.model <- model.matrix(~ luad.counts$samples$group)
write.csv(luad.exp.model, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/all.groups.luad.exp.model.3cm.2yr.csv")
  
luad.counts <- estimateDisp.DGEList(luad.counts, rowsum.filter = 10000, design = luad.exp.model) #  remove very low expression genes with filter here


luad.fit<- glmFit.DGEList(y = luad.counts, design = luad.exp.model)
luad.lrt<- glmLRT(luad.fit)

q.val <- p.adjust(luad.lrt$table$PValue, method = "BH")
gene <- rownames(luad.lrt$table)
luad.lrt$table <- cbind(luad.lrt$table, q.val, gene)
 
write.table(luad.lrt$table, file = "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/luad.diff.analysis.all.groups.3cm.2yr.csv", sep = "\t", quote = F)
luad.lrt$table <- plyr::arrange(df = luad.lrt$table, q.val)
significant.at.0.01 <- luad.lrt$table[which(luad.lrt$table$q.val<=0.01), ]
write.table(significant.at.0.01, file = "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.at.01.all.groups.3cm.2yr.csv", sep = "\t", quote = F)

# Overlap options- pick correct save file for clinical threshold used above

sig.others.df <- intersect(significant.at.01.n0df.v.progressed$X, significant.at.01.n0df.v.ns$X) #  from all disease free 
sig.others.n <- intersect(significant.at.01.n0df.v.ns$X, significant.at.01.progressed.v.ns$X) #  from all Ns if desired
write.csv(sig.others.df, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.n0.analyses.3cm.1yr.csv")
write.csv(sig.others.df, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.n0.analyses.3cm.2yr.csv")

sig.all <- intersect(significant.at.01.progressed.v.ns$X, intersect(significant.at.01.n0df.v.progressed$X, significant.at.01.n0df.v.ns$X))
write.csv(sig.all, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.analyses.3cm.1yr.csv")
write.csv(sig.all, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.analyses.3cm.2yr.csv")


# Plots

install.packages("plotly")
library(shiny)
library(ggplot2)
library(plotly)
library(circlize)
library(heatmap3)
library(colorRamps) #if using

# Volcano

#If needed to read in data, pick correct yr for clinical threshold!

significant.at.0.01 <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.at.01.all.groups.3cm.1yr.csv", row.names=1)
significant.at.0.01 <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.at.01.all.groups.3cm.2yr.csv", row.names=1)

sig.others.df <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.n0.analyses.3cm.1yr.csv", row.names=1)
sig.others.df <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.n0.analyses.3cm.2yr.csv", row.names=1)

sig.all <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.analyses.3cm.1yr.csv", row.names=1)
sig.all <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/significant.between.all.analyses.3cm.2yr.csv", row.names=1)

get.files <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.1yr.csv", row.names = 1)
get.files <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/sorted.groups.filenames.2yr.csv", row.names = 1)

# adjust titles as needed when selecting other data sets 

workingset <- significant.at.0.01 #  Make sure this remains sorted by Q
workingset$q.val <- log10(workingset$q.val)*-1

index <- na.omit(match(sig.others.df$X,workingset$gene))

highlight<-na.omit(match(sig.all$X,workingset$gene))

top.q <- as.character(workingset$gene[1:10])

objectplot <- ggplot(data = workingset, aes(x = `logFC`, y = `q.val`)) +
  geom_point(data = filter(workingset[-unlist(index), ]), aes(color = "Significant for Node-Negative Progression Only"), size = 3, alpha = 0.36) +
  geom_point(data = filter(workingset[unlist(index), ]), aes(color = "Significant for Node-Negative Progression and Nodal Status"), size = 3, alpha = 0.45) +
  geom_point(data = filter(workingset[unlist(highlight), ]), aes(color = "Significant in all 2-Way Comparisons"), size = 3, alpha = 0.45) +
  geom_text(data = filter(workingset[1:10, ]), label = top.q, nudge_x = 0.5, nudge_y = 0.5, check_overlap = TRUE, size = 4) +
  scale_color_manual(values = c("Significant for Node-Negative Progression Only" = "mediumslateblue", "Significant for Node-Negative Progression and Nodal Status" = "green1", "Significant in all 2-Way Comparisons"="tomato")) +
  geom_rug(color = "steelblue", alpha = 0.1) +
  theme_bw(base_size = 15) +
  theme(panel.border = element_rect(fill = NA), panel.background = element_rect(fill = "white"), legend.position = "bottom", legend.direction = "vertical") +
  xlab("Node-Negative Progression: log FC") +
  ylab("-log10 Q Value (BH)") +
  ggtitle("Single-Gene Variance for Node-Negative Progression", subtitle = "Lung Adenocarcinoma (TCGA-LUAD) Data, at minimum 24 months, q<0.01")

plot(objectplot)
  

# Heat Map

# Convert data to log2 rpkm
with.rpkm <- c(na.omit(match(gene.length$hgnc_symbol, rownames(luad.counts$counts))), na.omit(match(gene.length$ensembl_gene_id, rownames(luad.counts$counts))))
luad.rpkm <- luad.counts[with.rpkm, ] #  filter and gene length from above
length.order <- c(na.omit(match(rownames(luad.rpkm), gene.length$hgnc_symbol)), na.omit(match(rownames(luad.rpkm),gene.length$ensembl_gene_id)))
rpkm.bp <- gene.length$bp.length[length.order]
luad.rpkm <- edgeR::rpkm.DGEList(y = luad.rpkm, gene.length = rpkm.bp, log = TRUE, normalized.lib.sizes = TRUE) #  Using log2 Values for rpkm
write.csv(luad.rpkm, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/luad.rpkm.3cm.2yr.csv")

luad.rpkm <- luad.rpkm[na.omit(match(significant.at.0.01$gene,rownames(luad.rpkm))),] #  If filtering for significance from previous comparison

# Define groupswrite.csv(luad.rpkm, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/luad.rpkm.3cm.2yr.csv")

groups <- match(colnames(luad.rpkm), get.files$dash.proj.id)
groups <- as.character(get.files$group[groups])
groups <- gsub("T=1,N=0,Progressed", "orange1", x = groups)            
groups <- gsub("T=1,N=0,Disease Free", "skyblue1", x = groups)  
groups <- gsub("T=1,N>0", "mediumpurple1", x = groups)  

heatmap3::heatmap3(x = luad.rpkm,col = colorspace::diverge_hcl(200, power = 0.8, c = 350, h = c(219, 19), l = 90), balanceColor = TRUE, ColSideColors = groups, showColDendro = FALSE, scale = "row",
                   Colv = NA ,ColSideLabs = "Group", RowSideLabs = "Genes at q<0.01", ColSideAnn = c("T=1,N=0, Progression", "T=1,N=0, No Progression","T=1, Node-Positive"), showRowDendro = TRUE,
                  verbose = TRUE, margins = c(7,6), cexRow = 0.19, cexCol = 0.3, ylab = "Genes at q<0.01 Between All Groups", cex = 3)
#Alternate colors
heatmap3::heatmap3(x = luad.rpkm,col=colorRamps::blue2red(50), balanceColor = TRUE, ColSideColors = groups, showColDendro = TRUE, scale = "row",
                   Colv = NA ,ColSideLabs = "Group", RowSideLabs = "Genes at q<0.01", ColSideAnn = c("T=1,N=0, Progression", "T=1,N=0, No Progression","T=1, Node-Positive"), showRowDendro = TRUE,
                   verbose = TRUE, margins = c(7,6), cexRow = 0.19, cexCol = 0.3, ylab = "Genes at q<0.01 Between All Groups", cex = 3)

title(main = "TCGA LUAD Scaled Log2 RPKM", sub = "TNM Staging at Minimum 24 Months", outer = FALSE, xlab = "Subject (TCGA Barcode)")
       
legend("topleft", legend = c("T=1, N=0, No Progression", "T=1, N=0, Progression", "T=1, Node-Positive"),  
       fill = c("skyblue1", "orange1", "mediumpurple1"), border=FALSE, bty="n", x.intersp = 0.5, cex=0.6, horiz = TRUE, inset = 0.1)


# Manta Ray Plots

workingset <- t(luad.rpkm)

groups <- match(rownames(workingset),get.files$dash.proj.id)
groups <- as.matrix(get.files$group[groups])
groups[which(groups[,1] == "T=1,N=0,Disease Free"),] <- "T=1,N=0,No Progression"

color <- as.matrix(groups)
color[which(color[,1] == "T=1,N=0,Progressed"),1] <- "orange1"  
color[which(color[,1] == "T=1,N=0,No Progression"),] <- "skyblue1"
color[which(color[,1] == "T=1,N>0"),] <- "mediumpurple1" 

workingset <- cbind(workingset,groups,color)
workingset <- `colnames<-`(workingset, c(colnames(workingset)[1:(length(colnames(workingset))-2)], "groups", "color"))
workingset <- as.data.frame(workingset)

ggplot(workingset, aes(x = groups, y = as.numeric(workingset$SERPINE2), fill = color)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "area") +
  theme(legend.position = "none") +
  ggplot2::ggtitle("SERPINE2 Expression in LUAD Dataset", subtitle = "TNM Staging Past 24 Months")+
  ggplot2::ylab("log 2 RPKM")+
  scale_fill_manual(values = c("mediumpurple1","orange1","skyblue" )) +
  geom_signif(comparisons = list(c("T=1,N=0,Progressed", "T=1,N=0,No Progression"), c("T=1,N>0", "T=1,N=0,No Progression")), map_signif_level = TRUE)

# Alternative, Jeff's pref of std box plot

qplot( x = groups , y = as.numeric(workingset$SERPINE2) , data = workingset, geom=c("boxplot","jitter") , aes(fill = groups), ylab = "Log2 Counts (RPKM)")+
  ggplot2::ggtitle("SERPINE2 Expression in LUAD Dataset", subtitle = "TNM Staging Past 24 Months") +
  scale_fill_manual(values = c("skyblue", "orange1", "mediumpurple1")) +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("T=1,N=0,Progressed", "T=1,N=0,No Progression"), c("T=1,N>0", "T=1,N=0,No Progression")), map_signif_level = TRUE)


# Diff Co-expression
source("https://bioconductor.org/biocLite.R")
biocLite("diffcoexp")
library(diffcoexp)

nonprogressed <- get.files$dash.proj.id[which(get.files$group == "T=1,N=0,Disease Free")] #  if using file read in, may need to use dec.proj.id
nonprogressed <- luad.rpkm[,match(nonprogressed, colnames(luad.rpkm))]

progressed <- get.files$dash.proj.id[which(get.files$group == "T=1,N=0,Progressed")] #  if using file read in, may need to use dec.proj.id
progressed <- luad.rpkm[,match(progressed, colnames(luad.rpkm))]

node.pos <- get.files$dash.proj.id[which(get.files$group == "T=1,N>0")] #  if using file read in, may need to use dec.proj.id
node.pos <- luad.rpkm[,match(node.pos, colnames(luad.rpkm))]

diff.coexp.n0.nonprog.vs.prog <- diffcoexp::coexpr(exprs.1 = nonprogressed, exprs.2 = progressed, r.method = "pearson", q.method = "BH")
diff.coexp.n0.nonprog.vs.prog <- plyr::arrange(diff.coexp.n0.nonprog.vs.prog, q.diffcor)
write.csv(diff.coexp.n0.nonprog.vs.prog, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/top.tag.diff.coexp.prog.v.nonprog.luad.rpkm.3cm.2yr.csv")

diff.coexp.n0.nonprog.vs.npos <- diffcoexp::coexpr(exprs.1 = nonprogressed, exprs.2 = node.pos, r.method = "pearson", q.method = "BH")
diff.coexp.n0.nonprog.vs.npos<-plyr::arrange(diff.coexp.n0.nonprog.vs.npos, q.diffcor)
write.csv(diff.coexp.n0.nonprog.vs.npos, "C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/diff.coexp.npos.v.nonprog.luad.rpkm.3cm.2yr.csv")

#GSEA
source("https://bioconductor.org/biocLite.R")
biocLite("SeqGSEA")




