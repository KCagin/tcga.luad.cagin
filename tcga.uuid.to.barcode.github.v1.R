#Borrowed from https://www.biostars.org/p/306400/
#Finds TCGA IDs from UUIDs

setwd("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicDataCommons")

library(GenomicDataCommons)
library(magrittr)
manifest <- read.csv("C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/MANIFEST_key.csv", header = TRUE, row.names=1)

file_uuids <- manifest$id

TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files() %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

res = TCGAtranslateID(file_uuids)
head(res)

write.csv(res,"C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/tcga.id.index.for.manifest.csv")

#Add to Manifest Key for convenience

res$file_id<-as.matrix(res$file_id)
manifest$id<-as.matrix(manifest$id)

file.index<-match(manifest$id, res$file_id)
res<-res[file.index,]

manifest$ProjID<-res$submitter_id
manifest$ProjID<-as.matrix(manifest$ProjID)
manifest$Txt_File<-manifest$filename
manifest<-`colnames<-`(manifest, c(colnames(manifest[,1:5]),"proj.id", "txt.file"))

write.csv(manifest,"C:/Users/kelly/Desktop/Biorepository/data_sets/LUAD Clinical/RNAseq LUAD 08.20.18/Read Counts/Revised 3cm Groups/manifest.key.csv")
