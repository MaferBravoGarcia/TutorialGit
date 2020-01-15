#MaferBravoGarcía
#Ubicación
setwd("~/Microbiología/6sem/Genómica/BowtieAln")
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
#Libreria
library(ShortRead)
library(rtracklayer)
#Archivo
alnFile = "clean_reads.bwt"
aln = readAligned(alnFile, type="Bowtie")
aln
#Sread
sread(aln)
head(chromosome(aln))
table(strand(aln))
#Coexpresión: Que caigan en el mismo sitio
alnCov = as(coverage(aln),"RangedData")
alnCov
##Que al menos hayan caido 5 lecturas: Alta covertura
highCov = alnCov[alnCov$score >= 5,]
highCov
#Para exportar
export(highCov, "highCov.bed", format="bedGraph")
