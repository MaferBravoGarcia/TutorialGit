#MaríaFernanda
#Ubicación
setwd("~/Microbiología/6sem/Genómica/Anotación de alineamiento")
#Librerias
library(Biobase)
library(IRanges)
library(ShortRead)
library(GenomicRanges)
library(biomaRt)
#Ejemplito de la práctica 
x = IRanges(start=c(11,35,40), end=c(20,50,63))
x
##2 exones se traslapan por lo tanto se reduce a 2 exones 
exons = reduce(x)
##Lecturas obtenidad del alineamiento de rango 20
reads = IRanges(start=c(1,21,30,50,80), width=20)
reads
##En qué exones caen en las lecturas
countOverlaps(exons, reads)

##Ejercicio con Humanitos 
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Los atributos que queremos obtener
attribs = c("ensembl_transcript_id", "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end")
# Obtener la información de la bases de datos
transInfo = getBM(attributes=attribs, filters="ensembl_transcript_id", values="ENST00000337138", mart=mart)
transInfo

#Cromosoma 14
coords = read.table("human_chr14_annotation.tab", header=TRUE, sep="\t", check.names=FALSE)
head(coords)

##
exonGRs = GRanges(
seqnames = Rle(factor(coords$chromosome_name)),
ranges   = IRanges(start=coords$exon_chrom_start, end=coords$exon_chrom_end),
strand   = Rle(factor(coords$strand)),
biotype  = as.character(coords$transcript_biotype),
id       = as.character(coords$ensembl_transcript_id)
)

exonGRs

##Ejercicio
#Cómo le harían para quedarse exclusivamente con las anotaciones de "protein_coding"?
#y solamente aquellas anotaciones de la cadena "-"?
####No me salió 8887?

Prot=exonGRs[exonGRs$biotype == "Protrein_coding",]
Menos= Prot[strand(Prot)=="-",]

NegProt=exonGRs[exonGRs$biotype == "Protrein_coding" & strand(exonGRs)=="-",]

#Ubicación
setwd("~/Microbiología/6sem/Genómica/BowtieAln/")
#
aln = readAligned("clean_reads.bwt", type="Bowtie")
aln
# Conertir el aln a la base de datos 
alnGRs = as(aln, "GRanges")
alnGRs
#Caracteristicas del alineamiento vs a qué exones pega de la base de datos
values(exonGRs)$counts = countOverlaps(exonGRs, alnGRs)
values(exonGRs)
#
typeCounts = aggregate(values(exonGRs)$counts, by=list("biotype"=values(exonGRs)$biotype), sum)
head(typeCounts)
#
transCounts = aggregate(values(exonGRs)$counts, by=list("id"=values(exonGRs)$id), sum)
head(transCounts)
##Gráfica resumen 
minCount = 1000
typeCountsHigh = typeCounts[typeCounts$x > minCount,]
typeCountsHigh = typeCountsHigh[order(typeCountsHigh$x),]
typeCountsHigh = rbind(data.frame("biotype"="other",
                                  "x"=sum(typeCounts$x[typeCounts$x <= minCount])),
                       typeCountsHigh)

pie(typeCountsHigh$x, labels=typeCountsHigh$biotype,
    main="Number of aligned reads per biotype")