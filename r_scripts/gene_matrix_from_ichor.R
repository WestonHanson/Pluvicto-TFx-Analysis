# if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#   BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# }

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(data.table)
library(GenomicRanges)
library(dplyr)
options(stringsAsFactors=F)
library(GenomicRanges)
#source("/home/unix/gavinha/software/code/git/TitanCNA/R/utils.R")
source("/Volumes/ha_g/projects/Collaborations/mskAncestry/WGS/titan_latest/TitanCNA/R/utils.R")
args <- commandArgs(TRUE)
options(stringsAsFactors=F, scipen=999)

setwd('/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/test')


geneList <- '/Volumes/ha_g/projects/Collaborations/Lung_Lavage_Nair/CNV/gene_freq/test_doce/GRCh38.protein.coding.genes_unique.txt'
centromeres <- fread("/Volumes/ha_g/projects/Collaborations/Lung_Lavage_Nair/CNV/ichor_run/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt")
filterLen <- 1000
method <- 'common' #{'common','severity'}; default is most 'common' state
headerType <- 'Gene' #{'Gene','chrPosn'}


## Triplets
ichorFilePaths <- '/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/data-tables/gene_matrix_from_ichor/radium223_ichor_file_list.csv'
idFile <- '/Volumes/ha_g/projects/Collaborations/Lung_Lavage_Nair/CNV/gene_freq_lavage/matched/input_plasma/map_file.txt'
outroot <- '/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/data-tables/gene_matrix_from_ichor/radium223_ichor'

# seq.info <- Seqinfo(genome = "hg38")
# seq.info <- keepStandardChromosomes(seqinfo, pruning.mode = "coarse")
# seq.info <- seqinfo[paste0("chr", c(1:22, "X"))]

seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
seq.info <- seqinfo_hg38[paste0("chr", c(1:22, "X", "Y","M"))]

seqinfo_hg38 <- Seqinfo(genome = "hg38")
saveRDS(seqinfo_hg38, file = "seqinfo_hg38.rds")
seqinfo_hg38_loaded <- readRDS("seqinfo_hg38.rds")

filterSnp <- 0

#chrs <- as.character(c(1:22,"X"))
chrs <- paste0("chr", c(1:22, "X"))

## load in ulp files and params ##
df <- read.csv(ichorFilePaths)
paths <- df$curated_solution_path
files <- data.frame(
  basename = sub('.cna.seg', '', basename(paths)),
  path = paths,
  stringsAsFactors = FALSE
)
files <- as.matrix(files)

## load list of ids if specified ##
# if (idFile != "0"){
# 	idList <- read.delim(idFile, header=F, as.is=TRUE)[, 1]
# 	idToUse <- NULL
# 	for (i in 1:length(idList)){
# 		ind <- which(grepl(idList[i], names(files)))
# 		if (length(ind) > 0){
# 			idToUse <- c(idToUse, names(files)[ind])
# 		}
# 	}
# 	files <- files[idToUse]
# }
# 

## load in annotation file ##
genes <- fread(geneList)
#colnames(genes)[1:6] <- c("EnsgId","Gene","Chr","Start","Stop","Strand")#,"Band","TxnCount","Description","Status")
colnames(genes)[4:6] <- c("Chr","Start","End")
#genes[, Chr := gsub("chr", "", Chr)]
genes$Chr = paste0("chr",genes$Chr)
# genes[Chr==23, Chr := "X"]
# genes[Chr==24, Chr := "Y"]
# genes[Chr==25, Chr := "MT"]
genes[, chrPosn := paste0(Chr,":", Start, "-", End)]
if (!headerType %in% colnames(genes)){
   genes[[headerType]] <- genes$chrPosn
}
genes[, Chr := factor(Chr, levels = chrs)]
genes <- genes[!is.na(Chr)]

stateKey <- array(0:14); 
names(stateKey) <- c('HOMD','HETD','NEUT','GAIN','AMP','HLAMP', paste0(rep("HLAMP", 8), 2:10))
severity <- array(c(5,4,1,2,4,5,rep(5,9))); 
names(severity) <- c('HOMD','HETD','NEUT','GAIN','AMP','HLAMP', paste0(rep("HLAMP", 8), 2:10))

save.image(paste(outroot,".RData",sep=""))
## input: lohHits = loh rows that overlap region of interest
# start = start coordinate of region of interest
# end = end coordinate of region of interest
getOverlapLength <- function(cnHits, start, end){
	coords <- cbind(cnHits[, c("start","end")], as.numeric(start), as.numeric(end))
	coordsSort <- t(apply(coords, 1, sort))
	dist <- coordsSort[, 3] - coordsSort[, 2] + 1
	return(dist)
}


#Input: states is an array of state names (e.g. (DLOH,NLOH,...,))
#Uses global variable "severity"
getMostSevereState <- function(states){	
	severityValue <- 0
	severeState <- states[i]
	for (i in states){
		if (severity[i] > severityValue){
			severeState <- i
			severityValue <- severity[i]
		}
	}
	return(severeState)
}

# output the matrix to file
#Input: Matrix to output; output file name
writeMatrixToFile <- function(mat,outfile){
	outMat <- cbind(rownames(mat),mat)
	if (!is.null(colnames(outMat))){
		colnames(outMat)[1] <- "Sample"
	}
	write.table(outMat,file=outfile,row.names=F,col.names=T,quote=F,na="NaN",sep="\t")
}

#ids <- gsub(".seg.txt","",basename(files))
ids <- files[,1]

numSamples <- nrow(files)
numGenes <- nrow(genes)
geneCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(files[,1], genes[, get(headerType)]))
geneLogRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(files[,1], genes[, get(headerType)]))
geneLenmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(files[,1], genes[, get(headerType)])) 
geneSCmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(files[,1], genes[, get(headerType)])) 


for (i in 1:numSamples){
	caseId <- ids[i] #gsub(".cna.seg","",basename(files[i,2]))
	cat("Analyzing sample:\t",caseId,"\n")
	cat("Reading file ",files[i,1],"\n")
	#files[i, 3] = gsub(".seg.txt",".cna.seg",files[i, 2])
	#files[i, 3] <- gsub("\\.seg.txt$", ".cna.seg", files[i, 2])
	#CN
	cn <- fread(sub('fh/fast/', '/Volumes/', files[i, 2]))
	cn <- cn[chr %in% chrs]
	#cn = data.table(cn)
  colnames(cn)[c(1,2,3)] <- c("Chromosome","Start","End")
  cn$chr <- as.character(cn$Chromosome)
	cn <- extendSegments(cn, removeCentromeres = TRUE, centromeres = centromeres, 
											extendToTelomeres = FALSE, seqInfo = seqinfo_hg38_loaded)
	cn$Chromosome <- cn$chr
	cn[, Length.bp := End - Start + 1]
	cn$chr <- NULL
	## filter by length threshold ##
	cn <- cn[Length.bp>=filterLen, ]

  geneCN <- getOverlap(x=genes, y=cn, type="within", colToReturn=paste0(caseId, ".copy.number"))
	geneLogR <- getOverlap(x=genes, y=cn, type="within", colToReturn=paste0(caseId,".logR"))
  geneLen <- getOverlap(x=genes, y=cn, type="within", colToReturn="Length.bp")
  geneSC <- getOverlap(x=genes, y=cn, type="within", colToReturn=paste0(caseId,".subclone.status"))
  
	#build matrices
	geneCNmat[caseId,] <- geneCN #call matrix
	geneLogRmat[caseId,] <- geneLogR   #logR matrix
  geneLenmat[caseId, ] <- geneLen
	
}	
save(geneCNmat, geneLogRmat, geneLenmat, geneList, geneSCmat, chrs, file = paste(outroot,".RData",sep=""))

geneCNmat[is.na(geneCNmat)] <- NaN
geneLogRmat[is.na(geneLogRmat)] <- NaN
geneLenmat[is.na(geneLenmat)] <- NaN
geneSCmat[is.na(geneSCmat) <- NaN]

# output the call matrix
outfile <- paste(outroot,"_geneCN.txt",sep="")
writeMatrixToFile(geneCNmat,outfile)

# output the logR matrix
outfile <- paste(outroot,"_geneLogR.txt",sep="")
writeMatrixToFile(geneLogRmat,outfile)

# output the Length matrix
outfile <- paste(outroot,"_geneLength.txt",sep="")
writeMatrixToFile(geneLenmat,outfile)

# output the subclonal status matrix
outfile <- paste(outroot,"_geneSCstatus.txt",sep="")
writeMatrixToFile(geneSCmat,outfile)


