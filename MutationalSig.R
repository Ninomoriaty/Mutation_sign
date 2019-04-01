
# packages
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
library(deconstructSigs)

maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
basic_sample <- sample.mut.ref

# read .maf file
maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
# get mutationalSigs-related  infomation
dat.sample <- data.frame(maf_input[,ncol(maf_input)])
dat.chr <- data.frame(as.character(maf_input[,2]), stringsAsFactors=FALSE)
dat.chr[,1] <- paste("chr", dat.chr[,1], sep="")
dat.pos <- maf_input[,3]
dat.ref <- maf_input[,7]
dat.alt <- maf_input[,9]
mut.sig.ref <- data.frame(dat.sample, dat.chr, dat.pos, dat.ref, dat.alt)
colnames(mut.sig.ref) <- c("Sample", "chr", "pos", "ref", "alt")

# deconstructSigs
sigs.input <- mut.to.sigs.input(mut.ref = mut.sig.ref, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

sigs.which = whichSignatures(tumor.ref = sigs.input, 
                       signatures.ref = signatures.nature2013, 
                       sample.id = "311252-S",
                       contexts.needed = TRUE)

plotSignatures(sigs.which, sub = 'example')





