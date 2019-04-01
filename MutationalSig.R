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
branch = c("311252-TC1", "311252-TP1", "311252-TC2", "311252-TP2-2", "311252-V", "311252-S")

# main function
Mutational_Sigs_branch <- function(maf_file, branch){
  # generate branch name
  branch_name <- paste(branch, collapse = "+")
  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get mutationalSigs-related  infomation
  dat.sample <- data.frame(as.character(maf_input[,ncol(maf_input)]), stringsAsFactors=FALSE)
  dat.sample[which(dat.sample[,1] %in% branch), 1] <- branch_name
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
  
  sigs.which <- whichSignatures(tumor.ref = sigs.input, 
                                signatures.ref = signatures.nature2013, 
                                sample.id = branch_name,
                                contexts.needed = TRUE)
  
  sigs.max <- colnames(sigs.which[["weights"]][which.max(sigs.which[["weights"]])])
  sigs.max
}

# basic_sample <- sample.mut.ref

plotSignatures(sigs.which, sub = 'example')




