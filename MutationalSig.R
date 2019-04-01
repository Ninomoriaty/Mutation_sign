# packages
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
library(deconstructSigs)
library(plyr)

maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"

branch11 = data.frame(c("311252-S"))
branch10 = data.frame(c("311252-TP2-2"))
branch9 = data.frame(c("311252-TC2"))
branch8 = data.frame(c("311252-TP1"))
branch7 = data.frame(c("311252-TC1"))
branch6 = data.frame(c("311252-V"))
branch5 = data.frame(c("311252-TP2-2", "311252-S"))
branch4 = data.frame(c("311252-TP1", "311252-TC2"))
branch3 = data.frame(c("311252-TP1", "311252-TC2", "311252-V"))
branch2 = data.frame(c("311252-TP1", "311252-TC2", "311252-TP2-2", "311252-V", "311252-S"))
branch1 = data.frame(c("311252-TC1", "311252-TP1", "311252-TC2", "311252-TP2-2", "311252-V", "311252-S"))
branches =  rbind.fill(branch1, branch2, branch3, branch4, branch5, branch6, branch7, branch8, branch9, branch10, branch11)

# main function
# Usage: Mutational_Sigs_branch(maf_file, samples_vector)
Mutational_sigs_tree <- function(maf_file, branches){

  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get mutationalSigs-related  infomation
  dat.sample <- data.frame(as.character(maf_input[,ncol(maf_input)]), stringsAsFactors=FALSE)
  dat.chr <- data.frame(as.character(maf_input[,2]), stringsAsFactors=FALSE)
  dat.chr[,1] <- paste("chr", dat.chr[,1], sep="")
  dat.pos <- maf_input[,3]
  dat.ref <- maf_input[,7]
  dat.alt <- maf_input[,9]
  mut.sig.ref <- data.frame(dat.sample, dat.chr, dat.pos, dat.ref, dat.alt)
  colnames(mut.sig.ref) <- c("Sample", "chr", "pos", "ref", "alt")
  

}

# Confirm sets of mutation
Mutation_sets <- function(mut.sig.ref, branches){
  # generate branch name
  branch_name <- paste(branch, collapse = "+")
  dat.sample[which(dat.sample[,1] %in% branch), 1] <- branch_name
}
  
# Weight mutational Signature of each branch
Mutational_sigs_branch <- function(mut.sets.sig.ref, branch){
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



