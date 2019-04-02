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
branch_file = "/home/ninomoriaty/Nutstore Files/Nutstore/VAF_plot_beta/Mutation_sign/311252.NJtree.edges.txt"

# main function
# Usage: Mutational_Sigs_branch(maf_file, samples_vector)
Mutational_sigs_tree <- function(maf_file, branches){
  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get mutationalSigs-related  infomation
  dat.sample <- data.frame(as.character(maf_input[,ncol(maf_input)]), stringsAsFactors=FALSE)
  dat.chr <- data.frame(as.character(maf_input[,2]), stringsAsFactors=FALSE)
  dat.chr[,1] <- paste("chr", dat.chr[,1], sep="")
  dat.pos.start <- maf_input[,3]
  dat.pos.end <- maf_input[,4]
  dat.ref <- maf_input[,7]
  dat.alt <- maf_input[,9]
  dat.num <- 1:length(dat.alt)
  mut.sig.ref <- data.frame(dat.num, dat.sample, dat.chr, dat.pos.start, dat.pos.end, dat.ref, dat.alt)
  colnames(mut.sig.ref) <- c("ID", "Sample", "chr", "pos", "pos_end", "ref", "alt")
  
  patientID = strsplit(as.character(maf_input$Tumor_Sample_Barcode[1]), "-")[[1]][1]
  ID_prefix = paste(" ", patientID, "-", sep = "")
  
  branch_input <- gsub("\xa1\xc9", ID_prefix, readLines(branch_file))
  branches <- strsplit(as.character(paste(patientID, "-", branch_input, sep = "")), split=" ")
  
  # generate sets of different branches
  for (branch_counter in length(branches):1){
    branch <- Filter(Negate(is.na), branches[[branch_counter]])
    mut.branch <- mut.sig.ref[which(mut.sig.ref$Sample %in% branch), ]
    for (tsb in branch){
      mut.tsb <- mut.sig.ref[which(mut.sig.ref$Sample %in% tsb), ]
      mut.branch <- match_df(mut.branch, mut.tsb, on = c("chr", "pos", "pos_end", "ref", "alt"))
    }
    branch_name <- paste(branch, collapse = "+")
    mut.sig.ref[which(mut.sig.ref[,1] %in% mut.branch[,1]), 2] <- branch_name
    Mutational_sigs_branch(mut.sig.ref, branch_name)
  }
}

  
# Weight mutational Signature of each branch
Mutational_sigs_branch <- function(mut.branch, branch_name){
  # deconstructSigs
  sigs.input <- mut.to.sigs.input(mut.ref = mut.branch, 
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
  print(paste(branch_name, ":", sigs.max))
}


# Confirm sets of mutation

# Mutation_sets <- function(mut.sig.ref, branches){
#   # generate branch name
#   for (branch_counter in 1:length(branches)){
#     branch <- Filter(Negate(is.na), branches[,branch_counter])
#     mut.branch <- mut.sig.ref[which(mut.sig.ref$Sample %in% branch), ]
#     for (tsb in branch){
#       mut.tsb <- mut.sig.ref[which(mut.sig.ref$Sample %in% tsb), ]
#       mut.branch <- match_df(mut.branch, mut.tsb, on = c("chr", "pos", "pos_end", "ref", "alt"))
#     }
#     branch_name <- paste(branch, collapse = "+")
#     Mutational_sigs_branch(mut.branch, branch_name)
#     mut.sig.ref[which(mut.sig.ref[,1] %in% mut.branch[,1]), 2] <- branch_name
#   }
#   mut.sig.ref
# }

# test branches

# branch11 = data.frame(c("311252-S"))
# branch10 = data.frame(c("311252-TP2-2"))
# branch9 = data.frame(c("311252-TC2"))
# branch8 = data.frame(c("311252-TP1"))
# branch7 = data.frame(c("311252-TC1"))
# branch6 = data.frame(c("311252-V"))
# branch5 = data.frame(c("311252-TP2-2", "311252-S"))
# branch4 = data.frame(c("311252-TP1", "311252-TC2"))
# branch3 = data.frame(c("311252-TP1", "311252-TC2", "311252-V"))
# branch2 = data.frame(c("311252-TP1", "311252-TC2", "311252-TP2-2", "311252-V", "311252-S"))
# branch1 = data.frame(c("311252-TC1", "311252-TP1", "311252-TC2", "311252-TP2-2", "311252-V", "311252-S"))
# branches =  rbind.fill(branch1, branch2, branch3, branch4, branch5, branch6, branch7, branch8, branch9, branch10, branch11)

