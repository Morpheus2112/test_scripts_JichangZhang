library(Biostrings)
library(ggplot2)
library(ggseqlogo)

###### first question ######
###### 1. read reference fasta, bed file and amplicon fasta file ########
ref <- readDNAStringSet("genome.fa")
a1 <- readDNAStringSet("amplicon_1.fa")
a2 <- readDNAStringSet("amplicon_2.fa")
a3 <- readDNAStringSet("amplicon_3.fa")
pos <- 
  as.data.frame(read.table("amplicon_coordinates.bed",header = FALSE, sep="\t",
                           stringsAsFactors=FALSE, quote=""))

###### 2. select reference seqs ######

ref1 <- toString(subseq(ref, start = pos[1,2], end = pos[1,3]))
ref2 <- toString(subseq(ref, start = pos[2,2], end = pos[2,3]))
ref3 <- toString(subseq(ref, start = pos[3,2], end = pos[3,3]))

###### 3. calculate consensus matrix for amplicons ###### 
report <- data.frame(gene_name =c("amplicon1","amplicon2","amplicon3"),
                     mutation = NA,
                     freq = NA,
                     number = NA,
                     ref_seq = NA,
                     amp_seq = NA)

tab1 <- consensusMatrix(a1, as.prob=TRUE)
which(tab1<1 & tab1>0, arr.ind= TRUE)
tab1[,58]
report$mutation[1] <- "C->T"
report$freq[1] <- tab1["T",58]
report$number[1] <- tab1["T",58]*200
report$ref_seq[1] <- ref1
substr(ref1, 58, 58) <- "T"
report$amp_seq[1] <- ref1

tab2 <- consensusMatrix(a2, as.prob=TRUE)
which(tab2<1 & tab2>0, arr.ind= TRUE)
tab2[,4]
report$mutation[2] <- "C->T"
report$freq[2] <- tab2["T",4]
report$number[2] <- tab2["T",4]*200
report$ref_seq[2] <- ref2
substr(ref2, 4, 4) <- "T"
report$amp_seq[2] <- ref2

tab3 <- consensusMatrix(a3, as.prob=TRUE)
which(tab3<1 & tab3>0, arr.ind= TRUE)
tab3[,10]
report$mutation[3] <- "T->C"
report$freq[3] <- tab3["C",10]
report$number[3] <- tab3["C",10]*200
report$ref_seq[3] <- ref3
substr(ref3, 10, 10) <- "C"
report$amp_seq[3] <- ref3





###### check mutation and position are correct######
mutations <- function(str1, str2) {
  str1vec <- unlist(strsplit(str1, ""))
  str2vec <- unlist(strsplit(str2, ""))
  iMut <- (1:length(str1vec))[str1vec != str2vec]
  return(paste0(str1vec[iMut], iMut, str2vec[iMut]))
}

mutations(report$ref_seq[1],report$amp_seq[1])
mutations(report$ref_seq[2],report$amp_seq[2])
mutations(report$ref_seq[3],report$amp_seq[3])

###### same result to a csv file ######
write.csv(report, "report1.csv")
###### use sequence logo to illustrate the mutation 
ggseqlogo( as.data.frame(a1))
ggseqlogo( as.data.frame(a2))
ggseqlogo( as.data.frame(a3))



###### second question ######
a3 <- readDNAStringSet("amplicon_3.test_2.fa")
ref3 <- toString(subseq(ref, start = pos[3,2], end = pos[3,3]))
ggseqlogo( as.data.frame(a3))
# check logo graph, find two mutations, 10 and 55
report2 <- data.frame(gene_name =c("amplicon3","amplicon3"),
                     mutation = NA,
                     freq = NA,
                     number = NA,
                     ref_seq = NA,
                     amp_seq = NA)


tab3 <- consensusMatrix(a3, as.prob=TRUE)
which(tab3<1 & tab3>0, arr.ind= TRUE)

tab3[,10]
report2$mutation[1] <- "T->C"
report2$freq[1] <- tab3["C",10]
report2$number[1] <- tab3["C",10]*200

tab3[,55]
report2$mutation[2] <- "G->A"
report2$freq[2] <- tab3["A",55]
report2$number[2] <- tab3["A",55]*200

report2$ref_seq[1] <- ref3
report2$ref_seq[2] <- ref3
substr(ref3, 10, 10) <- "C"
substr(ref3, 55, 55) <- "A"
report2$amp_seq[1] <- ref3
report2$amp_seq[2] <- ref3

##### check mutation recorded correct ######
mutations(report2$ref_seq[1],report2$amp_seq[1])

### export csv #####
write.csv(report2, "report2.csv")


