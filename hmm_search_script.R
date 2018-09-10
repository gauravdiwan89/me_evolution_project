library(phytools)
library(seqinr)

final_tree <- read.tree("20161115_pruned_tree_cutoff_0_03.nwk")
int_file <- read.csv("20170222_int_file_tree_based.csv",row.names = 1)
trna_file <- read.csv("20161122_int_file_trnagcn_complete.csv",row.names = 1)

enzymes <- readLines("enzymes.faa")
enzymes <- grep(">",enzymes,value = T)
enzymes <- sapply(1:length(enzymes),function(x) strsplit(enzymes[x],split=" ")[[1]][2])
enzymes <- sapply(1:length(enzymes),function(x) paste(c(tolower(strsplit(enzymes[x],split="")[[1]][1]),strsplit(enzymes[x],split="")[[1]][2:length(strsplit(enzymes[x],split="")[[1]])]),collapse = ""))
enz_seq <- read.fasta("enzymes.faa",seqonly = T,as.string = T,strip.desc = T)

###Homology detection

##first download all available reviewed sequences for the enzymes of interest and name the file as uniprot-EnzName.fasta

pa_matrix <- matrix(nrow = nrow(int_file),ncol = length(enzymes))
time_per_enz <- NULL
par(mfrow=c(2,4),mai=rep(0.46,4),mgp=c(2,0.5,0)) ##these are the e-value cutoff plots. change mfrow acc to no of enz analysed
noofspecies_eval <- NULL
for(k in 1:length(enzymes)) {
  print(enzymes[k])
  s_time <- proc.time()
  
  ##jackhmmer of some enz: some enzyme have only one entry of the enz seq in UniProt. These will require jackhmmer. name the seq file in this case EnzName_jackhmmer.fasta
  
  if(file.exists(paste0(enzymes[k],"_jackhmmer.fasta"))) {
    system(paste0("jackhmmer --cpu 7 --noali --tblout /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],"_simple.out /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],"_jackhmmer.fasta /media/gaurav/Data/Idea/Ortho_detection/all_1093species_comb_genomes.fasta"))
  } else {
    ##first align all the seqs from UniProt
    system(paste0("t_coffee /media/gaurav/Data/Idea/Ortho_detection/uniprot-",enzymes[k],".fasta -mode quickaln >/media/gaurav/Data/Idea/Ortho_detection/uniprot-",enzymes[k],".aln"))
    ##convert .aln file to .stockholm for hmmer's sake
    system(paste0("t_coffee -other_pg seq_reformat -in /media/gaurav/Data/Idea/Ortho_detection/uniprot-",enzymes[k],".aln -output stockholm_aln >/media/gaurav/Data/Idea/Ortho_detection/uniprot-",enzymes[k],".stockholm"))
    ##build the HMM profile of the enzymes
    system(paste0("hmmbuild /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],".hmm /media/gaurav/Data/Idea/Ortho_detection/uniprot-",enzymes[k],".stockholm"))
    ##search for the enzymes in the combined fasta file (make this combined fasta file of all faa files of orgs of interest and change the name of this file in the next line)
    system(paste0("hmmsearch --cpu 7 --tblout /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],"_simple.out /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],".hmm /media/gaurav/Data/Idea/Ortho_detection/all_1093species_comb_genomes.fasta"))
  }

  
  test <- as.data.frame(read.table(paste0(enzymes[k],"_simple.out"),header = F,sep = "",skip = 3,fill = T,strip.white = T))
  
  ##determine the evalue cutoff based on the noofspecies vs bitscore cutoff curve
  
  V3_find <- as.character(test$V3[1])
  acc_find <- as.character(test$V1[which(test$V3==V3_find)]) ## cleanup of hmm output file
  writeLines(acc_find,paste0(enzymes[k],"_acc_find.txt"))
  
  ##next I look for the accessions of the hits in the individual faa files in the folder of interest. make sure the individual faa files are in the folder mentioned in the line below
  system(paste0("LC_ALL=C fgrep -F -f /media/gaurav/Data/Idea/Ortho_detection/",enzymes[k],"_acc_find.txt /media/gaurav/Data/Idea/Ortho_detection/*.faa >",enzymes[k],"_hits.txt"))
  
  where_hits <- readLines(paste0(enzymes[k],"_hits.txt"))
  
  hits_files <- sapply(1:length(where_hits),function(x) strsplit(where_hits[x],split=":>")[[1]][1])
  hits_acc <- sapply(1:length(where_hits),function(x) strsplit(where_hits[x],split=":>")[[1]][2])
  hits_acc <- sapply(1:length(hits_acc),function(x) strsplit(hits_acc[x],split=" ")[[1]][1])
  
  ##find range of e-values to check
  e_values <- as.numeric(as.character(test$V5[match(hits_acc,test$V1)]))
  e_limits <- range(na.omit(round(log10(e_values))))
  if("-Inf" %in% e_limits) {
    e_cuts <- 10^(0:-100)
  } else e_cuts <- 10^(e_limits[1]:e_limits[2])
  
  ## find no of species that arise after setting cutoffs
  noofspecies <- NULL
  for(i in 1:length(e_cuts)) {
    cutoff_hits <- hits_acc[which(e_values<e_cuts[i])]
    cutoff_files <- hits_files[which(e_values<e_cuts[i])]
    which_file_has_what <- split(cutoff_hits,as.factor(cutoff_files))
    noofspecies[i] <- length(which_file_has_what)
  }
  diff_spp <- diff(noofspecies)
  
  rr <- rle(diff_spp)
  
  cutoff <- sum(rr$lengths[1:((which(rr$values==0)[which(rr$lengths[which(rr$values==0)]==max(rr$lengths[which(rr$values==0)]))][1]))]) ##choosing the most stable plateau. if the plateau is in not very clear or wrongly chosen, you can set a manual cutoff. Note that the cutoff is the index number of the e-value vector
  
  plot(e_cuts,noofspecies,xlim=c(max(e_cuts),min(e_cuts)),log="x",main=enzymes[k],xlab="E-value cutoff",ylab="Number of Species")
  abline(v=e_cuts[cutoff],col=2,lwd=2)
  text(max(e_cuts)/100000,noofspecies[1]+100,labels = e_cuts[cutoff],cex=1)
  noofspecies_eval[k] <- noofspecies[cutoff]
  
  cutoff_hits <- hits_acc[which(e_values<e_cuts[cutoff])]
  cutoff_files <- hits_files[which(e_values<e_cuts[cutoff])]
  which_file_has_what <- split(cutoff_hits,as.factor(cutoff_files))
  
  test <- unlist(which_file_has_what)
  writeLines(test,con=paste0(enzymes[k],"_ortholog_acc.txt"))
  
  a <- sapply(1:length(which_file_has_what),function(x) strsplit(names(which_file_has_what)[x],split="detection/")[[1]][2])
  presence_acc <- sapply(1:length(a),function(x) strsplit(a[x],split="[.]")[[1]][1])
  
  presence_id <- sapply(1:length(presence_acc),function(x) as.character(int_file$ncbi_name[grep(presence_acc[x],as.vector(int_file$ncbi_acc))[1]]))
  writeLines(presence_id,con = paste0(enzymes[k],"_orgs_with_presence_evalue_cutoff.txt"))
  
  these_bact_have <- match(presence_id,int_file$ncbi_name)
  
  enz_pa <- rep(0,nrow(int_file))
  enz_pa[these_bact_have] <- 1
  
  pa_matrix[,k] <- enz_pa
  
  time_taken <- proc.time()-s_time
  time_per_enz[k] <- time_taken[3]
}

colnames(pa_matrix) <- enzymes
rownames(pa_matrix) <- as.vector(int_file$phylotree_tip_label)
write.csv(pa_matrix,file = "20180126_int_file_enzymes_pa_evalue_cutoff.csv")
