library(phytools)
library(seqinr)

data_folder_path <- "/media/gaurav/Data/Idea/Ortho_detection/"

enzymes <- readLines("enzymes.txt")
# enzymes <- grep(">",enzymes,value = T)
# enzymes <- sapply(1:length(enzymes),function(x) strsplit(enzymes[x],split=" ")[[1]][2])
# enzymes <- sapply(1:length(enzymes),function(x) paste(c(tolower(strsplit(enzymes[x],split="")[[1]][1]),strsplit(enzymes[x],split="")[[1]][2:length(strsplit(enzymes[x],split="")[[1]])]),collapse = ""))
# enz_seq <- read.fasta("enzymes.faa",seqonly = T,as.string = T,strip.desc = T)

###Homology detection

##first download all available reviewed sequences for the enzymes of interest and name the file as uniprot-EnzName.fasta

pa_matrix <- matrix(nrow = nrow(int_file),ncol = length(enzymes))
time_per_enz <- NULL
plot_rows <- round(length(enzymes)/4)
plot_cols <- round(length(enzymes)/plot_rows)
if(plot_rows*plot_cols<length(enzymes)) plot_rows <- plot_rows+1
par(mfrow=c(plot_rows,plot_cols),mai=rep(0.46,4),mgp=c(2,0.5,0)) ##these are the e-value cutoff plots.
noofspecies_eval <- NULL
for(k in 1:length(enzymes)) {
  print(enzymes[k])
  s_time <- proc.time()
  
  ##jackhmmer of some enz: some enzyme have only one entry of the enz seq in UniProt. These will require jackhmmer. name the seq file in this case uniprot-EnzName_jackhmmer.fasta
  
  if(file.exists(paste0(enzymes[k],"_jackhmmer.fasta"))) {
    system(paste0("jackhmmer --cpu 7 --noali --tblout ",enzymes[k],"_simple.out uniprot-",enzymes[k],"_jackhmmer.fasta ",data_folder_path,"all_1093species_comb_genomes.fasta"))
  } else {
    ##first align all the seqs from UniProt
    system(paste0("t_coffee ","uniprot-",enzymes[k],".fasta -mode quickaln >","uniprot-",enzymes[k],".aln"))
    ##convert .aln file to .stockholm for hmmer's sake
    system(paste0("t_coffee -other_pg seq_reformat -in uniprot-",enzymes[k],".aln -output stockholm_aln >uniprot-",enzymes[k],".stockholm"))
    ##build the HMM profile of the enzymes
    system(paste0("hmmbuild ",enzymes[k],".hmm uniprot-",enzymes[k],".stockholm"))
    ##search for the enzymes in the combined fasta file (make this combined fasta file of all faa files of orgs of interest and change the name of this file in the next line)
    system(paste0("hmmsearch --cpu 7 --tblout ",enzymes[k],"_simple.out ",enzymes[k],".hmm ",data_folder_path,"all_1093species_comb_genomes.fasta"))
  }

  
  test <- as.data.frame(read.table(paste0(enzymes[k],"_simple.out"),header = F,sep = "",skip = 3,fill = T,strip.white = T))
  
  ##determine the evalue cutoff based on the noofspecies vs bitscore cutoff curve
  
  V3_find <- as.character(test$V3[1])
  acc_find <- as.character(test$V1[which(test$V3==V3_find)]) ## cleanup of hmm output file
  writeLines(acc_find,paste0(enzymes[k],"_acc_find.txt"))
  
  ##next I look for the accessions of the hits in the individual faa files in the folder of interest. make sure the individual faa files are in the folder mentioned in the line below
  system(paste0("LC_ALL=C fgrep -F -f ",enzymes[k],"_acc_find.txt ",data_folder_path,"*.faa >",enzymes[k],"_hits.txt"))
  
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
  
  cutoff <- sum(rr$lengths[1:((which(rr$values==0)[which(rr$lengths[which(rr$values==0)]==max(rr$lengths[which(rr$values==0)]))][1]))]) ##choosing the most stable plateau. if the plateau is not very clear or wrongly chosen, you can set a manual cutoff. Note that the cutoff is the index number of the e-value vector
  
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
write.csv(pa_matrix,file = "enzymes_pa_matrix.csv")

##phylo_heatmap
final_tree <- read.tree("20161115_pruned_tree_cutoff_0_03.nwk")
int_file <- read.csv("20170222_int_file_tree_based.csv",row.names = 1)

rownames(pa_matrix) <- as.vector(int_file$ncbi_name)
final_tree$tip.label <- as.vector(int_file$ncbi_name)

families <- as.vector(int_file$families)

families_split <- split(1:length(final_tree$tip.label),families)
all_col<-toupper(c("#87003e","#ff6854","#00a653","#ff7166","#82ea9a","#b50066","#004285","#d6adff","#8be986","#dcd877","#34086a","#518300","#cc58c6","#ab9ae2","#ff778b","#d13393","#245900","#625100","#019a60","#773300","#b26393","#76b22e","#df6e25","#b7b824","#d27491","#df2f68","#9f1e00","#01ebc1","#1b4abb","#60006c","#0194f4","#e4c37b","#ff9beb"))
col_by_family <- NULL
for(i in 1:length(families_split)) {
  col_by_family[families_split[[i]]] <- all_col[i]
}

col_table <- data.frame(names(families_split),all_col)
colnames(col_table)[1] <- "families"

pdf(file = "phylo_heatmap.pdf",width = 15,height = 20)
phylo.heatmap(tree = final_tree,X = pa_matrix,fsize = c(0.1,2,0.1),colors = c(0,gray(0.6)),legend=F,split=c(0.6,0.4),lwd=0.5,mar=c(5,2,2,2))
tiplabels(tip = 1:length(final_tree$tip.label),pch = 15,col=col_by_family)
# nodelabels(node = c(1409,1430,1300,1685,1573,1309,1410,1097),pch=1,col=1,lwd=3,cex=4)
# nodelabels(node = c(1779,1687),pch=19,col=1,lwd=1,cex=4)
legend(-0.5,-0.003,legend = names(families_split),fill = all_col,ncol = 7,cex = 0.6,bty = "n",border = "white")
dev.off()
