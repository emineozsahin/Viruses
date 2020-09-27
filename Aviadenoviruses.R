##+++++++++++++++++++++++++++++++++++++++++++++++++++
#Emine Ozsahin
#Date: February 19, 2020
setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/Aviadenoviruses")

##+++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries ----
##+++++++++++++++++++++++++++++++++++++++++++++++++++

#Data obtained from NCBI
#This search was made on December 7th 2019, at 5:17pm (the run for this search)

# remotes::install_github("npcooley/FindHomology")
# install.packages("remotes")
# #BiocManager::install(c("DECIPHER"))
suppressMessages(library(DECIPHER))
install.packages("DECIPHER")
library(DECIPHER)
library(rentrez)
library(stringr)
library(Biostrings)
library(dplyr)
library(seqinr)
library(FindHomology)
library(ade4) #newick
library(ape)
library(dendextend)
library(stats)
#install.packages("genoPlotR")
library(genoPlotR)
library(remotes)


##+++++++++++++++++++++++++++++++++++++++++++++++++++
# Data Retrieve ----
##+++++++++++++++++++++++++++++++++++++++++++++++++++

#Aviadenoviruses <- entrez_search(db = "nuccore", term  = "Aviadenovirus[All Fields]) AND refseq[filter]")

# names(Aviadenoviruses)

# Aviadeno_complete_fetch <- entrez_fetch(db = "nuccore", id = Aviadenoviruses$ids, rettype = "fasta")

# write(Aviadeno_complete_fetch, "Aviadenoviruses_complete_genome_fetch.fasta", sep = "\n")

# stringSet_complete <- readDNAStringSet("Aviadenoviruses_complete_genome_fetch.fasta") 
# 
# summary(stringSet_complete)
# 
# names(stringSet_complete)
# 
# df_Aviadeno_complate <- data.frame(Original_Title = names(stringSet_complete), Sequence = paste(stringSet_complete)) 
# 
# df_Aviadeno_complate$Original_Title
# 
# df_Aviadeno_complate$Accession_Number <- word (df_Aviadeno_complate$Original_Title, 1L)
# 
# df_Aviadeno_complate$Species_names <- word (df_Aviadeno_complate$Original_Title, 2L, 4L)
# 
# df_Aviadeno_complate$Species_names <- str_replace(df_Aviadeno_complate$Species_names, "Fowl adenovirus ", "FAdV-")
# 
# df_Aviadeno_complate$Species_names <- str_replace(df_Aviadeno_complate$Species_names, ",", "")
# 
# df_Aviadeno_complate <- df_Aviadeno_complate[, c("Accession_Number", "Species_names", "Sequence", "Original_Title")]
# 
# #I could have made the search with accession numbers but I did not know when I made this seach that I could also use fttp adresses. To make the samples that I used same I filtered samples with their accession names. 
# 
# df_Aviadeno_complate <- df_Aviadeno_complate %>%
#   filter(str_detect(Species_names, "FAdV")) %>%
#   filter(Accession_Number == "NC_001720.1" | 
#          Accession_Number == "NC_021221.1" | 
#          Accession_Number == "NC_015323.1"| 
#          Accession_Number == "NC_000899.1" | 
#          Accession_Number == "NC_038332.1")
# 
# df_Aviadeno_complate$Accession_Number
# 
# df_Aviadeno_complate$Species_names
# #"FAdV-6" "FAdV-5" "FAdV-C" "FAdV-D" "FAdV-A"
# 
# df_Aviadeno_complate$Species_names <- df_Aviadeno_complate$Species_names%>%
#   str_replace_all(c("FAdV-6" = "FAdV-E", "FAdV-5" = "FAdV-B"))
# #"FAdV-E" "FAdV-B" "FAdV-C" "FAdV-D" "FAdV-A"
# 
# df_Aviadeno_complate$Species_names
# 
# #There are only one sample for each species as I filtered them with their accession numbers. But they were not in order, therefore I sampled them to make them in order. 
# 
# df_Aviadeno_complate <- df_Aviadeno_complate %>%
#   group_by(Species_names) %>%
#   sample_n(1)
# 
# dim(df_Aviadeno_complate) #5 4
# 
# hist(str_length(df_Aviadeno_complate$Sequence), main = "Distribution of complete genome lengths of FAdVs", ylab = "Frequency", xlab = "Base pairs")
# 
# dev.copy2pdf(file = "histogram.png")
# 
# 
# 
# write.fasta(as.list(df_Aviadeno_complate$Sequence), df_Aviadeno_complate$Species_names, "FAdVs_modified_complete_genome.fasta")

# Later I realized that I could have also use ftp adresses for syntheny alignment of the genomes as I made character vectors of ftp addresses for genomic GFF files (GeneCallAdds), genomic fasta files (GenomeAdds) and a vector of abbreviated names for the genomes (GenomeIDs) I am working with. I reached to the ftp adresses from https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get the gff and fasta files ----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GeneCallAdds_Gff <-  c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/105/GCF_000845105.1_ViralProj14522/GCF_000845105.1_ViralProj14522_genomic.gff.gz",
                       "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/907/375/GCF_000907375.1_ViralProj203280/GCF_000907375.1_ViralProj203280_genomic.gff.gz",
                       "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/890/915/GCF_000890915.1_ViralProj65223/GCF_000890915.1_ViralProj65223_genomic.gff.gz",
                       "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/844/285/GCF_000844285.1_ViralProj14523/GCF_000844285.1_ViralProj14523_genomic.gff.gz",
                       "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/817/995/GCF_002817995.1_ASM281799v1/GCF_002817995.1_ASM281799v1_genomic.gff.gz")

GenomeAdds_Fasta <- c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/105/GCF_000845105.1_ViralProj14522/GCF_000845105.1_ViralProj14522_genomic.fna.gz",
                      "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/907/375/GCF_000907375.1_ViralProj203280/GCF_000907375.1_ViralProj203280_genomic.fna.gz",
                      "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/890/915/GCF_000890915.1_ViralProj65223/GCF_000890915.1_ViralProj65223_genomic.fna.gz",
                      "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/844/285/GCF_000844285.1_ViralProj14523/GCF_000844285.1_ViralProj14523_genomic.fna.gz",
                      "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/817/995/GCF_002817995.1_ASM281799v1/GCF_002817995.1_ASM281799v1_genomic.fna.gz")

GenomeIDs <- c("FAdV-A", "FAdV-B", "FAdV-C", "FAdV-D", "FAdV-E")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Main Analysis Genome comparison (Database) ---- 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#I used library DECIPHER fistly.

#I made databases to find the homologous reagions of the genomes. 
#I used two ways to make databases either from 
# i) fasta file that I organised them or 
# ii) using ftp adresses of fasta files of the genomes.  
# The funcions Seq2DB made a database from sequences 
# as a database is necerrary for the function FindSynteny 
# which finds homologous regions between pairs of genomes, 
# which can then be aligned using AlignSynteny.
# Database from fttp adresses
GeneCalls <- FindHomology::GFFParser(GFFAddress = GeneCallAdds_Gff,Verbose = TRUE)

length(GeneCalls)

names(GeneCalls) <- GenomeIDs

# I did not store this database in my computer because I won't use it for this study. So, I assign it to a temperory file to display the database.   
# DBPath <- tempfile()
# 
# ?dbConnect
# DBConn <- dbConnect(SQLite(),DBPath)
# 
# for (i in seq_along(GenomeAdds_Fasta)) { 
#   Seqs2DB(seqs = GenomeAdds_Fasta[i], 
#           type = "FASTA", 
#           dbFile = DBConn,
#           identifier =  as.character(i), 
#           tblName = "Seqs",verbose = FALSE)  
#   }
# 
# 
# dbDisconnect(DBConn)
# 
# BrowseDB(DBPath)
# 
# SyntenyObject_from_ftp <- FindSynteny(dbFile = DBPath, verbose = TRUE)  
# 
# length(SyntenyObject_from_ftp)
# 
# unlink(DBPath)

# Database from fasta file
long_seqs <- readDNAStringSet(file.path("/Users/emineozsahin/Documents/R_worksheets_winter2020/Aviadenoviruses/FAdVs_modified_complete_genome.fasta"))

names(long_seqs) #"FAdV-A" "FAdV-B" "FAdV-C" "FAdV-D" "FAdV-E"

#Following code provide to write the database to the harddrive. 

Seqs2DB(seqs = long_seqs, type = "XStringSet", dbFile = "db_avia", identifier  = names(long_seqs), tblName = "Seqs", replaceTbl = TRUE) 

# To be able to see the database, displayed my data base from fasta file.There were difference for the descriptions because I organised the fasta file or I need to give a parsameter to the functionhwhen I made the database.
db_path_fasta = "/Users/emineozsahin/Documents/R_worksheets_winter2020/Aviadenoviruses/db_avia"

DBConn_file <- dbConnect(SQLite(), db_path_fasta)

BrowseDB(db_path_fasta)

dbDisconnect(DBConn_file)

SyntenyObject_from_fasta_file <- FindSynteny("db_avia") 

# As if I use database from fasta file I do not need to indicate genome names seperately to display on the plots I continue to use the SyntenyObject_from_fasta_file. 

#Visualisation of Synteny object. 

pairs(SyntenyObject_from_fasta_file,
      bounds = TRUE,
      boxBlocks = TRUE,
      gap = 0.5,
      line.main = 3,
      cex.labels = NULL,
      font.labels = 1,
      main = "Synteny of Genomes")

dev.copy2pdf(file = "Synteny_of_Genomes.png")

pairs(SyntenyObject_from_fasta_file[3:4, 3:4], boxBlocks = TRUE)

dev.copy2pdf(file = "Synteny_of_C_D.png")


#
plot(SyntenyObject_from_fasta_file,
     colorBy = 1,
     colorRamp = colorRampPalette(c("#FCF9EE", "#FFF272",
                                    "#FFAC28", "#EC5931",
                                    "#EC354D", "#0D0887")),
     barColor = "#CCCCCC",
     horizontal = TRUE,
     cex.labels = NULL,
     width = 0.7,
     scaleBar = TRUE)

dev.copy2pdf(file = "Plot_Genomes.png")

plot(SyntenyObject_from_fasta_file,"neighbor")
dev.copy2pdf(file = "Neighbor_Genomes.png")

plot(SyntenyObject_from_fasta_file[3:4, 3:4],"neighbor")
dev.copy2pdf(file = "Neighbor_C_D.png")

plot(SyntenyObject_from_fasta_file[3:4, 3:4])
dev.copy2pdf(file = "Synteny_C_D.png")

plot(SyntenyObject_from_fasta_file[3:4, 3:4], "frequency")
dev.copy2pdf(file = "Frequency_C_D.png")


# I combined Gene Calls variable that I obtained from fttp adresses above and SyntenyObject_from_fasta_file to obtain a matrix of what. 

MatrixObject <- NucleotideOverlap(SyntenyObject = SyntenyObject_from_fasta_file, GeneCalls = GeneCalls)

length(MatrixObject)

Homologs <- Catalog(MatrixObject, Verbose = TRUE)
length(Homologs)

hist(sapply(Homologs,function(x) nrow(x)),
            main = "Size of Agglomerations",
            ylab = "Number of Sets",
            xlab = "Numer of Gene Pairs",
            breaks = 30L, col = "black")

dev.copy2pdf(file = "Size_of_Agglomerations.png")


MaxRows <- max(sapply(Homologs,function(x) nrow(x))) #10

CoreSet <- which(sapply(Homologs, function(x) nrow(x)) == MaxRows)

# I assign a db_path to use later for CoreAligner function   



CoreGenome <- CoreAligner(Homologs[CoreSet],
                          PATH = db_path_fasta,
                          GeneCalls = GeneCalls,
                          Verbose = TRUE, 
                          IDs = GenomeIDs)

length(CoreGenome) #5

row.names(CoreGenome)

names(CoreGenome)

length(Homologs[CoreSet])

CoreDist <- DistanceMatrix(myXStringSet = CoreGenome,
                           verbose = FALSE,
                           correction = "Jukes-Cantor")

length(CoreDist) # 25
row.names(CoreDist) <- GenomeIDs



row.names(CoreDist) <- GenomeIDs

CoreDend <- IdClusters(myDistMatrix = CoreDist,
                       myXStringSet = CoreGenome,
                       method = "NJ",
                       verbose = FALSE,
                       showPlot = TRUE,
                       type = "dendrogram")

dev.copy2pdf(file = "tree_core_genomes.png")


PanGenomeMatrix <- LogicalPan(HomologList = Homologs,
                              GeneCalls = GeneCalls,
                              Verbose = TRUE,
                              Plot = FALSE)

row.names(PanGenomeMatrix) <- GenomeIDs

length(PanGenomeMatrix)

image(t(PanGenomeMatrix), col = c("white", "blue"), main = "Presence Absence")

dev.copy2pdf(file = "presence_absence.png")

PanGenome <- dist(PanGenomeMatrix,
                  method = "binary")


PanDend <- IdClusters(myDistMatrix = PanGenome,
                      method = "NJ",
                      type = "dendrogram",
                      showPlot = TRUE,
                      verbose = FALSE)

dev.copy2pdf(file = "tree_pan_genomes.png")


pan_core_dend <- dendlist(PanDend, CoreDend)

pan_core_dend <- untangle(pan_core_dend, method = "step2side")

tanglegram(pan_core_dend, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, margin_inner=6,lwd=2, main_left="Pan Genome", main_right="Core Genome")

dev.copy2pdf(file = "tree_pan_core_ comparison_genomes.png")

#Null hypothesis is that they have no relation or are not similar at all, the laternaive is that they in some way relate to one another. So this test shows that they are very similar if the p value is significantly low. Which agrees with the tanglegram
mantel.test(as.matrix(CoreDist), as.matrix(PanGenome))

CoreGenes <- matrix(data = NA_integer_,
                    ncol = length(GeneCalls),
                    nrow = length(CoreSet))

CoreAnnotations <- matrix(data = NA_character_,
                          ncol = length(GeneCalls),
                          nrow = length(CoreSet))


for (FAdV in seq_len(ncol(CoreGenes))) {
  for (block in seq_len(nrow(CoreGenes))) {
    CoreGenes[block, FAdV] <- unique(Homologs[CoreSet[block]][[1]][, FAdV])[!is.na(unique(Homologs[CoreSet[block]][[1]][, FAdV]))]
    CoreAnnotations[block, FAdV] <- GeneCalls[[FAdV]][CoreGenes[block, FAdV], "Annotation"]
  }
}

class(CoreAnnotations)

CoreAnnotations <- t(CoreAnnotations)

?t

data <- as.matrix(CoreAnnotations)

write.csv(CoreAnnotations, "Core_Annotations.cvs", )


CoreAnnotations <- data.frame(CoreAnnotations, stringsAsFactors = FALSE)
dim(CoreAnnotations)


row.names(CoreAnnotations) <- GenomeIDs

CoreAnnotations[, c(1L, 3L)]


CoreAnnotations[, 2L]

CoreAnnotations[, c(7L, 6L)]


#===================================================================
# Blast ----
#===================================================================

# I used BLAST alignment on the website of https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome  by using two accession number at a time and result alignment was downloaded by selecting the option Hit Table (text) .

# The accession numbers are NC_001720.1 (FAdV-A), NC_021221.1 (FAdV-B), NC_015323.1 (FAdV-C), NC_000899.1 (FAdVD) and NC_038332.1 (FAdV-E) which I obtained from at the beginning of the project. 

#To be able to write a tree in The Newick tree format I used Pan Genome dendogram with following codes.

my_tree <- as.phylo(PanDend) 

write.tree(phy=my_tree, file="exported_pan_tree.newick")

tree_blast <- newick2phylog("(FAdV-B:0.140701322,(FAdV-E:0.1402078048,((FAdV-D:0.0907618755,FAdV-C:0.0907618755):0.01020563044,FAdV-A:0.1009675059):0.03924029882):0.0004935172457);")

names <- c("FAdV-B", "FAdV-E",  "FAdV-D", "FAdV-C", "FAdV-A")
names(tree_blast$leaves) <- names
names(tree_blast$leaves)

tree_newick$leaves

# I downloaded the GFF files when I linked to the ftp server. 
FAdV_B <- try(read_dna_seg_from_file("FAdV-B.gb"))
FAdV_E <- try(read_dna_seg_from_file("FAdV-E.gb")) #gp20 hexon
FAdV_D <- try(read_dna_seg_from_file("FAdV-D.gb")) #gp15 hexon
FAdV_C <- try(read_dna_seg_from_file("FAdV-C.gb"))
FAdV_A <- try(read_dna_seg_from_file("FAdV-A.gb")) 

FAdV_B

dna_segs=list(FAdV_B, FAdV_E, FAdV_D, FAdV_C, FAdV_A)
names(dna_segs)
names(dna_segs) <- names
names(dna_segs)

# The order of the names are matter, I assign them as in order of blast tree. 

dna_segs[[1]]$end[2]
dna_segs[[2]]$name
dna_segs[[3]]$start
dna_segs[[3]]$length
dna_segs[[3]]$product

#dna_segs[[3]]$fill <- "green"

length(dna_segs[[3]]$fill) #29
#dna_segs[[3]]$fill[1:5] <- "red"
#dna_segs[[3]]$fill[24:29] <- "red"

#length(dna_segs[[4]]$fill) #46
#dna_segs[[4]]$fill[1:5] <- "red"
#dna_segs[[4]]$fill[40:46] <- "red"

# # I gave colors 
# comparisonsB_E$col <- apply_color_scheme(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),"grey")
# comparisonsE_D$col <- apply_color_scheme(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),"grey")
# comparisonsD_C$col <- apply_color_scheme(c(1, 2, 3, 4, 5),"grey")
# comparisonsC_A$col <- apply_color_scheme(c(1, 2, 3, 4),"grey")



dna_segs[[3]]$name <- gsub("(.+[A-Z])", "", dna_segs[[3]]$name)
dna_segs[[4]]$name <- gsub("(.+[A-Z]-[0-9]_)", "", dna_segs[[4]]$name)
dna_segs[[5]]$name <- gsub("(^[^L][A-Z][a-z][A-Z]+)", "", dna_segs[[5]]$name)
dna_segs[[1]]$name <- gsub("([A-Z][0-9]+_)", "", dna_segs[[1]]$name)
dna_segs[[2]]$name <- gsub("([A-Z][0-9][A-Z][0-9]+_)", "", dna_segs[[2]]$name)
length(dna_segs)
length(tree_newick$leaves)

# I downloaded the alignment file by selecting the option Hit Table (text).
comparisonsB_E <- try(read_comparison_from_blast("BLAST2-AlignmentB-E.txt"))
comparisonsE_D <- try(read_comparison_from_blast("BLAST2-AlignmentE_D.txt"))
comparisonsD_C <- try(read_comparison_from_blast("BLAST2-AlignmentD_C.txt"))
comparisonsC_A <- try(read_comparison_from_blast("BLAST2-AlignmentC_A.txt"))
comparisons = list(comparisonsB_E, comparisonsE_D, comparisonsD_C, comparisonsC_A)

comparisons[[1]]$e_value
comparisons[[1]]$gaps
comparisons[[1]]$aln_len
comparisons[[1]]$mism
comparisons[[1]]$direction
comparisons[[1]]$aln_len

seq_check <- df_Aviadeno_complate %>%
  group_by(Species_names) %>%
  summarise(len = nchar(as.character(Sequence))) %>%
  arrange(desc(len))

# Genome size are that I used has average of 45000bp
#xlims <- list(c(1,45781), c(1,45667), c(1,45063), c(1,43810), c(1,43804))
xlims <- list(c(1,46000), c(1,46000), c(1,46000), c(1,46000), c(1,46000))


main <- "Comparison of homologous segments in group of Fowl adenovirus genomes"

annots <- lapply(dna_segs, function(x){
    mid <- middle(x)
    annot <- genoPlotR::annotation(x1=mid, text=x$name, rot= 30)
    idx <- grep("g", annot$text, perl=TRUE)
    annot[idx[idx %% 5 == 0],] }) 

plot_gene_map(dna_segs=dna_segs, 
              comparisons=comparisons,
              dna_seg_scale=TRUE,
              annotations=annots,
              tree_scale = TRUE,
              scale = FALSE,
              main=main,
              gene_type="side_blocks",
              tree_width=1.5,
              annotation_height=1,
              xlims=xlims,
              limit_to_longest_dna_seg=FALSE,
              tree = tree_blast)

dev.copy2pdf(file = "NCBI_comparison_genomes.png")

#mauve alignment---- 
#I used mauve to align the sequences. When I connect the ftp server, I downloaded the gff files of the genomes. 

backbone <- read_mauve_backbone("Mauve_alignment_complate_genomes.backbone") 

backbone$dna_segs
backbone$comparisons

#ref:which of the dna segments will be the reference, i.e. which one will have its blocks in order , A numeric. filter_low If larger than 0, all blocks smaller that this number will be filtered out. Defaults to 0.

?read_mauve_backbone
length(backbone$dna_segs) # 5
backbone$dna_segs[1]
names(backbone$dna_segs)
labels_mauve <- c("FAdV-C",  "FAdV-A", "FAdV-E", "FAdV-D", "FAdV-B")
names(backbone$dna_segs) <- labels_mauve
names(backbone$dna_segs)

length(backbone$comparisons) #4
backbone$comparisons[1]


tree_mauve <- newick2phylog("((FAdV-C:0.280125,FAdV-A:0.13843):0.069215,((FAdV-E:0.180872,FAdV-D:0.229406):0.037601,FAdV-B:0.243013):0.0311133);")

names2 <- c("FAdV-C","FAdV-A", "FAdV-E", "FAdV-D", "FAdV-B")
  
names(tree_mauve$leaves) <- names2

names(tree_mauve$leaves)

tree_mauve$leaves

tree_mauve$nodes

for (i in 1:length(backbone$comparisons)){
    cmp <- backbone$comparisons[[i]]
    backbone$comparisons[[i]]$length <- abs(cmp$end1 - cmp$start1) + abs(cmp$end2 - cmp$start2)
}


plot_gene_map(dna_segs=backbone$dna_segs,
              comparisons=backbone$comparisons,
              global_color_scheme=c("length", "increasing", "red_blue", 0.7),
              override_color_schemes=TRUE,
              dna_seg_scale = TRUE,
              tree = tree_mauve,
              annotations=annots, tree_scale = TRUE, scale = FALSE, gene_type="side_blocks")



dev.copy2pdf(file = "mauve_FAdV_complete genome.png")

#--------------------------
#Extrassss  continue to explore how to use that blocks----
# library(genoPlotR) 
# 
# alignment <- AlignSynteny(SyntenyObject_from_fasta_file, "db")
# 
# class(alignment)
# 
# head(alignment)
# 
# names(alignment)
# 
# writeXStringSet(DNAStringSet(alignment), "alignment_out.fa")
# 
# blocksA_B <- alignment[[1]]
# blocksB_C <- alignment[[5]]
# blocksC_D <- alignment[[8]] # There are 15 blocks
# blocksD_E <- alignment[[10]]
# blocksD_E <- alignment[[10]]
# 
# 
# blocksA_B <- unlist(alignment[[1]])
# blocksA_C <- unlist(alignment[[2]])
# blocksA_D <- unlist(alignment[[3]]) # There are 15 blocks
# blocksA_E <- unlist(alignment[[4]])
# 
# fadvA_vs_fadvB <- try(read_dna_seg_from_fasta("genome_blocksA_B_out.fa"))
# writeXStringSet(blocksA_B, "genome_blocksA_B_out.fa")
# 
# fadv1_vs_fadv5 <- try(read_comparison_from_blast("genome_blocks5_4_out.fa"))
# 
# AlignProfiles(DNAStringSet(alignment))
# 
# 
#' # plot_gene_map(dna_segs=list(fadv1, fadv5), comparisons=list(blocks5_C), xlims=xlims, main= "not yet", gene_type="side_blocks", dna_seg_scale=TRUE, scale=FALSE)
#' 
#' 
#' #--------------------------
#' #BiocManager::install("GenomicAlignments")  # look at this package 
#' #library(GenomicAlignments)
#' 
#' #' #''''''
#' #' library(rBLAST)
#' #' install.packages("ncbi-blast+")
#' #' #db1 <- DNAStringSet("db")
#' #' seq <- readRNAStringSet("/Users/emineozsahin/Documents/R worksheets/BINF*6210/Assignment-5/Org_Aviadenovirus_complete_genome_fetch.fasta")
#' #' 
#' #' bl <- blast(db="./db1")
#' #' ?blast()
