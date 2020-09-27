##=========================================================*
#**********************************************************#
# Emine Ozsahin
# September 2020
# whispovirus (Nimaviridae)
# 15 complete genome aligned (MAFFT) on terminal

##=========================================================*
#**********************************************************#

# Libraries

getwd()

setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/Coronavirus_BVES")

##=========================================================*
#**********************************************************#
# Main package used in this study
library(phangorn)

# Not all these libraries were used 
library(seqinr)
library(cluster)
library(dendextend)
library(factoextra)
library("gplots")
library("pheatmap")
library("d3heatmap")
library(colorspace)

library(snpStats)
library(factoextra)
library(ggbiplot)
library(clustertend)
library(clues)
library(dendextend)
library(colorspace)
library(cluster)
library(tidyverse)
library(phylotools)
library(gplots)
library(factoextra)

#install.packages("circlize")
library(circlize) # you have to try this package 
library(ape)
##=========================================================*
#**********************************************************#

# Alignment

##=========================================================*
#**********************************************************#
# Data load in
aln <- read.alignment("wispoviruses.aln", format = "fasta")
mat <- as.matrix.alignment(aln)
mat[1:15,30000:30050]
mat[1, ]
df_mat <- as.data.frame(mat)
df_mat[1:15,1:10]
rownames(df_mat)

seq1 <- as.character(df_mat[1,])
seq1[1:100]
length(seq1) #MF768985.1 
class(seq1)
names(seq1) <- "MF768985.1" 
length(which(seq1 != "-"))
length(which(seq1 == "-"))

seq2<- as.character(df_mat[2,]); length(seq2) #AF332093.3
names(seq2) <- "AF332093.3" 

seq3<- as.character(df_mat[3,]) #MN840357.1
names(seq3) <- "MN840357.1" 
seq4<- as.character(df_mat[4,]) #KU216744.2
names(seq4) <- "KU216744.2" 

seq5<- as.character(df_mat[5,]) #MH090824.1
names(seq5) <- "MH090824.1" 
seq6<- as.character(df_mat[6,])#AF369029.2
names(seq6) <- "AF369029.2" 
seq7<- as.character(df_mat[7,])#KT995472.1
names(seq7) <- "KT995472.1" 
seq8<- as.character(df_mat[8,])#KT995471.1
names(seq8) <- "KT995471.1" 
seq9<- as.character(df_mat[9,])#KT995470.1
names(seq9) <- "KT995470.1" 
seq10<- as.character(df_mat[10,]) #MH663976.1
names(seq10) <- "MH663976.1" 
seq11 <- as.character(df_mat[11,])#KX686117.1
names(seq11) <- "KX686117.1" 
seq12<- as.character(df_mat[12,])#MG702567.1
names(seq12) <- "MG702567.1" 
seq13 <- as.character(df_mat[13,])#KY827813.1
names(seq13) <- "KY827813.1" 
seq14<- as.character(df_mat[14,])#AF440570.1
names(seq14) <- "AF440570.1" 
seq15<- as.character(df_mat[15,])#JX515788.1
names(seq15) <- "JX515788.1" 

seqlist <- list(noquote(paste0("seq", c(1:15))))

sq <- list(seq1,  seq2,  seq3 , seq4 , seq5,  seq6 , seq7,  seq8,  seq9,  seq10, seq11, seq12, seq13, seq14, seq15)
nms <- c("MF768985.1", "AF332093.3", "MN840357.1", "KU216744.2", "MH090824.1", "AF369029.2", "KT995472.1", "KT995471.1", "KT995470.1", "MH663976.1", "KX686117.1", "MG702567.1", "KY827813.1", "AF440570.1", "JX515788.1")
names(sq) <- nms
genome_seq <- lapply(sq, function(x) length(which(x != "-")))
lapply(sq, function(x) length(which(x != "-")))

for (z in 1:15) {  name <- names(sq[z])
print(name)
}

names(sq[1])
final <- list()
length(sq[[1]])
reference[1]
       
c <- data.frame() 
rfn = names(sq[15])
reference = sq[[15]]
for (z in 1:15) { 
  rn<-names(sq[z])
  record <- sq[[z]]
  pc=0
  for (i in 1:dim(df_mat)[2]) { 
    if (reference[i] == record[i] & (record[i] != "-" & reference[i] != "-")) {pc=pc+1}}
  temp=data.frame(ref_name=rfn, pcount=pc, record_name=rn)
  c <- rbind(c, temp)
  }
  
write.csv(c, "wispo_match_nucl_sep17.csv")  
  
# distances
d <- dist.alignment(aln, matrix = "identity")

# convert the distances to a dataframe 
z <- apply(as.data.frame(as.matrix(d)), 2, as.numeric)
dim(z)
rownames(z) <- colnames(z)
head(z)

#distance scale : scale which represents the number of differences between sequences 
#(e.g. 0.1 means 10 % differences between two sequences)
# Therefore I can calculate the similarity by substraction distance from 1
apply(z, 2, function(x) 1-x)


#######################################################################################
#Rooted phylogenetic tree refers to a phylogenetic tree that shows the ancestry relationship, 
# while unrooted phylogenetic tree refers to the phylogenetic tree that only shows the relatedness of organisms.

## set.seeds are important as every time different figure (slightly different) is produced
set.seeds()
fdir <- system.file("/Users/emineozsahin/Documents/R_worksheets/Coronavirus_BVES", package = "phangorn")
wispoviruses <- read.phyDat(file.path("/Users/emineozsahin/Documents/R_worksheets/Coronavirus_BVES/wispoviruses.aln"),format="fasta") 

nms <- c("WSSV-AU", "WSSV-CN", "WSSV-CN-95-DFPE", "WSSV-MEX2008", "WSSV-EC-15098",
         "WSSV-TH", "WSSV-CN01", "WSSV-CN03", "WSSV-CN02","WSSV-CNPc02", "WSSV-CN-Pc",
         "WSSV-IN-AP4RU", "WSSV-CN04", "WSSV-TW", "WSSV-K-LV1")

names(wispoviruses) <- nms

      #distances
dm  <- dist.ml(wispoviruses)

treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

plot(treeUPGMA, main="UPGMA")
plot(treeNJ, main="NJ")
plot(treeNJ, "unrooted", main="NJ")
plotBS(treeNJ, bs, "unrooted", main="NJ")
add.scale.bar(cex = 1, font = 3, col = "black") 

parsimony(treeUPGMA, wispoviruses) #3795
parsimony(treeNJ, wispoviruses) #3779

treePars  <- optim.parsimony(treeUPGMA, wispoviruses)
plot(treePars)
add.scale.bar()

treePars  <- optim.parsimony(treeNJ, wispoviruses)
plot(treePars)
add.scale.bar()

treeRatchet  <- pratchet(wispoviruses, trace = 0)
plot(treeRatchet)
add.scale.bar()

#Final p-score 3711 after  7 nni operations 
parsimony(c(treePars, treeRatchet), wispoviruses) #3711 3738

mtJC = modelTest(wispoviruses)
?modelTest

envJC <- attr(mtJC, "env")

mtJC$Model[which.min(mtJC$AICc)] #"GTR+G+I"

fitJC <- eval(get("GTR+G+I", env), env)

?pml
fitJC = pml(treeNJ, data=wispoviruses)
methods(class="pml")
fitJC  <- optim.pml(fitJC, TRUE)
logLik(fitJC)

      # Bootstrap 
bsJC = bootstrap.pml(fitJC, bs=100, optNni=TRUE)

png("wispovirus_nj.png", units="in", width=10, height=10, res=300)
plotBS(midpoint(fitt$tree), bsJC, p = 50, type="p", bs.adj = c(1, -0.3), cex = 1)
add.scale.bar(cex = 1, col = "black", x = 0.0000001, y = 0.7) 
dev.off()

plotBS(midpoint(fitt$tree), bsJC, p = 50, type="unrooted", bs.adj = c(1, -0.3), cex = 1)

?consensusNet
cnet <- consensusNet(bsJC, p=0.6)
plot(cnet, "2D", show.edge.label=TRUE)
add.scale.bar(x = -598.3842, y = 420.0865, cex = 0.5, font = 2, col = "black", length = 0.005) 


locator()
?add.scale.bar
cnet$Nnode
##################### 


####################
set.seed(1233)
names(wispoviruses) <- nms
names(wispoviruses) <- gsub("WSSV-", "", fit$tree$tip.label)
dmF81 = dist.ml(wispoviruses, "F81")
dmF81


df81 <- apply(as.data.frame(as.matrix(dmF81)), 2, as.numeric)
dim(df81)
rownames(df81) <- nms
colnames(df81) <- nms
head(df81)

#distance scale : scale which represents the number of differences between sequences 
#(e.g. 0.1 means 10 % differences between two sequences)
# Therefore I can calculate the similarity by substraction distance from 1
smlrt <- apply(df81, 2, function(x) (1-x)*100)
rownames(smlrt) <- nms
#write.csv(smlrt, "wispo_similarity.csv")


treeNJ81  <- NJ(dmF81)
?dist.ml

plot(treeNJ81, "unrooted")

add.scale.bar(cex = 1, col = "black")

# 2. alternative: modelTest
mt81 <- modelTest(wispoviruses, tree=treeNJ81, multicore=TRUE)
mt81[order(mt81$AICc),]

# choose best model from the table according to AICc
bestmodel <- mt$Model[which.min(mt$AICc)]
bestmodel#"GTR+G+I"
env81 = attr(mt81, "env")
#fitStart = eval(get("GTR+G+I", env), env)
# or let R search the table
fitStart = eval(get(bestmodel, env81), env81)
?eval
# equivalent to:   fitStart = eval(get("GTR+G+I", env), env)
fit = optim.pml(fitStart, rearrangement = "stochastic",
                optGamma=TRUE, optInv=TRUE, model="GTR")

?optim.pml

bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)

png("wispovirus_Maximum.png", units="in", width=6, height=6, res=600)
plotBS(midpoint(fit$tree), bs, p = 50, type="p", bs.adj = c(1, -0.3), cex = 1)
add.scale.bar(cex = 1, col = "black", x = 0.0000001, y = 0.7) 
dev.off()
# I used locator() to find x and y coordinates.
locator()

?add.scale.bar
png("temp.png", units="in", width=6, height=7, res=600)
plotBS(fit$tree, bs, p = 50, type="unrooted",cex = 0.5, bs.adj = c(1.1, 0), edge.width=1.5, show.tip.label=TRUE, node.depth = 1, bs.col ="black") #
add.scale.bar(cex = 0.5, col = "black")
dev.off()
?optim.pml
?plot.phylo
?plotBS

png("temp.png", units="in", width=5, height=5, res=300)
plotBS(fit$tree, bs, p = 50, type="unrooted",cex = 1, bs.adj = c(1, 0))
add.scale.bar(cex = 1, col = "black")
dev.off()

plotBS(midpoint(fit$tree), bs, p = 50, type="unrooted",cex = 0.5)

plot(fit$tree, type="unrooted",cex = 1)

midpoint(fit$tree)
?midpoint
length(fit$tree)
fit$tree[1]

fit

cnet <- consensusNet(bs, p=0.5)
?consensusNet
plot(cnet, "2D", show.edge.label=TRUE)
add.scale.bar(cex = 0.5, font = 2, col = "black") 
original_edge_lenght=(fit$tree$edge.length)
fit$tree$edge.length <- original_edge_lenght
plotBS(midpoint(fit$tree), bs, p = 50, type="unrooted", edge.width = 2)

?plotBS
fit$tree$edge.length

cnet <- consensusNet(bs, p=0.5)
?consensusNet
plot(cnet, "2D", show.edge.label=TRUE)

plot(cnet, edge.label = cnet$edge.labels, show.edge.label = T, "2D", col.edge.label = "blue", cex=.7, bs.adj = c(2, 2))
?plot

plotBS(midpoint(fit$tree), bs, p = 50, type="unrooted",cex = 1, )
add.scale.bar(cex = 1, col = "black") 

####################################################################
ngj <- nj(d)

install.packages("phangorn")

plot(as.phylo(ngj), type = "unrooted", cex = 0.6, no.margin = FALSE)
add.scale.bar(cex = 1, col = "black") 
               
hc <- hclust(d)
dhc <- as.dendrogram(hc)
library(ggplot2)
library(ggdendro)
ddata <- dendro_data(dhc, type = "rectangle")
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),colour = "dark grey", alpha = 250) + 
  scale_y_reverse(expand = c(0.2, 0))+ 
  coord_flip() + 
  theme_dendro() +
  geom_text(data = label(ddata), 
            aes(x = x, y = y, label = label), vjust = 0.5, hjust = -0.1, size = 3)

p

ggplot() +
  geom_segment(data = segment(ddata), 
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(ddata), 
            aes(x = x, y = y, label = label, hjust = 0), size = 3) +
  coord_flip() +
  theme_dendro() +
  scale_y_reverse(expand = c(0.2, 0))





ggdendrogram(hc, rotate = TRUE , size = 3)
?ggdendrogram





