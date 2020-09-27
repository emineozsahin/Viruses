##=========================================================*
#**********************************************************#
# Emine Ozsahin
# September 2020
# whispovirus (Nimaviridae)
# 18 complete genome aligned (MAFFT) on terminal

##=========================================================*
#**********************************************************#

# Libraries

getwd()

setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/Coronavirus_BVES")

##=========================================================*
#**********************************************************#

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

names(genome_seq)


for (i in genome_seq) {print(i)}
for (i in names(genome_seq)) {print(i)}

lapply(sq, function(x) length(which(x != "-")))

for (z in 1:15) {  name <- names(sq[z])
print(name)
}

names(sq[1])

pcount = 0 

final <- list()

length(sq[[1]])

 #

reference[1]
#c <- data.frame() #ref_name="AC899999", pcount=1, record_name="AC899999999"

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
  
 



#pcount = 0 
for (i in 1:dim(df_mat)[2]) {
  if (seq1[i] == seq2[i] & (seq1[i] != "-" & seq2[i] != "-")) {pcount=pcount+1} }




for (i in 1:dim(df_mat)[2]) {
  if (df_mat[1,][i] == df_mat[3,][i]) {pcount=pcount+1}}  

for (i in 1:dim(df_mat)[2]) {
  if (df_mat[1,][i] == df_mat[3,][i]) {pcount=pcount+1}} 







c

?range
dim(mat)

length(aln)
aln$com
sequinr
# Tidy the isolate names
aln$nam
length(aln$nam) #15
length(unique(aln$nam)) #15

# distances
d <- dist.alignment(aln, matrix = "identity")

?dist.alignment

class(d)

#Perfect heatmap + cluster for individuals

#Tuesday 9am 

?write.csv
#jallen1@ugcloud.ca

# convert the distances to a dataframe 
z <- apply(as.data.frame(as.matrix(d)), 2, as.numeric)
dim(z)
rownames(z) <- colnames(z)
head(z)

#distance scale : scale which represents the number of differences between sequences 
#(e.g. 0.1 means 10 % differences between two sequences)
# Therefore I can calculate the similarity by substraction distance from 1
apply(z, 2, function(x) 1-x)


# Some options which I did not use
#my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 147)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(-1,0,length=100),  # for red
#seq(0.01,0.8,length=100),           # for yellow
#seq(0.81,1,length=100))             # for green

# siply plot I don't like the shape 
#heatmap.2(z)

# This is better but there is an issue with names they do not match
heatmap.2(z,
          #cellnote = mat_data,  # same data set for cell labels
          #main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9))     # widens margins around plot
#col=my_palette,       # use on color palette defined earlier
#breaks=col_breaks,    # enable color transition at specified limits
#dendrogram="row")     # only draw a row dendrogram
#Colv="NA")            # turn off column clustering


#Pretty heat maps
png("pretty.heat.whispo.png", units="in", width=10, height=10, res=300)
pheatmap(z, fontsize_row = 3, fontsize_col = 3, cellwidth = 2, cellheight = 2) #cutree_cols  = 10, 
dev.off()

#Interactive heat maps, !!!! does not have a legend for color 
png("interactive.png", units="in", width=20, height=20, res=750)
d3heatmap(scale(z), colors = "RdYlBu",
          k_row = 4, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)
dev.off()

#Enhanced heat maps, it is good no problem
png("enhanced_heat.whispo.png", units="in", width=10, height=10, res=300)
heatmap.2(z, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
dev.off()

# dissimilarity plot
png("dissimilarity.whispo.png", units="in", width=10, height=10, res=300)
fviz_dist(d, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"), lab_size = 4)
dev.off()

#######################################################################################
#Rooted phylogenetic tree refers to a phylogenetic tree that shows the ancestry relationship, while unrooted phylogenetic tree refers to the phylogenetic tree that only shows the relatedness of organisms.
# Simply Dendogram
# PLOT() FUNCTION FOR PHYLO: https://www.rdocumentation.org/packages/ape/versions/5.2/topics/plot.phylo
#

## set.seeds are important as every time different figure (slightly different) is produced
set.seeds()
library(phangorn)
fdir <- system.file("/Users/emineozsahin/Documents/R_worksheets/Coronavirus_BVES", package = "phangorn")
wispoviruses <- read.phyDat(file.path("/Users/emineozsahin/Documents/R_worksheets/Coronavirus_BVES/wispoviruses.aln"),format="fasta") 

wispoviruses[14]


nms <- c("WSSV-AU", "WSSV-CN", "WSSV-CN-95-DFPE", "WSSV-MEX2008", "WSSV-EC-15098",
         "WSSV-TH", "WSSV-CN01", "WSSV-CN03", "WSSV-CN02","WSSV-CNPc02CN", "WSSV-CN-Pc",
         "WSSV-IN-AP4RU", "WSSV-CN04", "WSSV-TW", "WSSV-K-LV1")

for (i in nms) {print(i)}
names(wispoviruses) <- nms
dm  <- dist.ml(wispoviruses)
?dist.ml

dm

df <- apply(as.data.frame(as.matrix(dm)), 2, as.numeric)
dim(df)
rownames(df) <- nms
colnames(df) <- nms
head(df)

#distance scale : scale which represents the number of differences between sequences 
#(e.g. 0.1 means 10 % differences between two sequences)
# Therefore I can calculate the similarity by substraction distance from 1
apply(df, 2, function(x) 1-x)
write.csv(df, "wispo_similarity.csv")

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

# 2. alternative: preper with modelTest
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

# fit$tree$edge.length/10
# 
# fit$tree$edge.length  <- fit$tree$edge.length/10
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
?nj
ngj <- nj(d)

install.packages("phangorn")

plot(as.phylo(ngj), type = "unrooted", cex = 0.6, no.margin = FALSE)
add.scale.bar(cex = 1, col = "black") 

?plot

?hclust
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

# Fix the label size !!!!!!!!!!!!!!!!!!!!
plot(dhc)

# Horizontal plot
nodePar <- list(lab.cex = 1, pch = c(NA, 19), 
                cex = 1, col = "grey")

# Fix the label size !!!!!!!!!!!!!!!
plot(dhc, xlab = "", nodePar = nodePar, horiz = TRUE)

#change the labes to numbers then color
regions1 <- labels(dhc)%>%
  str_replace_all(c("BC.+" = "BC", "MB.+" = "MB", "NB.+" = "NB", "NL.+" = "NL", "NS.+" = "NS",
                    "ON.+" = "ON", "PE.+" = "PE", "SK.+" = "SK"))

regions_num <- labels(dhc)%>%
  str_replace_all(c("BC.+" = "1", "MB.+" = "2", "NB.+" = "3", "NL.+" = "4", "NS.+" = "5",
                    "ON.+" = "6", "PE.+" = "7", "SK.+" = "8"))
###
#labels(dhc) <- regions[order.dendrogram(dhc)]

length(labels(dhc))

labels_colors(dhc) <- rainbow_hcl(length(unique(regions1)))[as.numeric(regions_num)]
head(labels(dhc))
head(labels_colors(dhc))


png("circlize.png", units="in", width=11, height=11, res=300)
circlize_dendrogram(dhc)
legend("topleft", legend = unique(regions1), fill = c("#E495A5", "#64B5D6", "#D2A277", "#ABB065", "#D995CF", "#39BEB1", "#ACA4E2", "#72BB83"))
dev.off()

# Horizontal plot----
png("horizontal.png", units="in", width=5, height=5, res=300)
nodePar <- list(lab.cex = 0.5, pch = c(NA, 19), col = "blue", labels_cex= 0.5)
dhc %>% set("labels_cex", 0.3) %>% plot(xlab = "Height", nodePar = nodePar, horiz = TRUE)
legend("topleft", legend = unique(regions1), fill = c("#E495A5", "#64B5D6", "#D2A277", "#ABB065", "#D995CF", "#39BEB1", "#ACA4E2", "#72BB83"),
       xjust = -1, yjust = 0, x.intersp = 1, y.intersp = 1,)

dev.off()
?legend
#"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"

##=========================================================*
#**********************************************************#

# VCF file 

##=========================================================*
#**********************************************************#
vcf <-read.table("gisaid_covid19.vcf.vcf", stringsAsFactors = FALSE)

dim(vcf) #185 156

vcf_id <-readLines("gisaid_covid19.vcf.vcf")

vcf_id <- vcf_id[-(grep("#CHROM", vcf_id)+ 1):-(length(vcf_id))]

ids <- unlist(strsplit(vcf_id[length(vcf_id)],"\t"))
a <- gsub("/2020.+", "", ids)
b <- gsub("hCoV-19/Canada/", "", a)
b
names(vcf) <- b

length(ids) #156

length(names(vcf))

png("SNP_pos.png", units="in", width=10, height=10, res=300)
plot(vcf$POS, cumsum(vcf$POS), xlab = "SNP positions", ylab = "Cumulative summary of SNP positions")
dev.off()

length(vcf$POS)

png("SNP_freq.png", units="in", width=10, height=10, res=300)
plot(vcf$POS, frequency$freq, xlab = "SNP positions", ylab = "SNP frequency")
dev.off()


library(BioCircos)

# Chromosomes on which the points should be displayed
points_chromosomes = c(1:3) 
# Coordinates on which the points should be displayed
points_coordinates = c(vcf$POS) 
# Values associated with each point, used as radial coordinate 
#   on a scale going to minRadius for the lowest value to maxRadius for the highest value
points_values = 0:length(vcf$POS)

?BioCircosSNPTrack

tracklist = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates, 
                              points_values, colors = c("tomato2", "darkblue"))  #minRadius = 0.5, maxRadius = 0.9

# Background are always placed below other tracks
tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", 
                                                 minRadius = 0.5, maxRadius = 0.9,
                                                 borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")  

BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 18, genomeLabelDy = 0)




library("circ")
?circos.genomicPoints()
df <- data.frame(a=vcf$`#CHROM`, b=vcf$POS, c=vcf$POS)
circos.initialize()
circos.par("track.height" = 0.1)
circos.initialize(factors = df$a, x = (df$b) +10 )

circos.track(factors = df$a, y = df$b,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

circos.genomicPoints(df, data.frame(ref=vcf$INFO, 
                                    alt=rep(0, length(vcf$POS)), 
                                    c=rep("+", length(vcf$POS),
                                          d=vcf$POS, e=vcf$POS)))

vcf <- rbind(names(vcf), vcf)

vcf <- vcf[,-c(1:9)]

vcf_ord <- vcf[, order(vcf[1,])] 

dim(vcf_ord) #186 147

vcf_ord[1,]
vcf_ord[,1]
vcf_ord[2,]

vcf_ord[1,]
class(vcf_ord)
class(vcf_ord[1,1])
vcf_ord <- apply(vcf_ord, 2, as.numeric)

vcf_ord <-vcf_ord[-1,]

class(vcf_ord[1,1])
dim(vcf_ord) #185 147

rownames(vcf_ord) # NULL
rownames(vcf_ord) <- 1:185 # snps
colnames(vcf_ord) # individuals

class(vcf_ord)

Variant <- t(vcf_ord)
class(Variant[1,1])
rownames(Variant) # individuals
colnames(Variant) # variants
dim(Variant)#147 185

tail(Variant)

length(Variant[is.na(Variant)]) #0
Variant <- Variant[,apply(Variant, 2, function(x) length(unique(x))) != 1]
dim(Variant) #147 185

Variant[1,]

length(which(colSums(is.na(as.matrix(Variant))) > 0)) 

f <- data.frame()
for (i in 1:dim(Variant)[2]) { num <- length(grep(0, Variant[,i])) 
temp <- data.frame(num)
f <- rbind(f, temp)
}
head(f)


frequency <- as.data.frame(apply(f, 1, function(x) x=x/dim(Variant)[1]))
colnames(frequency) <- "freq"
head(frequency)

hist(frequency$freq, main = "Allele frequency spectrum", 
     xlab="Frequency in Population", 
     ylab = "Proportion of Variants", 
     col = "grey")

plot(frequency$freq)
boxplot(frequency$freq)

which(lapply(apply(Variant, 2, table), length) == 1) # 0 

snps <- new("SnpMatrix", Variant)
snps
dim(snps)
colnames(snps) 

# Check the final matrix statistics
summary(snps)

# Retrieve the column statistics
snpsum <- col.summary(snps)
summary(snpsum)
length(snpsum$z.HWE^2) ##185

# Check the multiple allelle frequescies (MAF) and HWE 
par(mfrow = c(1, 2))
hist(snpsum$MAF)
hist(snpsum$z.HWE)
par(mfrow = c(1, 1))

#Convert z values of HWE to p values 
pvalue = pnorm(-abs(snpsum$z.HWE))
pvalue

# Filter the snps based on HWE and MFA 
use <- snpsum$MAF > 0.01 & pvalue < 0.05
sum(use) 

rownames(Variant)

dim(Variant)
pca <- prcomp(Variant, center = TRUE, scale = FALSE)

summary(pca)

qqnorm(pca$x[,1]); qqline(pca$x[,1], col = "steelblue", lwd = 2)

hist(pca$x[,1], main = "Histogram of quantiles", xlab= "Quantiles")

plot(summary(pca)$importance[3,]*100, xlab = "Number of components", ylab = "Explained %")

barplot((pca[[1]]^2)[1:8],
        main= "Scree Plot", 
        names.arg=1:8, 
        xlab="PCs", 
        ylab="variance")

regions <- rownames(pca$x)%>%
  str_replace_all(c("BC.+" = "BC", "MB.+" = "MB", "NB.+" = "NB", "NL.+" = "NL", "NS.+" = "NS",
                    "ON.+" = "ON", "PE.+" = "PE", "SK.+" = "SK"))

unique(regions)
length(regions)
length(rownames(pca$x))

png("ggbiplot.png", units="in", width=5, height=5, res=300)
ggbiplot(pca, var.axes=FALSE, groups = regions) +
  scale_colour_manual(name="Regions",values= c("forest green", "red4", "dark blue", "orange", "purple", "red", "light blue", "yellow")) +
  theme_minimal() +
  theme(legend.direction = 'horizontal', legend.position = 'top') 
dev.off()

png("SNPs_distance.heat.png", units="in", width=7.5, height=7.5, res=300)
fviz_dist(dist(Variant), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()


hopkins(Variant, 5, byrow = F, header = F)

length(rownames(Variant))
length(regions)
unique(regions)
length(which(regions == "BC")) #59
length(which(regions == "MB")) #7
length(which(regions == "NB")) #4
length(which(regions == "NL"))# 1
length(which(regions == "NS")) #4
length(which(regions == "ON"))#66
length(which(regions == "PE"))#1
length(which(regions == "SK"))#5
num <- c(rep(1, 59), rep(2, 7), rep(3, 4), rep(4, 1), rep(5, 4), rep(6, 66), rep(7, 1), rep(8, 5))
length(num)

Variant_num <- Variant

length(rownames(Variant_num))

rownames(Variant_num) <- num

# Hierarchical Clustering
hc <- hclust(dist(Variant_num), method = "complete")
dend_hc <- as.dendrogram(hc)

# Agglomerative Nesting (Hierarchical Clustering)
ag <- agnes(dist(Variant_num), method = "complete")
dend_ag <- as.dendrogram(ag)

### Divisive Clustering
dv <- diana(dist(Variant_num))
dend_dv <- as.dendrogram(dv)

# obtain the dendogram leaf names as numbered group names to check the strenght of cluster measures 

hclust <- num[order.dendrogram(dend_hc)]
hagnes <- num[order.dendrogram(dend_ag)]
hdiana <- num[order.dendrogram(dend_dv)]

hcm <- adjustedRand(hclust, num)
ham <- adjustedRand(hagnes, num)
hdm <- adjustedRand(hdiana, num)

compare <- data.frame(hclust = hcm, hagnes=ham, hdiana=hdm)
compare

plot(dend_hc)

nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")


par(mfrow = c(2, 2))
dev.set(which = dev.next())
dev.cur()
dev.off(i) 

#variants
png("SNPs_dend.heat.png", units="in", width=10, height=10, res=300)
plot(dend_hc, xlab = "Height", nodePar = nodePar, horiz = TRUE)

nodePar <- list(lab.cex = 0.5, pch = c(NA, 19), col = "blue", labels_cex= 0.5)
dend_hc %>% set("labels_cex", 0.5) %>% plot(xlab = "Height", nodePar = nodePar, horiz = TRUE)
legend("topleft", legend = unique(regions), fill = c("#E495A5", "#64B5D6", "#D2A277", "#ABB065", "#D995CF", "#39BEB1", "#ACA4E2", "#72BB83"),
       xjust = -1, yjust = 0, x.intersp = 1, y.intersp = 1,)

dev.off()


#alignment
plot(dhc, xlab = "Height", nodePar = nodePar, horiz = TRUE)
par(mfrow = c(1, 1))


# Create dendrogram plot
library(ggdendro)
dendro.plot <- ggdendrogram(data = dend_hc, rotate = TRUE) + 
  theme(axis.text.y = element_text(size = 6))
dendro.plot
?fviz_dist

dist.plot <- fviz_dist(dist(scale(Variant)), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dist.plot

# All together
grid.newpage()
print(dist.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.92))

##*********************************************
# The standard R code for computing hierarchical clustering looks like this:
# Compute dissimilarity matrix
res.dist <- dist(as.data.frame(Variant_num), method = "euclidean") #The standard R code for computing hierarchical clustering looks like this:

# Compute hierarchical clustering
res.hc <- hclust(res.dist, method = "ward.D2")#The standard R code for computing hierarchical clustering looks like this:

# Visualize
plot(res.hc, cex = 0.5)#The standard R code for computing hierarchical clustering looks like this:
##*********************************************


# Correlation-based distance method
res.dist <- get_dist(as.data.frame(Variant), method = "pearson")
head(round(as.matrix(res.dist), 2))[, 1:6]

# Visualize the dissimilarity matrix
png("SNPs_dist.heat2.png", units="in", width=10, height=10, res=300)
fviz_dist(res.dist, lab_size = 5) #In the plot, similar objects are close to one another. Red color corresponds to small distance and blue color indicates big distance between observation.
dev.off()

# Enhanced k-means clustering
?eclust # pay attantion to the hc_metric and FUNcluster (here kmeans)
res.km <- eclust(as.data.frame(Variant), "kmeans", nstart = 25, hc_metric = "euclidean") # In the following R code, we’ll show some examples for enhanced k-means clustering and hierarchical clustering. Note that the same analysis can be done for PAM, CLARA, FANNY, AGNES and DIANA.
res.km
fviz_cluster(res.km)

# Gap statistic plot
fviz_gap_stat(res.km$gap_stat)

# Silhouette plot
fviz_silhouette(res.km)

# Optimal number of clusters using gap statistics
res.km$nbclust

# Print result
res.km

# Enhanced hierarchical clustering
res.hc <- eclust(as.data.frame(Variant), "hclust") # compute hclust

res.hc$labels

png("SNPs_eclust_dend.png", units="in", width=10, height=10, res=300)
fviz_dend(res.hc, rect = TRUE, cex = 0.5, main = "") # dendrogam
dev.off()

fviz_silhouette(res.hc) # silhouette plot

fviz_cluster(res.hc) # scatter plot

eclust(as.data.frame(Variant), "kmeans", k = 8)

##*********************************************
#PUBMED SEARCH
##*********************************************
l <- (read.table("pubmed_result.txt", sep = "\t"))

l$V1

strsplit(l$V1, ".")

?strsplit

l1 <- as.data.frame(unlist(strsplit(as.character(l$V1), "\t")))

l1 <- as.data.frame(strsplit(as.character(l$V1), "/n"))

l1

summary(l$V1)


# Heatmap
heatmap.plot <- ggplot(data = otter.long, aes(x = variable, y = accession)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")


c = CircosPlot(G,node_grouping="class",
               node_color="class",
               node_order="class",
               node_labels="class",
               group_label_position="middle",
               nodeprops ={"radius":5},
               node_label_layout="rotation",
               group_label_color=True, figsize = (25))
c.draw()
plt.tight_layout(rect=(0.5, 0.5, 0.5, 0.5))
plt.show()


