##=========================================================*
#**********************************************************#
# Emine Ozsahin
# April 27 2020
# SARS-CoV2 Canada sequences were downloaded from GISIAD (https://www.gisaid.org/)
# 147 SARS-CoV2 genomes were aligned using MAFFT in compute canada graham
# alignment converted to the VCF file (SNP-sites) in my computer 
# alignment and vcf were used for cluster analysis seperately.
# vcf is also used for PCA 
# try factor analysis with vcf
# try regression if you have population data

##=========================================================*
#**********************************************************#

# Libraries

setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/Coronavirus_BVES")

##=========================================================*
#**********************************************************#

# Not all of these libraries necessary
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
aln <- read.alignment("gisaid_covid19.aln", format = "fasta")
# Tidy the isolate names
a <- gsub("/2020.+", "", aln$nam)
b <- gsub("hCoV-19/Canada/", "", a)
b
aln$nam <- b
length(aln$nam) #147
length(unique(aln$nam)) #146

which(duplicated(aln$nam)) #34
aln$nam[34]

?write.csv()
write.csv(aln$nam, file = "SARS-CoV2.csv")

# double names cause circularized dendogram names labeled 1-lenght(rownames)       
rm.sequence.fasta("gisaid_covid19.aln", "unique_gisaid_covid19.aln", to.rm = c("hCoV-19/Canada/ON-VIDO-01/2020|EPI_ISL_413015|2020-01-23"))

aln <- read.alignment("unique_gisaid_covid19.aln", format = "fasta")
aln$nam
# Tidy the isolate names
a <- gsub("/2020.+", "", aln$nam)
b <- gsub("hCoV-19/Canada/", "", a)
b
aln$nam <- b
length(aln$nam) #146
length(unique(aln$nam)) #146
?dist.alignment

aln[[1]]

# distances
d <- dist.alignment(aln, matrix=c("similarity"))
head(as.data.frame(d))
class(d)

# heatmap + cluster for individuals

# convert the distances to a dataframe 
z <- apply(as.data.frame(as.matrix(d)), 2, as.numeric)
dim(z)
rownames(z) <- colnames(z)


# 
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
png("pretty.heat.png", units="in", width=10, height=10, res=300)
pheatmap(z, fontsize_row = 3, fontsize_col = 3, cellwidth = 2, cellheight = 2) #cutree_cols  = 10, 
dev.off()

#Interactive heat maps, !!!! does not have a legend for color 
png("interactive.png", units="in", width=10, height=10, res=300)
d3heatmap(scale(z), colors = "RdYlBu",
          k_row = 4, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)
dev.off()
 
#Enhanced heat maps
png("enhanced_heat.png", units="in", width=10, height=10, res=300)
heatmap.2(z, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
dev.off()

# dissimilarity plot
png("dissimilarity.png", units="in", width=10, height=10, res=300)
fviz_dist(d, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"), lab_size = 4)
dev.off()

# Simply Dendogram
hc <- hclust(d)
dhc <- as.dendrogram(hc)

# 
plot(dhc)

# Horizontal plot
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")

# 
plot(dhc, xlab = "Height", nodePar = nodePar, horiz = TRUE)

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


