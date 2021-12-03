
# SNMRS -------------------------------------------------------------------

rm(list=ls())

snmrs=function(x,y){
  sum=0;
  xm=0;
  ym=0;
  xnm=0;
  ynm=0;
  xnmm=0;
  ynmm=0;
  width=length(x)
  for (i in 1:width){
    xm=xm+x[i];
    ym=ym+y[i];
  }
  xm=xm/width;
  ym=ym/width;
  for (i in 1:width){
    sum=sum+abs(((x[i]-xm)-(y[i]-ym)))
    xnm=xnm+abs(x[i]-xm)
    ynm=ynm+abs(y[i]-ym)
  }
  diff<-abs(xnm-ynm)
  min<-min(xnm,ynm)
  ss<-(sum-diff)
  z<-1-((sum-diff)/(2*min))
  return(z)
}
snmrs_corr_upd = function (dataset) {
  SNMRS_MAT=matrix(data=NA,nrow = nrow(dataset), ncol = nrow(dataset));
  rownames(SNMRS_MAT)<-rownames(dataset)
  colnames(SNMRS_MAT)<-rownames(dataset)
  l <- length ( dataset[,1] )
  for (i in 1:l) {
    for (j in 1:l) {
      SNMRS_MAT[i,j] <- snmrs(t(dataset[i,]), t(dataset[j,]))
    }
  }
  return(SNMRS_MAT)
}

#####################SNMRS for iris data################################
BiocManager::install("NbClust")
library(NbClust)
data <- iris[, -5]
dim(data)


dist<-round(snmrs_corr_upd(data),4)
diag(dist)<-FALSE
dist<-1-dist

#internal validity - Average hierarchical clustering method
diss_matrix<-as.dist(dist,diag = FALSE,upper=FALSE)
nb<-NbClust(data, diss = diss_matrix, distance = NULL, min.nc = 2,max.nc = 20, method = "average", index = "alllong")
nb$All.index
setwd("F:/phd/1. SNMRS/iris")
write.csv(nb$All.index,"snmrs_iris.csv")
#internal validity - kmeans clustering method
diss_matrix<-as.dist(dist)
nb<-NbClust(data, diss = diss_matrix, distance = NULL, min.nc = 2,max.nc = 20, method = "kmeans", index = "alllong")
nb$All.index
setwd("F:/phd/1. SNMRS/iris/kmeans/")
write.csv(nb$All.index,"snmrs_iris_kmeans.csv")



# SNMRS for Yeast Sporulation data ----------------------------------------



library("dynamicTreeCut")
library("WGCNA")


setwd("F:/phd/1. SNMRS/YeastS/")
gexp <- read.csv("yeast_sporulation.csv")
rownames(gexp)<-make.names(gexp[,1], unique = TRUE)
gexp<-gexp[,-1]
dim(gexp)
x <- gexp[complete.cases(gexp), ]
dim(x)
gexp<-round(snmrs_corr_upd(x),4)
gexp<-gexp[,-1]
gexp[gexp<0.8]<-0
dist<-1-gexp
diag(dist)<-FALSE
s<-NbClust(data=NULL, diss = diss_matrix, distance = NULL, min.nc = 2,max.nc = 2, method = "average", index = "silhouette")
s$Best.nc

#Clustering
adjacency<-as.dist(dist,diag = FALSE,upper=FALSE)
geneTree = hclust(as.dist(adjacency), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on SNMRS dissimilarity",
     labels = FALSE, hang = 0.04);
minModuleSize = 30;

dynamicMods = cutreeDynamic(dendro = geneTree, distM = adjacency,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#Extraction of module
unique(dynamicColors)
table(dynamicColors)##########see the size of clusters

