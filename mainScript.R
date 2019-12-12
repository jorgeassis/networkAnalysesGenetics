
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

library(PopGenReport)
library(reshape2)
library(igraph)
library(mmod)
library(RColorBrewer)
library(hierfstat)
library(doParallel)

nClusters <- 6

distinctColors <- function(n) {
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  return((col_vector[1:n]))
  
}

setwd("/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera")
data <- "/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Genetic/All Genetic Data 6L Genepop Clean Final Genetix.gtx"

# ------------------------

data <- read.genetix(data)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

sharedDistance <- pairwise.propShared(data)
comb <- melt(as.matrix(sharedDistance))

D <- pairwise_D(data)
comb <- melt(as.matrix(D))
comb[,3] <- 1 - comb[,3]

Fst <- pairwise.fst(data)
comb <- melt(as.matrix(Fst))
comb[,3] <- 1 - comb[,3] 

# ---------------------------------------------------------
# ---------------------------------------------------------

clustering.method <- "leading.eigenvector.community" # leading.eigenvector.community walktrap.community

comb <- comb[ comb$Var1 != comb$Var2 ,]
comb <- as.data.frame( comb[ sort(comb[,"value"] , decreasing = TRUE, index.return =TRUE)$ix , c("Var1","Var2","value")] )
comb[comb[,3] < 0 ,3] <- 0
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
V(graph.obj)$name <- as.numeric(gsub("Pop","",names(V(graph.obj))))

graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight == 0))
graph.obj <- simplify(graph.obj)

# ------------------------
# ------------------------

removeVertex <- c(74)
clusterVertex <- 1

which(membership.graph[sorting] == clusterVertex)
subgraph.vector <- which(membership.graph[sorting] == clusterVertex)
subgraph.vector <- subgraph.vector[! subgraph.vector %in% removeVertex]
graph.obj <- induced.subgraph(graph.obj, sapply(subgraph.vector,function(x){ which(as.numeric(V(graph.obj)$name) == x )}))

# ------------------------
# ------------------------

pruned.graph.obj <- graph.obj
modularityTest.i <- 0

repeat {
  pruned.graph.obj.t <- delete.edges(pruned.graph.obj, which.min(E(pruned.graph.obj)$weight) )
  if( length(unique((clusters(pruned.graph.obj.t)$membership))) > 1 ) { break }
  pruned.graph.obj <- pruned.graph.obj.t
  modularityTest.i <- modularityTest.i + 1
}

edgeWeights <- sort(E(graph.obj)$weight,decreasing = F)[1:modularityTest.i]

Cluster <- makeCluster( nClusters )
registerDoParallel(Cluster)
modularityTest <- foreach(w=edgeWeights, .verbose=FALSE, .combine=c, .packages=c("igraph")) %dopar% {
  pruned.graph.obj.t <- delete.edges(graph.obj, which( E(graph.obj)$weight <= w) )
  return(modularity(pruned.graph.obj.t, get(clustering.method)(pruned.graph.obj.t)$membership))
}
stopCluster(Cluster) ; rm(Cluster)
plot(1:length(modularityTest),modularityTest)

pruned.graph <- delete.edges(graph.obj, which( E(graph.obj)$weight <= max(edgeWeights[1:which.max(modularityTest)]) ) )
modularity(pruned.graph, get(clustering.method)(pruned.graph)$membership)

# ------------------------

# pruned.graph <- graph.obj
clustering.graph <- get(clustering.method)(pruned.graph)
membership.graph <- clustering.graph$membership
sorting <- sort(as.numeric(V(pruned.graph)$name),index.return=T)$ix
membership.graph[sorting]


# Size 10:10
# https://www.r-graph-gallery.com/248-igraph-plotting-parameters.html
# https://medialab.github.io/iwanthue/

netColors <- c("#aaaea5","#3e413a")
membershipLevel1 <- netColors[membership.graph]
set.seed(6)
plot(pruned.graph, vertex.size=7, vertex.color=membershipLevel1,layout=layout.fruchterman.reingold , vertex.label=NA )
plot(pruned.graph, vertex.size=8, vertex.color=membershipLevel1,layout=layout.fruchterman.reingold , vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label=as.numeric(gsub("Pop","",names(V(pruned.graph)))) )

netColors <- distinctColors(length(unique(membership.graph)))
membershipLevel2 <- netColors[membership.graph]
set.seed(6)
plot(pruned.graph, vertex.size=7, vertex.color=membershipLevel2,layout=layout.fruchterman.reingold , vertex.label=NA )
plot(pruned.graph, vertex.size=7, vertex.color=membershipLevel2,layout=layout.fruchterman.reingold , vertex.label=as.numeric(gsub("Pop","",names(V(pruned.graph)))) )

# Number of clusters
length(unique(membership.graph))

sort(as.numeric(V(pruned.graph.obj)$name[membership.graph == 1]))
sort(as.numeric(V(pruned.graph.obj)$name[membership.graph == 2]))
sort(as.numeric(V(pruned.graph.obj)$name[membership.graph == 3]))
sort(as.numeric(V(pruned.graph.obj)$name[membership.graph == 4]))

vertexCentrality <- sort(eigen_centrality(graph.obj)$vector)
sort(vertexCentrality)
vertexCentrality[vertexCentrality > as.numeric(quantile(vertexCentrality,prob=0.95))]

# -------------

n.cells <- length(get(clustering.method)(pruned.graph.obj)$membership)
unique.members <- sort(unique(get(clustering.method)(pruned.graph.obj)$membership))

permutat <- 9999
mods <- sapply(1:permutat, function(i){
  modularity( pruned.graph.obj , sample( unique.members , n.cells , replace = TRUE) )
})

signif(sum( mods > modularity(pruned.graph.obj,get(clustering.method)(pruned.graph.obj)$membership) ) / permutat,digits=nchar(permutat))

# -------------
