doculects <- read.delim("../data/doculects.csv"
					, comment.char = "#"
					, row.names = 1
					, quote = ""
					)

read.PAD <- function(file) {
	a <- read.table(file
			, header = FALSE
			, quote = ""
			, comment.char = "#"
			, sep = "\t")
	
	a[,2] <- gsub(".", "", a[,2], fixed = T)
	
	getAnno <- a[,1] == ":ANN"
	annotation <- a[getAnno,]
	standard <- as.matrix(a[a[,2]=="STANDARD",-c(1,2)])[1,]
	standard <- tolower(standard)
	
	a <- a[!getAnno,]
	a <- as.matrix(a[!duplicated(a[,2]),])
	rownames(a) <- a[,2]
	a <- a[,-c(1,2)]
	
	a <- t(sapply(as.character(doculects$SHORT_NAME), function(x){
		if (sum(rownames(a) == x) == 0) {
			rep(NA, times = dim(a)[2])
		} else {
			a[x,]
		}
	}))
	
	colnames(a) <- standard	
	return(a)
	
}

setwd("../alignments/merged")
files <- list.files()
all <- sapply(files,read.PAD)
setwd("../../scripts")

cols <- unlist(sapply(all,colnames))
all <- do.call(cbind,all)
rows <- rownames(all)
words <- gsub("_.+\\.msa\\.V\\d+","",names(cols))
nr <- gsub(".+\\.V(\\d+)$","\\1",names(cols),perl=T)
nr <- as.numeric(nr)-2

#=====

library(stringi)

remove.diacritics <- function(x) {
	stri_replace_all_regex(stri_trans_nfd(x),"\\p{Diacritic}","")
}
all <- apply(all,2,remove.diacritics)
rownames(all) <- rows

all <- gsub("→", "", all)
all <- gsub("\\+", "", all)
all <- gsub("˻", "", all)
all <- gsub("ʿ", "", all)
all <- gsub("ʾ", "", all)

#=====

missing <- apply(all,2,function(x){sum(na.omit(x)=="-")+sum(is.na(x))})
all <- all[,missing<100]
cols <- cols[missing<100]
words <- words[missing<100]
nr <- nr[missing<100]

#=====

library(qlcMatrix)
library(igraph)

all2 <- all
all2[all2=="-"] <- NA
sim <- sim.obs(t(all2))
# sim@x[sim@x<0] <- 0
# sim <- drop0(sim)
# plot(hclust(as.dist(1-sim)),labels=rownames(sim),cex=.3)

g <- graph.adjacency(sim, mode="undirected", weighted = T)
cl <- infomap.community(g)
# sapply(1:max(cl$membership),function(x){table(cols[cl$membership == x])})

#=====
library(seriation)

source("plot.levels.R")

draw.cluster <- function(cluster, ncol = 6, order = "R2E") {

	sel <- cl$membership == cluster
	data <- all[,sel]
	colnames(data) <- paste(words, nr, cols)[sel]

	sim.cols <- sim.obs(t(data))
	sim.rows <- sim.obs(data)
	
	if (order == "pca") {
		
		# PCA on similarities (aka "correspondence analysis")
		order.cols <- order(prcomp(as.matrix(sim)[sel,sel])$x[,1])
		order.rows <- order(prcomp(as.matrix(sim.rows))$x[,1])
		
	} else if (order == "varimax") {
		
		# varimax rotations for ordering of correspondences
		order.cols <- order(varimax(
						prcomp(as.matrix(sim.cols))$x)$loadings[,1])
		order.rows <- order(varimax(
						prcomp(as.matrix(sim.sel))$x)$loadings[,1])
						
	} else if (order == "mds") {
		
		# classic MDS
		order.cols <- order(cmdscale(as.matrix(1-sim.cols))[,1])
		order.rows <- order(cmdscale(as.matrix(1-sim.rows))[,1])
		
	} else {
		# seriation: "R2E" works nice
		library(seriation)
		order.cols <- get_order(seriate(
						as.dist(as.matrix(sim.cols))
						, method = order
						))
		order.rows <- get_order(seriate(
						as.dist(as.matrix(sim.rows))
						, method = order
						))
	}
	
	# for visualitation of ordering differences use:
	# pimage(as.matrix(sim.cols)[order.cols,order.cols])
	# pimage(as.matrix(sim.rows)[order.rows,order.rows])

	plot.levels(data[order.rows,order.cols]
				, col = rainbow(ncol)
				, cex.axis=0.3
				, cex.legend = 0.7
				)
				
}

draw.cluster(1,10)


# === reordering of clusters


# split "ts" from "t/d" cluster
cl$membership[cl$membership==3][c(15,23,35,50,70,78,79,80,81)] <- 15
draw.cluster(3,order = "R2E")
draw.cluster(15,order = "R2E")
