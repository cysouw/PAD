setwd("/Users/cysouw/Documents/Github/Alignment/PAD/scripts/")

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

nr_words <- rep(1:length(all),sapply(all,function(x){dim(x)[2]}))
cols <- unlist(sapply(all,colnames))
all <- do.call(cbind,all)
rows <- rownames(all)
words <- gsub("_.+\\.msa\\.V\\d+","",names(cols))
nr_cols <- gsub(".+\\.V(\\d+)$","\\1",names(cols),perl=T)
nr_cols <- as.numeric(nr_cols)-2

#=====

# library(stringi)

remove.diacritics <- function(x) {
	stringi::stri_replace_all_regex(
		stringi::stri_trans_nfd(x)
		, "\\p{Mn}"
		, ""
		)
}
all <- apply(all,2,remove.diacritics)
rownames(all) <- rows

all <- gsub("→", "", all)
all <- gsub("\\+", "", all)
all <- gsub("˻", "", all)
all <- gsub("ʿ", "", all)
all <- gsub("ʾ", "", all)

# correct cedilla

all <- gsub("c", "ç", all)

#=====

missing <- apply(all,2,function(x){sum(na.omit(x)=="-")+sum(is.na(x))})
ignore <- 120
all <- all[,missing < ignore]
cols <- cols[missing < ignore]
words <- words[missing < ignore]
nr_words <- nr_words[missing < ignore]
nr_cols <- nr_cols[missing < ignore]

colnames(all) <- paste(1:ncol(all), words, nr_cols, cols)

#=====

library(qlcMatrix)

# tensor lect x correspondence x character
# library(slam)
# source("/Users/cysouw/Documents/Github/R/slam/R/unfold.R")

# turn "slam"-matrix into Matrix
# as.Matrix <- function(x) {
#  sparseMatrix(i=x$i,j=x$j,x=x$v,dims=c(x$nrow, x$ncol))
# }

# coor <- expand.grid(1:nrow(all), 1:ncol(all))
# letters <- as.factor(all)
# coor <- cbind(coor,as.numeric(letters))
# coor <- coor[!is.na(coor[,3]),]

# coor <- as.matrix(coor)
# attr(coor,"dimnames") <- NULL

# ten <- simple_sparse_array(
#		as.matrix(coor)
#		, v = rep(1, nrow(coor))
#		, dim = c( nrow(all), ncol(all), length(levels(letters)) )
#		, dimnames = list( rownames(all), colnames(all), levels(letters))
#		)

# library(reshape2)
long <- function(DF) {
	reshape2::melt(as.matrix(DF)
		 , varnames = c("observation", "attribute")
		 , na.rm = TRUE
		 )
}

A <- Array(long(all))

# rollup

chr <- as.Matrix( as.simple_triplet_matrix( rollup(A, 1, FUN=sum, DROP=T) ) )

sim.chars <- cosSparse(chr)
sort(sim.chars[2,],decreasing=T)[1:10]

#==========

library(qlcMatrix)

all2 <- all
all2[all2=="-"] <- NA
sim <- sim.obs(t(all2))

# context ???
# does not really improve things

last <- diff(nr_words)
first <- head(c(1,last),-1)
U <- bandSparse(n = nrow(sim), k = 1) * 1

sim_before <- t(U) %*% sim %*% U
sim_before[first == 1, first == 1] <- 1

sim_after <- U %*% sim %*% t(U)
sim_after[last == 1, last == 1] <- 1

sim_context <- sim + sim_before/3 + sim_after/3

# =========
# villages

sim.loc <- sim.obs(all)
#plot(hclust(as.dist(max(sim.loc)-sim.loc),method="ward.D2"),cex=.7)

library(apcluster)
clusters <- apcluster(sim.loc)
hier <- aggExCluster(sim.loc,clusters)

cl.vec <- c()
for (i in 1:length(clusters@clusters)){cl.vec[clusters[[i]]] <- i}

library(qlcVisualize)
library(mapdata)

loc <- doculects[,c(1,2)]

map("worldHires","Germany",fill=T,col="grey90")
lmap(loc,cl.vec,draw=1:max(cl.vec),levels=c(0.25),size=.5,col="rainbow", position="topleft",add=T)

#plot(hier)

# ===

#library(igraph)

# g <- igraph::graph.adjacency(sim_context, mode="undirected", weighted = T)
# cl <- igraph::infomap.community(g)$membership

# sapply(1:max(na.omit(cl)),function(x){table(cols[cl== x])})

#=====

# library(apcluster)
# colnames(sim) <- cols
# tmp <- apcluster::apcluster(as.matrix(sim))
# cl <- c()
# for (i in 1:length(tmp@clusters)){cl[tmp[[i]]] <- i}

#====

# library(fpc)
# p <- pamk(as.dist(1-sim),krange=20:40,critout=T)$pamobject

library(cluster)
p <- pam(as.dist(1-sim),30)

# this is necessary because of reordering
cl <- p$clustering
sil <- silhouette(cl,as.dist(1-sim_context)) # this is an ERROR!!!
cl[sil[,"sil_width"] < 0.01] <- NA

# tmp <- sapply(1:max(na.omit(cl)),function(x){table(cols[which(cl==x)])})
# n <- sapply(tmp,function(x){paste(names(x),collapse="/")})

# sapply(1:max(na.omit(cl)),function(x){table(cols[cl== x])})

# ===

library(qlcVisualize)
library(seriation)

draw.cluster <- function(cluster
						, ncol = 6
						, order = "R2E"
						, member = cl
						, control = NULL
						, ...
						) {

	sel <- which(member ==  cluster)
	data <- all[,sel]
	colnames(data) <- paste(1:ncol(all), words, nr_cols, cols)[sel]

	limage(data
		, col = rainbow(ncol)
		, order =  order
		, cex.axis = 0.3
		, cex.legend = 0.7
		, cex.remaining = 0.2
		, ...
		)

}

seg <- function(h, lwd = 1) {
	segments(rep(0,length(h))
			, rep(h, length(h))
			, rep(183, length(h))
			, rep(h, length(h))
			, lwd = lwd
			)
}

# =====
m <- cl

draw.cluster(1,10,"MDS_angle") # kurz a

draw.cluster(2) # w/v
	seg(25) # ben
	m[c(254,323,2,561)] <- 31 # b/_en
	draw.cluster(2,4,"MDS_angle",m=m)
	seg(21) # be, but note "Mittwoch"!
	m[c(416,328,602)] <- 32 # b/_e

draw.cluster(3) # m
	seg(c(5,7,9))
	m[c(573,3,324,255,563,623)] <- 33 # n/be_
	m[c(169,623)] <- 34 # n/fe_
	m[c(664,711)] <- 35 # Dativ m
	draw.cluster(3,m=m)

draw.cluster(4,8,"MDS") # t/d final

draw.cluster(5,8) # kurz a
	seg(5)
	m[c(144,579,136,462,235)] <- 36 # a/t_g

draw.cluster(6) # ach
	seg(6)
	m[c(580,463,236,137,145)] <- 37 # g/ta_

draw.cluster(7,5) # t/d, mostly "st", but not consistently
	seg(c(4,31))
	m[c(670,23,384,93)] <- 38 # t/zero
	m[c(146,122,119,138,49,598,161,131,596,629,150,153,159,170,585,593,588,593,588,126,129,234,476,461,578,683,248,707,409)] <- 56 # d
	draw.cluster(7,5,m=m)

draw.cluster(8,12,"MDS_angle") # e/a ???

draw.cluster(9) # p/b
	seg(20)
	m[c(30,9)] <- 39 # pf

draw.cluster(10)
	seg(c(22,56)) # en, er, other schwa
	seg(c(23,24,65,66)) # these seem to be "er" cases
	m[c(103, 690, 298, 274, 493, 157, 622, 398, 287, 615, 54, 422, 489, 431, 562, 226, 645, 31, 528, 10, 62)] <- 40 # schwa/_[nl]
	m[c(94, 175, 66, 125, 592, 141, 684, 308, 410, 603, 82, 118, 477, 671, 337, 652, 618, 202, 503, 625, 368, 543, 348, 217)] <- 41 # schwa/_r

	draw.cluster(10,6,"MDS_angle",m)
	seg(20)
	m[c(283,165,250,277,271,318,311,264,290,257,165,178)] <- 42 # schwa/[gb]_ (Prefix)

	draw.cluster(10,4,"MDS",m=m)
	draw.cluster(40,4,"MDS",m=m)
	draw.cluster(41,4,"MDS",m=m)
	draw.cluster(42,4,"MDS",m=m)

draw.cluster(11,5) # l

draw.cluster(12,12) # f/v
	seg(10)
	m[c(492,397,502,168,191,621,40,549,149,13)] <- 43 # f/p Note voicing!!!
	draw.cluster(12,m=m)
	draw.cluster(43,m=m)
	m[c(240,200,719)] <- 44 # f final

draw.cluster(13,11,"MDS_angle") # s
	seg(c(15,17,25,33))
	seg(c(4,6,7,11))
	m[c(644,488)] <- 45 # chs
	m[c(375,195,469,418)] <- 46 # s vs. zero (possible two groups: final s loss and s before a t)
	m[c(341,473,121,660,72,651,81,336)] <- 47 # s vs. t
	m[c(134,20,306,213,541,378)] <- 48 # s vs sch, typicaly before t
	draw.cluster(13,6,"MDS",m=m)

draw.cluster(14) # n/[g ch]e_
	seg(6)
	m[c(27,52)] <- 49 # ng

draw.cluster(15,6) # e/ö ???

draw.cluster(16) # o
	seg(4)
	m[116] <- 50 # brUder

draw.cluster(17) # och vs. k
	seg(1)
	m[42] <- 51 # auGenblick

draw.cluster(18) # u/o ???

draw.cluster(19) # i

draw.cluster(20) # au

draw.cluster(21) # r
	seg(14)
	m[c(635,76,148,322,163)] <- 52 # r/V_C

draw.cluster(22) # n, Lot's of differences, but no clear groups
	seg(21)

draw.cluster(23,8) # ei
	seg(11)
	m[c(412, 346, 216, 485)] <- 53 # eu

draw.cluster(24) # g initial
	seg(c(10,13,17))
	seg(c(15,16),0.3)
	m[c(310,282,317,270,276,186)] <- 54 # g/_e[sf]
	m[c(256,263,289)] <- 55 # g/_e[bk]
	m[c(225,421,77)] <- 51 # g/_e(n)

draw.cluster(25,o="MDS_angle") # sch

draw.cluster(26,o="MDS") # ich

draw.cluster(27) # o

draw.cluster(28) # ts

draw.cluster(29) # kh

draw.cluster(30) # h


# ======
# manual decision on difficult cases (only consonants)

# which(is.na(cl))

get.best <- function(corr,range=10,levels=.3){
	print(names(cl)[corr])
	table <- m[names(sort(sim[,corr],decreasing=T)[1:range])]
	print(table)
	print(sort(table(table),decreasing=T))
	tmp <- m[corr]
	m[corr] <<- 100
	plot.corr(100,draw=6,levels=levels,main=names(cl)[corr],size=0.6)
	m[corr] <<- tmp
}

#get.best(7)
m[7] <- 4

#get.best(38)
m[38] <- 17

#get.best(47)
m[47] <- 29

#get.best(53)
m[53] <- 29

#get.best(59)
#get.best(95)
m[95] <- 52

#get.best(99)
m[99] <- 12

#get.best(117)
m[117] <- 38

#get.best(294)
m[294] <- 4

#get.best(344)
m[344] <- 4

#get.best(357,20)
m[357] <- 47

#get.best(391,20)
m[391] <- 4

#get.best(413,20)
m[413] <- 38

#get.best(427,20)
m[427] <- 4

#get.best(500)
m[500] <- 39

#get.best(504)
m[504] <- 39

#get.best(511)
m[511] <- 4

#get.best(546)
m[546] <- 45

#get.best(584)
m[584] <- 29

#get.best(697)
m[697] <- 4

#get.best(701)
m[701] <- 4

# sound change models

library(corHMM)
tree <- nj(as.dist(1-sim.loc))

which(m==47)
tmp <- all[,na.omit(which(m==47))]
(sort(table(tmp))->tab)

# 43
tmp <- gsub("ː","",tmp)
tmp <- gsub("ʰ","",tmp)
tmp <- gsub("h","",tmp)
tmp <- gsub("ʱ","",tmp)
tmp <- gsub("ʼ","",tmp)
tmp <- gsub("β","v",tmp)
tmp <- gsub("ɸ","f",tmp)
tmp <- gsub("ff","f",tmp)
tmp <- gsub("fv","f",tmp)

# 47
tmp <- gsub("ː","",tmp)
tmp <- gsub("ʰ","",tmp)
tmp <- gsub("h","",tmp)
tmp <- gsub("ʱ","",tmp)
tmp <- gsub("ʼ","",tmp)
tmp <- gsub("ʔ","t",tmp)
tmp <- gsub("ð","d",tmp)
tmp <- gsub("θ","t",tmp)
tmp <- gsub("ss","s",tmp)
tmp <- gsub("zs","z",tmp)
tmp <- gsub("zz","z",tmp)
tmp <- gsub("sz","s",tmp)
tmp <- gsub("tt","t",tmp)
tmp <- gsub("ɾ","d",tmp)
tmp <- gsub("ɹ","d",tmp)
tmp <- gsub("ɻ","d",tmp)
tmp <- gsub("ʒ","z",tmp)
tmp <- gsub("ʈ","t",tmp)

(sort(table(tmp))->tab)
(rem <- names(tab[tab<100]))

for (i in rem) {
	tmp[tmp==i] <- NA
}

sort(table(tmp))
tmp <- as.matrix(tmp)
tmp[is.na(tmp)] <- "?"
tmp <- cbind(location=rownames(tmp),tmp)

models <- sapply(1:(ncol(tmp)-1),function(x){rayDISC(tree,tmp,charnum=x,model="ARD",node.states="marginal")},simplify=F)

# ====

# library(raster)
# deu <- getData("GADM", country = "DEU", level = 0)
# plot(deu, col = "grey90")


loc <- doculects[,c(1,2)]
library(mapdata)

plot.corr <- function(c, levels, main = "", ...) {

	tmp <- A[ , which(m == c),]
	tmp <- rollup(tmp, 2, FUN = sum, DROP = T)
	tmp <- as.matrix(as.simple_triplet_matrix(tmp))
	tmp <- tmp[ , colSums(tmp) > 0]

	if (length(which(m==c))==1) {
		cex=.7
	} else {
		cex=.1
	}

	map("worldHires","Germany",fill=T,col="grey90")
	lmap(loc, tmp, cex = cex, add = T, position =  "topleft", levels = levels, ...)

	text(14.5, 48, cex = 0.5, labels=paste("Level =",100*levels, "%"))
	title(main = main)
}

# ====

plot.corr(6, draw = 4, levels = 0.30, size = 0.7, main = "Correspondences “x”")
plot.corr(37, draw = 5, levels = 0.20, size = 0.6, main = "Correspondences “(ta)g”")
plot.corr(17, draw = 6, levels = 0.20, size = 0.6, main = "Correspondences “k/x”")
plot.corr(26, draw = 6, levels = 0.20, size = 0.6, main = "Correspondences “ç”")
plot.corr(29, draw = 4, levels = 0.10, size = 0.7, main = "Correspondences “k”")
plot.corr(24, draw = 4, levels = 0.15, size = 0.7, main = "Correspondences “k/g”")
plot.corr(45, draw = 4, levels = 0.25, size = 0.7, main = "Correspondences “ks”")
plot.corr(51, draw = 6, levels = 0.20, size = 0.6, main = "Correspondences “g loss”")
plot.corr(54, draw = 5, levels = 0.25, size = 0.6, main = "Correspondences “ge prefix”")
plot.corr(55, draw = 5, levels = 0.25, size = 0.6, main = "Correspondences “ge prefix”")


plot.corr(2, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “β/v”")
plot.corr(31, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “b before en”")
plot.corr(32, draw = 6, levels = 0.15, size = 0.5, main = "Correspondences “β/p”")


plot.corr(4, draw = 4, levels = 0.15, size = 0.7, main =  "Correspondences “t final”")
plot.corr(7, draw = 3, levels = 0.30, size = 0.7, main =  "Correspondences “t”")
plot.corr(38, draw = 4, levels = 0.25, size = 0.7, main = "Correspondences “t loss”")
plot.corr(56, draw = 4, levels = 0.25, size = 0.7, main = "Correspondences “t/d”")
plot.corr(28, draw = 4, levels = 0.25, size = 0.7, main = "Correspondences “ts”")


plot.corr(9, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “p/b”")
plot.corr(12, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “f”")
plot.corr(43, draw = 4, levels = 0.20, size = 0.7, main = "Correspondences “f/p”")
plot.corr(44, draw = 7, levels = 0.15, size = 0.5, main = "Correspondences “f final”")
plot.corr(39, draw = 5, levels = 0.20, size = 0.6, main = "Correspondences “pf”")


plot.corr(13, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “s/z”")
plot.corr(45, draw = 5, levels = 0.25, size = 0.6, main = "Correspondences “ks”")
plot.corr(46, draw = 4, levels = 0.20, size = 0.7, main = "Correspondences “s loss”")
plot.corr(47, draw = 5, levels = 0.15, size = 0.6, main = "Correspondences “s/t”")
plot.corr(48, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “s before t”")
plot.corr(25, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “sch”")


plot.corr(21, draw = 5, levels = 0.30, size = 0.6, main = "Correspondences “r”")
plot.corr(52, draw = 5, levels = 0.25, size = 0.6, main = "Correspondences “r loss”")
plot.corr(30, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “h”")
plot.corr(11, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “l”")


plot.corr(3, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “m”")
plot.corr(33, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “n after be”")
plot.corr(34, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “n after fe”")
plot.corr(35, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “Dativ m”")
plot.corr(22, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “n”")
plot.corr(14, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “n loss”")
plot.corr(49, draw = 3, levels = 0.30, size = 0.7, main = "Correspondences “ŋ”")


plot.corr(20, draw = 6, levels = 0.17, size = 0.6, main = "Correspondences “au”")
plot.corr(23, draw = 5, levels = 0.20, size = 0.6, main = "Correspondences “ei”")
plot.corr(53, draw = 8, levels = 0.15, size = 0.5, main = "Correspondences “eu”")


plot.corr(10, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “schwa”")
plot.corr(40, draw = 3, levels = 0.25, size = 0.7, main = "Correspondences “schwa before n/l”")
plot.corr(41, draw = 4, levels = 0.20, size = 0.7, main = "Correspondences “schwa before r”")
plot.corr(42, draw = 4, levels = 0.25, size = 0.7, main = "Correspondences “schwa in prefix”")


# monophtonghs

plot.corr(1, draw = 5, levels = 0.15)
plot.corr(5, draw = 5, levels = 0.15)
plot.corr(8, draw = 3, levels = 0.15)
plot.corr(15, draw = 4, levels = 0.30)
plot.corr(16, draw = 5, levels = 0.15)
plot.corr(18, draw = 5, levels = 0.15)
plot.corr(19, draw = 5, levels = 0.15)
plot.corr(27, draw = 5, levels = 0.15)




