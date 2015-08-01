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
ignore <- 120
all <- all[,missing < ignore]
cols <- cols[missing < ignore]
words <- words[missing < ignore]
nr_words <- nr_words[missing < ignore]
nr_cols <- nr_cols[missing < ignore]

colnames(all) <- paste(1:ncol(all), words, nr_cols, cols)

#=====

# tensor lect x correspondence x character

library(qlcMatrix)
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

library(reshape2)

long <- function(DF) {
	melt(as.matrix(DF)
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

# sim@x[sim@x<0] <- 0
# sim <- drop0(sim)
# plot(hclust(as.dist(1-sim)),labels=rownames(sim),cex=.3)

# ===

#library(igraph)
#g <- graph.adjacency(sim, mode="undirected", weighted = T)
#cl <- infomap.community(g)$membership
#sapply(1:max(cl),function(x){table(cols[cl== x])})

#=====

#library(apcluster)
#colnames(sim) <- cols
#tmp <- apcluster(as.matrix(sim))
#cl <- c()
#for (i in 1:length(tmp@clusters)){cl[tmp[[i]]] <- i}

#====

#library(fpc)
#p <- pamk(as.dist(1-sim),krange=10:30,critout=T)$pamobject

library(cluster)
p <- pam(as.dist(1-sim),24)

# this is necessary because of reordering
cl <- p$clustering
sil <- silhouette(cl,as.dist(1-sim))
cl[sil[,"sil_width"] < 0.01] <- NA

tmp <- sapply(1:max(na.omit(cl)),function(x){table(cols[which(cl==x)])})
n <- sapply(tmp,function(x){paste(names(x),collapse="/")})

# ===
library(qlcVisualize)
library(seriation)

draw.cluster <- function(cluster
						, ncol = 6
						, order = "R2E"
						, member = cl
						, control = NULL
						, screen.width = 12
						) {
	sel <- which(member ==  cluster)
	data <- all[,sel]
	colnames(data) <- paste(1:ncol(all), words, nr_cols, cols)[sel]

	
#	quartz(width=screen.width,height=2+(screen.width-4)*ncol(data)/nrow(data))

	limage(data
		, col = rainbow(ncol)
		, order =  order
		, show.remaining = F
		, cex.axis=0.3
		, cex.legend = 0.7
		, cex.remaining = .2
		)

}

seg <- function(h, lwd) {
	segments(rep(0,length(h))
			, rep(h, length(h))
			, rep(183, length(h))
			, rep(h, length(h))
			, lwd = lwd
			)
}


draw.cluster(14,10,"R2E") # u
draw.cluster(1,12,"MDS_angle") # o
draw.cluster(5,12,"MDS_angle") # a
draw.cluster(7,12,"MDS") # e
draw.cluster(13,12,"R2E") # i

draw.cluster(9,11) # schwa
	seg(c(22,56),1) # en, er, other schwa
	seg(c(23,24,65,66),0.5) # these seem to be "er" cases

draw.cluster(16,10) # au
	seg(c(3),1) # au/u, au/o
draw.cluster(19,10) # ei/eu
	seg(c(4,9),1) # ei/i, eu, ei/e

draw.cluster(8,4,"MDS_angle") # bp
	seg(4,1) # b/p, p/pf
draw.cluster(11) # fv
	seg(10,1) # p/f, f

draw.cluster(4) # dt
	seg(c(20,51,59),1) # t/th (final) vs t/0 vs d/t vs t
draw.cluster(23,6,"MDS") # t>ts

draw.cluster(15,8,"MDS_angle") # k
	seg(c(9,12),1) # inital, final, ach
draw.cluster(20) # g
	seg(c(10,13,17),1) # ges/gef-, geb/gek-, -ge/gen, initial
	seg(c(15,16),0.5)
draw.cluster(22,9,"R2E") # ich
draw.cluster(6,12) # ach
	seg(7,1) # ach vs. tag

draw.cluster(12,10) # sz
	seg(c(5,7,23,27,36),1) # ks, t>s, s
draw.cluster(21,6,"MDS") # sch

draw.cluster(3) # m
	seg(7,1) # ben
draw.cluster(18,5,"MDS") # n
	seg(c(7,45),1) # ken/gen/chen , nt/nd/schn/n-syllable initial
draw.cluster(10,6,"MDS_angle") # l
	seg(c(1,18,37,38),1) # initial/gl/fl/bl, final/ld/lb/lt
draw.cluster(17,6,"MDS_angle") # r
	seg(c(6),1) # initial vs. final
draw.cluster(24,3,"MDS") # h
draw.cluster(2,6,"MDS_angle") # vw
	seg(25,1) # ben

# === manual

m <- cl

#------
# schwa

# schwa + N or L
25 -> m[c(103, 690, 298, 274, 493, 157, 622, 398, 287, 615, 54, 422, 489, 431, 562, 226, 645, 31, 528, 10, 62)]

# schwa + R
26 -> m[c(94, 175, 66, 125, 592, 141, 684, 308, 410, 603, 82, 118, 477, 671, 337, 652, 618, 202, 503, 625, 368, 543, 348, 217)]

draw.cluster(9,6,"MDS_angle",m)

# GB + schwa prefix
27 -> m[c(283,165,250,277,271,318,311,264,290,257,165,178)]

# 675 wievIEl wrongly categorized
13 -> m[675]

draw.cluster(9,6,"MDS",m)
draw.cluster(25,6,"MDS",m)
draw.cluster(26,8,"MDS_angle",m)
draw.cluster(27,8,"MDS_angle",m)

#------
# au

# au/o
28 -> m[c(41, 620, 327)]

draw.cluster(16,12,"MDS",m)
draw.cluster(28,12,"MDS_angle",m)

#----
# ei

# ei/e
29 -> m[c(548,356,220,174)]

# eu
30 -> m[c(412,482,346,485,216)]

draw.cluster(19,14,"MDS_angle",m)
draw.cluster(29,12,"MDS_angle",m)
draw.cluster(30,12,"MDS_angle",m)

#----
# b/p

# pf
31 -> m[c(30,9,504,500)]

draw.cluster(8,4,"MDS",m)
draw.cluster(31,4,"MDS_angle",m)

#-----
# f

# f/p
32 -> m[c(13,149,549,40,621,191,168,502,397,492)]

draw.cluster(11,4,"MDS",m)
draw.cluster(32,6,"MDS_angle",m)

#------
# t/d

# t-d
33 -> m[c(511,294,697,113,427,574,701,4,391,262,89,206,269,595,173)]

# t-ø
34 -> m[c(210,466,413,23,384,670,93,59)]

draw.cluster(33,6,"MDS",m)
draw.cluster(34,6,"MDS_angle",m)























# ===

m <- cl

pf <- c(30,9,500,504)
m[pf] <- 30
draw.cluster(8,4,"MDS",member=m)

b_en <- c(254,323,2,561)
m[b_en] <- 31
draw.cluster(2,6,"MDS_angle",member=m)

dativ <- c(711,664)
m[dativ] <- 32

be_n <- c(573, 3, 255, 324, 563, 623, 43)
m[be_n] <- 18
draw.cluster(18,7,"R2E",member=m)


chs <- c(644, 488)
std <- c(336, 81, 651, 72, 660, 121, 341, 473, 357, 418, 195, 469, 375)

m[chs] <- 40
m[std] <- 41
draw.cluster(12,10,"MDS",m=m)
draw.cluster(41,8,"MDS_angle",m=m)

dt <- c(511,294,697,113,427,574,701,4,391,262,89,206,269,595,173)
d0 <- c(210,466,413,23,384,670,93)
m[dt] <- 50
m[d0] <- 51
draw.cluster(4,6,"MDS_angle",m) # dt

#===
kgx <- cl
kgx[kgx==15] <- 6
kgx[kgx==20] <- 6
kgx[kgx==22] <- 6
draw.cluster(6,8,"R2E",kgx)

# ===

#colnames(all) <- paste(words, nr_cols, cols)
#quartz(width=12,height=48)
#limage(all, rainbow(60), order="GW", show=T, cex.axis=.3,cex.remaining=.2,font="CharisSIL")

# ===




draw.cluster(1,8,"R2E") # schwa
draw.cluster(2,6,"mds") # t/d
draw.cluster(3,6,"R2E") # n
draw.cluster(4,10,"R2E") # i/e
draw.cluster(5,10,"R2E") # o/u/a
draw.cluster(6,8,"R2E") # k/x/g
draw.cluster(7,4,"R2E")	# l
draw.cluster(8,8,"R2E") # f/p/v
draw.cluster(9,10,"R2E") # s/z
draw.cluster(10,6,"mds") # b/v ~ w
draw.cluster(11,4,"mds") # p/b/pf
draw.cluster(12,4,"mds") # sch
draw.cluster(13,8,"R2E") # r
draw.cluster(14,4,"mds") # m
draw.cluster(15,4,"mds") # h

#=======

seg <- function(h, lwd) {
	segments(rep(0,length(h))
			, h
			, rep(183, length(h))
			, h
			, lwd = lwd
			)
}

# == fpb ===

fpb <- cl$membership
fpb[fpb==11] <- 8
draw.cluster(8,8,"R2E",fpb)

seg(c(4,25,49),2)
seg(c(2,24,47),0.5)

text(-4,54,"*p")
text(-4,35,"*f")
text(-4,15,"*b")
text(-4,2,"*pf")

# === fpb++ with *w ===

fpb[fpb==10] <- 8
draw.cluster(8,10,"R2E",fpb)

# === schwa ===

draw.cluster(1,8,"R2E")

seg(c(10,32,55),2)
seg(c(4,5,6,7,51,65,71,72),0.5)

text(-4,74,"*ə",cex=.6)
text(-4,60,"*ə/g_",cex=.6)
text(-4,42,"*ə/_[nl]",cex=.6)
text(-4,22,"*ə/_r",cex=.6)
text(-4,3,"*ə",cex=.6)

# === tdz ====

td <- cl$membership
z <- which(colnames(all)=="z")
td[z] <- 20

draw.cluster(20,8,"mds",td)

draw.cluster(2,8,"R2E",td)

seg(c(17,60),2)
seg(c(1,21,27,29,49),0.5)

text(-4,70,"*t/_#",cex=.5)
text(-4,54,"*d/V_V",cex=.5)
text(-4,35,"*d",cex=.5)
text(-4,8,"*t/[sx]_",cex=.5)
text(-4,10,"*-t",cex=.5)

# === sz ===

sz <- cl$membership
sz[sz == 12] <- 9
draw.cluster(9,10,"R2E",sz)
