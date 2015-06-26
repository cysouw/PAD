plot.levels <- function(m, col, cex.axis = 1, cex.legend = 1) {

	aspect <- dim(m)[2]/dim(m)[1]
#	dev.new(width=10,height=10*aspect)
	plot.new()
	plot.window(xlim=c(0,dim(m)[1])
				, ylim=c(0,dim(m)[2])
#				, asp = 1 + aspect
				)
				
	axis(1
		, at = c(1:dim(m)[1])-0.5
		, labels = rownames(m)
		, tick = FALSE
		, las = 2
		, cex.axis = cex.axis
		, mgp = c(3,0,0)
		)
	axis(2
		, at = c(1:dim(m)[2])-0.5
		, labels = colnames(m)
		, tick = FALSE
		, las = 2
		, cex.axis = cex.axis
		, mgp = c(3,0,0)
		)

	points(which(is.na(m), arr.ind = T)-0.5
			, pch = 20
			, col = "grey"
			, cex = .2
			)

	if (is.list(col)) {
		levs <- names(col)
		col <- unlist(col)
	} else {
		levs <- names(sort(table(as.vector(m)),decreasing = T))
		if (length(levs) > length(col)) {
			levs <- levs[1:length(col)]
		} else {
			col <- col[1:length(levs)]
		}		
		names(col) <- levs
	}

# trick from http://stackoverflow.com/questions/15627674/efficiency-of-drawing-rectangles-on-image-matrix-in-r

	cuts <- function(x) {
	  n <- length(x) %/% 4
	  map <- rep(c(rep(TRUE,4),FALSE), n)
	  result <- rep(NA, n*5)
	  result[map] <- x
	  result
	}

	for (i in levs) {
		todo <- which(m == i, arr.ind=T)
		todoX <- cbind(todo[,1]-1, todo[,1]-1, todo[,1], todo[,1])
		todoY <- cbind(todo[,2]-1, todo[,2], todo[,2], todo[,2]-1)
		polygon(  cuts(t(todoX))
				, cuts(t(todoY))
				, col = col[i]
				, border = NA
				)
	}

	par(family = "CharisSIL")
	legend(x = dim(m)[1]
		 , y = dim(m)[2]
		 , legend = c(levs, "other", "NA")
		 , xpd = TRUE
		 , pch = c(rep(15, times =  length(levs)), 0, 20)
		 , col = c(col, "white", "grey")
		 , bty = "n"
		 , cex = cex.legend
		 , ncol = 1
		 )
	par(family = "")
	
}