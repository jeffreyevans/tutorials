#### Example of using focal on two rasters to evaluate the 
####   mutual information content and correlation

library(raster)
  setwd("C:/evans/TMP/mvfocal")
  r <- stack("lai.year.tif")
    r <- r[[1:2]]

w <- c(5,5)
out <- writeStart(r[[1]], "test.tif", overwrite=TRUE)
  for( rl in 1:nrow(r) ) { 
    v <- getValuesFocal(r[[1:2]], row=rl, nrows=1, ngb = w, array = FALSE)
      mi <- rep(NA,nrow(v[[1]]))
      mi.var <- rep(NA,nrow(v[[1]]))
        for(i in 1:nrow(v[[1]]) ) {
		  xy <- na.omit( data.frame(x=v[[1]][i,],y=v[[2]][i,]) )	
            if( nrow(xy) > 3 ) {
              x <- as.vector(xy$x)
			  y <- as.vector(xy$y) 
		    } else {  
		      x <- NA
              y <- NA
            }				    
      n <- length(x) - 1
    if(n > 1) {
         #mi[i] <- stats::median(FNN::mutinfo(x, y, k = n, direct=FALSE))
		 mi[i] <- stats::cor(x,y)
         #mi.var[i] <- stats::var(FNN::mutinfo(x, y, k = n, direct=FALSE))
    } else {
         mi[i] <- NA
         #mi.var[i] <- NA
       }
     }
    out <- writeValues(out, mi, start=rl)
  }
writeStop(out)



#### Example of using focal on multiple rasters to evaluate the 
####   dissimilarity using MDS 
library(raster)
  setwd("C:/evans/TMP/mvfocal")
  r <- stack("lai.year.tif")

w <- c(5,5)
fidx <- median(1:(w[1]*w[2]))
lai.mds <- writeStart(r[[1]], "lai_mds.tif", overwrite=TRUE)
  for( rl in 1:nrow(r) ) {
    n <- w[1] * w[2]  
    v <- getValuesFocal(r, row=rl, nrows=1, ngb = w, array = FALSE)
    diss <- rep(NA,nrow(v[[1]]))
        for(i in 1:nrow(v[[1]]) ) {
		  xy <- data.frame(id=1:length(v[[1]][1,]))
		    for(j in 1:length(v)) { xy[,j+1] <- v[[j]][i,] }
			  xy <- xy[,2:ncol(xy)]
			  na.idx <- unique(which(is.na(xy), arr.ind=TRUE)[,1])
			  na.idx.col <- unique(which(is.na(xy), arr.ind=TRUE)[,2])	  
		    if(length(na.idx.col) == ncol(xy)) { xy <- xy[,2:ncol(xy)][-na.idx,] }
          if( nrow(xy) == n ) {
            diss[i] <- cmdscale(dist(xy), k=1)[fidx]
		  } else if( nrow(xy) > 4) {
		    diss[i] <- median(cmdscale(dist(xy), k=1))
		  } else {  
		    diss[i] <- NA
          }				    
        }
    lai.mds <- writeValues(lai.mds, diss, start=rl)
  }
writeStop(lai.mds)

if(fidx %in% na.idx) xy[fidx,] <- rep(0,ncol(xy))



######################################################
# Regression example
library(raster)
library(gstat)                                         
library(sp)

setwd("C:/evans/TMP")
out.raster <- "lm.tif"

data(meuse)                                            
data(meuse.grid)                                       
coordinates(meuse) <- ~x + y                           
coordinates(meuse.grid) <- ~x + y                      

v1 <- variogram(log(copper) ~ 1, meuse)                  
x1 <- fit.variogram(v1, vgm(1, "Sph", 800, 1))           
copper <- krige(copper ~ 1, meuse, meuse.grid, x1, nmax = 30)
gridded(copper) <- TRUE                                      
copper@data = as.data.frame(copper@data[,-2])

v2 <- variogram(log(elev) ~ 1, meuse)                  
x2 <- fit.variogram(v2, vgm(.1, "Sph", 1000, .6))        
elev <- krige(elev ~ 1, meuse, meuse.grid, x2, nmax = 30)
gridded(elev) <- TRUE    
elev@data <- as.data.frame(elev@data[,-2])
elev@data[,1] <- elev@data[,1]

r <- stack( raster(copper), raster(elev) )

w <- c(3,3)
lapse <- r[[1]]
  lapse[] <- NA
  for( rl in 1:nrow(r) ) { 
    v <- getValuesFocal(r[[1:2]], row=rl, nrows=1, ngb = w, array = FALSE)
      fit <- rep(NA,nrow(v[[1]]))
        for(i in 1:nrow(v[[1]]) ) {
		  xy <- na.omit( data.frame(x=v[[1]][i,],y=v[[2]][i,]) )	
          if( nrow(xy) > 4 ) {
            fit[i] <- coefficients(lm(as.numeric(xy$y) ~ as.numeric(xy$x)))[2]
	          if(is.null(fit)) fit = 1
          } else {
            fit[i] <- NA 
         }
     }
    lapse[i] <- fit
  }


