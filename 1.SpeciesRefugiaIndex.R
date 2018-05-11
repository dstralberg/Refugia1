library(raster) #raster manipulation
library(rgdal) #spatial data library
library(yaImpute) # for ann function
library(gdata) # for nobs function

#ref = output directory
#species = SDM directory
#ecor = study area extent raster
#sample = sample SDM raster

#Set working directory and obtain list of species with SDMs
setwd(species)
spec <- list.dirs(species)
sample <- sample*0 + 1

#Calculate probability from fat-tailed negative exponential distribution for a given alpha value and distance.
fattail <- function(x, alpha, c) {
	left <- (c/(2*alpha*gamma(1/c)))
	right <- exp(-1*((abs(x/alpha))^c))
	result <- left*right
	return(right)
}

#Calculate mean dispersal distance for a given fat-tailed distribution
ftmean <- function(alpha, c) {
		result <-(alpha*gamma(2/c))/gamma(1/c)
		return(result)
		}

#Calculate mean prevelance and generate thresholded binary distribution raster layers
prevtab <- data.frame("spec" = gsub(paste(species,"/",sep=""),"",spec), "prev" = 0)
for (i in 2:length(spec)){
  velstack<-emptystack
  p <- list.files(spec[i],pattern="current.tif$")
  pres <- raster(paste(spec[i],"\\",p,sep=""))
  presc <- crop(pres,ecor)
  prescm <- mask(presc,ecor)
  prev <- cellStats(prescm,'mean')
  prevtab[i,2] <- prev
  present<-as.data.frame(rasterToPoints(prescm))
  names(present)[3] <- "prob"
  presprev <- prescm
  values(presprev) <- ifelse(values(presprev) > prev, 1, 0) 
  writeRaster(presprev,filename=paste(spec[i],"presprev",sep=""),format="GTiff", overwrite=TRUE)
  f <- list.files(spec[i],pattern="future.tif$")
  for (j in 1:4) {
    fut <- raster(paste(spec[i],"/",f[j],sep=""))
    futc <- crop(fut,ecor)
    futcm <- mask(futc,ecor)
    future<-as.data.frame(rasterToPoints(futcm))
    names(future)[3] <- "prob"	
    futprev <- futcm
    values(futprev) <- ifelse(values(futprev) > prev, 1, 0) 
    writeRaster(futprev,filename=paste(gsub(paste(species,"/",sep=""),"",f[j]),"futprev",sep=""),format="GTiff", overwrite=TRUE)
  }
}

#Create empty stack to hold species-specific refugia layers
emptystack<-sample
emptystack<-addLayer(emptystack,emptystack)
emptystack<-dropLayer(emptystack,c(1,2))

#Calculate species-specific refugia layers based on backward climate velocity
for (i in 2:length(spec)){
  refstack<-emptystack
  futprevstack<-emptystack
  p <- list.files(spec[i],pattern="presprev.tif$")
  presprev <- raster(paste(spec[i],"//",p,sep=""))
  present<-as.data.frame(rasterToPoints(presprev))
  names(present)[3] <- "prev"
  f <- list.files(spec[i],pattern="futprev.tif")
  for (j in 1:4) {
    futprev <- raster(paste(spec[i],"\\",f[j],sep=""))
    future <- as.data.frame(rasterToPoints(futprev))
    names(future)[3] <- "prev"	
    p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$prev)
    f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$prev)
    p.xy2<-p.xy[p.xy[,4]>0.1,1:3,drop=FALSE]
    f.xy2<-f.xy[f.xy[,4]>0.1,1:3,drop=FALSE]
    if(nrow(f.xy)>0){
      d.ann <- as.data.frame(ann(
        as.matrix(p.xy2[,-1,drop=FALSE]),
        as.matrix(f.xy2[,-1,drop=FALSE]),
        k=1, verbose=F)$knnIndexDist)
      d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
      names(d1b) <- c("ID","X","Y","bvel")
    } else {
      print(spec[i])
    }
    f.xy <- as.data.frame(f.xy)
    names(f.xy) <- c("ID","X","Y","Pres")
    d1b<-merge(f.xy,d1b,by=c("ID","X","Y"),all.x=TRUE)
    d1b$fat <- fattail(d1b$bvel, 8333.3335, 0.5) # 8333.335 = alpha value resulting in mean of 50 km / century or 500 m/year
    sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
    sppref[is.na(sppref[])] <- 0
    refstack<-addLayer(refstack,sppref)
    futprevstack<-addLayer(futprevstack,futprev)
  }
  futprevmean <- overlay(futprevstack, fun='mean')
  writeRaster(futprevmean,filename=paste(gsub(paste(species,"/",sep=""),"",spec[i]),"futprev",sep=""),format="GTiff",overwrite=TRUE)
  refmean <- overlay(refstack, fun='mean')
  writeRaster(refmean,filename=paste(ref,gsub(paste(species,"/",sep=""),"",spec[i]),"refugia",sep=""),format="GTiff",overwrite=TRUE)
}


	


