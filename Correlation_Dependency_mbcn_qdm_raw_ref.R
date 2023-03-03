rm(list=ls())
library(raster)
library(ppcor)
library(rasterVis)
library(sp)
library(rnaturalearth)
library("RColorBrewer")
library(spatialEco)
library(rasterVis)
library(latticeExtra)
library(colorRamps)
library("RColorBrewer")
require(colorRamps)
cols=colorRampPalette(matlab.like2(255))


wp<-ne_countries(continent = c('africa','asia','europe','oceania','north america','south america'), returnclass = "sp") #plots some continents
crs(wp)='+proj=longlat +datum=WGS84'

basin<-shapefile("V:/water_balance/Basins/Major_Basins_of_the_World.shp")
af<-shapefile("V:/water_balance/Basins/hybas_af_lev01-12_v1c/hybas_af_lev02_v1c.shp")
af$HYBAS_ID[[8]] <- "Lake \nChad"
b1<-subset(basin,basin$NAME==c("Amazon"));b2<-subset(basin,basin$NAME==c("Parana"))
b3<-subset(basin,basin$NAME==c("Mississippi"));b4<-subset(basin,basin$NAME==c("Mackenzi"))
b5<-subset(basin,basin$NAME==c("Nile"));b6<-subset(basin,basin$NAME==c("Volta"))
b7<-subset(basin,basin$NAME==c("Orange"));b8<-subset(basin,basin$NAME==c("Danube"))
#b9<-subset(basin,basin$NAME==c("Dnieper"));
b10<-subset(basin,basin$NAME==c("Murray-Darling"))

b11<-subset(basin,basin$NAME==c("Yangtze"));b12<-subset(basin,basin$NAME==c("Lena"))
b13<-subset(basin,basin$NAME==c("Kolyma"));b14<-subset(basin,basin$NAME==c("Volga"))
b15<-subset(basin,basin$NAME==c("Ob"))
#b17<-subset(basin,basin$NAME==c("Hwang Ho"))
b18<-subset(basin,basin$NAME==c("Zambezi"))
b19<-subset(basin,basin$NAME==c("Amur"))#;b20<-subset(basin,basin$NAME==c("Colorado"))
b21<-subset(basin,basin$NAME==c("Yukon"))#;b22<-subset(basin,basin$NAME==c("Zaire"))

basin<-rbind(b13,b2,b3,b4,b5,b6,b7,b8,b10,b11,b12,b14,b15,b18,b19,b21,b1)
k2<-data.frame(basin$NAME[c(-16,-18)]) ##removing repeated Amazon label

k1<-data.frame(1:nrow(k2),k2)
names(k1)<-c("ID","Names")
basin1<-basin

basin1$NAME[c(-16,-18)] <- k1$ID
b16<-subset(af,af$HYBAS_ID==c("Lake \nChad")) ##eatracting lake chad from anothe r shapefile
b161<-b16
b161$HYBAS_ID <- nrow(k2)+1
wp12<-ne_coastline(returnclass = "sp") 
crs(wp12)='+proj=longlat +datum=WGS84'
wp<-shapefile("H:/mask_files_shp_file/World_Continents/World_Continents.shp")
crs(wp)='+proj=longlat +datum=WGS84'




#########################mbcn
#################################################### world monthly


ENS_H_qdm_wp_pr_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/qdm/ENS_HIST.nc',varname='pr'),wp)
ENS_H_qdm_wp_mrro_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/qdm/ENS_HIST.nc',varname='mrro'),wp)    
ENS_H_qdm_wp_et_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/qdm/ENS_HIST.nc',varname='et'),wp)
ENS_H_qdm_wp_WC_monthly<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/qdm/ENS_WB_HIST.nc'),wp)

ENS_H_raw_wp_pr_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/raw/ENS_HIST.nc',varname='pr'),wp)
ENS_H_raw_wp_mrro_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/raw/ENS_HIST.nc',varname='mrro'),wp)    
ENS_H_raw_wp_et_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/raw/ENS_HIST.nc',varname='et'),wp)
ENS_H_raw_wp_WC_monthly<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/raw/ENS_WB_HIST.nc'),wp)

ENS_H_mbcn_wp_pr_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/ENS_HIST.nc',varname='pr'),wp)
ENS_H_mbcn_wp_mrro_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/ENS_HIST.nc',varname='mrro'),wp)    
ENS_H_mbcn_wp_et_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/ENS_HIST.nc',varname='et'),wp)
ENS_H_mbcn_wp_WC_monthly<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_HIST.nc'),wp)


Ref_H_mbcn_wp_pr_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/Ref_HIST.nc',varname='pr'),wp)
Ref_H_mbcn_wp_mrro_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/Ref_HIST.nc',varname='mrro'),wp)    
Ref_H_mbcn_wp_et_monthly<-mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/for_BC/mbcn/Ref_HIST.nc',varname='et'),wp)
Ref_H_mbcn_wp_WC_monthly<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST.nc'),wp)


###############correlation with thermal stress


gridcorts <- function(rasterstack, method, type=c("corel","pval","both")){
  # Values for (layers, ncell, ncol, nrow, method, crs, extent) come straight from the input raster stack
  # e.g. nlayers(rasterstack), ncell(rasterstack)... etc.
  print(paste("Start Gridcorts:",Sys.time()))
  print("Loading parameters")
  layers=nlayers(rasterstack);ncell=ncell(rasterstack);
  ncol=ncol(rasterstack);nrow=nrow(rasterstack);crs=crs(rasterstack);
  extent=extent(rasterstack);pb = txtProgressBar(min = 0, max = ncell, initial = 0)
  print("Done loading parameters")
  mtrx <- as.matrix(rasterstack,ncol=layers)
  empt <- matrix(nrow=ncell, ncol=2)
  print("Initiating loop operation")
  if (type == "corel"){
    for (i in 1:ncell){
      setTxtProgressBar(pb,i)
      if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
        empt[i,1] <- NA 
      } else 
        if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
          empt[i,1] <- NA 
        } else 
          empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate)
    }
    print("Creating empty raster")
    corel <- raster(nrows=nrow,ncols=ncol,crs=crs)
    extent(corel) <- extent
    print("Populating correlation raster")
    values(corel) <- empt[,1]
    print(paste("Ending Gridcorts on",Sys.time()))
    corel
  } 
  else
    if (type == "pval"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,2] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,2] <- NA 
          } else 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
      }
      pval <- raster(nrows=nrow,ncols=ncol,crs=crs)
      extent(pval) <- extent
      print("Populating significance raster")
      values(pval) <- empt[,2]
      print(paste("Ending Gridcorts on",Sys.time()))
      pval
    }
  else
    if (type == "both"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,] <- NA 
          } else {
            empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate) 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
          }
      }
      c <- raster(nrows=nrow,ncols=ncol,crs=crs)
      p <- raster(nrows=nrow,ncols=ncol,crs=crs)
      print("Populating raster brick")
      values(c) <- empt[,1]
      values(p) <- empt[,2]
      brk <- brick(c,p)
      extent(brk) <- extent
      names(brk) <- c("Correlation","Pvalue")
      print(paste("Ending Gridcorts on",Sys.time()))
      brk
    }
}

###########correlation

ENS_H_qdm_WC_pr_monthly<-gridcorts(stack(ENS_H_qdm_wp_WC_monthly,ENS_H_qdm_wp_pr_monthly),method = "pearson", type="both")
ENS_H_qdm_WC_mrro_monthly<-gridcorts(stack(ENS_H_qdm_wp_WC_monthly,ENS_H_qdm_wp_mrro_monthly),method = "pearson", type="both")    
ENS_H_qdm_WC_et_monthly<-gridcorts(stack(ENS_H_qdm_wp_WC_monthly,ENS_H_qdm_wp_et_monthly),method = "pearson", type="both")

ENS_H_mbcn_WC_pr_monthly<-gridcorts(stack(ENS_H_mbcn_wp_WC_monthly,ENS_H_mbcn_wp_pr_monthly),method = "pearson", type="both")
ENS_H_mbcn_WC_mrro_monthly<-gridcorts(stack(ENS_H_mbcn_wp_WC_monthly,ENS_H_mbcn_wp_mrro_monthly),method = "pearson", type="both")    
ENS_H_mbcn_WC_et_monthly<-gridcorts(stack(ENS_H_mbcn_wp_WC_monthly,ENS_H_mbcn_wp_et_monthly),method = "pearson", type="both")

ENS_H_raw_WC_pr_monthly<-gridcorts(stack(ENS_H_raw_wp_WC_monthly,ENS_H_raw_wp_pr_monthly),method = "pearson", type="both")
ENS_H_raw_WC_mrro_monthly<-gridcorts(stack(ENS_H_raw_wp_WC_monthly,ENS_H_raw_wp_mrro_monthly),method = "pearson", type="both")    
ENS_H_raw_WC_et_monthly<-gridcorts(stack(ENS_H_raw_wp_WC_monthly,ENS_H_raw_wp_et_monthly),method = "pearson", type="both")

Ref_H_mbcn_WC_pr_monthly<-gridcorts(stack(Ref_H_mbcn_wp_WC_monthly,Ref_H_mbcn_wp_pr_monthly),method = "pearson", type="both")
Ref_H_mbcn_WC_mrro_monthly<-gridcorts(stack(Ref_H_mbcn_wp_WC_monthly,Ref_H_mbcn_wp_mrro_monthly),method = "pearson", type="both")    
Ref_H_mbcn_WC_et_monthly<-gridcorts(stack(Ref_H_mbcn_wp_WC_monthly,Ref_H_mbcn_wp_et_monthly),method = "pearson", type="both")
############significant
ENS_H_qdm_WC_pr_monthly$Pvalue[ENS_H_qdm_WC_pr_monthly$Pvalue >= 0.01] <- NA
ENS_H_qdm_WC_pr_monthly.mask <- rasterToPoints(ENS_H_qdm_WC_pr_monthly$Pvalue, spatial=TRUE)##masking out significant trends as dots
ENS_H_qdm_WC_mrro_monthly$Pvalue[ENS_H_qdm_WC_mrro_monthly$Pvalue >= 0.01] <- NA
ENS_H_qdm_WC_mrro_monthly.mask <- rasterToPoints(ENS_H_qdm_WC_mrro_monthly$Pvalue, spatial=TRUE)
ENS_H_qdm_WC_et_monthly$Pvalue[ENS_H_qdm_WC_et_monthly$Pvalue >= 0.01] <- NA
ENS_H_qdm_WC_et_monthly.mask <- rasterToPoints(ENS_H_qdm_WC_et_monthly$Pvalue, spatial=TRUE)

ENS_H_mbcn_WC_pr_monthly$Pvalue[ENS_H_mbcn_WC_pr_monthly$Pvalue >= 0.01] <- NA
ENS_H_mbcn_WC_pr_monthly.mask <- rasterToPoints(ENS_H_mbcn_WC_pr_monthly$Pvalue, spatial=TRUE)
ENS_H_mbcn_WC_mrro_monthly$Pvalue[ENS_H_mbcn_WC_mrro_monthly$Pvalue >= 0.01] <- NA
ENS_H_mbcn_WC_mrro_monthly.mask <- rasterToPoints(ENS_H_mbcn_WC_mrro_monthly$Pvalue, spatial=TRUE)
ENS_H_mbcn_WC_et_monthly$Pvalue[ENS_H_mbcn_WC_et_monthly$Pvalue >= 0.01] <- NA
ENS_H_mbcn_WC_et_monthly.mask <- rasterToPoints(ENS_H_mbcn_WC_et_monthly$Pvalue, spatial=TRUE)

ENS_H_raw_WC_pr_monthly$Pvalue[ENS_H_raw_WC_pr_monthly$Pvalue >= 0.01] <- NA
ENS_H_raw_WC_pr_monthly.mask <- rasterToPoints(ENS_H_raw_WC_pr_monthly$Pvalue, spatial=TRUE)
ENS_H_raw_WC_mrro_monthly$Pvalue[ENS_H_raw_WC_mrro_monthly$Pvalue >= 0.01] <- NA
ENS_H_raw_WC_mrro_monthly.mask <- rasterToPoints(ENS_H_raw_WC_mrro_monthly$Pvalue, spatial=TRUE)
ENS_H_raw_WC_et_monthly$Pvalue[ENS_H_raw_WC_et_monthly$Pvalue >= 0.01] <- NA
ENS_H_raw_WC_et_monthly.mask <- rasterToPoints(ENS_H_raw_WC_et_monthly$Pvalue, spatial=TRUE)

Ref_H_mbcn_WC_pr_monthly$Pvalue[Ref_H_mbcn_WC_pr_monthly$Pvalue >= 0.01] <- NA
Ref_H_mbcn_WC_pr_monthly.mask <- rasterToPoints(Ref_H_mbcn_WC_pr_monthly$Pvalue, spatial=TRUE)
Ref_H_mbcn_WC_mrro_monthly$Pvalue[Ref_H_mbcn_WC_mrro_monthly$Pvalue >= 0.01] <- NA
Ref_H_mbcn_WC_mrro_monthly.mask <- rasterToPoints(Ref_H_mbcn_WC_mrro_monthly$Pvalue, spatial=TRUE)
Ref_H_mbcn_WC_et_monthly$Pvalue[Ref_H_mbcn_WC_et_monthly$Pvalue >= 0.01] <- NA
Ref_H_mbcn_WC_et_monthly.mask <- rasterToPoints(Ref_H_mbcn_WC_et_monthly$Pvalue, spatial=TRUE)




brk <- seq(min(c(ENS_H_qdm_WC_pr_monthly@data@min,ENS_H_qdm_WC_mrro_monthly@data@min,ENS_H_qdm_WC_et_monthly@data@min,ENS_H_mbcn_WC_pr_monthly@data@min,ENS_H_mbcn_WC_mrro_monthly@data@min,ENS_H_mbcn_WC_et_monthly@data@min,
                 ENS_H_raw_WC_pr_monthly@data@min,ENS_H_raw_WC_mrro_monthly@data@min)),
           max(c(ENS_H_qdm_WC_pr_monthly@data@max,ENS_H_qdm_WC_mrro_monthly@data@max,ENS_H_qdm_WC_et_monthly@data@max,ENS_H_mbcn_WC_pr_monthly@data@max,ENS_H_mbcn_WC_mrro_monthly@data@max,ENS_H_mbcn_WC_et_monthly@data@max,
                 ENS_H_raw_WC_pr_monthly@data@max,ENS_H_raw_WC_mrro_monthly@data@max)),length.out=100)# where the colors change

lab <- round(seq(min(c(ENS_H_qdm_WC_pr_monthly@data@min,ENS_H_qdm_WC_mrro_monthly@data@min,ENS_H_qdm_WC_et_monthly@data@min,ENS_H_mbcn_WC_pr_monthly@data@min,ENS_H_mbcn_WC_mrro_monthly@data@min,ENS_H_mbcn_WC_et_monthly@data@min,
                       ENS_H_raw_WC_pr_monthly@data@min,ENS_H_raw_WC_mrro_monthly@data@min)),
                 max(c(ENS_H_qdm_WC_pr_monthly@data@max,ENS_H_qdm_WC_mrro_monthly@data@max,ENS_H_qdm_WC_et_monthly@data@max,ENS_H_mbcn_WC_pr_monthly@data@max,ENS_H_mbcn_WC_mrro_monthly@data@max,ENS_H_mbcn_WC_et_monthly@data@max,
                       ENS_H_raw_WC_pr_monthly@data@max,ENS_H_raw_WC_mrro_monthly@data@max)),by=0.1),2)


myColorkey<- list(at=brk,space="right", title = list("Correlation ",cex = 0.8,fontface = 7,col = 'black'),
                  title.gpar = list(" ",cex = 0.001,font = 0.5, col = 'white'),#this is stubborn, so i reduced the font because i dont want to see it on the map
                  
                  ## where the colors change
                  labels=list(at=lab,rot=0,font=6,fontface=7,cex=0.8), #fontface="bold"## where to print labels
                  axis.line=list(col='black'),#legend outline (including its ticks)
                  height=1,width=1.6, #height and width of legend
                  tri.lower = TRUE, tri.upper = TRUE) #adds triangular arrows at both the top and end of legend




A<-levelplot(ENS_H_qdm_WC_pr_monthly$Correlation,colorkey=F,#layout=c(3,3), # create a 4x4 layout for the data
             col.regions=cols, # add a color ramp
             par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1,ylab.key.padding=10.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10,xlab.key.padding=0.1,ylab.key.padding=10.1),#adjust the map positioning#parSettings, # suppress axes and legend outline
             scales=list(x=list(draw=T), y=list(draw=FALSE),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"),
             ylim=c(-48, 76),
             #names.attr=list('(a) Reference'), #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(sp.points(ENS_H_qdm_WC_pr_monthly.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+  #when combinining the plots, the titles are missing. This adds text to lattice plot
  layer(ltext(x=-120,y=-51.8,labels = '(a) QDM',cex=0.9,col="black",fontface = 7)) 
update(A, aspect=0.5) #

B<-levelplot(ENS_H_qdm_WC_mrro_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_qdm_WC_mrro_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))
update(B, aspect=0.5)


C<-levelplot(ENS_H_qdm_WC_et_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_qdm_WC_et_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))
update(C, aspect=0.5)


D<-levelplot(ENS_H_mbcn_WC_pr_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_mbcn_WC_pr_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+
  layer(ltext(x=-120,y=-51.8,labels = '(b) MBCN',cex=0.9,col="black",fontface = 7))
update(D, aspect=0.5)

E<-levelplot(ENS_H_mbcn_WC_mrro_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_mbcn_WC_mrro_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))
update(E, aspect=0.5)



F1<-levelplot(ENS_H_mbcn_WC_et_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_mbcn_WC_et_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=75,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))
update(F1, aspect=0.5)



G<-levelplot(ENS_H_raw_WC_pr_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_raw_WC_pr_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+
  layer(ltext(x=-130,y=-51.8,labels = '(c) Raw CMIP6',cex=0.9,col="black",fontface = 7))
update(G, aspect=0.5)


H<-levelplot(ENS_H_raw_WC_mrro_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_raw_WC_mrro_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))
update(H, aspect=0.5)


I<-levelplot(ENS_H_raw_WC_et_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(ENS_H_raw_WC_et_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=75,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))
update(I, aspect=0.5)

J<-levelplot(Ref_H_mbcn_WC_pr_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(Ref_H_mbcn_WC_pr_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+
  layer(ltext(x=-130,y=-51.8,labels = '(d) Reference',cex=0.9,col="black",fontface = 7))
update(J, aspect=0.5)


K<-levelplot(Ref_H_mbcn_WC_mrro_monthly$Correlation,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(Ref_H_mbcn_WC_mrro_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))
update(K, aspect=0.5)


L<-levelplot(Ref_H_mbcn_WC_et_monthly$Correlation,col.regions=cols,colorkey=myColorkey,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             #names.attr=rasterNames1, #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+layer(sp.polygons(wp, lwd=0.5))+
  layer(sp.points(Ref_H_mbcn_WC_et_monthly.mask, pch=1, cex=0.1, alpha=0.8, col="black")) +
  layer(ltext(x=75,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))
update(L, aspect=0.5)


cObj1<- c(A,B,C,D,E,F1,G,H,I,J,K,L,  layout = c(3,4), merge.legends = T,x.same = F, y.same = F)

update(cObj1,fontface = 7,scales = list(draw=T,y = list(rot = 90)),xlab=list("Longitude",cex=1,fontface = 7),ylab=list("Latitude",cex=1,fontface = 7), 
       main=list('',fontface = 7))





###########################RMSE for Deprendency structure

#####RMSE
ENS_H_mbcn_WC_pr_RMSE<- sqrt(mean(ENS_H_mbcn_WC_pr_monthly$Correlation - Ref_H_mbcn_WC_pr_monthly$Correlation)^2)
ENS_H_mbcn_WC_mrro_RMSE<- sqrt(mean(ENS_H_mbcn_WC_mrro_monthly$Correlation - Ref_H_mbcn_WC_mrro_monthly$Correlation)^2)
ENS_H_mbcn_WC_et_RMSE<- sqrt(mean(ENS_H_mbcn_WC_et_monthly$Correlation - Ref_H_mbcn_WC_et_monthly$Correlation)^2)

ENS_H_qdm_WC_pr_RMSE<- sqrt(mean(ENS_H_qdm_WC_pr_monthly$Correlation - Ref_H_mbcn_WC_pr_monthly$Correlation)^2)
ENS_H_qdm_WC_mrro_RMSE<- sqrt(mean(ENS_H_qdm_WC_mrro_monthly$Correlation - Ref_H_mbcn_WC_mrro_monthly$Correlation)^2)
ENS_H_qdm_WC_et_RMSE<- sqrt(mean(ENS_H_qdm_WC_et_monthly$Correlation - Ref_H_mbcn_WC_et_monthly$Correlation)^2)

ENS_H_raw_WC_pr_RMSE<- sqrt(mean(ENS_H_raw_WC_pr_monthly$Correlation - Ref_H_mbcn_WC_pr_monthly$Correlation)^2)
ENS_H_raw_WC_mrro_RMSE<- sqrt(mean(ENS_H_raw_WC_mrro_monthly$Correlation - Ref_H_mbcn_WC_mrro_monthly$Correlation)^2)
ENS_H_raw_WC_et_RMSE<- sqrt(mean(ENS_H_raw_WC_et_monthly$Correlation - Ref_H_mbcn_WC_et_monthly$Correlation)^2)




brk2 <- seq(min(c(ENS_H_mbcn_WC_pr_RMSE@data@min,ENS_H_mbcn_WC_mrro_RMSE@data@min,ENS_H_mbcn_WC_et_RMSE@data@min,
                  ENS_H_qdm_WC_pr_RMSE@data@min,ENS_H_qdm_WC_mrro_RMSE@data@min,ENS_H_qdm_WC_et_RMSE@data@min,
                  ENS_H_raw_WC_pr_RMSE@data@min,ENS_H_raw_WC_mrro_RMSE@data@min,ENS_H_raw_WC_et_RMSE@data@min)), 
            max(c(ENS_H_mbcn_WC_pr_RMSE@data@max,ENS_H_mbcn_WC_mrro_RMSE@data@max,ENS_H_mbcn_WC_et_RMSE@data@max,
                  ENS_H_qdm_WC_pr_RMSE@data@max,ENS_H_qdm_WC_mrro_RMSE@data@max,ENS_H_qdm_WC_et_RMSE@data@max,
                  ENS_H_raw_WC_pr_RMSE@data@max,ENS_H_raw_WC_mrro_RMSE@data@max,ENS_H_raw_WC_et_RMSE@data@max)), by=0.05) # where the colors change

lab2<-round(seq(min(c(ENS_H_mbcn_WC_pr_RMSE@data@min,ENS_H_mbcn_WC_mrro_RMSE@data@min,ENS_H_mbcn_WC_et_RMSE@data@min,
                      ENS_H_qdm_WC_pr_RMSE@data@min,ENS_H_qdm_WC_mrro_RMSE@data@min,ENS_H_qdm_WC_et_RMSE@data@min,
                      ENS_H_raw_WC_pr_RMSE@data@min,ENS_H_raw_WC_mrro_RMSE@data@min,ENS_H_raw_WC_et_RMSE@data@min)), 
                max(c(ENS_H_mbcn_WC_pr_RMSE@data@max,ENS_H_mbcn_WC_mrro_RMSE@data@max,ENS_H_mbcn_WC_et_RMSE@data@max,
                      ENS_H_qdm_WC_pr_RMSE@data@max,ENS_H_qdm_WC_mrro_RMSE@data@max,ENS_H_qdm_WC_et_RMSE@data@max,
                      ENS_H_raw_WC_pr_RMSE@data@max,ENS_H_raw_WC_mrro_RMSE@data@max,ENS_H_raw_WC_et_RMSE@data@max)),by=0.1),2) # where to print labels


myColorkey2 <- list(at=brk2,space="right", title = list("Joint dependency",cex = 0.8,fontface = 7, col = 'black'),
                    title.gpar = list(" ",cex = 0.001,font = 0.5, col = 'white'),#this is stubborn, so i reduced the font because i dont want to see it on the map
                    ## where the colors change
                    labels=list(at=lab2,rot=0,font=6,fontface=7,cex=0.8), #fontface="bold"## where to print labels
                    axis.line=list(col='black'),#legend outline (including its ticks)
                    height=1,width=1.6, #height and width of legend
                    tri.lower = TRUE, tri.upper = TRUE) #adds triangular arrows at both the top and end of legend


A1<-levelplot(ENS_H_qdm_WC_pr_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
             col.regions=cols, # add a color ramp
             par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"),
             ylim=c(-48, 76),
             #names.attr='(b) Historical', #names on top of each plot
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
             ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+  #when combinining the plots, the titles are missing. This adds text to lattice plot
  layer(ltext(x=-120,y=-51.8,labels = '(a) QDM',cex=0.9,col="black",fontface = 7)) 
update(A1, aspect=0.5) 

B1<-levelplot(ENS_H_qdm_WC_mrro_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(B1, aspect=0.5) 

C1<-levelplot(ENS_H_qdm_WC_et_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(C1, aspect=0.5) 


D1<-levelplot(ENS_H_mbcn_WC_pr_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
 # layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+  #when combinining the plots, the titles are missing. This adds text to lattice plot
  layer(ltext(x=-120,y=-51.8,labels = '(b) MBCN',cex=0.9,col="black",fontface = 7)) 
update(D1, aspect=0.5) 

E1<-levelplot(ENS_H_mbcn_WC_mrro_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(E1, aspect=0.5) 

F11<-levelplot(ENS_H_mbcn_WC_et_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(F11, aspect=0.5) 





G1<-levelplot(ENS_H_raw_WC_pr_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(i) Water Conservation - Precipitation',cex=0.9,col="black",fontface = 7))+  #when combinining the plots, the titles are missing. This adds text to lattice plot
  layer(ltext(x=-120,y=-51.8,labels = '(c) Raw CMIP6',cex=0.9,col="black",fontface = 7)) 
update(G1, aspect=0.5) 

H1<-levelplot(ENS_H_raw_WC_mrro_RMSE,at=brk2,colorkey=F,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(ii) Water Conservation - Runoff',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(H1, aspect=0.5) 

I1<-levelplot(ENS_H_raw_WC_et_RMSE,at=brk2,colorkey=myColorkey2,panel = panel.levelplot.raster,#layout=c(3,3), # create a 4x4 layout for the data
              col.regions=cols, # add a color ramp
              par.settings=list(panel.background=list(col="white"),axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                                strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
              layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, # suppress axes and legend outline
              scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
              main=list('',cex=1,col="black"),
              ylim=c(-48, 76),
              #names.attr='(b) Historical', #names on top of each plot
              margin=F, auto.key=F,maxpixels=2e5,xlab=list('',fontface='bold'),
              ylab=list('', rot=90, fontface='bold'))+#layer(sp.polygons(basin, lwd=2,col = "red"))+
  layer(sp.polygons(wp, lwd=1,col = "black"))+#layer(sp.polygons(b16, lwd=2,col = "red"))+
  #layer(sp.text(coordinates(basin1), txt = basin1$NAME,cex=0.8,font=7,col="black"))+layer(sp.text(coordinates(b161), txt = b161$HYBAS_ID,cex=0.8,font=7,col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(iii) Water Conservation - Evapotranspiration',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(I1, aspect=0.5) 




cObj2<- c(A1,B1,C1,D1,E1,F11,G1,H1,I1,  layout = c(3,3), merge.legends = T,x.same = F, y.same = F)

update(cObj2,fontface = 7,scales = list(draw=T,y = list(rot = 90)),xlab=list("Longitude",cex=1,fontface = 7),ylab=list("Latitude",cex=1,fontface = 7), 
       main=list('',fontface = 7))






