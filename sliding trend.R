rm(list=ls())
library(xts)
library(hydroTSM)
library(zoo)

library(raster)
library(rasterVis)
library(sp)
library(rnaturalearth)
library("RColorBrewer")
library(spatialEco)
library(rasterVis)
library(latticeExtra)

require(colorRamps)
cols=colorRampPalette(matlab.like2(255))



AFR<-ne_countries(continent = "africa")
AS<-ne_countries(continent = "asia")
EU<-ne_countries(continent = "europe")
OC<-ne_countries(continent = "oceania")
NA1<-ne_countries(continent = "north america")
SA<-ne_countries(continent = "south america")
wp<-ne_countries(continent = c('africa','asia','europe','oceania','north america','south america'), returnclass = "sp") #plots some continents

crs(SA)=crs(NA1)=crs(OC)=crs(EU)=crs(AS)=crs(AFR)=crs(wp)='+proj=longlat +datum=WGS84'

year_H <- seq(as.Date('1959-01-01'), as.Date('2014-12-31'), by='year')
year_F <- seq(as.Date('2045-01-01'), as.Date('2100-12-31'), by='year')
###############################
######################################

###############################################Annual##############################

Ref_mbcn_wp_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),wp),mean),order.by=year_H)
ENS_370_mbcn_wp_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),wp),mean),order.by=year_F)
ENS_585_mbcn_wp_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),wp),mean),order.by=year_F)

Ref_mbcn_AFR_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),AFR),mean),order.by=year_H)
ENS_370_mbcn_AFR_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),AFR),mean),order.by=year_F)
ENS_585_mbcn_AFR_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),AFR),mean),order.by=year_F)

Ref_mbcn_AS_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),AS),mean),order.by=year_H)
ENS_370_mbcn_AS_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),AS),mean),order.by=year_F)
ENS_585_mbcn_AS_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),AS),mean),order.by=year_F)

Ref_mbcn_EU_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),EU),mean),order.by=year_H)
ENS_370_mbcn_EU_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),EU),mean),order.by=year_F)
ENS_585_mbcn_EU_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),EU),mean),order.by=year_F)

Ref_mbcn_OC_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),OC),mean),order.by=year_H)
ENS_370_mbcn_OC_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),OC),mean),order.by=year_F)
ENS_585_mbcn_OC_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),OC),mean),order.by=year_F)

Ref_mbcn_NA1_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),NA1),mean),order.by=year_H)
ENS_370_mbcn_NA1_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),NA1),mean),order.by=year_F)
ENS_585_mbcn_NA1_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),NA1),mean),order.by=year_F)

Ref_mbcn_SA_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),SA),mean),order.by=year_H)
ENS_370_mbcn_SA_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),SA),mean),order.by=year_F)
ENS_585_mbcn_SA_WC_annual<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),SA),mean),order.by=year_F)



#################################################################
##########################################summer
#########################################

Ref_mbcn_wp_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),wp),mean),order.by=year_H)
ENS_370_mbcn_wp_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),wp),mean),order.by=year_F)
ENS_585_mbcn_wp_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),wp),mean),order.by=year_F)

Ref_mbcn_AFR_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),AFR),mean),order.by=year_H)
ENS_370_mbcn_AFR_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),AFR),mean),order.by=year_F)
ENS_585_mbcn_AFR_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),AFR),mean),order.by=year_F)

Ref_mbcn_AS_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),AS),mean),order.by=year_H)
ENS_370_mbcn_AS_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),AS),mean),order.by=year_F)
ENS_585_mbcn_AS_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),AS),mean),order.by=year_F)

Ref_mbcn_EU_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),EU),mean),order.by=year_H)
ENS_370_mbcn_EU_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),EU),mean),order.by=year_F)
ENS_585_mbcn_EU_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),EU),mean),order.by=year_F)

Ref_mbcn_OC_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),OC),mean),order.by=year_H)
ENS_370_mbcn_OC_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),OC),mean),order.by=year_F)
ENS_585_mbcn_OC_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),OC),mean),order.by=year_F)

Ref_mbcn_NA1_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),NA1),mean),order.by=year_H)
ENS_370_mbcn_NA1_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),NA1),mean),order.by=year_F)
ENS_585_mbcn_NA1_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),NA1),mean),order.by=year_F)

Ref_mbcn_SA_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),SA),mean),order.by=year_H)
ENS_370_mbcn_SA_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),SA),mean),order.by=year_F)
ENS_585_mbcn_SA_WC_summer<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),SA),mean),order.by=year_F)
#################################################################
##########################################Winter
#########################################

Ref_mbcn_wp_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),wp),mean),order.by=year_H)
ENS_370_mbcn_wp_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),wp),mean),order.by=year_F)
ENS_585_mbcn_wp_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),wp),mean),order.by=year_F)

Ref_mbcn_AFR_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),AFR),mean),order.by=year_H)
ENS_370_mbcn_AFR_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),AFR),mean),order.by=year_F)
ENS_585_mbcn_AFR_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),AFR),mean),order.by=year_F)

Ref_mbcn_AS_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),AS),mean),order.by=year_H)
ENS_370_mbcn_AS_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),AS),mean),order.by=year_F)
ENS_585_mbcn_AS_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),AS),mean),order.by=year_F)

Ref_mbcn_EU_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),EU),mean),order.by=year_H)
ENS_370_mbcn_EU_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),EU),mean),order.by=year_F)
ENS_585_mbcn_EU_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),EU),mean),order.by=year_F)

Ref_mbcn_OC_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),OC),mean),order.by=year_H)
ENS_370_mbcn_OC_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),OC),mean),order.by=year_F)
ENS_585_mbcn_OC_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),OC),mean),order.by=year_F)

Ref_mbcn_NA1_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),NA1),mean),order.by=year_H)
ENS_370_mbcn_NA1_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),NA1),mean),order.by=year_F)
ENS_585_mbcn_NA1_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),NA1),mean),order.by=year_F)

Ref_mbcn_SA_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),SA),mean),order.by=year_H)
ENS_370_mbcn_SA_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),SA),mean),order.by=year_F)
ENS_585_mbcn_SA_WC_winter<- xts(cellStats(mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),SA),mean),order.by=year_F)


#####start from here


#data_H_annual<-xts(Ref_mbcn_wp_WC_annual[,2],Ref_mbcn_AFR_WC_annual[,2],Ref_mbcn_AS_WC_annual[,2],Ref_mbcn_EU_WC_annual[,2],
   #                  Ref_mbcn_SA_WC_annual[,2],Ref_mbcn_NA1_WC_annual[,2],Ref_mbcn_OC_WC_annual[,2])

#names(data_H_annual)<-c("Date","World","Africa","Asia","Europe","South_America","North_America","Oceania")


#data_H_annual.Date=as.Date(data_H_annual$Date, "%Y-%m-%d") #setting the format of my date as it was in the csv file


#data_H_annual.xts=xts(data_H_annual[,2:8],order.by=data_H_annual.Date) #treating my csv as a matrix of all the variables in my csv file


###################################################
##complaining about git_hub key
#use Sys.unsetenv("GITHUB_PAT") to remove it

#remotes::install_github("metno/esd", auth_token = NULL)

library(esd)
library(ggplotify)
library(ggplot2)
source("G:/Host/script/water_balance/sliding trend_source_file_WB.R")
#par(mfrow=c(2,2),order.by=year_F)
#source("G:/Host/script/Homogenisation_Jhdrol_R_script/sliding_trend_ppt_source.R")
##CALCULATE SLIDING TRENDS and convert to ggplot
D<-ggplotify::as.ggplot(function()vist(Ref_mbcn_wp_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                             lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white");D

E<-ggplotify::as.ggplot(function()vist(Ref_mbcn_AFR_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                            lwd=2,unitlabel = "mm",new=F, varlabel = "\n(a) Annual - Historical - Trend    \n[ii] Africa"));E


F1<-ggplotify::as.ggplot(function()vist(Ref_mbcn_SA_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                              lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white");F1



A<-as.ggplot(function()vist(ENS_370_mbcn_wp_WC_annual,minlen = 10,pmax=0.1,vmax=20,
                                  lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white");A


#+annotate("text",x=.7, y=.9, label= "man") ##this adds label man on the plot

B<-ggplotify::as.ggplot(function()vist(ENS_370_mbcn_AFR_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                             lwd=2,unitlabel = "mm",new=F, varlabel = "\n(b) Annual - SSP 370 - Trend     \n[ii] Africa"))



C<-ggplotify::as.ggplot(function()vist(ENS_370_mbcn_SA_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                             lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white")




##########################
F12<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_wp_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                               lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white")

F13<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_AFR_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                               lwd=2,unitlabel = "mm",new=F, varlabel =  "\n(c) Annual - SSP 585 - Trend     \n[ii] Africa"))

F14<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_SA_WC_annual,minlen = 10,pmax=0.05,vmax=20,
                                               lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white")

gridExtra::grid.arrange(D,E,F1,A,B,C,F12,F13,F14,ncol=3)
###########################################


###############################summer
D<-ggplotify::as.ggplot(function()vist(Ref_mbcn_wp_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                       lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white");D

E<-ggplotify::as.ggplot(function()vist(Ref_mbcn_AFR_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                       lwd=2,unitlabel = "mm",new=F, varlabel = "\n(a) Winter - Historical - Trend    \n[ii] Africa"));E


F1<-ggplotify::as.ggplot(function()vist(Ref_mbcn_SA_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                        lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white");F1



A<-as.ggplot(function()vist(ENS_370_mbcn_wp_WC_summer,minlen = 10,pmax=0.1,vmax=20,
                            lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white");A


#+annotate("text",x=.7, y=.9, label= "man") ##this adds label man on the plot

B<-ggplotify::as.ggplot(function()vist(ENS_370_mbcn_AFR_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                       lwd=2,unitlabel = "mm",new=F, varlabel = "\n(b) Winter - SSP 370 - Trend     \n[ii] Africa"))



C<-ggplotify::as.ggplot(function()vist(ENS_370_mbcn_SA_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                       lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white")




##########################
F12<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_wp_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                         lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[i] World"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), #this plots a box on the plot to cover the unwanted legend
           ymax=c(.88,96),alpha=1, color="white",fill="white")

F13<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_AFR_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                         lwd=2,unitlabel = "mm",new=F, varlabel =  "\n(c) Winter - SSP 585 - Trend     \n[ii] Africa"))

F14<-ggplotify::as.ggplot(function()vist(ENS_585_mbcn_SA_WC_summer,minlen = 10,pmax=0.05,vmax=20,
                                         lwd=2,unitlabel = "mm",new=F, varlabel = "\n      \n[iii] South America"))+
  annotate("rect",xmin = c(.7,.9),xmax = c(.91,0.95),ymin=c(.5,.7), 
           ymax=c(.88,96),alpha=1, color="white",fill="white")

gridExtra::grid.arrange(D,E,F1,A,B,C,F12,F13,F14,ncol=3)






#####################Summer

