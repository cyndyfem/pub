rm(list=ls())
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


wp<-ne_countries(continent = c('africa','asia','europe','oceania','north america','south america'), returnclass = "sp") #plots some continents
crs(wp)='+proj=longlat +datum=WGS84'

###############################
######################################
###############################################Annual##############################

Ref_mbcn_wp_WC_annual<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_annual.nc'),wp)
ENS_H_mbcn_wp_WC_annual<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_HIST_annual.nc'),wp)
ENS_370_mbcn_wp_WC_annual<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_annual.nc'),wp)
ENS_585_mbcn_wp_WC_annual<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_annual.nc'),wp)
ENS_H_raw_wp_WC_annual<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/raw/ENS_WB_HIST_annual.nc'),wp)

#################################################################
##########################################summer
#########################################

Ref_mbcn_wp_WC_summer<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_summer.nc'),wp)
ENS_H_mbcn_wp_WC_summer<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_HIST_summer.nc'),wp)
ENS_370_mbcn_wp_WC_summer<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_summer.nc'),wp)
ENS_585_mbcn_wp_WC_summer<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_summer.nc'),wp)
ENS_H_raw_wp_WC_summer<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/raw/ENS_WB_HIST_summer.nc'),wp)

#################################################################
##########################################Winter
#########################################

Ref_mbcn_wp_WC_winter<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/Ref_WB_HIST_winter.nc'),wp)
ENS_H_mbcn_wp_WC_winter<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_HIST_winter.nc'),wp)
ENS_370_mbcn_wp_WC_winter<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp370_mbcn_winter.nc'),wp)
ENS_585_mbcn_wp_WC_winter<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/mbcn/ENS_WB_ssp585_mbcn_winter.nc'),wp)
ENS_H_raw_wp_WC_winter<- mask(brick('V:/Runoff/for_BC/corrected_BC/runoff_bc/outputs_WB/raw/ENS_WB_HIST_winter.nc'),wp)

#Climate Change Signal
#annual
ENS_CC_370_WC_annual<-ENS_370_mbcn_wp_WC_annual-ENS_H_mbcn_wp_WC_annual
ENS_CC_585_WC_annual<-ENS_585_mbcn_wp_WC_annual-ENS_H_mbcn_wp_WC_annual

#winter
ENS_CC_370_WC_winter<-ENS_370_mbcn_wp_WC_winter-ENS_H_mbcn_wp_WC_winter
ENS_CC_585_WC_winter<-ENS_585_mbcn_wp_WC_winter-ENS_H_mbcn_wp_WC_winter
#summer
ENS_CC_370_WC_summer<-ENS_370_mbcn_wp_WC_summer-ENS_H_mbcn_wp_WC_summer
ENS_CC_585_WC_summer<-ENS_585_mbcn_wp_WC_summer-ENS_H_mbcn_wp_WC_summer


############Trend
Ref_mbcn_wp_WC_annual_trnd <- raster.kendall(Ref_mbcn_wp_WC_annual,  p.value = TRUE,confidence = TRUE)
Ref_mbcn_wp_WC_annual_trnd$slope <-Ref_mbcn_wp_WC_annual_trnd$slope*10
names(Ref_mbcn_wp_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
Ref_mbcn_wp_WC_annual_trnd$p.value[Ref_mbcn_wp_WC_annual_trnd$p.value >= 0.05] <- NA
Ref_mbcn_wp_WC_annual_trnd.mask <- rasterToPoints(Ref_mbcn_wp_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#Ref_mbcn_wp_WC_annual_trnd$slope <- mask(Ref_mbcn_wp_WC_annual_trnd$slope, Ref_mbcn_wp_WC_annual_trnd$p.value) ###Plot only significant trend



ENS_H_mbcn_wp_WC_annual_trnd <- raster.kendall(ENS_H_mbcn_wp_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_H_mbcn_wp_WC_annual_trnd$slope <-ENS_H_mbcn_wp_WC_annual_trnd$slope*10
names(ENS_H_mbcn_wp_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_mbcn_wp_WC_annual_trnd$p.value[ENS_H_mbcn_wp_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_H_mbcn_wp_WC_annual_trnd.mask <- rasterToPoints(ENS_H_mbcn_wp_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_mbcn_wp_WC_annual_trnd$slope <- mask(ENS_H_mbcn_wp_WC_annual_trnd$slope, ENS_H_mbcn_wp_WC_annual_trnd$p.value) ###Plot only significant trend


ENS_H_raw_wp_WC_annual_trnd <- raster.kendall(ENS_H_raw_wp_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_H_raw_wp_WC_annual_trnd$slope <-ENS_H_raw_wp_WC_annual_trnd$slope*10
names(ENS_H_raw_wp_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_raw_wp_WC_annual_trnd$p.value[ENS_H_raw_wp_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_H_raw_wp_WC_annual_trnd.mask <- rasterToPoints(ENS_H_raw_wp_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_raw_wp_WC_annual_trnd$slope <- mask(ENS_H_raw_wp_WC_annual_trnd$slope, ENS_H_raw_wp_WC_annual_trnd$p.value) ###Plot only significant trend


ENS_370_mbcn_wp_WC_annual_trnd <- raster.kendall(ENS_370_mbcn_wp_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_370_mbcn_wp_WC_annual_trnd$slope <-ENS_370_mbcn_wp_WC_annual_trnd$slope*10
names(ENS_370_mbcn_wp_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_370_mbcn_wp_WC_annual_trnd$p.value[ENS_370_mbcn_wp_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_370_mbcn_wp_WC_annual_trnd.mask <- rasterToPoints(ENS_370_mbcn_wp_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_370_mbcn_wp_WC_annual_trnd$slope <- mask(ENS_370_mbcn_wp_WC_annual_trnd$slope, ENS_370_mbcn_wp_WC_annual_trnd$p.value) ###Plot only significant trend


ENS_585_mbcn_wp_WC_annual_trnd <- raster.kendall(ENS_585_mbcn_wp_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_585_mbcn_wp_WC_annual_trnd$slope <-ENS_585_mbcn_wp_WC_annual_trnd$slope*10
names(ENS_585_mbcn_wp_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_585_mbcn_wp_WC_annual_trnd$p.value[ENS_585_mbcn_wp_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_585_mbcn_wp_WC_annual_trnd.mask <- rasterToPoints(ENS_585_mbcn_wp_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_585_mbcn_wp_WC_annual_trnd$slope <- mask(ENS_585_mbcn_wp_WC_annual_trnd$slope, ENS_585_mbcn_wp_WC_annual_trnd$p.value) ###Plot only significant trend

Ref_mbcn_wp_WC_summer_trnd <- raster.kendall(Ref_mbcn_wp_WC_summer,  p.value = TRUE,confidence = TRUE)
Ref_mbcn_wp_WC_summer_trnd$slope <-Ref_mbcn_wp_WC_summer_trnd$slope*10
names(Ref_mbcn_wp_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
Ref_mbcn_wp_WC_summer_trnd$p.value[Ref_mbcn_wp_WC_summer_trnd$p.value >= 0.05] <- NA
Ref_mbcn_wp_WC_summer_trnd.mask <- rasterToPoints(Ref_mbcn_wp_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#Ref_mbcn_wp_WC_summer_trnd$slope <- mask(Ref_mbcn_wp_WC_summer_trnd$slope, Ref_mbcn_wp_WC_summer_trnd$p.value) ###Plot only significant trend


ENS_H_mbcn_wp_WC_summer_trnd <- raster.kendall(ENS_H_mbcn_wp_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_H_mbcn_wp_WC_summer_trnd$slope <-ENS_H_mbcn_wp_WC_summer_trnd$slope*10
names(ENS_H_mbcn_wp_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_mbcn_wp_WC_summer_trnd$p.value[ENS_H_mbcn_wp_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_H_mbcn_wp_WC_summer_trnd.mask <- rasterToPoints(ENS_H_mbcn_wp_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_mbcn_wp_WC_summer_trnd$slope <- mask(ENS_H_mbcn_wp_WC_summer_trnd$slope, ENS_H_mbcn_wp_WC_summer_trnd$p.value) ###Plot only significant trend

ENS_H_raw_wp_WC_summer_trnd <- raster.kendall(ENS_H_raw_wp_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_H_raw_wp_WC_summer_trnd$slope <-ENS_H_raw_wp_WC_summer_trnd$slope*10
names(ENS_H_raw_wp_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_raw_wp_WC_summer_trnd$p.value[ENS_H_raw_wp_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_H_raw_wp_WC_summer_trnd.mask <- rasterToPoints(ENS_H_raw_wp_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_raw_wp_WC_summer_trnd$slope <- mask(ENS_H_raw_wp_WC_summer_trnd$slope, ENS_H_raw_wp_WC_summer_trnd$p.value) ###Plot only significant trend

ENS_370_mbcn_wp_WC_summer_trnd <- raster.kendall(ENS_370_mbcn_wp_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_370_mbcn_wp_WC_summer_trnd$slope <-ENS_370_mbcn_wp_WC_summer_trnd$slope*10
names(ENS_370_mbcn_wp_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_370_mbcn_wp_WC_summer_trnd$p.value[ENS_370_mbcn_wp_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_370_mbcn_wp_WC_summer_trnd.mask <- rasterToPoints(ENS_370_mbcn_wp_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_370_mbcn_wp_WC_summer_trnd$slope <- mask(ENS_370_mbcn_wp_WC_summer_trnd$slope, ENS_370_mbcn_wp_WC_summer_trnd$p.value) ###Plot only significant trend


ENS_585_mbcn_wp_WC_summer_trnd <- raster.kendall(ENS_585_mbcn_wp_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_585_mbcn_wp_WC_summer_trnd$slope <-ENS_585_mbcn_wp_WC_summer_trnd$slope*10
names(ENS_585_mbcn_wp_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_585_mbcn_wp_WC_summer_trnd$p.value[ENS_585_mbcn_wp_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_585_mbcn_wp_WC_summer_trnd.mask <- rasterToPoints(ENS_585_mbcn_wp_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_585_mbcn_wp_WC_summer_trnd$slope <- mask(ENS_585_mbcn_wp_WC_summer_trnd$slope, ENS_585_mbcn_wp_WC_summer_trnd$p.value) ###Plot only significant trend


Ref_mbcn_wp_WC_winter_trnd <- raster.kendall(Ref_mbcn_wp_WC_winter,  p.value = TRUE,confidence = TRUE)
Ref_mbcn_wp_WC_winter_trnd$slope <-Ref_mbcn_wp_WC_winter_trnd$slope*10
names(Ref_mbcn_wp_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
Ref_mbcn_wp_WC_winter_trnd$p.value[Ref_mbcn_wp_WC_winter_trnd$p.value >= 0.05] <- NA
Ref_mbcn_wp_WC_winter_trnd.mask <- rasterToPoints(Ref_mbcn_wp_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#Ref_mbcn_wp_WC_winter_trnd$slope <- mask(Ref_mbcn_wp_WC_winter_trnd$slope, Ref_mbcn_wp_WC_winter_trnd$p.value) ###Plot only significant trend


ENS_H_mbcn_wp_WC_winter_trnd <- raster.kendall(ENS_H_mbcn_wp_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_H_mbcn_wp_WC_winter_trnd$slope <-ENS_H_mbcn_wp_WC_winter_trnd$slope*10
names(ENS_H_mbcn_wp_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_mbcn_wp_WC_winter_trnd$p.value[ENS_H_mbcn_wp_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_H_mbcn_wp_WC_winter_trnd.mask <- rasterToPoints(ENS_H_mbcn_wp_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_mbcn_wp_WC_winter_trnd$slope <- mask(ENS_H_mbcn_wp_WC_winter_trnd$slope, ENS_H_mbcn_wp_WC_winter_trnd$p.value) ###Plot only significant trend

ENS_H_raw_wp_WC_winter_trnd <- raster.kendall(ENS_H_raw_wp_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_H_raw_wp_WC_winter_trnd$slope <-ENS_H_raw_wp_WC_winter_trnd$slope*10
names(ENS_H_raw_wp_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_H_raw_wp_WC_winter_trnd$p.value[ENS_H_raw_wp_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_H_raw_wp_WC_winter_trnd.mask <- rasterToPoints(ENS_H_raw_wp_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_H_raw_wp_WC_winter_trnd$slope <- mask(ENS_H_raw_wp_WC_winter_trnd$slope, ENS_H_raw_wp_WC_winter_trnd$p.value) ###Plot only significant trend


ENS_370_mbcn_wp_WC_winter_trnd <- raster.kendall(ENS_370_mbcn_wp_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_370_mbcn_wp_WC_winter_trnd$slope <-ENS_370_mbcn_wp_WC_winter_trnd$slope*10
names(ENS_370_mbcn_wp_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_370_mbcn_wp_WC_winter_trnd$p.value[ENS_370_mbcn_wp_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_370_mbcn_wp_WC_winter_trnd.mask <- rasterToPoints(ENS_370_mbcn_wp_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_370_mbcn_wp_WC_winter_trnd$slope <- mask(ENS_370_mbcn_wp_WC_winter_trnd$slope, ENS_370_mbcn_wp_WC_winter_trnd$p.value) ###Plot only significant trend


ENS_585_mbcn_wp_WC_winter_trnd <- raster.kendall(ENS_585_mbcn_wp_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_585_mbcn_wp_WC_winter_trnd$slope <-ENS_585_mbcn_wp_WC_winter_trnd$slope*10
names(ENS_585_mbcn_wp_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_585_mbcn_wp_WC_winter_trnd$p.value[ENS_585_mbcn_wp_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_585_mbcn_wp_WC_winter_trnd.mask <- rasterToPoints(ENS_585_mbcn_wp_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_585_mbcn_wp_WC_winter_trnd$slope <- mask(ENS_585_mbcn_wp_WC_winter_trnd$slope, ENS_585_mbcn_wp_WC_winter_trnd$p.value) ###Plot only significant trend


ENS_CC_370_WC_annual_trnd <- raster.kendall(ENS_CC_370_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_CC_370_WC_annual_trnd$slope <-ENS_CC_370_WC_annual_trnd$slope*10
names(ENS_CC_370_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_370_WC_annual_trnd$p.value[ENS_CC_370_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_CC_370_WC_annual_trnd.mask <- rasterToPoints(ENS_CC_370_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_370_WC_annual_trnd$slope <- mask(ENS_CC_370_WC_annual_trnd$slope, ENS_CC_370_WC_annual_trnd$p.value) ###Plot only significant trend


ENS_CC_585_WC_annual_trnd <- raster.kendall(ENS_CC_585_WC_annual,  p.value = TRUE,confidence = TRUE)
ENS_CC_585_WC_annual_trnd$slope <-ENS_CC_585_WC_annual_trnd$slope*10
names(ENS_CC_585_WC_annual_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_585_WC_annual_trnd$p.value[ENS_CC_585_WC_annual_trnd$p.value >= 0.05] <- NA
ENS_CC_585_WC_annual_trnd.mask <- rasterToPoints(ENS_CC_585_WC_annual_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_585_WC_annual_trnd$slope <- mask(ENS_CC_585_WC_annual_trnd$slope, ENS_CC_585_WC_annual_trnd$p.value) ###Plot only significant trend



ENS_CC_370_WC_summer_trnd <- raster.kendall(ENS_CC_370_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_CC_370_WC_summer_trnd$slope <-ENS_CC_370_WC_summer_trnd$slope*10
names(ENS_CC_370_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_370_WC_summer_trnd$p.value[ENS_CC_370_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_CC_370_WC_summer_trnd.mask <- rasterToPoints(ENS_CC_370_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_370_WC_summer_trnd$slope <- mask(ENS_CC_370_WC_summer_trnd$slope, ENS_CC_370_WC_summer_trnd$p.value) ###Plot only significant trend


ENS_CC_585_WC_summer_trnd <- raster.kendall(ENS_CC_585_WC_summer,  p.value = TRUE,confidence = TRUE)
ENS_CC_585_WC_summer_trnd$slope <-ENS_CC_585_WC_summer_trnd$slope*10
names(ENS_CC_585_WC_summer_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_585_WC_summer_trnd$p.value[ENS_CC_585_WC_summer_trnd$p.value >= 0.05] <- NA
ENS_CC_585_WC_summer_trnd.mask <- rasterToPoints(ENS_CC_585_WC_summer_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_585_WC_summer_trnd$slope <- mask(ENS_CC_585_WC_summer_trnd$slope, ENS_CC_585_WC_summer_trnd$p.value) ###Plot only significant trend



ENS_CC_370_WC_winter_trnd <- raster.kendall(ENS_CC_370_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_CC_370_WC_winter_trnd$slope <-ENS_CC_370_WC_winter_trnd$slope*10
names(ENS_CC_370_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_370_WC_winter_trnd$p.value[ENS_CC_370_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_CC_370_WC_winter_trnd.mask <- rasterToPoints(ENS_CC_370_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_370_WC_winter_trnd$slope <- mask(ENS_CC_370_WC_winter_trnd$slope, ENS_CC_370_WC_winter_trnd$p.value) ###Plot only significant trend


ENS_CC_585_WC_winter_trnd <- raster.kendall(ENS_CC_585_WC_winter,  p.value = TRUE,confidence = TRUE)
ENS_CC_585_WC_winter_trnd$slope <-ENS_CC_585_WC_winter_trnd$slope*10
names(ENS_CC_585_WC_winter_trnd) <- c("slope", "p.value", "LCI", "UCI")
ENS_CC_585_WC_winter_trnd$p.value[ENS_CC_585_WC_winter_trnd$p.value >= 0.05] <- NA
ENS_CC_585_WC_winter_trnd.mask <- rasterToPoints(ENS_CC_585_WC_winter_trnd$p.value, spatial=TRUE)##masking out significant trends as dots
#ENS_CC_585_WC_winter_trnd$slope <- mask(ENS_CC_585_WC_winter_trnd$slope, ENS_CC_585_WC_winter_trnd$p.value) ###Plot only significant trend



trend_WB_annual<- stack(Ref_mbcn_wp_WC_annual_trnd$slope,
                        ENS_H_raw_wp_WC_annual_trnd$slope,      
                 ENS_H_mbcn_wp_WC_annual_trnd$slope,
                 ENS_370_mbcn_wp_WC_annual_trnd$slope,
                 ENS_585_mbcn_wp_WC_annual_trnd$slope)


trend_WB_summer<- stack(Ref_mbcn_wp_WC_summer_trnd$slope,
                        ENS_H_raw_wp_WC_summer_trnd$slope,
                        ENS_H_mbcn_wp_WC_summer_trnd$slope,
                        ENS_370_mbcn_wp_WC_summer_trnd$slope,
                        ENS_585_mbcn_wp_WC_summer_trnd$slope)
  
trend_WB_winter<- stack(
  Ref_mbcn_wp_WC_winter_trnd$slope,
  ENS_H_raw_wp_WC_winter_trnd$slope,
  ENS_H_mbcn_wp_WC_winter_trnd$slope,
  ENS_370_mbcn_wp_WC_winter_trnd$slope,
  ENS_585_mbcn_wp_WC_winter_trnd$slope)


WB_CC<-stack(ENS_CC_370_WC_annual_trnd$slope,ENS_CC_585_WC_annual_trnd$slope,
            ENS_CC_370_WC_summer_trnd$slope,ENS_CC_585_WC_summer_trnd$slope,
            ENS_CC_370_WC_winter_trnd$slope,ENS_CC_585_WC_winter_trnd$slope)

########### start from here!!!!!!!!!!!







rasterNames_trend_WB_annual  <- c("(a) Reference - Annual ","(b) Raw CMIP6","(c) Ens - Historical - Annual","(d) Ens - SSP370 - Annual","(e) Ens - SSP585 - Annual" )

rasterNames_trend_WB_summer  <- c("(a) Reference - Summer","(b) Ens - Historical - Summer","(c) Ens - SSP370 - Summer"," (d) Ens - SSP585 - Summer" )

rasterNames_trend_WB_winter  <- c( "(a) Reference - Winter","(b) Ens - Historical - Winter","(c) Ens - SSP370 - Winter"," (d) Ens - SSP585 - Winter" )



rasterNames_trend_CC  <- c("(a) SSP585 - Historical - Annual", "(b) SSP370 - Historical - Annual",
                           "(c) SSP585 - Historical - Summer", "(d) SSP370 - Historical - Summer",
                           "(e) SSP585 - Historical - Winter", "(f) SSP370 - Historical - Winter")


brk_annual <- round(seq(min(c(ENS_585_mbcn_wp_WC_annual_trnd@data@min,
                              ENS_370_mbcn_wp_WC_annual_trnd@data@min,
                              ENS_H_mbcn_wp_WC_annual_trnd@data@min,
                              Ref_mbcn_wp_WC_annual_trnd@data@min)),
                 max(c(ENS_585_mbcn_wp_WC_annual_trnd@data@max,
                       ENS_370_mbcn_wp_WC_annual_trnd@data@max,
                       ENS_H_mbcn_wp_WC_annual_trnd@data@max,
                       Ref_mbcn_wp_WC_annual_trnd@data@max)),by=48),0);brk_annual# where the colors change


brk_summer <- round(seq(min(c(ENS_585_mbcn_wp_WC_summer_trnd@data@min,
                              ENS_370_mbcn_wp_WC_summer_trnd@data@min,
                              ENS_H_mbcn_wp_WC_summer_trnd@data@min,
                              Ref_mbcn_wp_WC_summer_trnd@data@min)),
                        max(c(ENS_585_mbcn_wp_WC_summer_trnd@data@max,
                              ENS_370_mbcn_wp_WC_summer_trnd@data@max,
                              ENS_H_mbcn_wp_WC_summer_trnd@data@max,
                              Ref_mbcn_wp_WC_summer_trnd@data@max)),by=22),0)


brk_winter <- round(seq(min(c(ENS_585_mbcn_wp_WC_winter_trnd@data@min,
                              ENS_370_mbcn_wp_WC_winter_trnd@data@min,
                              ENS_H_mbcn_wp_WC_winter_trnd@data@min,
                              Ref_mbcn_wp_WC_winter_trnd@data@min)),
                        max(c(ENS_585_mbcn_wp_WC_winter_trnd@data@max,
                              ENS_370_mbcn_wp_WC_winter_trnd@data@max,
                              ENS_H_mbcn_wp_WC_winter_trnd@data@max,
                              Ref_mbcn_wp_WC_winter_trnd@data@max)),by=48),0)

brk_winter_CC <- round(seq(min(c(ENS_CC_585_WC_winter_trnd@data@min,
                              ENS_CC_370_WC_winter_trnd@data@min )),
                        max(c(ENS_CC_585_WC_winter_trnd@data@max,
                              ENS_CC_370_WC_winter_trnd@data@max)),by=48),0);brk_winter_CC

brk_summer_CC <- round(seq(min(c(ENS_CC_585_WC_summer_trnd@data@min,
                                 ENS_CC_370_WC_summer_trnd@data@min )),
                           max(c(ENS_CC_585_WC_summer_trnd@data@max,
                                 ENS_CC_370_WC_summer_trnd@data@max)),by=21),0);brk_summer_CC

brk_annual_CC <- round(seq(min(c(ENS_CC_585_WC_annual_trnd@data@min,
                                 ENS_CC_370_WC_annual_trnd@data@min )),
                           max(c(ENS_CC_585_WC_annual_trnd@data@max,
                                 ENS_CC_370_WC_annual_trnd@data@max)),by=48),0);brk_annual_CC


lab_annual <- round(seq(min(c(ENS_585_mbcn_wp_WC_annual_trnd@data@min,
                              ENS_370_mbcn_wp_WC_annual_trnd@data@min,
                              ENS_H_mbcn_wp_WC_annual_trnd@data@min,
                              Ref_mbcn_wp_WC_annual_trnd@data@min)),
                        max(c(ENS_585_mbcn_wp_WC_annual_trnd@data@max,
                              ENS_370_mbcn_wp_WC_annual_trnd@data@max,
                              ENS_H_mbcn_wp_WC_annual_trnd@data@max,
                              Ref_mbcn_wp_WC_annual_trnd@data@max)),by=48),0)# where the colors change


lab_summer <- round(seq(min(c(ENS_585_mbcn_wp_WC_summer_trnd@data@min,
                              ENS_370_mbcn_wp_WC_summer_trnd@data@min,
                              ENS_H_mbcn_wp_WC_summer_trnd@data@min,
                              Ref_mbcn_wp_WC_summer_trnd@data@min)),
                        max(c(ENS_585_mbcn_wp_WC_summer_trnd@data@max,
                              ENS_370_mbcn_wp_WC_summer_trnd@data@max,
                              ENS_H_mbcn_wp_WC_summer_trnd@data@max,
                              Ref_mbcn_wp_WC_summer_trnd@data@max)),by=22),0)


lab_winter <- round(seq(min(c(ENS_585_mbcn_wp_WC_winter_trnd@data@min,
                              ENS_370_mbcn_wp_WC_winter_trnd@data@min,
                              ENS_H_mbcn_wp_WC_winter_trnd@data@min,
                              Ref_mbcn_wp_WC_winter_trnd@data@min)),
                        max(c(ENS_585_mbcn_wp_WC_winter_trnd@data@max,
                              ENS_370_mbcn_wp_WC_winter_trnd@data@max,
                              ENS_H_mbcn_wp_WC_winter_trnd@data@max,
                              Ref_mbcn_wp_WC_winter_trnd@data@max)),by=48),0)

lab_winter_CC <- round(seq(min(c(ENS_CC_585_WC_winter_trnd@data@min,
                                 ENS_CC_370_WC_winter_trnd@data@min )),
                           max(c(ENS_CC_585_WC_winter_trnd@data@max,
                                 ENS_CC_370_WC_winter_trnd@data@max)),by=48),0);lab_winter_CC

lab_summer_CC <- round(seq(min(c(ENS_CC_585_WC_summer_trnd@data@min,
                                 ENS_CC_370_WC_summer_trnd@data@min )),
                           max(c(ENS_CC_585_WC_summer_trnd@data@max,
                                 ENS_CC_370_WC_summer_trnd@data@max)),by=21),0);lab_summer_CC

lab_annual_CC <- round(seq(min(c(ENS_CC_585_WC_annual_trnd@data@min,
                                 ENS_CC_370_WC_annual_trnd@data@min )),
                           max(c(ENS_CC_585_WC_annual_trnd@data@max,
                                 ENS_CC_370_WC_annual_trnd@data@max)),by=48),0);lab_annual_CC



myColorkey_annual<- list(at=brk_annual,space="right", title = list("Water Conservation trend [ mm /decade] ",cex = 0.8,fontface = 7,col = 'black'),
                  title.gpar = list(" ",cex = 0.001,font = 0.5, col = 'white'),#this is stubborn, so i reduced the font because i dont want to see it on the map
                  
                  ## where the colors change
                  labels=list(at=lab_annual,rot=0,font=6,fontface=7,cex=0.8), #fontface="bold"## where to print labels
                  axis.line=list(col='black'),#legend outline (including its ticks)
                  height=1,width=1.6, #height and width of legend
                  tri.lower = TRUE, tri.upper = TRUE) #adds triangular arrows at both the top and end of legend

require(colorRamps)
cols=colorRampPalette(matlab.like2(255))
#cols<-colorRampPalette(rev(c("navyblue","blue","steelblue1","slategray1","red","red3")))
cols<-colorRampPalette(rev(c("navyblue","blue","white","pink","red")))
A<-levelplot(Ref_mbcn_wp_WC_annual_trnd$slope,col.regions=cols,colorkey=F,
          par.settings=list(panel.background=list(col="white"),
                            axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                            strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
          layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
          scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
          main=list('',cex=1,col="black"), ylim=c(-48, 76),
          margin=F, auto.key=F,maxpixels=2e5,xlab=list('Longitude',fontface='bold'),
          ylab=list('Latitude', rot=90, fontface='bold'))+
  layer(sp.polygons(wp, lwd=0.6,col="black"))+
  layer(sp.points(Ref_mbcn_wp_WC_annual_trnd.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(a) Reference - Annual',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(A, aspect=0.5) 

B<-levelplot(ENS_H_raw_wp_WC_annual_trnd$slope,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('Longitude',fontface='bold'),
             ylab=list('Latitude', rot=90, fontface='bold'))+
  layer(sp.polygons(wp, lwd=0.6,col="black"))+
  layer(sp.points(ENS_H_raw_wp_WC_annual_trnd.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(b) Raw CMIP6 Historical - Annual',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(B, aspect=0.5) 


C<-levelplot(ENS_H_mbcn_wp_WC_annual_trnd$slope,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('Longitude',fontface='bold'),
             ylab=list('Latitude', rot=90, fontface='bold'))+
  layer(sp.polygons(wp, lwd=0.6,col="black"))+
  layer(sp.points(ENS_H_mbcn_wp_WC_annual_trnd.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(c) MBCN Historical - Annual',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(C, aspect=0.5)

D<-levelplot(ENS_370_mbcn_wp_WC_annual_trnd$slope,col.regions=cols,colorkey=F,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('Longitude',fontface='bold'),
             ylab=list('Latitude', rot=90, fontface='bold'))+
  layer(sp.polygons(wp, lwd=0.6,col="black"))+
  layer(sp.points(ENS_370_mbcn_wp_WC_annual_trnd.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(d) MBCN SSP370 - Annual',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(D, aspect=0.5)

E<-levelplot(ENS_585_mbcn_wp_WC_annual_trnd$slope,col.regions=cols,colorkey=myColorkey_annual,
             par.settings=list(panel.background=list(col="white"),
                               axis.line=list(lwd=0.8,col="blue"), #background colour, axis line 
                               strip.border=list(lwd=1.1),layout.heights=list(xlab.key.padding=0.1)), #border line, adjusts layout ,xlab.key.padding adjust the closeness of the horizontal legend to the main plot
             layout.widths=list(left.padding=4,right.padding=10),#adjust the map positioning#parSettings, 
             scales=list(x=list(draw=T), y=list(draw=T),relation="free",tick.number=6,alternating=2,cex=0.6,col="black"), # remove axes labels & ticks if set to false
             main=list('',cex=1,col="black"), ylim=c(-48, 76),
             margin=F, auto.key=F,maxpixels=2e5,xlab=list('Longitude',fontface='bold'),
             ylab=list('Latitude', rot=90, fontface='bold'))+
  layer(sp.polygons(wp, lwd=0.6,col="black"))+
  layer(sp.points(ENS_585_mbcn_wp_WC_annual_trnd.mask, pch=1, cex=0.1, alpha=1, col="black"))+
  layer(ltext(x=80,y=-51.8,labels = '(e) MBCN SSP585 - Annual',cex=0.9,col="black",fontface = 7))  #when combinining the plots, the titles are missing. This adds text to lattice plot
update(E, aspect=0.5)

cObj1<- c(A,B,C,D,E,  layout = c(2,3), merge.legends = T,x.same = F, y.same = F)

update(cObj1,fontface = 7,scales = list(draw=T,y = list(rot = 90)),xlab=list("Longitude",cex=1,fontface = 7),ylab=list("Latitude",cex=1,fontface = 7), 
       main=list('',fontface = 7)) 







trend_WB_summer<- stack(Ref_mbcn_wp_WC_summer_trnd$slope,
                        ENS_H_raw_wp_WC_summer_trnd$slope,
                        ENS_H_mbcn_wp_WC_summer_trnd$slope,
                        ENS_370_mbcn_wp_WC_summer_trnd$slope,
                        ENS_585_mbcn_wp_WC_summer_trnd$slope)

trend_WB_winter<- stack(
  Ref_mbcn_wp_WC_winter_trnd$slope,
  ENS_H_raw_wp_WC_winter_trnd$slope,
  ENS_H_mbcn_wp_WC_winter_trnd$slope,
  ENS_370_mbcn_wp_WC_winter_trnd$slope,
  ENS_585_mbcn_wp_WC_winter_trnd$slope)


WB_CC<-stack(ENS_CC_370_WC_annual_trnd$slope,ENS_CC_585_WC_annual_trnd$slope,
             ENS_CC_370_WC_summer_trnd$slope,ENS_CC_585_WC_summer_trnd$slope,
             ENS_CC_370_WC_winter_trnd$slope,ENS_CC_585_WC_winter_trnd$slope)
