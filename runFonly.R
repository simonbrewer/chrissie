## Only runs the fprtran code
library(rgdal)
library(raster)
library(RColorBrewer)
source("helpers.R")

###############################################################################
## Lake border
cllake = readOGR("./chrissie_lake/chrissie_lake.shp")

###############################################################################
## Climate data
load("testdata.RData")

###############################################################################
## Lake border
cllake = readOGR("./chrissie_lake/chrissie_lake.shp")

###############################################################################
## MODEL SETUP
## Parameters
delt = 60*60 ## Time step (s)
deltu = 24 ## Number of time steps to run model for
bpf = 0 ## Proportion of runoff to put in baseflow [0-1]
effvol = 0.3 

## Base files
dem.r = raster("./dem/chrissie_dem.nc")
bas.r = raster("./dem/chrissie_bas.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")
ldd.r = crop(ldd.r, extent(dem.r))
mflac.r = raster("./dem/cllake_mflac.nc")
iout.r = raster("./dem/cllake_iout.nc")
jout.r = raster("./dem/cllake_jout.nc")

## Clip lake basin
mask.r = bas.r == 18872 ## 18872 is the basin corresponding to Chrissie lake

###############################################################################
## Grid sizes for outout
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]

nyrs = 1
ndays = 365
sim.out = rhydra(gridx, gridy, nyrs, ndays, startyear=1, 
                 res=30, converg=1, laket=0, spin=1,
                 dem=dem, mask=mask, area=cella, rivdir=ldd, mflac=mflac,
                 outnewi=iout, outnewj=jout, basin=bas, 
                 prcpi=pre, evapi=evp, runin=sro, drainin=bro) 
