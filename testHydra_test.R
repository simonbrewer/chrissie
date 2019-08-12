###############################################################################
## Runs version of hydra model modified from original code
##
## Ref: 
## M.T. Coe (1998). A linked global model of terrestrial hydrologic processes:
## Simulation of modern rivers, lakes, and wetlands. JGR, 103, D8, 8885-8899
##
## M.T. Coe (2000) Modeling terrestrial hydrological systems at the continental 
## scale: Testing the accuracy of an atmospheric GCM. J.Clim., 13, 686-704
###############################################################################

###############################################################################
## Libraries
library(rgdal)
library(raster)
library(RColorBrewer)
source("helpers.R")

###############################################################################
## Lake border
# cllake = readOGR("./chrissie_lake/chrissie_lake.shp")

###############################################################################
## MODEL SETUP
## Parameters
delt = 60*60 ## Time step (s)
deltu = 24 ## Number of time steps to run model for
bpf = 0 ## Proportion of runoff to put in baseflow [0-1]
effvol = 0.3 

## Base files
oldw <- getOption("warn")
options(warn = -1)

dem.r = raster("./dem/test_dem.nc")
myext = extent(dem.r)
bas.r = raster("./dem/test_bas.nc")
ldd.r = raster("./dem/test_ldd.nc")
mflac.r = raster("./dem/test_mflac.nc")
iout.r = raster("./dem/test_iout.nc")
jout.r = raster("./dem/test_jout.nc")

options(warn = oldw)

# dem.r = raster("./dem/chrissie_dem.nc")
# bas.r = raster("./dem/chrissie_bas.nc")
# ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")
# ldd.r = crop(ldd.r, extent(dem.r))
# mflac.r = raster("./dem/cllake_mflac.nc")
# iout.r = raster("./dem/cllake_iout.nc")
# jout.r = raster("./dem/cllake_jout.nc")

## Replace nulls
dem.r[is.na(dem.r)] <- 0
bas.r[is.na(bas.r)] <- 0
ldd.r[is.na(ldd.r)] <- 0
mflac.r[is.na(mflac.r)] <- 0
iout.r[is.na(iout.r)] <- 0
jout.r[is.na(jout.r)] <- 0

## Clip lake basin
# mask.r = bas.r == 18872 ## 18872 is the basin corresponding to Chrissie lake
mask.r = dem.r > 0

###############################################################################
## Estimate cell areas
area.r = area(dem.r) * 1e6

###############################################################################
## Grid sizes for outout
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]

###############################################################################
## Forcing data
dpre.stk = crop(brick("./inputs/dpre_cl.nc"), myext) * 2
dpre.stk[is.na(dpre.stk)] <- 0
dpet.stk = crop(brick("./inputs/dpet_cl.nc"), myext)
dpet.stk[is.na(dpet.stk)] <- 0
devp.stk = crop(brick("./inputs/devp_cl.nc"), myext)
devp.stk[is.na(devp.stk)] <- 0
dcn.stk = crop(brick("./inputs/dcn_cl.nc"), myext)
dcn.stk[is.na(dcn.stk)] <- 0

###############################################################################
## Calculate runoff
dro.stk = (dpre.stk + dcn.stk) - dpet.stk

## Set constant values
dro.stk[] <- 1
dpre.stk[] <- 1
devp.stk[] <- 0
area.r [] <- 1000

## Plot 'em
dpre.reg = cellStats(dpre.stk, mean)
dpet.reg = cellStats(dpet.stk, mean)
devp.reg = cellStats(devp.stk, mean)
dro.reg = cellStats(dro.stk, mean)

plot(dpre.reg, type='l', ylim=c(0,7))
lines(dpet.reg, col=2)
lines(devp.reg, col=3)
lines(dro.reg, col=4)
###############################################################################
## Water storage rasters
wse.r = wvl.r = war.r = setValues(dem.r, 0)

###############################################################################
## Convert to matrices
dem = as.matrix(dem.r)
ldd = as.matrix(ldd.r)
bas = as.matrix(bas.r)
mflac = as.matrix(mflac.r)
iout = as.matrix(iout.r)
jout = as.matrix(jout.r)
outelev = as.matrix(dem.r)  # Needs changing
mask = as.matrix(mask.r)*1
cella = as.matrix(area.r)
# celld = as.matrix(dist.r)
wvl = as.matrix(wvl.r)
wse = as.matrix(wse.r)
war = as.matrix(war.r)

pre = as.array(dpre.stk)
evp = as.array(devp.stk)
ro = as.array(clamp(dro.stk, lower=0, useValues=TRUE))
sro = ro * (1-bpf)
bro = ro * bpf

## Test with 5m everywhere
# wse = wse + 5

###############################################################################
cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)

nyrs = 1
nspin = 1
ndays = 365

sim.out = rhydra(gridx, gridy, nyrs, ndays, startyear=1, 
                 res=30, converg=1, laket=0, spin=nspin,
                 dem=dem, mask=mask, area=cella, rivdir=ldd, mflac=mflac,
                 outnewi=iout, outnewj=jout, basin=bas, 
                 prcpi=pre, evapi=evp, runin=sro, drainin=bro) 

tmp.r = setValues(dem.r, matrix(sim.out$outelv, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r - dem.r)

tmp.r = setValues(dem.r, matrix(sim.out$larea, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r)
stop()
tmp.r = setValues(dem.r, matrix(sim.out$volr, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r)
plot(cllake, add=TRUE)


