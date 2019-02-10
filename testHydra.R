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
cllake = readOGR("./chrissie_lake/chrissie_lake.shp")

###############################################################################
## MODEL SETUP
## Parameters
delt = 60*60 ## Time step (s)
deltu = 24 ## Number of time steps to run model for
bpf = 0 ## Proportion of runoff to put in baseflow [0-1]
effvol = 0.3 

## Base files
myext = extent(c(30,31,-27,-26))
dem.r = crop(raster("./dem/chrissie_dem.nc"), myext)
bas.r = crop(raster("./dem/chrissie_bas.nc"), myext)
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")
ldd.r = crop(ldd.r, extent(dem.r))
mflac.r = crop(raster("./dem/cllake_mflac.nc"), myext)
iout.r = crop(raster("./dem/cllake_iout.nc"), myext)
jout.r = crop(raster("./dem/cllake_jout.nc"), myext)

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
mask.r = bas.r == 18872 ## 18872 is the basin corresponding to Chrissie lake

###############################################################################
## Estimate cell areas
area.r = area(dem.r) * 1e6

###############################################################################
## Grid sizes for outout
gridx = dim(dem.r)[1]
gridy = dim(dem.r)[2]

###############################################################################
## Assign outlet cell
out.x = 20.364
out.y = -30.88
out.sp = SpatialPoints(cbind(out.x,out.y))
out.cell <- cellFromXY(dem.r, out.sp)
outlet.r = setValues(dem.r, 0)
outlet.r[out.cell] <- 1

###############################################################################
## Forcing data
dpre.stk = crop(brick("./inputs/dpre_cl.nc"), myext)
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
ndays = 365

sim.out = rhydra(gridx, gridy, nyrs, ndays, startyear=1, 
                 res=30, converg=0, laket=0, spin=50,
                 dem=dem, mask=mask, area=cella, rivdir=ldd, mflac=mflac,
                 outnewi=iout, outnewj=jout, basin=bas, 
                 prcpi=pre, evapi=evp, runin=sro, drainin=bro) 

tmp.r = setValues(dem.r, matrix(sim.out$outelv, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r - dem.r)
plot(cllake, add=TRUE)

tmp.r = setValues(dem.r, matrix(sim.out$larea, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r)
plot(cllake, add=TRUE)

tmp.r = setValues(dem.r, matrix(sim.out$volr, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r, col=RColorBrewer::brewer.pal(7, "Blues"))
plot(cllake, add=TRUE)


