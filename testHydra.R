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
bclake = readOGR("~/Documents/grassdata/hydrosheds/bclake/929m lake level.kml")

## Lake center 
lon = 20.0833
lat = -30.75
bccent = SpatialPoints(cbind(lon,lat))

###############################################################################
## MODEL SETUP
## Parameters
delt = 60*60 ## Time step (s)
deltu = 24 ## Number of time steps to run model for
bpf = 0 ## Proportion of runoff to put in baseflow [0-1]
effvol = 0.3 

## Base files
dem.r = raster("./dem/bclake_dem.nc")
bas.r = raster("./dem/bclake_bas.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")
ldd.r = crop(ldd.r, extent(dem.r))
mflac.r = raster("./dem/bclake_mflac.nc")
iout.r = raster("./dem/bclake_iout.nc")
jout.r = raster("./dem/bclake_jout.nc")

## Clip lake basin
mask.r = bas.r == 19238

###############################################################################
## Estimate cell areas
area.r = area(dem.r) * 1e6

## Distance between cells (set as constant)
## Needs to be calculated in fortran from coordinates
dist.r =  setValues(dem.r, 857) ## Approximately 857m cell centers

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
dpre.stk = brick("./inputs/dpre_bc.nc") * 5
dpet.stk = brick("./inputs/dpet_bc.nc")
devp.stk = brick("./inputs/devp_bc.nc")
dcn.stk = brick("./inputs/dcn_bc.nc")

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
mask = as.matrix(mask.r)
cella = as.matrix(area.r)
celld = as.matrix(dist.r)
wvl = as.matrix(wvl.r)
wse = as.matrix(wse.r)
war = as.matrix(war.r)

pre = as.array(dpre.stk)
evp = as.array(devp.stk)
ro = as.array(clamp(dpre.stk, lower=0, useValues=TRUE))
sro = ro * (1-bpf)
bro = ro * bpf

## Test with 5m everywhere
# wse = wse + 5

###############################################################################
cols <- colorRampPalette(brewer.pal(9,"Blues"))(100)

nyrs = 1
ndays = 365

sim.out = rhydra(gridx, gridy, nyrs, ndays, startyear=1, 
                 res=30, converg=1, laket=0, spin=1,
                 dem=dem, mask=mask, area=cella, rivdir=ldd, mflac=mflac,
                 outnewi=iout, outnewj=jout, basin=bas, 
                 prcpi=pre, evapi=evp, runin=sro, drainin=bro) 

tmp.r = setValues(dem.r, matrix(sim.out$outelv, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r - dem.r)
plot(bclake, add=TRUE)

tmp.r = setValues(dem.r, matrix(sim.out$larea, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r)
plot(bclake, add=TRUE)
stop()
tmp.r = setValues(dem.r, matrix(sim.out$volr, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(tmp.r)
plot(bclake, add=TRUE)


