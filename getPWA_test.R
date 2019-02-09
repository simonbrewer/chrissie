###############################################################################
## Potential water area estimates
##
###############################################################################

require(raster)
require(rgdal)
source("helpers.R")

dem.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_dem_30g.nc")
ldd.r = raster("~/Documents/grassdata/hydrosheds/bclake/af_ldd_30g.nc")

bc.r = raster("./dem/chrissie_dem.nc")
myext = extent(bc.r)
# myext = extent(c(19.9,20.3,-30.95,-30.6))
# myext = extent(c(17.5,18.5,-31.5,-31))
plot(crop(ldd.r, myext))

dem.r = crop(dem.r, myext)
ldd.r = crop(ldd.r, myext)
###############################################################################
## Lake border
cllake = readOGR("./chrissie_lake/chrissie_lake.shp")

plot(dem.r)
plot(cllake, add=TRUE)

#dem.r = crop(dem.r, cllake)
#ldd.r = crop(ldd.r, cllake)
## r.watershed directions
## "Provides the "aspect" for each cell measured CCW from East. 
## Multiplying positive values by 45 will give the direction in degrees 
## that the surface runoff will travel from that cell."
# 1: 45:  NE
# 2: 90:  N
# 3: 135: NW
# 4: 180: W
# 5: 225: SW
# 6: 270: S
# 7: 315: SE
# 8: 360: E

offx = c(1, 0, -1, -1, -1, 0, 1, 1) ## cols
offy = c(-1, -1, -1, 0, 1, 1, 1, 0) ## rows

## Buld mask
mask.r = setValues(dem.r, 1)
mask.r[is.na(dem.r)] <- 0
## Set NAs to -9999
dem.r[is.na(dem.r)] <- -9999
ldd.r[is.na(ldd.r)] <- -9999

gridx = nrow(dem.r)
gridy = ncol(dem.r)
dem = as.matrix(dem.r)
ldd = as.matrix(ldd.r)
mask = as.matrix(mask.r)

pwa.out = fpwa(gridx, gridy, dem, ldd, mask)

pwa.r = setValues(dem.r, matrix(pwa.out$pwa, 
                                  nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(pwa.r)
plot(cllake, add=TRUE)
writeRaster(pwa.r, "./dem/cllake_pwa.nc", format="netCDF", overwrite=TRUE)

drain.r = setValues(dem.r, matrix(pwa.out$drain, 
                                nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(drain.r)
plot(cllake, add=TRUE)
writeRaster(drain.r, "./dem/cllake_drain.nc", format="netCDF", overwrite=TRUE)

outelev.r = setValues(dem.r, matrix(pwa.out$outelev, 
                                  nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(outelev.r-dem.r)
plot(cllake, add=TRUE)
mflac.r = outelev.r*pwa.r
writeRaster(mflac.r, "./dem/cllake_mflac.nc", format="netCDF", overwrite=TRUE)

iout.r = setValues(dem.r, matrix(pwa.out$iout, 
                                 nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(iout.r)
plot(cllake, add=TRUE)
writeRaster(iout.r, "./dem/cllake_iout.nc", format="netCDF", overwrite=TRUE)

jout.r = setValues(dem.r, matrix(pwa.out$jout, 
                                 nrow=dim(dem.r)[1], ncol=dim(dem.r)[2]))
plot(jout.r)
plot(cllake, add=TRUE)
writeRaster(jout.r, "./dem/cllake_jout.nc", format="netCDF", overwrite=TRUE)

