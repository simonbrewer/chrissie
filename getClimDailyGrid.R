###############################################################################
## getClimDailyGrid.R
## 
## Basic script to extract data for single year from CRU and run through Splash
## and Morton's CRWE/CRLE across a set of grid points
##
## Will need reference DEM for extracting climate
###############################################################################

## Extract climate data and run splash
require(raster)
require(rgdal)
require(Evapotranspiration)
source("./helpers.R")

## Location
lon = 20.0833
lat = -30.75
yr = 1975
elv = 929
whc = 150.0
mdays = c(31,28,31,30,31,30,31,31,30,31,30,31)

## Read raster grid file
dem.r = raster("./dem/chrissie_dem.nc")
bas.r = raster("./dem/chrissie_bas.nc")
ncell.bc = ncell(dem.r)

## Lake center 
bccent = SpatialPoints(cbind(lon,lat))

## Lake polygon
bclake = readOGR("~/Documents/grassdata/hydrosheds/bclake/chrissie_lakes/chrissie_lakes.shp")
plot(dem.r)
plot(bclake, add=TRUE)

plot(bas.r == 18872) ## 18872 is the basin corresponding to Chrissie lake
plot(bclake, add=TRUE)

## Get monthly series
tmp.r = stack("./cru/cl/cl_cru_10min_tmp.nc")
dtr.r = stack("./cru/cl/cl_cru_10min_dtr.nc")
pre.r = stack("./cru/cl/cl_cru_10min_pre.nc")
rhm.r = stack("./cru/cl/cl_cru_10min_reh.nc")
sun.r = stack("./cru/cl/cl_cru_10min_sun.nc")
elv.r = raster("./cru/cl/cl_cru_10min_elv.nc")
ncell.cru = ncell(tmp.r)

## Output stacks
dcn.stk = devp.stk = dpet.stk = dpre.stk = 
  brick(nrows = nrow(tmp.r), ncols = ncol(tmp.r), 
        xmn = xmin(tmp.r), xmx = xmax(tmp.r), 
        ymn = ymin(tmp.r), ymx = ymax(tmp.r), nl = 365)
# dpre.stk = array(dim = c(nrow(dem.r) * ncol(dem.r), 365))
## Load constants for CRLE run
data(constants)

## Main loop
for (i in 1:ncell.cru) {
# for (i in 1:10) {
    ## Counter
  print(paste("Doing",i,"of",ncell.cru))
  
  ## Get current cell location
  curr.xy = xyFromCell(dem.r, i)
  
  ## Extract CRU climate
  # tmp.bc = extract(tmp.r, curr.xy) #- 273.15
  # dtr.bc = extract(dtr.r, curr.xy) 
  # tmn.bc = tmp.bc - (dtr.bc/2)
  # tmx.bc = tmp.bc + (dtr.bc/2)
  # pre.bc = extract(pre.r, curr.xy) #* 60*60*24
  # rhm.bc = extract(rhm.r, curr.xy) 
  # sun.bc = extract(sun.r, curr.xy) / 100
  # elv = extract(elv.r, curr.xy) * 1000

  tmp.bc = tmp.r[i] #- 273.15
  dtr.bc = dtr.r[i]
  tmn.bc = tmp.bc - (dtr.bc/2)
  tmx.bc = tmp.bc + (dtr.bc/2)
  pre.bc = pre.r[i] #* 60*60*24
  rhm.bc = rhm.r[i]
  sun.bc = sun.r[i] / 100
  elv = elv.r[i] * 1000
  
  if (!is.na(elv)) {
    
    ## Get daily values
    dtmp = daily(c(tmp.bc))$dly
    dtmn = daily(c(tmn.bc))$dly
    dtmx = daily(c(tmx.bc))$dly
    dpre = daily(c(pre.bc/mdays))$dly
    drhm = daily(c(rhm.bc))$dly
    dsun = daily(c(sun.bc))$dly
    
    ## sunp is listed as 
    ## sunp	sunshine percent of maximum possible (percent of daylength)
    
    ## Run SPLASH
    aetpet.df = splashf(dtmp,dpre,dsun,lat=lat,
                        yr=yr, elv=elv)
    
    ## Estimate Tdew
    ## From humidity.to.dewpoint from the weathermetrics package
    ## dewpoint <- (rh/100)^(1/8) * (112 + (0.9 * t)) - 112 + (0.1 * t)
    ## This is originally from the source code for the 
    ## US National Weather Service's online heat index calculator.
    
    ddew = (drhm/100)^(1/8) * (112 + (0.9 * dtmp)) - 112 + (0.1 * dtmp)
    
    ## Run Morton CRWE
    bcdates = seq.Date(as.Date("1975-01-01"), length.out = 365, by = 1)
    
    ## Using n for sunshine hours
    clim.df = data.frame(Year = as.numeric(format(bcdates, "%Y")),
                         Month = as.numeric(format(bcdates, "%m")),
                         Day = as.numeric(format(bcdates, "%d")),
                         Tmax = dtmx, Tmin = dtmn, 
                         n = aetpet.df$dsl, Tdew = ddew, Precip = dpre)
    
    clim.in = ReadInputs(varnames = c("Tmax","Tmin","Tdew","n"),
                         clim.df, 
                         stopmissing=c(10,10,3))
    
    ## Update constants for CRLE run
    constants$lat = curr.xy[1,2]
    constants$lat_rad = lat * pi / 180
    constants$Elev = elv
    constants$PA = sum(dpre)
    
    lake.out <- ET.MortonCRWE(clim.in, constants, ts="monthly",
                              est="shallow lake ET", solar="sunshine hours", Tdew= TRUE, 
                              alpha = NULL, message="no", save.csv="no")
    
    ## Output
    dpre.stk[i] <- dpre
    dpet.stk[i] <- aetpet.df$dpet
    dcn.stk[i] <- aetpet.df$dcn
    devp.stk[i] <- daily(lake.out$ET.MonthlyAve)$dly
    
  }
  
}

writeRaster(dpre.stk, "./inputs/dpre_cru.nc", format='CDF', overwrite = TRUE)
writeRaster(dpet.stk, "./inputs/dpet_cru.nc", format='CDF', overwrite = TRUE)
writeRaster(dcn.stk, "./inputs/dcn_cru.nc", format='CDF', overwrite = TRUE)
writeRaster(devp.stk, "./inputs/devp_cru.nc", format='CDF', overwrite = TRUE)
