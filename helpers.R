###############################################################################
## Helper functions for "runSplash.R"
dyn.load("./fortran/splash.so")
dyn.load("./fortran/getpwa.so")
dyn.load("./fortran/rhydra.so")
###############################################################################

###############################################################################
## aetpet: Function to calculate aet and pet
splashf <- function(dtemp,dprec,dsun,lat,yr,elv) {
  dyn.load('fortran/splash.so')
  retdata <- .Fortran("spin_up_sm",
                      yr = as.integer(yr),
                      lat = as.double(lat),
                      elv = as.double(elv),
                      pr = as.double(dprec),
                      tc = as.double(dtemp),
                      sf = as.double(dsun),
                      maet = double(12),
                      mpet = double(12),
                      mcn = double(12),
                      mro = double(12),
                      msm = double(12),
                      daet = double(365),
                      dpet = double(365),
                      dcn = double(365),
                      dro = double(365),
                      dsm = double(365),
                      sm = double(1),
                      ddl = double(365),
                      dsl = double(365))
  return(retdata)
}
###############################################################################

###############################################################################
# ## DAILY: Function to interpolate from monthly to daily
# ## Replace with 'approx'?
# 
daily <- function(mly) {
  retdata <- .Fortran("daily",
                      mly = as.double(mly),
                      dly = double(365))
  
  return(retdata)
}
###############################################################################

###############################################################################
## Potential water area code (fortran)
fpwa <- function(gridx, gridy, dem, ldd, mask) {
  
  pwaout = .Fortran("getpwa",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), 
                   ldd = as.integer(ldd),
                   mask = as.integer(mask),
                   pwa = integer(gridx*gridy),
                   drain = integer(gridx*gridy),
                   outelev = double(gridx*gridy),
                   iout = integer(gridx*gridy),
                   jout = integer(gridx*gridy)
                   )
  return(pwaout)
  
}
###############################################################################

###############################################################################
## Calls rhydra fortran routine
## Should be able to eliminate converg, laket
rhydra <- function(gridx, gridy, nyrs, ndays, startyear, 
                   res = 30, converg = 1, laket = 0, spin = 1,
                   dem, mask, area, rivdir, mflac,
                   outnewi, outnewj, basin, 
                   prcpi, evapi, runin, drainin) 
  {
  
  simcf = .Fortran("rhydra",
                   nc = as.integer(gridx), nr = as.integer(gridy),
                   nyrs = as.integer(nyrs), ndays = as.integer(ndays),
                   startyear = as.integer(startyear), 
                   res = as.integer(res),
                   converg = as.integer(converg), 
                   laket = as.integer(laket), 
                   spin = as.integer(spin),
                   dem = as.double(dem), 
                   mask = as.integer(mask), 
                   area = as.double(area), 
                   outdir = as.integer(rivdir), sillh = as.double(mflac),
                   outnewi = as.double(outnewi), outnewj = as.double(outnewj),
                   basin = as.double(basin), 
                   prcpi = as.double(prcpi), evapi = as.double(evapi),
                   runin = as.double(runin), drainin = as.double(drainin),
                   outelv = double(gridx*gridy), lakem = double(gridx*gridy),
                   lakevolm = double(12*(nyrs+spin)), lakevola = double(nyrs+spin), 
                   voll = double(gridx*gridy), volb = double(gridx*gridy), volr = double(gridx*gridy), 
                   tempdl = double(gridx*gridy), larea = double(gridx*gridy))
  return(simcf)
}
###############################################################################

###############################################################################
## binOutlet
## Function to find outlet cells
## Defined as cells next to the border mask, with no lower elevation neighbors
binOutlet <- function (x) {
  outlet = FALSE
  nna = length(is.na(x))
  if (length(nna) >= 1) { ## Are we next to an edge?
    if (!is.na(x[5])) { ## Is the cell an edge cell?
      if (which.min(x) == 5) {
        outlet = TRUE
      }
    } 
  } 
  return(outlet)
}

## Older version
# binOutlet <- function (x) {
#   outlet = FALSE
#   if (max(x) == 1e6) { ## Are we next to an edge?
#     if (x[5] != 1e6) { ## Is the cell an edge cell?
#       if (which.min(x) == 5) {
#         outlet = TRUE
#       }
#     } 
#   } 
# }
###############################################################################

###############################################################################
## Function to convert radians to degrees
rad2deg <- function(x) {
  x*180/pi
}
###############################################################################

###############################################################################
## Function to convert radians to degrees
deg2rad <- function(x) {
  x/180*pi
}
###############################################################################

###############################################################################
## Function to calculate great circle distances
gcDist <- function(lon1, lat1, lon2, lat2, r=6378) {
  ## Degree conversion
  lon1 <- lon1 * pi/180
  lon2 <- lon2 * pi/180
  lat1 <- lat1 * pi/180
  lat2 <- lat2 * pi/180
  ## Central angle
  ca <- acos((sin(lat1)*sin(lat2)) + 
               (cos(lat1)*cos(lat2) * cos(abs(lon1-lon2))))
  d <- r * ca
  return(d)
}
###############################################################################

###############################################################################
## Function to lower border next to outlets
modBorder <- function(outlet, dem, mask) {
  require(geosphere)
  oID <- Which(outlet==1, cells=TRUE)
  if (length(oID) > 0) {
    for (i in 1:length(oID)) {
      cen.crds = xyFromCell(mask, oID[i])
      ngb.rc = adjacent(mask, oID[i], directions = 8)
      ngb.crds = xyFromCell(mask, ngb.rc[,2])
      ngb.vals = extract(mask, ngb.rc[,2])
      ngb.crds <- ngb.crds[which(ngb.vals==1),]
      ngb.rc <- ngb.rc[which(ngb.vals==1),]
      ngb.dist = distCosine(cen.crds, ngb.crds)
      
      bID <- which.min(ngb.dist)
      dem[ngb.rc[bID,2]] <- dem[ngb.rc[bID,2]]*-1
      mask[ngb.rc[bID,2]] <- mask[ngb.rc[bID,2]]*-1
      
    }
    
  }
  return(list(dem=dem,mask=mask))
}
###############################################################################

###############################################################################
## Function to find pits in DEM
findPit <- function(x) {
  pit = FALSE
  if (which.min(x) == 5) {
    pit = TRUE
  }
  return(pit)
}
###############################################################################

############################################################################
# Matrix manipulation methods
#
# For simplicity we have avoided to create generic functions for 'flip' etc.
# and therefore we have to call the corresponding methods coupled to the
# 'matrix' class explicitly, i.e. flip.matrix().
############################################################################
# Flip matrix (upside-down)
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

# Rotate matrix 90 clockworks
rotate90.matrix <- function(x) {
  t(mirror.matrix(x))
}

# Rotate matrix 180 clockworks
rotate180.matrix <- function(x) { 
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}

# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

###############################################################################
## Older version based on original fortran code
## Calls rhydra fortran routine
# dyn.load("./src/rhydra.so")
# 
# rhydra <- function(nyrs, startyear, converg = 1, laket = 0, spin = 1,
#                    normal = 1, leap = 1, irrig = 1,
#                    outnewi, outnewj, basin, dem, rivdir, mflac,
#                    prcpi, evapi, runin, drainin,
#                    gridxf, gridyf) {
#   
#   simcf = .Fortran("rhydra",
#                    nyrs = as.integer(nyrs),
#                    startyear = as.integer(startyear), 
#                    converg = as.integer(converg), 
#                    laket = as.integer(laket), 
#                    spin = as.integer(spin),
#                    normal = as.integer(normal), 
#                    leap = as.integer(leap), 
#                    irrig = as.integer(irrig),
#                    outnewi = as.double(outnewi), outnewj = as.double(outnewj),
#                    basin = as.double(basin), dem = as.double(dem),
#                    outdir = as.double(rivdir), sillh = as.double(mflac),
#                    prcpi = as.double(prec), evapi = as.double(evap),
#                    runin = as.double(runoff), drainin = as.double(drain),
#                    outelv = double(gridxf*gridyf), lakem = double(gridxf*gridyf),
#                    lakevolm = double(12*(nyrs+spin)), lakevola = double(nyrs+spin))
#   return(simcf)
# }
###############################################################################

