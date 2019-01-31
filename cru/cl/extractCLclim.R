## Extract cru data for the region
## THis is the first step, then followed by calculating the climate indices
## Useful to make this larger than the regional extent

library(raster)

crudir = "~/Dropbox/Data/climate/cru_cl_2.00/"
filelist = list.files(crudir, pattern=".nc")

r = raster(paste0(crudir,"cru_10min_elv.nc"))

## CL extent
cl.r = raster("~/Dropbox/Data/devtools/chrissie/dem/chrissie_dem.nc")
# myext = extent(c(10,40,-35,-10))
# myext = extent(c(19,27.5,-33,-27.5))
myext = extent(cl.r)
myext = extent(c(29,34,-28,-25))

for (i in 1:length(filelist)) {
  varname = unlist(strsplit(unlist(strsplit(filelist[i], "_"))[3],"\\."))[1]
  outfile = paste0("cl_cru_10min_",varname,".nc")
  tmp.stk = stack(paste0(crudir,filelist[i]))
  tmp.stk.2 = crop(tmp.stk, myext)
  
  writeRaster(tmp.stk.2, filename = outfile, varname=varname, overwrite=TRUE)
}
