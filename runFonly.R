## Only runs the fprtran code
source("helpers.R")
nyrs = 1
sim.out = rhydra(gridx, gridy, nyrs, ndays, startyear=1, 
                 res=30, converg=1, laket=0, spin=1,
                 dem=dem, mask=mask, area=cella, rivdir=ldd, mflac=mflac,
                 outnewi=iout, outnewj=jout, basin=bas, 
                 prcpi=pre, evapi=evp, runin=sro, drainin=bro) 
