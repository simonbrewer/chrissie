!-------------------------------------------------------------------------------
! RHYDRA: Version of hydra hydrological model to be called from R
!
! Coe, M.T. 2000: Modeling terrestrial hydrological systems at the continental 
! scale: Testing the accuracy of an atmospheric GCM, 
! Journal of Climate, 13, 686-704.
!
! Ver. 0.1 Simple template to read in a DEM
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
      subroutine rhydra( nc, nr, nyrs, ndays, startyear,
     *                   res, converg, laket, spin,
     *                   dem, mask, area, outdir, sillh,
     *                   outnewi, outnewj, basin, 
     *                   prcpi, evapi, runin, drainin,
     *                   outelv, lakem, lakevolm, lakevola,
     *                   voll, volb, volr, tempdl, larea )

      ! Input variables
      integer nc,nr ! Grid dimensions
      ! Forcings
      integer mask(nc,nr) ! Binary mask
      integer nyrs,ndays
      integer res ! Grid resolution in arc seconds
      double precision dem(nc,nr),area(nc,nr),
     *     outdir(nc,nr),sillh(nc,nr),
     *     outnewi(nc,nr),outnewj(nc,nr),basin(nc,nr)
      ! Topographic variables
      double precision drainin(nc,nr,ndays),runin(nc,nr,ndays),
     *     prcpi(nc,nr,ndays),evapi(nc,nr,nmons) 

      ! Outputs
      double precision larea(nc,nr) ! lake area of cell (0/1)
      double precision laream(nc,nr) ! monthly average area
      double precision outelv(nc,nr) ! water surface elevation (m)
      double precision sflux(nc,nr) ! surface flux (?)
      double precision elevm(nc,nr),lakem(nc,nr),deptm(nc,nr) ! (?) 
      double precision voll(nc,nr) ! lake reservoir (m3)
      double precision volb(nc,nr) ! baseflow reservoir (m3)
      double precision volr(nc,nr) ! runoff reservoir (m3)
      double precision volt(nc,nr) ! total volume (m3) (?)
      double precision dvoll(nc,nr) ! volume change (m3)
      double precision tempdl(nc,nr),tempdr(nc,nr),temp(nc,nr)
      double precision fluxout(nc,nr),sfluxin(nc,nr)
      double precision areat(nc,nr),basin2(nc,nr)
      double precision sfluxout(nc,nr,12)
      
      ! Model internals
      double precision real circ,dy,dx,pi,rad,phi,delt,ic,io,ioo
      double precision grideps,dveps,gridif
      double precision timer,timed,timeg
      integer i,j,k,ii,jj,kk,k2,tmpdir,tmpdir2
      integer spin
      integer ioff(8),joff(8)
      integer ndaypm(12)

c--------------------------------------------------------------
      ! Parameters
      parameter (grideps = 1.e-10,dveps = 1.e-10)
c     parameter (grideps = 0.01,dveps = 1.e-06)
      parameter (circ=4.0024E+7, pi = 3.1415926536)
      !parameter (rad = pi/180., dy = circ*((5./60.)/360.))
      parameter (rad = pi/180.) !, dy = circ*((res/(60.*60.))/360.))
      parameter (dgper = 1/360)  !number of degrees/input grid cell

c--------------------------------------------------------------
c set values for searching all gridcells around a given gridcell
c 8 possible directions !! NEEDS UPDATING TO GRASS DIRECTIONS
c
      !data ioff /0,1,1,1,0,-1,-1,-1/ ! traditional way
      !data joff /-1,-1,0,1,1,1,0,-1/
      data ioff /-1,-1,-1,0,+1,+1,+1,0/ ! GRASS r.watershed
      data joff /+1,0,-1,-1,-1,0,+1,+1/ ! Offsets for movement

c--------------------------------------------------------------
c initialize variables and constants
c nyr is the number of years to run the model, delt is the timestep.
c
      ioo  = 2.5e-03 
      dy = circ*((res/(60.*60.))/360.)
c     delt = 1800. !shorter timestep for poleward locations
      delt = 3600. !good for tropics and mid-latitudes
      timed = 13.0e+05   !when seperating baseflow from surface runoff
      timer = 7200.    
c     timed = timer     !for IBIS runs, which already calculates residence time
      effref = 0.8      !average of channel and floodplain flow
      ndaypm = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      iday = 0
      imon = 1
      iyear = 1
      itime = 0
      ktstep = 0
      nspday = 24./(delt/3600.)
      volchk = 0.
      !write(*,*) "nspday",nspday
c--------------------------------------------------------------
c initialize output variables
c
      do 205 j = 1,nr
       do 206 i = 1,nc
        outelv(i,j) = dem(i,j)
        elevm(i,j) = 0.
        sflux(i,j) = 0.
        lakem(i,j) = 0.
        deptm(i,j) = 0.
c       vollm(i,j) = 0.
c       dvollm(i,j) = 0.
        voll(i,j) = 0.
        volb(i,j) = 0.
        volr(i,j) = 0.
        volt(i,j) = 0.
        dvoll(i,j) = 0.
        tempdl(i,j) = 0.
        tempdr(i,j) = 0.
        temp(i,j)  = 0.
        fluxout(i,j) = 0.
        sfluxin(i,j) = 0.
        areat(i,j) = 0.
        basin2(i,j) = 0.
        do 207 k = 1,12
         sfluxout(i,j,k) = 0.
207     continue
206    continue
205   continue
c
      do 210 j = 1,nr
       do 220 i = 1,nc
        if(laket .eq. 0)then
         larea(i,j) = 0.
        else
         larea(i,j) = min(larea(i,j),1.)
        endif
        if(sillh(i,j) .eq. 0.)then
         outnewi(i,j) = 0.
         outnewj(i,j) = 0.
        endif
        if(outnewi(i,j) .eq. 0.)sillh(i,j) = 0.
 220   continue
 210  continue
c
c
c--------------------------------------------------------------
c convert input from mm/day to m/s
c
      do 334 k = 1,ndays
       do 335 j = 1,nr
        do 336 i = 1,nc
c
         if (mask(i,j).eq.1) then
         if(runin(i,j,k) .ge. 1.e+20) then
          runin(i,j,k) = 0.
         endif
         if(drainin(i,j,k) .ge. 1.e+20) then 
          drainin(i,j,k) = 0.
         endif
         if(prcpi(i,j,k) .ge. 1.e+20) then 
          prcpi(i,j,k) = 0.
         endif
         if(evapi(i,j,k) .ge. 1.e+20) then 
          evapi(i,j,k) = 0.
         endif
         !if(i.eq.10.and.j.eq.12) then
         !        write(*,*) i,j,prcpi(i,j,k)
         !end if
c
         runin(i,j,k)   = max((runin(i,j,k))/0.864e+8,0.) !mm/d>m3/s
         drainin(i,j,k) = max((drainin(i,j,k))/0.864e+8,0.) !mm/d>m3/s
         prcpi(i,j,k)   = max((prcpi(i,j,k))/0.864e+8,0.) !mm/d>m3/s
         evapi(i,j,k)   = max((evapi(i,j,k))/0.864e+8,0.) !mm/d>m3/s
c
c Test to see response of open lake
c
c        evapi(i,j,k)   = 0.
c        runin(i,j,k)   = runin(i,j,k)*100.
         end if ! Mask check
c
 336    continue
 335   continue
 334  continue
cc
c--------------------------------------------------------------
c Calculate the total volume within each potential lake area (volt).
c First set the volume of the outlet location to initial value.
c volt is used to determine when the lake is full and outflow
c will occur to the river downstream of the sill.
c
      do 914 j = 1,nr
       do 915 i = 1,nc
c
        if (mask(i,j).eq.1) then
        if(sillh(i,j) .gt. 0.)then
         ii = outnewi(i,j)
         jj = outnewj(i,j)
         if ((ii .gt. 0).and. (ii .le. nc))then
          if((jj .gt. 0).and. (jj .le. nr))then
           volt(ii,jj) = volt(ii,jj) + 
     *        max(sillh(i,j)-dem(i,j),0.1)*area(i,j)
          endif
         endif
        endif
        end if ! Mask check
c
915    continue
914   continue
c
c--------------------------------------------------------------
c CONVERGENCE - set converg to 0. in .inf file
c
c  This convergence function is designed to approximately fill all
c potential lake basins in the region of simulation. It is useful
c in humid regions where all lake basins are full and you do not
c wish to wait for the model to fill them. It can take many hundreds
c of years to fill large lakes (such as Lake Michiga/Huron). Therefore,
c this saves time. It is only approximate, therefore, don't depend on the
c answer to be exactly right until the model has run for some time.
c 
c  if converg = 0. then set the outlet volume to the lake
c  maximum volume. This allows the model to start out with
c  all lakes full and discharge to the ocean without first 
c  having to fill large lake basins (Great Lakes).
c  This is only used if are interested in river discharge only.
c  If want to predict lakes must have Converg = 1
c
      if(converg .eq. 0)then
      write(*,*)
      write(*,*)'###################'
      write(*,*)'converge on'
      write(*,*)'###################'
c
      do 918 j = 1,nr
       do 919 i = 1,nc
         if (mask(i,j).eq.1) then
         ii = outnewi(i,j)
         jj = outnewj(i,j)
         !iii = min(max(outnewi(i,j)-(istart-1),0.),REAL(incf))
         !jjj = min(max(outnewj(i,j)-(jstart-1),0.),REAL(inrf))
         if(outnewi(i,j) .gt. 0.)then
           voll(ii,jj)  = max(max(dem(i,j)-sillh(i,j),0.1)*area(i,j),
     *     volt(ii,jj))
           larea(i,j)  = 1.
           outelv(i,j) = max(dem(i,j),sillh(i,j))  !10/7/98
           !areat(iii,jjj) = areat(iii,jjj) + area(j)   !10/7/98
         endif
         endif
919    continue
918   continue
c
      endif
c
c--------------------------------------------------------------
c FILL LOCATION
c This is used if you would like to specify the location that
c a lake will start filling from. This is useful for closed 
c basin lakes that have multiple sub-basins. If nothing is
c is changed here than the lake will begin to fill at the low
c point nearest to the sill location.
c
c Find the lake kernal location, the location that a lake will
c start to fill from. For now the kernal will be the first pit
c encountered. Later it can be modified to provide a more 
c physically realistic kernal location.
c
c set basin2 as a mask for lake kernal location. Currently am
c setting as the outlet location. Corrections to that will
c be made outside this loop
c
      do 916 j = i,nr
       do 917 i = 1,nc
        if(outnewi(i,j) .gt. 0.)then
          write(*,*) i,j,outnewi(i,j),outnewj(i,j)
          !! THIS ISN'T BEING TRIGGERED
         if((i .eq. outnewi(i,j)) .and. (j .eq. outnewj(i,j)))then
          basin2(i,j) = 1.
         endif
        endif
          !write(*,*) i,j,basin2(i,j)
917    continue
916   continue
c
c
c--------------------------------------------------------------
c
c Main loop in two parts 120 and 121 loops. 
c First loop (120) calculates the change in volume with time.
c Second loop (121) spreads volume in lake area
c
c     goto 441
c
      write(*,*)'to main loop'
c
c loop 130 is the total number of years that the model will be run.
c 
c This version is set for equilibrium runs so the number of years will
c be simply spinup plus 1 with repeated climate
c 
c Recommended spin up is between 30 and several 100 years although
c (To avoid long spin up, use the convergence function)
c
      do 130 iyear = 1,spin+1
c
c monthly loop
c
       k = 1 ! Doy counter
       do 131 imon = 1,12
c
c
 131   continue   ! end imonth loop
c
       write(*,*)'year = ',iyear
c
 130  continue    ! end iyear loop
c

c
!-------------------------------------------------------------------------------
      end
!-------------------------------------------------------------------------------

