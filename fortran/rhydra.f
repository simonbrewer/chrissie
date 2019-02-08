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
      integer outdir(nc,nr) ! Drainage direction
      integer nyrs,ndays
      integer res ! Grid resolution in arc seconds
      double precision dem(nc,nr),area(nc,nr),
     *     sillh(nc,nr),
     *     outnewi(nc,nr),outnewj(nc,nr),basin(nc,nr)
      ! Topographic variables
      double precision drainin(nc,nr,ndays),runin(nc,nr,ndays),
     *     prcpi(nc,nr,ndays),evapi(nc,nr,ndays) !! THIS LOOKS WRONG !! 

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
      integer edgecell
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
c
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
         !! THIS ISN'T BEING TRIGGERED
         if((i .eq. outnewi(i,j)) .and. (j .eq. outnewj(i,j)))then
          !write(*,*) i,j,outnewi(i,j),outnewj(i,j)
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
c re-initialize some variables each month.
c
       do j = 1,nr
        do i = 1,nc
         sfluxout(i,j,imon) = 0.
         laream(i,j) = 0.
        enddo
       enddo
c
c write a few diagnostics
c
       write(*,*)'begin mon ',imon
c
c start the daily loop
c
       do 132 iday = 1,ndaypm(imon)
c
c start the sub-daily timestep (one hour in this case but depends
c on how delt is set).
c
       do 133 kt = 1,int(nspday)
        ktstep = ktstep + 1 !NOT USED?
c
c start first spatial loop.
c
       do 120 j = 1,nr
        do 110 i = 1,nc
c
        if(mask(i,j) .eq. 1)then
        !write(*,*) i,j,kt,iday
c----------------------------------------------------
c Initialize hourly climate input variables.
c
          rin     = max(runin(i,j,k)*area(i,j),0.)  !surface runoff
          bin     = max(drainin(i,j,k)*area(i,j),0.) !sub-surface runoff
          bout    = 0. !sub-surface flow to river
          rout    = 0. !surface flow to river
          irrout  = 0. !irrigation rate (prescribed below)
          prcpl   = max(prcpi(i,j,k)*area(i,j),0.)  !precipitation rate
          evapl   = max(evapi(i,j,k)*area(i,j),0.)  !evaporation rate
c
c         if(i.eq.25.and.j.eq.25.and.kt.eq.1) then
c                 write(*,*) i,j,k,iday,runin(i,j,k),rin
c         end if
c --------------------------------------------------------------------
c IRRIGATION - irrigation withdrawals to go here
c --------------------------------------------------------------------
c
c --------------------------------------------------------------------
c RESERVOIRS         
c --------------------------------------------------------------------
c calculate volume in runoff reservoir for land or lake
c
         rout = volr(i,j)/timer
         volr(i,j) = max(volr(i,j) + (rin-rout)*delt,0.)
c         if(i.eq.25.and.j.eq.25.and.kt.eq.1) then
c                 write(*,*) i,j,k,kt,rin,rout,volr(i,j)
c         end if
c
c calculate volume in baseflow reservoir, land only
c
c        bout = 0.75*volb(ii,jj)/timed + 0.25*volb(ii,jj)/timeg
         bout = volb(i,j)/timed
         volb(i,j) = max(volb(i,j) + (bin-bout)*delt,0.)
c
c calculate volume in transport reservoir
c needs outlet locations (ii,jj)
c
         ii = outnewi(i,j)
         jj = outnewj(i,j)
c lake water balance for this timestep
c
         temp(i,j) = (rout+bout)*(1.-larea(i,j))
c
c land water balance for this timestep.
c
         tempdr(i,j) = ((rout+
     *           bout)*(1.-larea(i,j))-irrout
     *          + (sfluxin(i,j) - fluxout(i,j)))*delt
c
c subtract any evaporation from the lake from the outlet
c location. The outlet is the accountant for the entire lake.
c
         if(larea(i,j) .gt. 0.)then
          tempdl(ii,jj) = tempdl(ii,jj) 
     *          + ((prcpl-evapl)*larea(i,j))*delt
         !if(i.eq.25.and.j.eq.25.and.kt.eq.1) then
         !        write(*,*) "a",i,j,k,kt,tempdl(ii,jj),larea(i,j)
         !end if
         endif
c
c This should be accounted for in mask check
!         else                 !basin .ne. requested number
!          voll(ii,jj) = 0.
!         endif
c
        end if ! Mask check
c
 110    continue
 120   continue
c 
c----------------------------------------------------------------
c Calculate the change in volume (dvoll) as the sum of the P-E (tempdl)
c and the river flow (tempdr). Set the minimum value of dvoll.
c Add dvoll to existing reservoir volume (voll). Set variables
c to 0. for next timestep.
c
       do 123 j = 1,nr
        do 113 i = 1,nc
        if(mask(i,j) .eq. 1)then
         dvoll(i,j)  = tempdr(i,j) + tempdl(i,j)
         if(abs(dvoll(i,j))/delt/area(i,j) .lt. dveps) then
          dvoll(i,j) = 0.
         endif
         !if(i.eq.25.and.j.eq.25.and.kt.eq.1) then
         !        write(*,*) i,j,k,kt,dvoll(i,j)
         !end if
         voll(i,j)  = max(voll(i,j) + dvoll(i,j),0.)
         tempdl(i,j) = 0.
         tempdr(i,j) = 0.
        end if ! Mask check
113     continue
123    continue
c
c----------------------------------------------------------------
c Distribute the dvoll of this timestep roughly into the existing lake area
c It will be evened out in loop 121
c
       do 122 j = 1,nr
        do 112 i = 1,nc
        if(mask(i,j) .eq. 1)then
c
         !ii = i - (istart-1)
         !jj = j - (jstart-1)
         i2 = outnewi(i,j)
         j2 = outnewj(i,j)
c
         if(((i2 .gt. 0).and. (j2 .gt. 0)) !IF1 
     *          .and. (laket .eq. 0))then
c                             write(*,*) "IF1"
c
c if volume in basin > lake volume larea = 1. everywhere
c outelv = sill height
c
         if(voll(i2,j2) .ge. volt(i2,j2))then !IF2
                             !write(*,*) "IF2a"
           outelv(i,j) = max(outelv(i,j) + larea(i,j)*
     *            dvoll(i2,j2)/areat(i2,j2),dem(i,j))  !0.0 if larea = 0.
c           outelv(i,j) = sillh(i,j)  !simpler way of handling it
           larea(i,j) = max(min(outelv(i,j)-dem(i,j),1.),0.)
c
c if there is no volume then larea = 0.
c
         elseif(voll(i2,j2) .eq. 0.)then 
                             !write(*,*) "IF2a"
           larea(i,j)  = 0.
           outelv(i,j) = dem(i,j)
           areat(i2,j2) = 0.   !probably not necessary
c
         elseif((voll(i2,j2).gt.0.).and.
     *          (voll(i2,j2).lt.volt(i2,j2)))then
                             !write(*,*) "IF2a"
c
c if some lake already exists in closed basin distribute dvoll
c evenly to those existing cells 
c
         if(areat(i2,j2) .gt. 0.)then !IF3
                             !write(*,*) "IF3"
c
c set outelv if larea > 0. Add depth of water if positive
c or negative. 
c
           outelv(i,j) = max(outelv(i,j) + larea(i,j)*
     *            dvoll(i2,j2)/areat(i2,j2),dem(i,j))  !0.0 if larea = 0.
c
c If no existing lake; set larea of outlet location = 1.
c Ideally would like to choose a kernal location in a 
c realistic location within a lake. This would be stored
c in array basin2.
c
           else   !if(areat(i2,j2) .eq. 0.))then 
             if(basin2(i,j) .eq. 1.)then !?? basin2 !IF4
                             !write(*,*) "IF4"
               outelv(i,j) = max(outelv(i,j) +
     *             dvoll(i,j)/area(i2,j2),dem(i,j))
               larea(i,j) = max(min(outelv(i,j)-dem(i,j),1.),0.)
               areat(i2,j2) = max(area(i,j)*larea(i,j),0.)
c
             endif !IF4
           endif !IF3
c
         endif !IF2
       endif !IF1
c
       else                 ! masked cell
         voll(i,j) = 0.
       endif
c
       fluxout(i,j) = 0.
       sfluxin(i,j) = 0.
c
 112   continue
 122   continue
c
c--------------------------------------------------------------
c In this loop calculate the flux out of the cell, to which direction,
c and the sum of the fluxes into cells. Also calculate the
c water depth and lake area (larea) for each cell. 
c
       do 121 j = 1,nr
        do 111 i = 1,nc
c
        if(mask(i,j) .eq. 1) then
c
         !ii = i - (istart-1)
         !jj = j - (jstart-1)
         !iii = min(max(outnewi(i,j)-(istart-1),1.),REAL(incf))
         !jjj = min(max(outnewj(i,j)-(jstart-1),1.),REAL(inrf))
         ii = outnewi(i,j) ! outlet...
         jj = outnewj(i,j)
         !write(*,*) i,j,ii,jj,"here1",mask(i,j),outdir(i,j)
c
c use river directions computed in getpwd.f
c
c GRASS directions use negatives for edge/coast pour points
c Should track these for total outflow
         if (outdir(i,j).le.0) then
          edgecell = 1
         else 
          edgecell = 0
         end if
         tmpdir = abs( outdir(i,j) )
c
         j2 = j + joff(tmpdir)
         i2 = i + ioff(tmpdir)
         dx = (area(i,j)/dy+area(i2,j2)/dy)/2.
         dist = sqrt((dx*dx*abs(i2-i)*abs(i2-i))
     *        + (dy*dy*abs(j2-j)*abs(j2-j)))
         dist = max(dx,dist)
c
c set effective velocity of each cell. It is dependent on the
c lake volume of the cell. For large lakes the value is quite
c slow, for small lakes it approaches the reference value,
c effref. For non-lake points it is a function of the
c gradient as in Miller et al. 1994. This is fairly rough
c and could be improved with considerations of sinuosity or
c stream order for example.
c
c Invoke the top half of this if statement if you want the volume
c of the lake to impact the river discharge velocity. It is unique
c to each lake and should be studied before use
c
c       if((larea(i,j) .gt. 0.) .and.
c    *     (voll(ii,jj) .gt. 0.))then
c        vollocal = 1.*area(j)
c        volref = sqrt(vollocal/volt(iii,jjj))
c        effvel = min(effref,0.1*effref*volref)
c        effvel = min(effref,0.08*effref*volref)
c        effvel = max(effvel,1.0e-02)
c       else                         !non-lake
         ic = max(dem(i,j)-dem(i2,j2),1.)/dist
         io = ioo/(dx/dy)   !scale reference gradient to latitude
         effvel = effref*sqrt(ic/io)
         effvel = min(effvel,3.0)
c       endif
c calculate fluxout of each cell and send it as fluxin to
c either the cell downstream if it is not a lake downstream
c or to the outlet location. The fluxout of the cell which
c corresponds to the sill is calculated for only that water 
c volume in excess of the volume required to fill the lake.
c
c       effvel = 0.5 !alternatively could set velocity to a  constant
c
         fluxout(i,j) = max((voll(i,j)-volt(i,j))*
     *                  (effvel/dist),0.)
c
         fluxout(i,j) = max(min(fluxout(i,j),
     *           sfluxin(i,j) + temp(i,j) +
     *           ((voll(i,j)-volt(i,j))/(delt*2.))),0.)
c
c Truncate fluxout if too small for computation. 
c
         if(fluxout(i,j)/area(i,j) .lt. dveps) then 
          fluxout(i,j) = 0. 
         endif
         !i3 = i2-(istart-1)
         !j3 = j2-(jstart-1)
         if((i2.gt.0).and.(i2.le.nc).and.
     *      (j2.gt.0).and.(j2.le.nr)) then
          sfluxin(i2,j2) =
     *    sfluxin(i2,j2) + fluxout(i,j)
         endif
c
c--------------------------------------------------------------
c Adjust height of water column in each cell using the cellular
c automata. Distribute water height within a lake basin only
c This flattens the lake surface so that there are no hills or
c valleys due to differences in the local water budget.
c
        if((laket .eq. 0) .and. ((ii .ne. 0).and. (jj.ne.0)))then
c
         if((outelv(i,j) .gt. dem(i,j)) .and.
     *     (voll(ii,jj) .lt.
     *      volt(ii,jj)))then
c
         if((voll(ii,jj) .ge. 0.).and.
     *      (sillh(i,j) .gt. 0.))then
c
          tmpdir2 = 0
          tmph = outelv(i,j)
!!!! NEED TO FIGURE OUT THESE I3/J3 COORDINATES         
!!!! I THINK THIS IS THE COORDINATES OF THE 8 NEIGHBOR CELLS
           do k2 = 1,8
            !j3 = (j + joff(k2))-(jstart-1)
            !i3 = (i + ioff(k2))-(istart-1)
            j2 = j + joff(k2)
            i2 = i + ioff(k2)

            if((outelv(i2,j2) .lt. outelv(i,j)) .and.
     *        (sillh(i2,j2) .eq. sillh(i,j)))then
             if(outelv(i2,j2) .lt. tmph)then
              tmpdir2 = k2
              tmph = outelv(i2,j2)
             endif
            endif
           enddo
c
          if(tmpdir2 .ne. 0.)then
c
           !j3 = (j + joff(tmpdir2))-(jstart-1)
           !i3 = (i + ioff(tmpdir2))-(istart-1)
           j2 = j + joff(tmpdir2)
           i2 = i + ioff(tmpdir2)
c
           gridif = min(max((outelv(i,j) -
     *               outelv(i2,j2)),0.),
     *               max((outelv(i,j) - dem(i,j)),0.))
           if(abs(gridif) .lt. grideps) then
            gridif = 0.
           endif
           gridif = gridif*area(i,j)*0.5  !move only 0.5 of difference
c
           outelv(i2,j2) = outelv(i2,j2) +
     *           gridif/area(i2,j2)
           outelv(i,j) = outelv(i,j) -
     *           gridif/area(i,j)
c
          endif
c
         else
          larea(i,j) = 0.
          outelv(i,j) = dem(i,j)
          areat(ii,jj) = 0.
         endif                !voll(i2,j2) .ge. 0.
         endif                !outelv .gt. dem
c        endif
c
c Set lake area = either; 1 if depth is greater than 1 meter, to
c a % of the lake cell equal to the depth of water, or to 0. if
c depth = 0.  depth = outelv-grid
c Adjust the total lake area (areat) to represent current larea.
c Do this by first subtracting current larea from total then
c adding new larea to total. If larea = 0. area added = 0.
c
         if((outnewi(i,j) .gt. 0.) .and.
     *     (voll(ii,jj) .lt. volt(ii,jj)))then
c
          if((voll(ii,jj) .gt. 0.) .and.
     *       (sillh(i,j) .gt. 0.))then
c
          ! These are just the outlets 
!           ix = min(max(int(outnewi(i,j)-(istart-1)),1),incf)
!           jx = min(max(int(outnewj(i,j)-(jstart-1)),1),inrf)
           areat(ii,jj) = max(areat(ii,jj) - area(i,j)*
     *         larea(i,j),0.)
           larea(i,j) = max(min(outelv(i,j)-
     *         dem(i,j),1.),0.)
           areat(ii,jj) = max(areat(ii,jj) + area(i,j)*
     *         larea(i,j),0.) !previous
          else
           larea(i,j) = 0.
           outelv(i,j) = dem(i,j)
           areat(ii,jj) = 0.
          endif
         endif
c
       laream(i,j) = (laream(i,j)+((outelv(i,j)-dem(i,j))
     *                        /(ndaypm(imon)*nspday)))
       if(outelv(i,j)-dem(i,j) .lt. 0.1)then
        laream(i,j) = 0.
       endif
c
c     endif        !iday
      endif        !laket = 0
c
       sfluxout(i,j,imon) = (sfluxout(i,j,imon)+fluxout(i,j)
     *                        /(ndaypm(imon)*nspday))
c
c
        endif              !end mask loop
c
 111    continue
 121   continue          !end fluxout loop
c
c
 133   continue   !end hourly loop
       k = k + 1
 132   continue   !end daily loop
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

