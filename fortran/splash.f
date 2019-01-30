c  !////////////////////////////////////////////////////////////////
c  ! SPLASH
c  ! This module contains all functions and parameter values to run
c  ! the SPLASH model, given inputs
c  ! - temperature (deg C, daily values)
c  ! - precipitation (mm/day, daily values)
c  ! - sunshine fraction (unitless, daily values)
c  ! - latitude (deg N)
c  ! - elevation (m.a.s.l.)
c  !
c  ! VERSION: 1.0-r1
c  ! LAST UPDATED: 2016-05-27
c  !
c  ! Copyright (C) 2016 Prentice Lab
c  !
c  ! SPLASH is free software: you can redistribute it and/or modify it under
c  ! the terms of the GNU Lesser General Public License as published by
c  ! the Free Software Foundation, either version 2.1 of the License, or
c  ! (at your option) any later version.
c  !
c  ! SPLASH is distributed in the hope that it will be useful,
c  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
c  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  ! GNU Lesser General Public License for more details.
c  !
c  ! You should have received a copy of the GNU Lesser General Public License
c  ! along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
c  !
c  ! Citation:
c  ! T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
c  ! Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
c  ! led algorithms for simulating habitats (SPLASH): Robust indices of radiation,
c  ! evapotranspiration and plant-available moisture, Geoscientific Model
c  ! Development, 2016 (in progress)
c  !----------------------------------------------------------------
c  !----------------------------------------------------------------
c  ! Version modified to run from R by Simon Brewer 160614
c  !
c  ! Note that the Julian days have been replaced by a standard 365 
c  ! day year
c  !      
c  ! To Do
c  ! Soil data - whc
c  ! Elevation as non-constant
c  ! Snow pack
c  ! Irrigation amount (or other water input)      
c  !----------------------------------------------------------------

c******************************************************************************c
        subroutine spin_up_sm( yr, lat, elv, pr, tc, sf,
     >                         maet,mpet,mcn,mro,msm,
     >                         daet,dpet,dcn,dro,dsm,sm,
     >                         ddl, dsl )
                
        !----------------------------------------------------------------
        ! Spins up the daily soil moisture
        !----------------------------------------------------------------
        ! arguments
        integer yr    ! year AD
        double precision lat   ! latitude (degrees)
        double precision elv   ! altitude (m)
        double precision pr(365)   ! monthly or daily precip (mm)
        double precision tc(365)   ! mean monthly or daily temperature (deg C)
        double precision sf(365)   ! mean monthly or daily sunshine fraction (unitless)

        ! local variables
        integer spin_count         ! counter variable
        double precision start_sm  ! soil moisture in first day of year
        double precision end_sm    ! soil moisture in last day of year
        double precision diff_sm   ! difference in soil moisture between first and last day of year
        double precision ro    ! daily runoff (mm)
        double precision sm    ! soil moisture (=water content, mm)

        ! output variables
        double precision daet(365)
        double precision dpet(365)
        double precision dcn(365)
        double precision dro(365)
        double precision dsm(365)
        double precision ddl(365)
        double precision dsl(365)

        double precision maet(12)
        double precision mpet(12)
        double precision mcn(12)
        double precision mro(12)
        double precision msm(12)

         ! control variables
        spin_count = 1
        diff_sm = 9999.0

        ! initialize soil moisture and runoff
        sm = 0.0
        ro = 0.0

        do 10 while (diff_sm.gt.1.0e-4)

        ! Run one year
        !print*,'********************************************'
        !print*,'spinup year', spin_count
        !print*,'--------------------------------------------'
        start_sm = sm
        call run_one_year( lat, elv, yr, pr, tc, sf,
     >                     maet,mpet,mcn,mro,msm,
     >                     daet,dpet,dcn,dro,dsm,sm,
     >                     ddl,dsl )

        end_sm  = sm

        ! Check to see if 1 Jan soil moisture matches 31 Dec:
        diff_sm  = abs(end_sm - start_sm)

        !print*,'soil moisture in equil. when diff < 0.0001'
        !print*,'start ',start_sm
        !print*,'end ',end_sm
        !print*,'diff ',diff_sm

        ! Increase counter
        !print*,spin_count
        spin_count = spin_count + 1

 10     continue

        end ! spin_up_sm
c******************************************************************************c

c******************************************************************************c
        subroutine run_one_year( lat, elv, yr, pr, tc, sf,
     >                            maet,mpet,mcn,mro,msm,
     >                            daet,dpet,dcn,dro,dsm,sm,
     >                            ddl,dsl )
        implicit none
        !----------------------------------------------------------------
        ! Calculates daily and monthly quantities for one year
        !----------------------------------------------------------------
        ! Rewritten to only use 365 daily input data SB 160616
        !----------------------------------------------------------------
        ! arguments
        double precision lat   ! latitude (degrees)
        double precision elv   ! altitude (m)
        integer yr    ! year AD
        double precision pr(365)   ! monthly or daily precip (mm)
        double precision tc(365)   ! mean monthly or daily temperature (deg C)
        double precision sf(365)   ! mean monthly or daily sunshine fraction (unitless)
        double precision aet   ! daily actual evapotranspiration (mm)
        double precision pet   ! daily potential evapotranspiration (mm)
        double precision cn    ! daily condensation (mm)
        double precision ro    ! daily runoff (mm)
        double precision sm    ! soil moisture (=water content, mm)
        double precision dayl       ! day length (hrs)
        double precision sunl       ! sunlight length (hrs)

        ! local variables
        double precision use_tc ! mean monthly or daily temperature (deg C)
        double precision use_sf ! mean monthly or daily sunshine fraction (unitless)
        double precision use_pr ! monthly or daily precipitation (mm)
        double precision sw     ! evaporative supply rate (mm/h)

        integer doy             ! day of year
        integer moy,nmonth      ! month of year
        integer ndaymonth(12),jd1,jd2         ! number of days in month (formerly 'nm')
        integer idx                       ! day of year corresponding to yesterday
        integer dm                        ! day of month

        data (ndaymonth(moy),moy=1,12) 
     *   / 31,28,31,30,31,30,31,31,30,31,30,31          /
 

        ! output variables
        double precision daet(365)
        double precision dpet(365)
        double precision dcn(365)
        double precision dro(365)
        double precision dsm(365)
        double precision ddl(365)
        double precision dsl(365)

        double precision maet(12)
        double precision mpet(12)
        double precision mcn(12)
        double precision mro(12)
        double precision msm(12)

        ! Parameters
        parameter(nmonth = 12)           ! number of months in year
        
        ! Reset monthly totals (output variables)
        ! call initmonthly

        ! initialize DOY
        doy = 0

        ! Iterate through months
        do 10 moy=1,nmonth
        ! initialize monthly output variables
        maet(moy) = 0
        mpet(moy) = 0
        mcn(moy) = 0
        mro(moy) = 0
        msm(moy) = 0


        ! Calculate number of days in month
        ! Replaced with simple 365 day cycle using month lengths
        ! defined above: SB 160616
        !call get_julian_day(yr, moy, 1, jd1)
        !call get_julian_day(yr, moy+1, 1, jd2)
        !ndaymonth = jd2 - jd1
        ! Iterate through days in this month
        do 20 dm=1,ndaymonth(moy)

        ! Calculate the day of the year
        ! xxx is this '+1' really necessary here?
        !call get_julian_day(yr, 1, 1, jd1)
        !call get_julian_day(yr, moy, dm, jd2)
        !doy = (jd2 - jd1) + 1
        doy = doy + 1
        ! write(*,*) "DOY:",doy

        ! Select daily data
        use_tc   = tc(doy)
        use_sf   = sf(doy)
        use_pr   = pr(doy)

        call run_one_day( lat, elv, doy, yr, use_pr, use_tc, use_sf,
     >                    aet, pet, sm, cn, ro, dayl, sunl )

        !print*,'currentsm ',doy,sm
         ! Collect daily output variables
        ! call getout_daily( doy, moy )
        !write(*,*) "aet:",aet
        daet(doy) = aet
        dpet(doy) = pet
        dcn(doy) = cn
        dro(doy) = ro
        dsm(doy) = sm
        ddl(doy) = dayl
        dsl(doy) = sunl

        ! Collect monthly output variables
        maet(moy) = maet(moy) + aet
        mpet(moy) = mpet(moy) + pet
        mcn(moy) = mcn(moy) + cn
        mro(moy) = mro(moy) + ro
        msm(moy) = msm(moy) + sm

 20     continue

       
 10     continue

        end ! run_one_year
      
c******************************************************************************c

c******************************************************************************c
        subroutine run_one_day( lat, elv, doy, yr, pr, tc, sf,
     >                          aet, pet, sm, cn, ro, dayl, sunl )
        implicit none
        !----------------------------------------------------------------
        ! Calculates daily and monthly quantities for one year
        !----------------------------------------------------------------
        ! arguments
        double precision lat   ! latitude (degrees)
        double precision elv   ! altitude (m)
        integer doy   ! day of year
        integer yr    ! year AD
        double precision pr    ! monthly or daily precip (mm)
        double precision tc    ! mean monthly or daily temperature (deg C)
        double precision sf    ! mean monthly or daily sunshine fraction (unitless)
        double precision aet   ! daily actual evapotranspiration (mm)
        double precision pet   ! daily potential evapotranspiration (mm)
        double precision cn    ! daily condensation (mm)
        double precision ro    ! daily runoff (mm)
        double precision sm    ! soil moisture (=water content, mm)
        double precision sw    ! evaporative supply rate (mm/h)
        double precision dayl       ! day length (hrs)
        double precision sunl       ! sunlight length (hrs)
        double precision kCw, kWm

        ! Parameters
        parameter(kCw = 1.05)           ! supply constant, mm/hr (Federer, 1982)
        parameter(kWm = 150)            ! soil moisture capacity, mm (Cramer & Prentice, 1988)
        
        ! Calculate evaporative supply rate, mm/h
        sw = kCw * sm / kWm

        ! Calculate radiation and evaporation quantities
        call evap( lat, doy, elv, yr, sf, tc, sw, cn, aet, pet,
     >             dayl, sunl )

        ! Update soil moisture
        sm = sm + pr + cn - aet

        if (sm>kWm) then
                ! Bucket is full
                ! * set soil moisture to capacity
                ! * add remaining water to monthly runoff total
                ro = sm - kWm  ! xxx add in sofun
                sm = kWm

        elseif (sm<0) then
                ! Bucket is empty
                ! * set soil moisture to zero
                aet = aet + sm
                ro = 0.0
                sm = 0.0 ! xxx check in sofun

        else
                ! No runoff occurrs
                ro = 0.0

        endif
        !write(*,*) 'sm:',sm
        !write(*,*) 'ro:',ro

        end ! subroutine run_one_day

c******************************************************************************c

c******************************************************************************c
        subroutine evap( lat, doy, elv, yr, sf, tc, sw,
     >                   cn, aet, pet, dayl, sunl)
        implicit none
        !----------------------------------------------------------------
        ! This subroutine calculates daily radiation and evapotranspiration
        ! quantities
        ! - daily photosynthetic photon flux density (ppfd), mol/m^2
        ! - daily equilibrium evapotranspiration (eet), mm
        ! - daily potential evapotranspiration (pet), mm
        ! - daily actual evapotranspiration (aet), mm
        ! - daily condensation (wc), mm
        !-------------------------------------------------------------
        ! arguments
        double precision lat           ! latitude, degrees
        integer doy           ! day of the year (formerly 'n')
        double precision elv ! elevation, metres
        integer yr  ! year
        double precision sf  ! fraction of sunshine hours
        double precision tc  ! mean daily air temperature, C
        double precision sw  ! evaporative supply rate, mm/hr

        ! out_evap variables
        double precision ra         ! daily top-of-atmosphere solar irradiation (J/m2)
        double precision rn         ! daily net radiation (J/m2)
        double precision ppfd       ! daily photosynthetic photon flux density (mol/m^2)
        double precision eet        ! daily equilibrium evapotranspiration (mm)
        double precision pet        ! daily potential evapotranspiration (mm)
        double precision aet        ! daily actual evapotranspiration (mm)
        double precision cn         ! daily condensation (mm)
        double precision dayl       ! day length (hrs)
        double precision sunl       ! sunlight length (hrs)
        
        ! out_berger variables
        double precision nu
        double precision lambda

        integer ndayyear,jd1,jd2

        double precision my_rho
        double precision dr                 ! distance factor
        double precision delta              ! declination angle
        double precision ru                 ! variable substitute for u
        double precision rv                 ! variable substitute for v
        double precision hs                 ! sunset hour angle
        double precision tau_o, tau         ! transmittivity (unitless)
        double precision rnl                ! net longwave radiation (W/m^2)
        double precision rnn_d              ! nighttime net radiation (J/m^2)
        double precision rw                 ! variable substitute (W/m^2)
        double precision hn                 ! net radiation cross-over hour angle
        double precision s                  ! slope of saturation vap press temp curve, Pa/K
        double precision pw                 ! density of water, kg/m^3
        double precision lv                 ! enthalpy of vaporization, J/kg
        double precision pres               ! air pressure (Pa)
        double precision g                  ! psychrometric constant, Pa/K
        double precision econ               ! water-to-energy conversion, m^3/J
        double precision rx                 ! variable substitute (mm/hr)/(W/m^2)
        double precision hi, cos_hi         ! intersection hour angle, degrees

        !----------------------------------------------------------------
        ! model parameters
        !----------------------------------------------------------------
        double precision kA,kalb_sw,kalb_vis,kb,kc,kCw,kd,kfFEC
        double precision kG,kGsc,kL,kMa,kMv,kPo,kR,kTo,kWm,kw
        double precision pi  
        parameter(kA = 107)             ! constant for Rnl (Monteith & Unsworth, 1990)
        parameter(kalb_sw = 0.17)       ! shortwave albedo (Federer, 1968)
        parameter(kalb_vis = 0.03)      ! visible light albedo (Sellers, 1985)
        parameter(kb = 0.20)            ! constant for Rnl (Linacre, 1968)
        parameter(kc = 0.25)            ! cloudy transmittivity (Linacre, 1968)
        parameter(kCw = 1.05)           ! supply constant, mm/hr (Federer, 1982)
        parameter(kd = 0.50)            ! angular coefficient of transmittivity (Linacre, 1968)
        parameter(kfFEC = 2.04)         ! from flux to energy conversion, umol/J (Meek et al., 1984)
        parameter(kG = 9.80665)         ! gravitational acceleration, m/s^2 (Allen, 1973)
        parameter(kGsc = 1360.8)        ! solar constant, W/m^2 (Kopp & Lean, 2011)
        parameter(kL = 0.0065)          ! temperature lapse rate, K/m (Cavcar, 2000)
        parameter(kMa = 0.028963)       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
        parameter(kMv = 0.01802)        ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)
        parameter(kPo = 101325)         ! standard atmosphere, Pa (Allen, 1973)
        parameter(kR = 8.31447)         ! universal gas constant, J/mol/K (Moldover et al., 1988)
        parameter(kTo = 288.15)         ! base temperature, K (Berberan-Santos et al., 1997)
        parameter(kWm = 150)            ! soil moisture capacity, mm (Cramer & Prentice, 1988)
        parameter(kw = 0.26)            ! entrainment factor (Lhomme, 1997; Priestley & Taylor, 1972)
        parameter(pi = 3.14159)

        ! Paleoclimate variables:
        double precision ke, keps, komega
        parameter(ke = 0.0167) ! eccentricity for 2000 CE (Berger, 1978)
        parameter(keps = 23.44) ! obliquity for 2000 CE, degrees (Berger, 1978)
        parameter(komega = 283.0) ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
        
        ! Error handle and assign required public variables:
        !if (lat>90.0 .or. lat<(-90.0)) then
        !print*,"Latitude outside range of validity (-90 to 90)!"
        !stop
        !endif

        !if (doy<1 .or. doy>366) then
        !print*,"Day of year outside range of validity (1 to 366)!"
        !print*,"doy:",doy
        !stop
        !endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 1. Calculate number of days in year (ndayyear), days (take from argument)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !call get_julian_day(yr, 1, 1, jd1)
        !call get_julian_day(yr+1, 1, 1, jd2)
        !ndayyear = jd2 - jd1
        ndayyear = 365
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 2. Calculate heliocentric longitudes (nu and lambda), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Berger (1978)
        call berger_tls( doy, ndayyear, nu, lambda )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 3. Calculate distance factor (dr), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Berger et al. (1993)
        my_rho = (1.0 - ke**2)/(1.0 + ke * cos(nu*pi/180))
        dr = (1.0/my_rho)**2
        !write(*,*) my_rho,dr

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 4. Calculate declination angle (delta), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Woolf (1968)
        delta = asin( sin(lambda*pi/180) * sin(keps *pi/180) )
        delta = delta * 180/pi
        !write(*,*) "delta:",delta
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 5. Calculate variable substitutes (u and v), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = sin(delta*pi/180) * sin(lat*pi/180)
        rv = cos(delta*pi/180) * cos(lat*pi/180)
        !write(*,*) "ru/rv:",ru,rv

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 6. Calculate the sunset hour angle (hs), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Note: u/v == tan(delta)*tan(lat)
        ! Eq. 3.22, Stine & Geyer (2001)
        if ((ru/rv) >= 1.0) then
                ! Polar day (no sunset)
                hs = 180.0
        elseif ((ru/rv) <= -1.0) then
                ! Polar night (no sunrise)
                hs = 0.0
        else
                hs = acos(-1.0*ru/rv) *180/pi
        endif
        !write(*,*) "hs:",hs

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 6.5. Calculate the day length and sunlight hours
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Equation 31 from splash_doc.pdf
        ! d_s = (24 * h_s) / pi
        dayl = (24 * hs / 180 * pi) / pi
        !write(*,*) "ds:",dayl
        sunl = dayl * sf
        !write(*,*) "sh:",sunl,sf
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 7. Calculate daily extraterrestrial solar radiation / irradiation (out_evap%ra), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 1.10.3, Duffy & Beckman (1993)
        ra = (86400.0/pi)*kGsc*dr*(ru*(hs*pi/180) + rv * sin(hs*pi/180))
        !write(*,*) "ra:",ra

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 8. Calculate transmittivity (tau), unitless
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 11, Linacre (1968)
        tau_o = (kc + kd*sf)

        ! Eq. 2, Allen (1996)
        tau = tau_o*(1.0 + (2.67e-5)*elv)
        !write(*,*) "tau_o/tau:",tau_o,tau
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 9. Calculate daily photosynthetic photon flux density (out_evap%ppfd), mol/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd = (1.0e-6)*kfFEC*(1.0 - kalb_vis)*tau*ra
        !write(*,*) "ppfd:",ppfd

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !  10. Estimate net longwave radiation (rnl), W/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Eq. 11, Prentice et al. (1993); Eq. 5 and 6, Linacre (1968)
        rnl = (kb + (1.0 - kb)*sf)*(kA - tc)
        !write(*,*) "rnl:",rnl

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 11. Calculate variable substitute (rw), W/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rw = (1.0-kalb_sw)*tau*kGsc*dr
        !write(*,*) "rw:",rw

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 12. Calculate net radiation cross-over hour angle (hn), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ((rnl - rw*ru)/(rw*rv) >= 1.0) then
                ! Net radiation negative all day
                hn = 0.0
        else if ((rnl - rw*ru)/(rw*rv) <= -1.0) then
                ! Net radiation positive all day
                hn = 180.0
        else
                hn = ( acos((rnl - rw*ru)/(rw*rv)) ) * 180/pi
        endif
        !write(*,*) "hn:",hn

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 13. Calculate daytime net radiation (out_evap%rn), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rn = (86400.0/pi) * (hn*(pi/180.0)*(rw*ru - rnl) + 
     *       rw*rv*sin(hn*pi/180))
        !write(*,*) "rn:",rn

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 14. Calculate nighttime net radiation (rnn_d), J/m^2
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rnn_d = (86400.0/pi)*((rw*ru*(hs-hn)*pi/180) + 
     *          rw*rv*(sin(hs*pi/180)-sin(hn*pi/180)) + 
     *          rnl*(pi - 2.0*(hs*pi/180) + (hn*pi/180)))
        !write(*,*) "rnn_d:",rnn_d
 
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 15. Calculate water-to-energy conversion (econ), m^3/J
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Slope of saturation vap press temp curve, Pa/K
        call get_sat_slope( tc, s )
        !write(*,*) "s:",s
        
        ! Enthalpy of vaporization, J/kg
        call get_enthalpy_vap( tc, lv )
        !write(*,*) "lv:",lv

        ! Density of water, kg/m^3
        call elv2pres(elv,pres)
        !write(*,*) "pres:",pres
        call get_density_h2o( tc, pres, pw )
        !write(*,*) "pw:",pw

        ! Psychrometric constant, Pa/K
        call get_psychro(tc, pres, lv, g)
        !write(*,*) "g:",g

        econ = s/(lv*pw*(s + g))
        !write(*,*) "econ:",econ
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 16. Calculate daily condensation (out_evap%cn), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cn = 1000.0 * econ * abs(rnn_d)
        !write(*,*) "cn:",cn

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 17. Estimate daily equilibrium evapotranspiration (out_evap%eet), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eet = 1000.0 * econ * rn
        !write(*,*) "eet:",eet

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 18. Estimate daily potential evapotranspiration (out_evap%pet), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pet = (1.0+kw)*eet
        !write(*,*) "pet:",pet

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 19. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rx = 1000.0*3600.0*(1.0+kw)*econ
        !write(*,*) "rx:",rx
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 20. Calculate the intersection hour angle (hi), degrees
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv
        !write(*,*) "cos_hi:",cos_hi

        if (cos_hi >= 1.0) then
                ! Supply exceeds demand:
                hi = 0.0
        elseif (cos_hi <= -1.0) then
                ! Supply limits demand everywhere:
                hi = 180.0
        else
                hi = acos(cos_hi) *180/pi
        endif
        !write(*,*) "hi:",hi
    
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! 21. Estimate daily actual evapotranspiration (out_evap%aet), mm
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        aet = (24.0/pi)*((sw*hi*pi/180) + 
     *        rx*rw*rv*(sin(hn*pi/180) - 
     *        sin(hi*pi/180)) + 
     *        (rx*rw*ru - rx*rnl)*(hn - hi)*pi/180)
        !write(*,*) "aet:",aet

        end ! End evap
c******************************************************************************c

c******************************************************************************c
        subroutine dgcos( x, out_dgcos )
        implicit none
        !----------------------------------------------------------------
        ! Calculates the cosine of an angle given in degrees. Equal to
        ! 'dsin' in Python version.
        !----------------------------------------------------------------
        ! arguments
        double precision x  ! angle, degrees (0-360)
        double precision pi  

        ! function return value
        double precision out_dgcos ! cosine value of x when x is in degrees
        parameter(pi = 3.14159)

        out_dgcos = cos(x*pi/180.0)

        end ! dgcos
c******************************************************************************c


c******************************************************************************c
        subroutine dgsin( x, out_dgsin )
        implicit none
        !----------------------------------------------------------------
        ! Calculates the sinus of an angle given in degrees. Equal to
        ! 'dsin' in Python version.
        !----------------------------------------------------------------
        ! arguments
        double precision x  ! angle, degrees (0-360)
        double precision pi  

        ! function return value
        double precision out_dgsin ! sinus value of x when x is in degrees
        parameter(pi = 3.14159)

        out_dgsin = sin(x*pi/180.0)

        end ! dgsin
c******************************************************************************c

c******************************************************************************c
        subroutine berger_tls( day, ndayyear, nu, lambda )
        implicit none
        !----------------------------------------------------------------
        ! Returns true anomaly and true longitude for a given day
        ! Reference: Berger, A. L. (1978), Long term variations of daily
        ! insolation and quaternary climatic changes, J. Atmos. Sci., 35,
        ! 2362-2367.
        !----------------------------------------------------------------
        ! arguments
        integer day   ! day of the year
        integer ndayyear

        ! function return variable
        double precision nu   
        double precision lambda   

        ! local variables
        double precision anm, ranm, anv, ranv
        double precision dlamm                ! Mean longitude for day of year
        double precision xee, xec, xse        ! variable substitutes
        double precision xlam                 ! Mean longitude for vernal equinox
        double precision tmp1, tmp2, tmp3     ! variable substitutes
        double precision pi  

        parameter(pi = 3.14159)

        ! Paleoclimate variables:
        double precision ke, keps, komega
        parameter(ke = 0.0167) ! eccentricity for 2000 CE (Berger, 1978)
        parameter(keps = 23.44) ! obliquity for 2000 CE, degrees (Berger, 1978)
        parameter(komega = 283.0) ! longitude of perihelion for 2000 CE, degrees (Berger, 1978)
        
        ! Variable substitutes:
        xee = ke**2
        xec = ke**3
        xse = sqrt(1.0 - xee)

        ! Mean longitude for vernal equinox:
        tmp1 = (ke/2.0 + xec/8.0)*(1.0 + xse)*
     *          sin(komega*pi/180)
        ! tmp1 = (ke/2.0 + xec/8.0)*(1.0 + xse)*dgsin(komega)
        tmp2 = xee/4.0*(0.5 + xse)*
     *         sin(2.0*komega*pi/180)
        ! tmp2 = xee/4.0*(0.5 + xse)*dgsin(2.0*komega)
        tmp3 = xec/8.0*(1.0/3.0 + xse)*
     *         sin(3.0*komega*pi/180)
        ! tmp3 = xec/8.0*(1.0/3.0 + xse)*dgsin(3.0*komega)
        xlam = tmp1 - tmp2 + tmp3
        xlam = (2.0*xlam)*180/pi

        ! Mean longitude for day of year:
        dlamm = xlam + (day - 80.0)*(360.0/ndayyear)

        ! Mean anomaly:
        anm = dlamm - komega
        ranm = anm*pi/180

        ! True anomaly:
        ranv = (ranm + (2.0*ke - xec/4.0) * sin(ranm) + 
     *          5.0/4.0*xee * sin(2.0*ranm) + 13.0/12.0*xec * 
     *          sin(3.0*ranm) )
        anv = ranv*180/pi

        ! True longitude:
        lambda = anv + komega
        if (lambda < 0.0) then
                lambda = lambda + 360.0
        else if (lambda > 360.0) then
                lambda = lambda - 360.0
        endif

        ! True anomaly:
        nu = (lambda - komega)
        if (nu < 0.0) then
                nu = nu + 360.0
        endif
        !write(*,*)nu,lambda

        end ! End berger_tls

c******************************************************************************c
        subroutine get_julian_day( yr, moy, dom, out_julian_day )
        implicit none
        !----------------------------------------------------------------
        ! Converts Gregorian date (year, month, day) to Julian
        ! Ephemeris Day
        ! * valid for dates after -4712 January 1 (i.e., jde >= 0)
        ! Reference:  Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day,"
        ! Astronomical Algorithms
        !----------------------------------------------------------------
        ! arguments
        integer yr      ! year
        integer moy     ! month of year (i.e., 1--12)
        integer dom     ! day of month (i.e., 1--31)

        ! local variables
        integer my_yr
        integer my_moy
        integer a, b

        ! function return value
        integer out_julian_day                   ! Julian Ephemeris Day

        if (moy <= 2) then
        my_yr  = yr - 1
        my_moy = moy + 12
        else
        my_yr  = yr
        my_moy = moy
        endif

        a = int(real(my_yr)/real(100))
        b = 2 - a + int(real(a)/real(4))
        out_julian_day=int(real(int(365.25*real(my_yr+4716)))+
     *          real(int(30.6001*real(my_moy+1)))+real(dom)+
     *          real(b)-1524.5)

        end ! get_julian_day
c******************************************************************************c
        
c******************************************************************************c
        subroutine get_sat_slope(tc, out_sat_slope)
        implicit none
        !----------------------------------------------------------------
        ! Calculates the slope of the sat pressure temp curve, Pa/K
        ! Ref:  Eq. 13, Allen et al. (1998)
        !----------------------------------------------------------------
        ! arguments
        double precision tc ! air temperature, degrees C

        ! function return value
        double precision out_sat_slope  ! slope of the sat pressure temp curve, Pa/K

        out_sat_slope = (17.269)*(237.3)*(610.78)*
     *              (exp(tc*17.269/(tc + 237.3))/((tc + 237.3)**2))

        end ! get_sat_slope
c******************************************************************************c

c******************************************************************************c
        subroutine get_enthalpy_vap(tc, out_enthalpy_vap)
        implicit none
        !----------------------------------------------------------------
        ! Calculates the enthalpy of vaporization, J/kg
        ! Ref:  Eq. 8, Henderson-Sellers (1984)
        !----------------------------------------------------------------
        ! arguments
        double precision tc ! air temperature, degrees C

        ! function return value
        double precision out_enthalpy_vap ! enthalpy of vaporization, J/kg

        out_enthalpy_vap = 1.91846e6*((tc + 273.15)/
     *                     (tc + 273.15 - 33.91))**2

        end ! get_enthalpy_vap
c******************************************************************************c

c******************************************************************************c
        subroutine elv2pres(alt, press)
        implicit none
        !----------------------------------------------------------------
        ! Calculates atm. pressure for a given elevation
        ! Ref:  Allen et al. (1998)
        !----------------------------------------------------------------
        ! arguments
        double precision alt ! elevation above sea level, m

        ! function return value
        double precision press ! atm. pressure for a given elevation

        ! parameters
        double precision kG,kL,kMa,kPo,kR,kTo
        parameter(kG = 9.80665)         ! gravitational acceleration, m/s^2 (Allen, 1973)
        parameter(kL = 0.0065)          ! temperature lapse rate, K/m (Cavcar, 2000)
        parameter(kMa = 0.028963)       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
        parameter(kPo = 101325)         ! standard atmosphere, Pa (Allen, 1973)
        parameter(kR = 8.31447)         ! universal gas constant, J/mol/K (Moldover et al., 1988)
        parameter(kTo = 288.15)         ! base temperature, K (Berberan-Santos et al., 1997)
        
        press = kPo*(1.0 - kL*alt/kTo)**(kG*kMa/(kR*kL))

        end ! elv2pres
c******************************************************************************c
        
c******************************************************************************c
        subroutine get_density_h2o(tc, press, density_h2o)
        implicit none
        !----------------------------------------------------------------
        ! Calculates density of water at a given temperature and pressure
        ! Ref:  Chen et al. (1977)
        !----------------------------------------------------------------
        ! arguments
        double precision tc     ! air temperature (degrees C)
        double precision press  ! atmospheric pressure (Pa)

        ! local variables
        double precision po, ko, ca, cb
        double precision pbar               ! atmospheric pressure (bar)

        ! function return value
        double precision density_h2o  ! density of water, kg/m^3

        ! Calculate density at 1 atm:
        po = (0.99983952
     *        + 6.788260e-5  *tc
     *        - 9.08659e-6   *tc*tc
     *        + 1.022130e-7  *tc*tc*tc  
     *        - 1.35439e-9   *tc*tc*tc*tc
     *        + 1.471150e-11 *tc*tc*tc*tc*tc
     *        - 1.11663e-13  *tc*tc*tc*tc*tc*tc
     *        + 5.044070e-16 *tc*tc*tc*tc*tc*tc*tc
     *        - 1.00659e-18  *tc*tc*tc*tc*tc*tc*tc*tc)

        ! Calculate bulk modulus at 1 atm:
        ko = (19652.17
     *        + 148.1830   *tc
     *        - 2.29995    *tc*tc
     *        + 0.01281    *tc*tc*tc
     *        - 4.91564e-5 *tc*tc*tc*tc
     *        + 1.035530e-7*tc*tc*tc*tc*tc)

        ! Calculate temperature dependent coefficients:
        ca = (3.26138
     *        + 5.223e-4  *tc
     *        + 1.324e-4  *tc*tc
     *        - 7.655e-7  *tc*tc*tc
     *        + 8.584e-10 *tc*tc*tc*tc)

        cb = (7.2061e-5
     *        - 5.8948e-6  *tc
     *        + 8.69900e-8 *tc*tc
     *        - 1.0100e-9  *tc*tc*tc
     *        + 4.3220e-12 *tc*tc*tc*tc)

        ! Convert atmospheric pressure to bar (1 bar = 100000 Pa)
        pbar = (1.0e-5)*press

        density_h2o = 1000.0*po*(ko + ca*pbar + cb*pbar**2.0)/
     *                (ko + ca*pbar + cb*pbar**2.0 - pbar)

        end ! get_density_h2o
        
c******************************************************************************c
        subroutine get_psychro(tc, press, lv, psychro)
        implicit none
        !----------------------------------------------------------------
        ! Calculates the psychrometric constant for a given temperature 
        ! and pressure
        ! Ref:  Allen et al. (1998); Tsilingiris (2008)
        !----------------------------------------------------------------
        ! arguments
        double precision tc ! air temperature, degrees C
        double precision press  ! atmospheric pressure, Pa

        ! local variables
        double precision lv  ! latent heat of vaporization (J/kg)
        double precision cp

        ! function return value
        double precision psychro  ! psychrometric constant, Pa/K

        ! Parameters
        double precision kMa,kMv
        parameter(kMa = 0.028963)       ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
        parameter(kMv = 0.01802)        ! molecular weight of water vapor, kg/mol (Tsilingiris, 2008)

        ! Calculate the specific heat capacity of water, J/kg/K
        ! Eq. 47, Tsilingiris (2008)
        cp = 1.0e3*(1.0045714270
     *            + 2.050632750e-3  *tc
     *            - 1.631537093e-4  *tc*tc
     *            + 6.212300300e-6  *tc*tc*tc
     *            - 8.830478888e-8  *tc*tc*tc*tc
     *            + 5.071307038e-10 *tc*tc*tc*tc*tc)

        ! Calculate latent heat of vaporization, J/kg
        ! Passed from evap subroutine
        !lv = get_enthalpy_vap( tc )

        ! Calculate psychrometric constant, Pa/K
        ! Eq. 8, Allen et al. (1998)
        psychro = cp*kMa*press/(kMv*lv)

        end ! get_psychro
c******************************************************************************c

c******************************************************************************c
      subroutine daily(mly,dly)
      implicit none
      double precision mly(12),dly(365),midday(12),vinc
      integer im,id
      data (midday(im),im=1,12)/16., 44., 75.,105.,136.,166.,
     >                         197.,228.,258.,289.,319.,350./

      vinc=(mly(1)-mly(12))/31.0
      dly(350)=mly(12)
      do 100 id=351,365
         dly(id)=dly(id-1)+vinc
100      continue
      dly(1)=dly(365)+vinc
      do 101 id=2,15
         dly(id)=dly(id-1)+vinc
101      continue
      do 103 im=1,11
         vinc=(mly(im+1)-mly(im))/(midday(im+1)-midday(im))
         dly(int(midday(im)))=mly(im)
         do 104 id=int(midday(im))+1,int(midday(im+1))-1
            dly(id)=dly(id-1)+vinc
104         continue
103      continue
      return
      end   
c******************************************************************************c

