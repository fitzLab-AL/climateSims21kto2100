!calculate potential and actual evapotranspiration
!gfortran -O3 -o ET ET.f90 read_final.o write_final.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program ET0

  use netcdf
  use constants
  use read_final
  use write_final
!  use nc_write
  implicit none

  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1
  real, parameter :: albedo = 0.15, albedo_snow = 0.65, soilmax = 150.

  ! prcp, met, snow fall in mm/month
  ! radiation in MJ/m2/day

  ! prescribed:
  real :: prcp,T,R ! prcp, mean temperature, net surface radiation
  real :: Tmax,Tmin,P,u ! daily max/min Temperature, Pressure, wind speed
  real :: vap,dT ! vapor pressure (kPa), change in temperature (imonth + 1 minus imonth - 1)

  real :: rain,snowfall,snow ! rain, snowfall,snow on ground
  real :: melt, dsnow ! snow melt (includes parts from falling and accumulated snow)
  real :: Win ! liquid water input to system
  real :: PET, AET ! potential and actual evapotranspiration
  real :: Fm ! fraction  of precip that's rain
  real :: soil ! soil water
  real :: dsoil

  real :: rain_fraction, snow_melt, Penman_Monteith

  real, parameter :: days(12) = [31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]

  real, allocatable :: xout(:,:,:),et(:,:,:,:)
  integer*2, allocatable :: ix(:,:,:,:)
  real, parameter :: sf=0.1, ao=0.0
  real :: temp

  integer :: idecade,ilat,ilon,imonth,iyear,imp,imm,ivar

  do idecade=decade1,decade2
     print*, idecade

     call  get_final(missing,missing,missing,missing,idecade, &
          'prcp tmax tmin vap wind ps swd lwn swu')
     !     1    2    3    4   5    6  7   8   9
     if(.not.allocated(xout)) then
        allocate(xout(nlonc,nlatc,4),et(nlonc,nlatc,12,2),ix(nlonc,nlatc,12,2))
        xout = 0.0
        call ncopen(1,'/Users/djlorenz/Research/jack/final/CCSM/ET.nc','pet aet', &
             vunits='mm',vnames='potential evapotranspiration; actual evapotranspiration;', &
             lons=lonsc,lats=latsc,title='Downscaled TraCE CCSM runs. Monthly means by decade.', &
             deflate=1,sf=[sf],ao=[ao])
     endif

     ! if first decade, then spin up: 
     if(idecade==decade1) then
        
        do iyear=1,50
           do imonth=1,12
              imp = modulo(imonth,12) + 1
              imm = modulo(imonth-2,12) + 1
              
              do ilat=1,nlatc
                 do ilon=1,nlonc
                    
                    if(.not.any(xc(ilon,ilat,:,1)==missing)) then
                       
                       prcp = xc(ilon,ilat,1,imonth)
                       Tmax = xc(ilon,ilat,2,imonth)
                       Tmin = xc(ilon,ilat,3,imonth)
                       ! convert vap from hPa to kPa:
                       vap = 0.1*xc(ilon,ilat,4,imonth)
                       T = 0.5*(Tmin + Tmax)                 
                       u = xc(ilon,ilat,5,imonth)
                       P = xc(ilon,ilat,6,imonth)
                       
                       soil = xout(ilon,ilat,1)
                       snow = xout(ilon,ilat,2)

!                       if(snow==0.) then
                          ! 0.0864 converts from average W/m**2 to MJ/m**2/day
!                          R = 0.0864*((1.0 - albedo)*xc(ilon,ilat,7,imonth) + xc(ilon,ilat,8,imonth))
                          R = 0.0864*(xc(ilon,ilat,7,imonth) - xc(ilon,ilat,9,imonth) + xc(ilon,ilat,8,imonth))
!                       else
!                          ! 0.0864 converts from average W/m**2 to MJ/m**2/day
!                          R = 0.0864*((1.0 - albedo_snow)*xc(ilon,ilat,7,imonth) + xc(ilon,ilat,8,imonth))
!                       endif
                       
                       dT = 0.5*(xc(ilon,ilat,2,imp) + xc(ilon,ilat,3,imp) - &
                            xc(ilon,ilat,2,imm) - xc(ilon,ilat,3,imm))
                       
                       Fm = rain_fraction(T)
                       rain = Fm*prcp
                       snowfall = (1.0 - Fm)*prcp
                       
                       melt = snow_melt(T,R,snow,snowfall)
                       dsnow = snowfall - melt
                       
                       Win = rain + melt
                       PET = days(imonth)*Penman_Monteith(R,Tmax,Tmin,vap,P,u,dT) ! days = days in month
                       
                       if(Win >= PET) then
                          AET = PET
                          ! from Lutz et al. (2010)
                          dsoil = min(soilmax,(Win - PET) + soil) - soil
                       else
                          ! One way, from Lutz:
                          ! from Lutz et al. (2010) (note minus relative to Lutz)
                          !dsoil = soil*(exp(-(PET - Win)/soilmax) - 1.0) ! this is negative
                          !AET = min(PET,Win - dsoil)
                          ! Or second Way, more consistent soil water budget:
                          dsoil = max(soil*(exp(-(PET - Win)/soilmax) - 1.0),Win - PET) ! max means least negative
                          AET = Win - dsoil
                       endif
                       
                       !last step
                       snow = snow + dsnow
                       soil = soil + dsoil
                       
                       xout(ilon,ilat,1) = soil
                       xout(ilon,ilat,2) = snow
                       xout(ilon,ilat,3) = PET
                       xout(ilon,ilat,4) = AET
                       
                    else
                       xout(ilon,ilat,:) = missing
                    endif
                    
                 enddo
              enddo
              
           enddo
        enddo
        
     endif

     et = 0.0
     do iyear=1,10
        do imonth=1,12
           imp = modulo(imonth,12) + 1
           imm = modulo(imonth-2,12) + 1

           do ilat=1,nlatc
              do ilon=1,nlonc

                 if(.not.any(xc(ilon,ilat,:,1)==missing)) then

                    prcp = xc(ilon,ilat,1,imonth)
                    Tmax = xc(ilon,ilat,2,imonth)
                    Tmin = xc(ilon,ilat,3,imonth)
                    ! convert vap from hPa to kPa:
                    vap = 0.1*xc(ilon,ilat,4,imonth)
                    T = 0.5*(Tmin + Tmax)                 
                    u = xc(ilon,ilat,5,imonth)
                    P = xc(ilon,ilat,6,imonth)
                    
                    soil = xout(ilon,ilat,1)
                    snow = xout(ilon,ilat,2)
                    
!                    if(snow==0.) then
                       ! 0.0864 converts from average W/m**2 to MJ/m**2/day
!                       R = 0.0864*((1.0 - albedo)*xc(ilon,ilat,7,imonth) + xc(ilon,ilat,8,imonth))
                       R = 0.0864*(xc(ilon,ilat,7,imonth) - xc(ilon,ilat,9,imonth) + xc(ilon,ilat,8,imonth))
!                    else
!                       ! 0.0864 converts from average W/m**2 to MJ/m**2/day
!                       R = 0.0864*((1.0 - albedo_snow)*xc(ilon,ilat,7,imonth) + xc(ilon,ilat,8,imonth))
!                    endif

                    dT = 0.5*(xc(ilon,ilat,2,imp) + xc(ilon,ilat,3,imp) - &
                         xc(ilon,ilat,2,imm) - xc(ilon,ilat,3,imm))

                    Fm = rain_fraction(T)
                    rain = Fm*prcp
                    snowfall = (1.0 - Fm)*prcp
                    
                    melt = snow_melt(T,R,snow,snowfall)
                    dsnow = snowfall - melt
                    
                    Win = rain + melt
                    PET = days(imonth)*Penman_Monteith(R,Tmax,Tmin,vap,P,u,dT) ! days = days in month
                    
                    if(Win >= PET) then
                       AET = PET
                       ! from Lutz et al. (2010)
                       dsoil = min(soilmax,(Win - PET) + soil) - soil
                    else
                       ! One way, from Lutz:
                       ! from Lutz et al. (2010) (note minus relative to Lutz)
                       !dsoil = soil*(exp(-(PET - Win)/soilmax) - 1.0) ! this is negative
                       !AET = min(PET,Win - dsoil)
                       ! Or second Way, more consistent soil water budget:
                       dsoil = max(soil*(exp(-(PET - Win)/soilmax) - 1.0),Win - PET) ! max means least negative
                       AET = Win - dsoil
                    endif
                    
                    !last step
                    snow = snow + dsnow
                    soil = soil + dsoil
                    
                    xout(ilon,ilat,1) = soil
                    xout(ilon,ilat,2) = snow
                    xout(ilon,ilat,3) = PET
                    xout(ilon,ilat,4) = AET

                 else
                    xout(ilon,ilat,:) = missing
                 endif

              enddo
           enddo
           
           et(:,:,imonth,:) = merge(et(:,:,imonth,:) + xout(:,:,3:4),missing, &
                et(:,:,imonth,:)/=missing .and. xout(:,:,3:4)/=missing)

        enddo
     enddo

     et = merge(et/10.0,missing,et/=missing)

     do ivar=1,2
        do imonth=1,12
           do ilat=1,nlatc
              do ilon=1,nlonc
                 if(et(ilon,ilat,imonth,ivar)/=missing) then
                    temp = (et(ilon,ilat,imonth,ivar)-ao)/sf
                    if(abs(temp)>32767.) then
                       print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                       stop
                    endif
                    ix(ilon,ilat,imonth,ivar) = nint(temp)
                 else
                    ix(ilon,ilat,imonth,ivar) = missing_short
                 endif
              enddo
           enddo
        enddo
     enddo
     
     call ncwrite(1,im=ix)

  enddo

  call ncclose(1)

end program ET0

! potential evapotranspiration
! R = MJ/m2/day
! Tmax, Tmin = C   (2m observation) 
! P,vap = kPa
! u = m/s    (10m observation)
! dT = C,  dT is T at following month  minus  T in preceeding month  
function Penman_Monteith(R,Tmax,Tmin,vap,P,u,dT)

  implicit none

  real :: Penman_Monteith
  ! vap = vapor pressure in kPa
  real :: R,Tmax,Tmin,vap,P,u
  real :: T,dT,G

  real :: ra,rs,gamma,es,Delta,rho

  real, parameter :: h = 0.12
  real, parameter :: rau = 277.71759 ! aerodynamic resistance * u

  real, parameter :: Rd = 287., cp = 1004., Lv = 2.45e6
  real, parameter :: eps = 0.622

  T = 0.5*(Tmax + Tmin)

  ra = rau/u
  rs = surface_resistance(T)

  gamma = psychrometric_constant(P)
  es = 0.5*(esat(Tmax) + esat(Tmin)) - vap
  Delta = slope_of_esat(T)
  rho = P/(Rd*(T+273.15))

  ! R - G in MJ/m2/day
  ! rho*cp in kPa/C (assuming P in kPa)
  ! es/ra in kPa*m/s = kN*m / m**2 / s => multiply by 86400./1000. = 86.4

  G = 0.07*dT

  Penman_Monteith = (Delta*(R - G) + 86.4*rho*cp*es/ra)/(Delta + gamma*(1.0 + rs/ra))

  Penman_Monteith = max(0.408*Penman_Monteith,0.) ! convert from MJ/m2/day to mm/day

contains 

  function slope_of_esat(T)  ! kPa/C

    real :: slope_of_esat,T

    slope_of_esat = 4098.*(0.6108*exp(17.27*T/(T + 237.3)))/(T + 237.3)**2

  end function slope_of_esat

  ! from Allen et al (1998)
  function esat(T) ! kPa,   also: T in Celsius

    real :: esat,T

    esat = 0.6108*exp(17.27*T/(T + 237.3))

  end function esat

  ! from Allen et al (1998)
  function psychrometric_constant(P) ! kPa/C

    real :: psychrometric_constant,P

    psychrometric_constant = cp*P/(eps*Lv) ! units = (units of P)/Celsius

  end function psychrometric_constant

  ! from Dobrowski et al (2012)
  function surface_resistance(T) ! s/m

    real :: surface_resistance,T

    real :: ks

    real, parameter :: rl = 100.0
    real, parameter :: LAIactive = 1.44

    real, parameter :: Tl = -10.0, To = 5.0, Th = 100.0
    real, parameter :: ksmin = 0.01

    real, parameter :: b4 = (Th - To)/(Th - Tl)
    real, parameter :: b3 = ((To - Tl)*(Th - To))**(-b4)

    if(T >= To) then
       ks = 1.0
    elseif(T <= Tl) then
       ks = ksmin
    else
       ks = b3*((T - Tl)*(Th - T))**b4
       ks = max(ks,ksmin)
    endif

    surface_resistance = rl/ks/LAIactive

  end function surface_resistance

end function Penman_Monteith

! from Dobrowski et al (2012)
function snow_melt(T,R,snow,snowfall) ! in mm/month

  implicit none

  real :: snow_melt
  real :: T,R,snow,snowfall

  real,  parameter :: b0 = -398., b1 = 81.7, b2 = 25.0

  snow_melt = b0 + b1*T + b2*R
  snow_melt = max(snow_melt,0.0)
  snow_melt = min(snow_melt,snow + snowfall)

end function snow_melt

! from Dobrowski et al (2012)
function rain_fraction(T)

  implicit none

  real :: rain_fraction
  real :: T

  real, parameter :: TL = -4.6, TH = 6.3

  if(T >= TH) then
     rain_fraction = 1.0
  elseif(T <= TL) then
     rain_fraction = 0.0
  else
     rain_fraction = (T - TL)/(TH - TL)
  endif

end function rain_fraction
