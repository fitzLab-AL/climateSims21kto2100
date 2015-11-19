!gfortran -O3 -o ET ET.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program temp0

  use netcdf
  use constants
  use read_data
  use nc_write
  implicit none

  integer, parameter :: year1 = 1950, year2 = 2005

  ! future scenarios:
  integer, parameter :: year11 = 2006, year22 = 2100
  character (len = 40) :: scens(2) = [character (len = 40) :: 'rcp45','rcp85']

  real, parameter :: soilmax = 150.

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
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  integer, parameter :: nmodel = 12
  character (len = 60) :: models(nmodel) = [character (len=60) :: &
       'ACCESS1-3','CanESM2','CESM1-CAM5','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-CM3', &
       'GISS-E2-R','HadGEM2-ES','inmcm4','IPSL-CM5A-MR','MIROC5','MRI-CGCM3']
  
  real, parameter :: days(12) = [31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]

  real, allocatable :: TT(:,:,:,:), ps(:,:)
  real, allocatable :: xout(:,:,:),xout_init(:,:,:),et(:,:,:,:)
  integer*2, allocatable :: ix(:,:,:)
  real, parameter :: sf=0.1, ao=0.0
  integer, parameter :: nspinup = 100 ! number of years to spinup soil moisture/snow
  
  character (len=400) :: fname
  real :: temp
  integer :: imodel,iscen,ilat,ilon,imonth,iyear,imm,imp,iym,iyp,mspin,ispin,ivar

  
  file_expr = '/Users/djlorenz/Research/jack/downscale/elevation.nc '
  vars = 'ps '
  call read_nc(missing,missing,.false.,missing,missing,missing,missing, &
       [all,all,1,2001],[all,all,1,2001],1,'month')
  ! CHECK x1 allocation **************************:
  allocate(ps(nlon,nlat),ix(nlon,nlat,2))
  allocate(et(nlon,nlat,2,12),xout(nlon,nlat,4),xout_init(nlon,nlat,4))
  ps = x1(:,:,1,1,1)
  deallocate(x1)

  xout = 0.0
  ! do models:
  do imodel = 1,nmodel
     print*, imodel,trim(models(imodel))

     ! read full temperature data for "dT"
     file_expr = '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
          trim(models(imodel))//'/temp.nc'
     vars = 'tmax tmin'
     call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
          [all,all,all,year1],[all,all,all,year2],1,'month')
     allocate(TT(nlon,nlat,12,nyear))
     TT = merge(0.5*sum(x1,dim=3),missing,x1(:,:,1,:,:)/=missing)
     deallocate(x1)

     fname = '/Users/djlorenz/Research/jack/final/cmip5/historical/'//trim(models(imodel))// &
          '/ET.nc'
     call ncopen(1,trim(fname),'pet aet', &
          vunits='mm',vnames='potential evapotranspiration; actual evapotranspiration;', &
          lons=lons,lats=lats,title='Downscaled CMIP5 models, monthly mean.', &
          deflate=1,sf=[sf],ao=[ao],syear=year1,tunit='M')

     mspin = nspinup
     do iyear=year1,year2
        
        file_expr = '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/prcp.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/temp.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/vap.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/wind.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/swd.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/lwn.nc '// &
             '/Users/djlorenz/Research/jack/final/cmip5/historical/'// &
             trim(models(imodel))//'/swu.nc '
        vars = 'prcp tmax tmin vap wind swd lwn swu'
        call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
             [all,all,all,iyear],[all,all,all,iyear],1,'month')

        do ispin = 1,mspin
           do imonth=1,12
              iym = iyear; iyp = iyear
              imm = imonth - 1
              if(imm < 1) then
                 imm = imm + 12
                 iym = iym - 1
                 if(iym < year1) then
                    iym = year1
                 endif
              endif
              imp = imonth + 1
              if(imp > 12) then
                 imp = imp - 12
                 iyp = iyp + 1
                 if(iyp > year2) then
                    iyp = year2
                 endif
              endif
              iym = iym - year1 + 1
              iyp = iyp - year1 + 1
              
              do ilat=1,nlat
                 do ilon=1,nlon
                    
                    if(.not.any(x1(ilon,ilat,:,1,1)==missing)) then
                       
                       prcp = x1(ilon,ilat,1,imonth,1)
                       Tmax = x1(ilon,ilat,2,imonth,1)
                       Tmin = x1(ilon,ilat,3,imonth,1)
                       ! convert vap from hPa to kPa:
                       vap = 0.1*x1(ilon,ilat,4,imonth,1)
                       T = 0.5*(Tmin + Tmax)                 
                       u = x1(ilon,ilat,5,imonth,1)
                       P = ps(ilon,ilat)
                       
                       soil = xout(ilon,ilat,1)
                       snow = xout(ilon,ilat,2)
                       
                       R = 0.0864*(x1(ilon,ilat,6,imonth,1) - x1(ilon,ilat,8,imonth,1) + x1(ilon,ilat,7,imonth,1))
                       
                       dT = TT(ilon,ilat,imp,iyp) - TT(ilon,ilat,imm,iym)
                       
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
              
              et(:,:,:,imonth) = xout(:,:,3:4)
              
           enddo
           
        enddo
        mspin = 1


        do imonth = 1,12
           do ivar=1,2
              do ilat=1,nlat
                 do ilon=1,nlon
                    if(et(ilon,ilat,ivar,imonth)/=missing) then
                       temp = (et(ilon,ilat,ivar,imonth)-ao)/sf
                       if(abs(temp)>32767.) then
                          print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                          stop
                       endif
                       ix(ilon,ilat,ivar) = nint(temp)
                    else
                       ix(ilon,ilat,ivar) = missing_short
                    endif
                 enddo
              enddo
           enddo
           call ncwrite(1,i1=ix)
        enddo

     enddo
     xout_init = xout

     
     call ncclose(1)
     deallocate(TT)

     ! future scenarios:
     do iscen = 1,2
        
        ! read full temperature data for "dT"
        file_expr = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/temp.nc'
        vars = 'tmax tmin'
        call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
             [all,all,all,year11],[all,all,all,year22],1,'month')
        allocate(TT(nlon,nlat,12,nyear))
        TT = merge(0.5*sum(x1,dim=3),missing,x1(:,:,1,:,:)/=missing)
        deallocate(x1)
        
        fname = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
             trim(models(imodel))//'/ET.nc'
        call ncopen(1,trim(fname),'pet aet', &
             vunits='mm',vnames='potential evapotranspiration; actual evapotranspiration;', &
             lons=lons,lats=lats,title='Downscaled CMIP5 models, monthly mean.', &
             deflate=1,sf=[sf],ao=[ao],syear=year11,tunit='M')
        
        xout = xout_init
        do iyear=year11,year22
           
           file_expr = '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/prcp.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/temp.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/vap.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/wind.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/swd.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/lwn.nc '// &
                '/Users/djlorenz/Research/jack/final/cmip5/'//trim(scens(iscen))//'/'// &
                trim(models(imodel))//'/swu.nc '
                vars = 'prcp tmax tmin vap wind swd lwn swu'
           call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
                [all,all,all,iyear],[all,all,all,iyear],1,'month')
           
           do imonth=1,12
              iym = iyear; iyp = iyear
              imm = imonth - 1
              if(imm < 1) then
                 imm = imm + 12
                 iym = iym - 1
                 if(iym < year11) then
                    iym = year11
                 endif
              endif
              imp = imonth + 1
              if(imp > 12) then
                 imp = imp - 12
                 iyp = iyp + 1
                 if(iyp > year22) then
                    iyp = year22
                 endif
              endif
              iym = iym - year11 + 1
              iyp = iyp - year11 + 1
              
              do ilat=1,nlat
                 do ilon=1,nlon
                    
                    if(.not.any(x1(ilon,ilat,:,1,1)==missing)) then
                       
                       prcp = x1(ilon,ilat,1,imonth,1)
                       Tmax = x1(ilon,ilat,2,imonth,1)
                       Tmin = x1(ilon,ilat,3,imonth,1)
                       ! convert vap from hPa to kPa:
                       vap = 0.1*x1(ilon,ilat,4,imonth,1)
                       T = 0.5*(Tmin + Tmax)                 
                       u = x1(ilon,ilat,5,imonth,1)
                       P = ps(ilon,ilat)
                       
                       soil = xout(ilon,ilat,1)
                       snow = xout(ilon,ilat,2)
                       
                       R = 0.0864*(x1(ilon,ilat,6,imonth,1) - x1(ilon,ilat,8,imonth,1) + x1(ilon,ilat,7,imonth,1))
                       
                       dT = TT(ilon,ilat,imp,iyp) - TT(ilon,ilat,imm,iym)
                       
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
              
              et(:,:,:,imonth) = xout(:,:,3:4)
              
           enddo
           
           do imonth = 1,12
              do ivar=1,2
                 do ilat=1,nlat
                    do ilon=1,nlon
                       if(et(ilon,ilat,ivar,imonth)/=missing) then
                          temp = (et(ilon,ilat,ivar,imonth)-ao)/sf
                          if(abs(temp)>32767.) then
                             print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                             stop
                          endif
                          ix(ilon,ilat,ivar) = nint(temp)
                       else
                          ix(ilon,ilat,ivar) = missing_short
                       endif
                    enddo
                 enddo
              enddo
              call ncwrite(1,i1=ix)
           enddo
           
        enddo
        
        call ncclose(1)
        deallocate(TT)
        
     enddo
     
  enddo
  
  
  
  
end program temp0

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
