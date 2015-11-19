!gfortran -O3 -o gdd_final gdd_final.f90 read_final.o write_final.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program gdd0

  use netcdf
  use constants
  use read_data
  use read_final
  use write_final
  implicit none

  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1

  real, parameter :: days(12) = [31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]

  real, allocatable :: tmstd(:,:,:),T(:,:,:)
  integer*2, allocatable :: ix(:,:,:,:)
  real, parameter :: sf=0.1, ao=0.0

  real :: temp,gdd,g,tmean(12),tastd(12),tstd(12)

  integer :: idecade,ilat,ilon,imonth,iyear,imp,imm,ivar,im,ip,i0
  real :: Lm,L0,Lp,am,ap

  file_expr = 'gdd_obs.nc '
  vars = 'tmstd '
  call read_nc(missing,missing,.false.,missing,missing,missing,missing, &
       [all,all,all,2001],[all,all,all,2001],1,'month')
  allocate(tmstd(nlon,nlat,12))
  tmstd = x1(:,:,1,:,1)
  deallocate(x1)

  do idecade=decade1,decade2
     print*, idecade

     call  get_final(missing,missing,missing,missing,idecade, &
          'tmax tmin')
     !     1    2   
     if(.not.allocated(T)) then
        allocate(T(nlonc,nlatc,12),ix(nlonc,nlatc,12,2))
        call ncopen(1,'/Users/djlorenz/Research/jack/final/CCSM/gdd.nc','gdd0 gdd5', &
             vunits='C',vnames='growing degree days base 0C; growing degree days base 5C;', &
             lons=lonsc,lats=latsc,title='Downscaled TraCE CCSM runs. Monthly means by decade.', &
             deflate=1,sf=[sf],ao=[ao])
     endif

     T = merge(0.5*(xc(:,:,1,:) + xc(:,:,2,:)),missing, &
          xc(:,:,1,:)/=missing .and. xc(:,:,2,:)/=missing)
     
     do ilat=1,nlat
        do ilon=1,nlon
           if(T(ilon,ilat,1)/=missing) then
              Tmean = T(ilon,ilat,:)
              do imonth=1,12
                 i0 = imonth
                 im = modulo(i0-2,12)+1
                 ip = modulo(i0,12)+1
                 L0 = days(i0); Lm = days(im); Lp = days(ip)
                 am = 2.0*(Tmean(i0) - Tmean(im))/(L0 + Lm)
                 ap = 2.0*(Tmean(ip) - Tmean(i0))/(Lp + L0)
                 tastd(imonth) = L0*sqrt((5.*ap**2 + 6.*ap*am + 5.*am**2)/192.)
              enddo
              tstd = sqrt(0.178211466801*days*tmstd(ilon,ilat,:)**2 + tastd**2)
              do imonth=1,12
                 g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),0.0)
                 temp = (g-ao)/sf
                 if(abs(temp)>32767.) then
                    print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                    stop
                 endif
                 ix(ilon,ilat,imonth,1) = nint(temp)
                 g = days(imonth)*gdd(Tmean(imonth),tstd(imonth),5.0)
                 temp = (g-ao)/sf
                 if(abs(temp)>32767.) then
                    print*, 'ERROR: conversion to 2 byte: out of range. ',temp
                    stop
                 endif
                 ix(ilon,ilat,imonth,2) = nint(temp)
              enddo
           else
              ix(ilon,ilat,:,:) = missing_short
           endif
        enddo
     enddo
   
     call ncwrite(1,im=ix)

  enddo
  
  call ncclose(1)

end program gdd0
   
function gdd(Tmean,sigma,T0)

  use constants
  implicit none

  real :: gdd
  real :: Tmean,sigma,T0

!  gdd = sigma/sqrt(2.*pi)*exp(-0.5*((Tmean-T0)/sigma)**2) + &
!       0.5*(Tmean - T0)*(1.0 + erf((Tmean - T0)/(sqrt(2.)*sigma)))
  gdd = sigma/sqrt(2.*pi)*exp(-0.5*((Tmean-T0)/sigma)**2) + &
       0.5*(Tmean - T0)*erfc((T0 - Tmean)/(sqrt(2.)*sigma))

end function gdd
