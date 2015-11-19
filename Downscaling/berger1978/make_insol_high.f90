!DiFFERENCE with make_insol.f90: high resolution & monthly
!make insolation dataset (latitude, month and decade)
!gfortran -O3 -o make_insol_high insol2.o make_insol_high.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program make_insol

  use netcdf
  use constants
  use nc_write
  implicit none

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1

  integer, parameter :: nlat=360
  real :: lats(nlat)
  real :: kyear
  integer, parameter :: days(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  real :: x(nlat,decade1:decade2,12)
  integer :: ix(nlat,12)

  integer :: ilat,idecade,imonth,day
  real*8 :: phi,t,ecc,pre,perh,xob,ww,dayl

  do ilat=1,nlat
     lats(ilat) = 180./nlat*(ilat - nlat/2 - 0.5) 
  enddo

  call ctlwrite('insol_high.ctl','i',lats=lats,nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=12)

  x = 0.0
  do idecade=decade1,decade2
     print*, idecade
     kyear = 0.01*(idecade + 0.5)
     t = dble(kyear)
     ix = 0
     do imonth=1,12
        do day=1,days(imonth)
           do ilat=1,nlat
              phi = dble(lats(ilat))
              ix(ilat,imonth) = ix(ilat,imonth) + 1
              call insol(phi,t,imonth,day,ecc,pre,perh,xob,ww,dayl)
              x(ilat,idecade,imonth) = x(ilat,idecade,imonth) + real(ww/86.4d0)
           enddo
        enddo
     enddo
     x(:,idecade,:) = x(:,idecade,:)/ix
  enddo
  open(1,file='insol_high.dat',access='stream',form='unformatted')
  write(1) x
  close(1)

end program make_insol
