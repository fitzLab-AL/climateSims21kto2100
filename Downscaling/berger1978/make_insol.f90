!make insolation dataset (latitude, season and decade)
!gfortran -O3 -o make_insol insol2.o read_CCSM.o make_insol.f90 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program make_insol

  use netcdf
  use constants
  use nc_write
  use read_CCSM
  implicit none

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1

  integer, parameter :: nlat=48
  real :: lats(nlat)
  real :: kyear
  integer, parameter :: days(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  real :: x(nlat,decade1:decade2,4)
  integer :: ix(nlat,4)

  integer :: ilat,idecade,isea,imonth,month,day
  real*8 :: phi,t,ecc,pre,perh,xob,ww,dayl

  ! get CCSM grid:
  call get_CCSM(missing,missing,missing,missing,0,0,'PS')
  if(nlatc /= nlat) then
     print*, 'ERROR: nlatc = ',nlatc
     stop
  endif

  lats = latsc

  call ctlwrite('insol.ctl','i',lats=lats,nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=4)

  x = 0.0
  do idecade=decade1,decade2
     print*, idecade
     kyear = 0.01*(idecade + 0.5)
     t = dble(kyear)
     ix = 0
     do isea=1,4
        do imonth=3*isea-3,3*isea-1
           month = modulo(imonth-1,12) + 1
           do day=1,days(month)
              do ilat=1,nlat
                 phi = dble(lats(ilat))
                 ix(ilat,isea) = ix(ilat,isea) + 1
                 call insol(phi,t,month,day,ecc,pre,perh,xob,ww,dayl)
                 x(ilat,idecade,isea) = x(ilat,idecade,isea) + real(ww/86.4d0)
              enddo
           enddo
        enddo
     enddo
     x(:,idecade,:) = x(:,idecade,:)/ix
  enddo
  open(1,file='insol.dat',access='stream',form='unformatted')
  write(1) x
  close(1)

end program make_insol
