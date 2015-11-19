!downscale RELATIVE downward shortwave radiation at surface
!gfortran -O3 -o swd swd.f90 read_CCSM.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program swd

  use netcdf
  use constants
  use nc_write
  use read_CCSM
  implicit none

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2200, decade2 = 3, ndecade = decade2-decade1+1

  ! reference climate decades:
  integer, parameter :: decade_ref1 = -7, decade_ref2 = 3, ndecade_ref = decade_ref2-decade_ref1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:,:),xave(:,:,:),y(:,:,:)

  integer, parameter :: nlatT = 48
  integer :: lat0
  real :: latsT(nlatT),insol(nlatT,ndecade,4)

  integer :: idecade,isea,jdecade,ilat,ilon

  ! get grid
  call get_CCSM(missing,missing,missing,missing,0,0,'FSDS')
  if(nlatT /= nlatc) then
     print*, 'ERROR: nlatc = ', nlatc
     stop
  endif
  latsT = latsc
  do ilat=1,nlatT
     if(latsT(ilat) >= minlat) exit
     lat0 = ilat
  enddo
  print*, latsT(lat0),latsT(lat0+1)

  open(1,file='../../berger1978/insol.dat',access='stream',form='unformatted')
  read(1) insol
  close(1)

  call get_CCSM(minlon,maxlon,minlat,maxlat,-7,3,'FSDS')
  do idecade = 1,ndecade_ref
     jdecade = idecade - ndecade_ref + ndecade
     print*, jdecade
     do isea=1,4
        do ilat=1,nlatc
           do ilon=1,nlonc
              xc(ilon,ilat,1,isea,idecade) = xc(ilon,ilat,1,isea,idecade)/ &
                   insol(ilat+lat0,jdecade,isea)
           enddo
        enddo
     enddo
  enddo
  call fill_miss(nlonc,nlatc,4*ndecade_ref,xc)

  allocate(x(nlonc,nlatc,4,ndecade_ref),xave(nlonc,nlatc,4),y(nlonc,nlatc,4))

  x = xc(:,:,1,:,:)
  xave = sum(x,dim=4)/ndecade_ref
  xave = merge(xave,missing,abs(xave)<1e30)
  deallocate(x)

  call get_CCSM(minlon,maxlon,minlat,maxlat,decade1,decade2,'FSDS')
  do idecade = 1,ndecade
     do isea=1,4
        do ilat=1,nlatc
           do ilon=1,nlonc
              xc(ilon,ilat,1,isea,idecade) = xc(ilon,ilat,1,isea,idecade)/ &
                   insol(ilat+lat0,idecade,isea)
           enddo
        enddo
     enddo
  enddo

  call ctlwrite('swd.ctl','swd',lons=lonsc,lats=latsc, &
       nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=4)
  open(1,file='swd.dat',access='direct',form='unformatted',recl=4*nlonc*nlatc)

  do idecade=1,ndecade
     print*, idecade
     call fill_miss(nlonc,nlatc,4,xc(:,:,1,:,idecade))
     y = xc(:,:,1,:,idecade)
     do isea=1,4
        y(:,:,isea) = merge(log(y(:,:,isea))/log(xave(:,:,isea)),missing, &
             y(:,:,isea)/=missing.and.xave(:,:,isea)/=missing)
     enddo
     do isea=1,4
        if(any(isnan(y(:,:,isea)))) then
           print*, 'NaN ',isea
        endif
        write(1,rec=ndecade*(isea-1)+idecade) y(:,:,isea)
     enddo
  enddo

  close(1)

end program swd

! assumes NOT periodic in lon
subroutine fill_miss(nlon,nlat,nsam,x)
  
  use constants

  integer :: nlon,nlat
  real :: x(nlon,nlat,nsam),x2(nlon,nlat)

  real :: xsum
  integer :: ilon,ilat,isam,jlon,jlat,isum

  do isam=1,nsam
     do ilat=1,nlat
        do ilon=1,nlon
           if(x(ilon,ilat,isam)==missing) then
              isum = 0; xsum = 0.0
              do jlat=max(ilat-1,1),min(ilat+1,nlat)
                 do jlon=max(ilon-1,1),min(ilon+1,nlon)
                    if(x(jlon,jlat,isam)/=missing) then
                       isum = isum + 1
                       xsum = xsum + x(jlon,jlat,isam)
                    endif
                 enddo
              enddo
              if(isum/=0) then
                 x2(ilon,ilat) = xsum/isum
              else
                 x2(ilon,ilat) = missing
              endif
           else
              x2(ilon,ilat) = x(ilon,ilat,isam)
           endif
        enddo
     enddo
     x(:,:,isam) = x2
  enddo

end subroutine fill_miss
