!downscale temperature
!gfortran -O3 -o temp temp.f90 read_CCSM.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program temp

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

  real, allocatable :: xave(:,:,:,:),y(:,:,:,:)

  integer :: idecade,isea,ivar

  call get_CCSM(minlon,maxlon,minlat,maxlat,-7,3,'TREFMXAV TREFMNAV')
  ! 2 for 2 variables:
  call fill_miss(nlonc,nlatc,2*4*ndecade_ref,xc)

  allocate(xave(nlonc,nlatc,2,4),y(nlonc,nlatc,2,4))

  xave = sum(xc,dim=5)/ndecade_ref
  xave = merge(xave,missing,abs(xave)<1e30)
  deallocate(xc)

  call get_CCSM(minlon,maxlon,minlat,maxlat,decade1,decade2,'TREFMXAV TREFMNAV')

  call ctlwrite('temp.ctl','tmax tmin',lons=lonsc,lats=latsc, &
       nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=4)
  open(1,file='temp.dat',access='direct',form='unformatted',recl=4*nlonc*nlatc)

  do idecade=1,ndecade
!     print*, idecade
     ! 2 for 2 variables:
     call fill_miss(nlonc,nlatc,2*4,xc(:,:,:,:,idecade))
     do isea=1,4
        do ivar=1,2
           y(:,:,ivar,isea) = merge(xc(:,:,ivar,isea,idecade)-xave(:,:,ivar,isea),missing, &
                xc(:,:,ivar,isea,idecade)/=missing.and.xave(:,:,ivar,isea)/=missing)
        enddo
     enddo
     do isea=1,4
        do ivar=1,2
           if(any(isnan(y(:,:,ivar,isea)))) then
              print*, 'NaN ',ivar,isea
           endif
           ! 2 for 2 variables:
           write(1,rec=2*ndecade*(isea-1)+ndecade*(ivar-1)+idecade) y(:,:,ivar,isea)
        enddo
     enddo
  enddo

  close(1)

end program temp

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
