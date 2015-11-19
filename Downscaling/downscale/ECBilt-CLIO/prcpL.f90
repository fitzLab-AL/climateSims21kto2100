!downscale precipitation
!gfortran -O3 -o prcpL prcpL.f90 read_ECBilt.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program prcp

  use netcdf
  use constants
  use read_data
  use nc_write
  use read_ECBilt
  implicit none

  integer, parameter :: year1 = 1901, year2 = 2011, nyear0 = (year2 - year1 + 1)/10

  ! 3 == 1981-1990, 0 == 1951-1960
  integer, parameter :: decade1 = -2100, decade2 = -1, ndecade = decade2-decade1+1

  ! reference climate decades:
  integer, parameter :: decade_ref1 = -11, decade_ref2 = -1, ndecade_ref = decade_ref2-decade_ref1+1
  
  real, parameter :: minlon = -173.0, maxlon = -48.0
  real, parameter :: minlat = 10.0, maxlat = 80.0

  real, allocatable :: x(:,:,:,:),y(:,:,:),xbox(:,:,:,:,:)
  real, allocatable :: xave(:,:,:),b(:,:,:,:),xcmin(:,:,:)
  real, parameter :: days(4) = [90.25,92.0,92.0,91.0]

  integer :: idecade,isea,ilat,ilon,jyear,iyear

  if(nyear0 /= ndecade_ref) then
     print*, 'ERROR: nyear0 /= ndecade_ref'
     print*, 'nyear0 = ',nyear0 
     print*, 'ndecade_ref = ',ndecade_ref
     stop
  endif

  ! sort "present" of climate model:
  call get_ECBilt(minlon,maxlon,minlat,maxlat,decade_ref1,decade_ref2,'prcp')

!  call fill_miss(nlonc,nlatc,4*ndecade_ref,xc)

  ! sort:
  allocate(xcmin(nlonc,nlatc,4))
  do isea=1,4
     do ilat=1,nlatc
        do ilon=1,nlonc
           call shell0(nyear0,xc(ilon,ilat,1,isea,:))
           xcmin(ilon,ilat,isea) = xc(ilon,ilat,1,isea,1)
        enddo
     enddo
  enddo
  
  ! sort "present" observations:
  allocate(xbox(nlonc,nlatc,1,4,nyear0))
  jyear = 0
  do iyear=year1,year2-9,10
     print*, iyear
     ! get prcp, Tmax Tmin, vap:
     file_expr = '/Users/djlorenz/Research/jack/data/prcp.nc '
     vars = 'prcp '
     call read_nc(minlon,maxlon,.false.,minlat,maxlat,missing,missing, &
          [all,all,all,iyear],[all,all,all,iyear+9],1,'month')
     if(.not.allocated(x)) allocate(x(nlon,nlat,1,4))
     jyear = jyear + 1
     x(:,:,:,1) = sum(31.0*x1(:,:,:,1,:) + 28.25*x1(:,:,:,2,:) + 31.0*x1(:,:,:,12,:),dim=4) &
          /days(1)/10
     x(:,:,:,2) = sum(31.0*x1(:,:,:,3,:) + 30.0*x1(:,:,:,4,:) + 31.0*x1(:,:,:,5,:),dim=4) &
          /days(2)/10
     x(:,:,:,3) = sum(30.0*x1(:,:,:,6,:) + 31.0*x1(:,:,:,7,:) + 31.0*x1(:,:,:,8,:),dim=4) &
          /days(3)/10
     x(:,:,:,4) = sum(30.0*x1(:,:,:,9,:) + 31.0*x1(:,:,:,10,:) + 30.0*x1(:,:,:,11,:),dim=4) &
          /days(4)/10
     x = merge(x,missing,abs(x)<=1e30)
     call box(lons,lats,nlon,nlat,lonsc,latsc,nlonc,nlatc,x,xbox(:,:,:,:,jyear),4,0)
  enddo

  ! sort:
  do isea=1,4
     do ilat=1,nlatc
        do ilon=1,nlonc
           call shell0(nyear0,xbox(ilon,ilat,1,isea,:))
        enddo
     enddo
  enddo
  
  ! find parameters for downscaling:
  allocate(xave(nlonc,nlatc,4),b(nlonc,nlatc,4,0:1))
  b = missing
  xave = sum(xbox(:,:,1,:,:),dim=4)/ndecade_ref
  do isea=1,4
     do ilat=1,nlatc
        do ilon=1,nlonc
           if(xbox(ilon,ilat,1,isea,1)/=missing) then
              call linear1(ndecade_ref,xc(ilon,ilat,1,isea,:),1,xbox(ilon,ilat,1,isea,:), &
                   b(ilon,ilat,isea,:))
           endif
        enddo
     enddo
  enddo

!  call ctlwrite('prcp_check.ctl','prcp',lons=lonsc,lats=latsc)
!  open(1,file='prcp_check.dat',access='stream',form='unformatted')
!  write(1) xave(:,:,1)
!  close(1)
!  print*, count(xave(:,:,1)==0.0); stop

  
  ! find downscale factors
  allocate(y(nlonc,nlatc,4))

  call get_ECBilt(minlon,maxlon,minlat,maxlat,decade1,decade2,'prcp')

  call ctlwrite('prcpL.ctl','prcp',lons=lonsc,lats=latsc, &
       nlev=ndecade,minlev=0.01*decade1,dlev=0.01,nens=4)
  open(1,file='prcpL.dat',access='direct',form='unformatted',recl=4*nlonc*nlatc)

  do idecade=1,ndecade
!     print*, idecade
!     call fill_miss(nlonc,nlatc,4,xc(:,:,1,:,idecade))
     do isea=1,4
        do ilat=1,nlatc
           do ilon=1,nlonc
              if(b(ilon,ilat,isea,0)/=missing) then
                 if(b(ilon,ilat,isea,0)>=0. .or. &
                      xc(ilon,ilat,1,isea,idecade) >= xcmin(ilon,ilat,isea)) then
                    y(ilon,ilat,isea) = b(ilon,ilat,isea,0) + &
                         b(ilon,ilat,isea,1)*xc(ilon,ilat,1,isea,idecade)
                 else
                    y(ilon,ilat,isea) = (b(ilon,ilat,isea,0) + &
                         b(ilon,ilat,isea,1)*xcmin(ilon,ilat,isea))/xcmin(ilon,ilat,isea)* &
                         xc(ilon,ilat,1,isea,idecade)
                 endif
              else
                 y(ilon,ilat,isea) = missing
              endif
           enddo
        enddo
     enddo
     do isea=1,4
        y(:,:,isea) = merge(merge(y(:,:,isea)/xave(:,:,isea),missing, &
             y(:,:,isea)/=missing.and.xave(:,:,isea)/=missing),1.0,xave(:,:,isea)/=0.)
     enddo
     ! SPECIAL: ridiculous increase seen at one point in summer, set equal to 
     !    average of surrounding grid points
     y(14,4,3) = 0.25*(y(13,4,3) + y(15,4,3) + y(14,3,3) + y(14,5,3))
     do isea=1,4
        if(any(isnan(y(:,:,isea)))) then
           print*, 'NaN ',isea
        endif
        write(1,rec=ndecade*(isea-1)+idecade) y(:,:,isea)
     enddo
  enddo

  close(1)

end program prcp

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
