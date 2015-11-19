!fix shorelines
!
! VARIABLES:
! tmax tmin gdd0 gdd5 pet aet prcp vap  swd  swu  lwn  (ps only one time step)
! temp.nc   gdd.nc    ET.nc
!
!gfortran -O3 -o shore shore.f90 read_final.o write_final.o -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -L/Users/djlorenz/Research/netcdf/lib -L/Users/djlorenz/lib -lnetcdff -lnetcdf  -llapack -lblas -lmylib
program shore

  use netcdf
  use constants
  use read_final
  use write_final
  implicit none

  real, parameter :: pow = -1., wlat = 0.75, wdis = 1.0 - wlat, fac = 1.0/sqrt(wdis)
  real, parameter :: Re = 6371.

  integer, parameter :: decade1 = -2100, decade2 = -1, ndecade = decade2-decade1+1
  integer, parameter :: nvar = 11
  !                               temp        gdd        ET    prcp  vap     sw      lwn 
  real, parameter :: sf(nvar) = [0.01,0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.1, 0.1, 0.01]
  real, parameter :: ao(nvar) = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0]

  real, allocatable :: land(:,:,:)
  integer*2, allocatable :: ix(:,:,:,:)

  real :: dis,dislat,dlat,dmax2
  real :: w,wsum,xsum,temp
  integer :: iland,idecade,imonth,ivar,ilat,ilon,jlat,jlon,irlat,irlon,j

  logical :: first = .true.
  logical :: done

  do idecade=decade1,decade2
     print*, idecade
     
     call  get_final(missing,missing,missing,missing,idecade, &
          'tmax tmin gdd0 gdd5 pet aet prcp vap swd swu lwn')

     if(first) then
        allocate(land(nlonc,nlatc,21),ix(nlonc,nlatc,12,nvar))
        open(1,file='/Users/djlorenz/Research/jack/shorelines/Raster/binary.dat', &
             access='stream',form='unformatted')
        read(1) land
        close(1)
        first = .false.

        dlat = (latsc(2) - latsc(1))*pi/180.

        call ncopen(1,'/Users/djlorenz/Research/jack/finalshore/ECBilt/temp.nc','tmax tmin', &
             vunits='C',vnames='mean maximum daily temperature; mean minimum daily temperature;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=sf(1:2),ao=ao(1:2))
        call ncopen(2,'/Users/djlorenz/Research/jack/finalshore/ECBilt/gdd.nc','gdd0 gdd5', &
             vunits='C',vnames='growing degree days base 0C; growing degree days base 5C;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=sf(3:4),ao=ao(3:4))
        call ncopen(3,'/Users/djlorenz/Research/jack/finalshore/ECBilt/ET.nc','pet aet', &
             vunits='mm',vnames='potential evapotranspiration; actual evapotranspiration;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=sf(5:6),ao=ao(5:6))
        call ncopen(4,'/Users/djlorenz/Research/jack/finalshore/ECBilt/prcp.nc','prcp ', &
             vunits='mm',vnames='mean monthly precipitation accumulation;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=[sf(7)],ao=[ao(7)])
        call ncopen(5,'/Users/djlorenz/Research/jack/finalshore/ECBilt/vap.nc','vap ', &
             vunits='hPa',vnames='mean monthly vapor pressure;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=[sf(8)],ao=[ao(8)])
        call ncopen(6,'/Users/djlorenz/Research/jack/finalshore/ECBilt/swd.nc','swd ', &
             vunits='Wm-2',vnames='mean monthly downward shortwave radiation at the surface;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=[sf(9)],ao=[ao(9)])
        call ncopen(7,'/Users/djlorenz/Research/jack/finalshore/ECBilt/swu.nc','swu ', &
             vunits='Wm-2',vnames='mean monthly upward shortwave radiation at the surface;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=[sf(10)],ao=[ao(10)])
        call ncopen(8,'/Users/djlorenz/Research/jack/finalshore/ECBilt/lwn.nc','lwn ', &
             vunits='Wm-2',vnames='mean monthly net longwave radiation at the surface;', &
             lons=lonsc,lats=latsc,title='Downscaled ECBilt runs. Monthly means by decade.', &
             deflate=1,sf=[sf(11)],ao=[ao(11)])
     endif

     iland = min(nint(-real(idecade)/100.),21)
!     print*, idecade,iland

     if(iland/=0) then
        do imonth=1,12
           do ivar=1,nvar
              do ilat=1,nlatc
                 do ilon=1,nlonc
                    if(land(ilon,ilat,iland)==0. .and. &
                         xc(ilon,ilat,ivar,imonth)==missing) then
                       done = .false.
                       wsum = 0.; xsum = 0.
                       do irlat=1,nlatc
                          dmax2 = ((real(irlat)+0.5)*dlat*Re)**2
                          ! reach for more grid points in longitude because
                          !  longitudes are closer together: 
                          irlon = ceiling(fac*real(irlat)/cos(pi*latsc(ilat)/180.))
                          do jlat=max(ilat-irlat,1),min(ilat+irlat,nlatc)
                             do jlon=max(ilon-irlon,1),min(ilon+irlon,nlonc)
                                if(xc(jlon,jlat,ivar,imonth)/=missing) then
                                   call distance(lonsc(ilon),latsc(ilat), &
                                        lonsc(jlon),latsc(jlat),dis,dislat)
                                   w = (wdis*dis**2 + wlat*dislat**2)
                                   if(w <= dmax2) then
                                      done = .true.
                                      w = w**pow
                                      wsum = wsum + w
                                      xsum = xsum + w*xc(jlon,jlat,ivar,imonth)
                                   endif
                                endif
                             enddo
                          enddo
                          if(done) then
                             temp = xsum/wsum
                             ix(ilon,ilat,imonth,ivar) = nint((temp - ao(ivar))/sf(ivar))
                             exit
                          endif
                       enddo
                    elseif(land(ilon,ilat,iland)/=0. .and. &
                         xc(ilon,ilat,ivar,imonth)/=missing) then
                       ix(ilon,ilat,imonth,ivar) = missing_short
                    else
                       temp = xc(ilon,ilat,ivar,imonth)
                       if(temp/=missing) then
                          ix(ilon,ilat,imonth,ivar) = nint((temp - ao(ivar))/sf(ivar))
                       else
                          ix(ilon,ilat,imonth,ivar) = missing_short
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
     else
        do imonth=1,12
           do ivar=1,nvar
              do ilat=1,nlatc
                 do ilon=1,nlonc
                    temp = xc(ilon,ilat,ivar,imonth)
                    if(temp/=missing) then
                       ix(ilon,ilat,imonth,ivar) = nint((temp - ao(ivar))/sf(ivar))
                    else
                       ix(ilon,ilat,imonth,ivar) = missing_short
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif
   
     call ncwrite(1,im=ix(:,:,:,1:2))
     call ncwrite(2,im=ix(:,:,:,3:4))
     call ncwrite(3,im=ix(:,:,:,5:6))
     call ncwrite(4,im=ix(:,:,:,7))
     call ncwrite(5,im=ix(:,:,:,8))
     call ncwrite(6,im=ix(:,:,:,9))
     call ncwrite(7,im=ix(:,:,:,10))
     call ncwrite(8,im=ix(:,:,:,11))

  enddo
     
  do j=1,8
     call ncclose(j)
  enddo

end program shore

subroutine distance(lon1,lat1,lon2,lat2,dis,dislat)

  implicit none

  real :: dis,dislat,lon1,lon2,lat1,lat2
  real*8 :: x

  real*8, parameter :: pi=3.14159265358979323846
  real*8, parameter :: rad=pi/180.d0
  real, parameter :: Re = 6371.

  real*8 :: c1,c2,s1,s2

  c1 = cos(rad*dble(lat1))
  c2 = cos(rad*dble(lat2))
  s1 = sin(rad*dble(lat1))
  s2 = sin(rad*dble(lat2))
  
  x = c1*c2*cos(rad*(lon2 - lon1)) + s1*s2
  dis = Re*real(acos(x))
  dislat = Re*real(rad*abs(dble(lat2)-dble(lat1)))

end subroutine distance
