!gfortran -O3 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -c read_CCSM.f90
module read_CCSM

  use netcdf
  implicit none

  private
  public get_CCSM
  
  integer, parameter :: nlon0=96, nlat0=48,nt0=2204
  real :: lons0(nlon0),lats0(nlat0),ts0r(nt0)
  integer :: ts0(nt0)
  integer, public :: nlonc,nlatc,ntc
  real, allocatable, public :: lonsc(:),latsc(:)
  real, allocatable, public :: xc(:,:,:,:,:) ! lon,lat,var,season,year

contains

  ! time from = -2200 to 3 (unit == decades)

  subroutine get_CCSM(lon11,lon22,lat1,lat2,time1,time2,vars)

    use constants
    use strings

    real, intent(in) :: lon11, lon22
    real, intent(in) :: lat1, lat2
    integer, intent(in) :: time1,time2
    character(len=*), intent(in) :: vars

    real:: lon1,lon2
    integer :: ncid(3)
    integer :: start(3),ncount(3)
    integer :: lon_start,lon_end,lat_start,lat_end,t_start,t_end
    integer :: dmin,dmax
    real :: test41,test42,test4

    integer :: nvar
    integer, parameter :: maxvar = 30 ! maximum number of variables
    character (len=1000) :: vars2
    character (len=20) :: vars0(maxvar)
    integer :: filen(maxvar),varn(maxvar)
    integer :: ivar,ilon,ilat,jj,it,isea


    if(lat1/=missing.and.lat2/=missing.and.lat2<lat1) then
       print*, 'ERROR: lat2 < lat1'
       stop
    endif
    lon1 = lon11
    lon2 = lon22
    if(lon1/=missing) then
       if(lon1<0.) then
          lon1 = lon1 + 360.
       elseif(lon1>=360.) then
          lon1 = lon1 - 360.
       endif
    endif
    if(lon2/=missing) then
       if(lon2<0.) then
          lon2 = lon2 + 360.
       elseif(lon1>=360.) then
          lon2 = lon2 - 360.
       endif
    endif
    
    ! deallocate arrays if necessary
    if(allocated(xc)) then
       deallocate(xc)
    endif
    if(allocated(lonsc)) then
       deallocate(lonsc)
    endif
    if(allocated(latsc)) then
       deallocate(latsc)
    endif

    vars2 = vars
    call parse(vars2,' ',vars0,nvar)
    if(nvar == 0) then
       print*, 'ERROR: get_CCSM: the string "vars" is empty'
       stop
    endif

    do ivar=1,nvar
       select case(vars0(ivar))
       case('FLNS')
          filen(ivar) = 1
          varn(ivar) = 1
       case('FSDS')
          filen(ivar) = 1
          varn(ivar) = 2
       case('TSMN')
          filen(ivar) = 1
          varn(ivar) = 3
       case('PRECT')
          filen(ivar) = 1
          varn(ivar) = 7
       case('AGDD0')
          filen(ivar) = 2
          varn(ivar) = 1
       case('AGDD5')
          filen(ivar) = 2
          varn(ivar) = 2
       case('FLDS')
          filen(ivar) = 2
          varn(ivar) = 3
       case('Q2M')
          filen(ivar) = 2
          varn(ivar) = 4
       case('TDA')
          filen(ivar) = 2
          varn(ivar) = 5
       case('TREFMNAV')
          filen(ivar) = 2
          varn(ivar) = 6
       case('TREFMXAV')
          filen(ivar) = 2
          varn(ivar) = 7
       case('WIND')
          filen(ivar) = 2
          varn(ivar) = 8
       case('PS')
          filen(ivar) = 3
          varn(ivar) = 1
       end select
    enddo
    
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.djf.Dave.nc',nf90_nowrite,ncid(1)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.clm2.ncrcat.djf.Dave.nc',nf90_nowrite,ncid(2)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.djf.PS.nc',nf90_nowrite,ncid(3)))
    call handle_err(nf90_get_var(ncid(1),5,lons0,[1],[nlon0]))
    call handle_err(nf90_get_var(ncid(1),4,lats0,[1],[nlat0]))
    call handle_err(nf90_get_var(ncid(1),6,ts0r,[1],[nt0]))
    ts0 = nint(100.*ts0r)
    
    ! if lons == missing then do all longitudes
    if(lon1==missing .or. lon2==missing) then
       dmin = 1; dmax = nlon0
    else
       ! make lons between 0 & 360
       do ilon=1,nlon0
          if(lons0(ilon)<0.) then
             lons0(ilon) = lons0(ilon) + 360.
          elseif(lons0(ilon)>=360.) then
             lons0(ilon) = lons0(ilon) - 360.
          endif
       enddo
       test42 = -huge(0.0)
       test41 = huge(0.0)
       dmax = 0; dmin = 0
       do ilon=1,nlon0
          if(lons0(ilon)<=lon2) then
             if(lons0(ilon)>test42) then
                test42 = lons0(ilon)
                dmax = ilon
             endif
          endif
          if(lons0(ilon)>=lon1) then
             if(lons0(ilon)<test41) then
                test41 = lons0(ilon)
                dmin = ilon
             endif
          endif
       enddo
       if(dmin==0 .or. dmax==0) then
          print*, 'ERROR: read_CCSM: longitude range outside dataset'
          stop
       endif
    endif
    if(dmax>=dmin) then
       nlonc = dmax - dmin + 1
    else
       nlonc = (nlon0 - dmin + 1) + dmax
    endif
    allocate(lonsc(nlonc))
    test4 = 0.
    if(lons0(dmin)>=180.) then
       test4 = -360.
    endif
    if(dmax>=dmin) then
       do ilon=dmin,dmax
          if(lons0(ilon)<180.) then
             lonsc(ilon - dmin + 1) = lons0(ilon)
          else
             lonsc(ilon - dmin + 1) = lons0(ilon) + test4
          endif
       enddo
    else
       do ilon=dmin,nlon0
          lonsc(ilon - dmin + 1) = lons0(ilon) - 360.
       enddo
       do ilon=1,dmax
          lonsc(nlon0 - dmin + 1 + ilon) = lons0(ilon)
       enddo
    endif
    lon_start=dmin
    lon_end=dmax
    
    ! latitudes
    if(lat1==missing .or. lat2==missing) then
       dmin = 1; dmax = nlat0
    else
       test42 = -huge(0.0)
       test41 = huge(0.0)
       dmax = 0; dmin = 0
       do ilat=1,nlat0
          if(lats0(ilat)<=lat2) then
             if(lats0(ilat)>test42) then
                test42 = lats0(ilat)
                dmax = ilat
             endif
          endif
          if(lats0(ilat)>=lat1) then
             if(lats0(ilat)<test41) then
                test41 = lats0(ilat)
                dmin = ilat
             endif
          endif
       enddo
       if(dmin==0 .or. dmax==0) then
          print*, 'ERROR: read_CCSM: latitude range outside dataset'
          stop
       endif
    endif
    if(dmax<dmin) then
       jj=dmin
       dmin=dmax
       dmax=jj
    endif
    nlatc = dmax - dmin + 1
    allocate(latsc(nlatc))
    latsc = lats0(dmin:dmax)
    lat_start=dmin
    lat_end=dmax
    
    ! times
    ntc = time2 - time1 + 1
    if(ntc < 1) then
       print*, 'ERROR: read_CCSM: time2 < time1'
       stop
    endif
    t_start = time1 + 2200 + 1
    if(t_start < 1 .or. t_start > nt0) then
       print*, 'ERROR: read_CCSM: time1 out of range'
       stop
    endif
    t_end = time2 + 2200 + 1
    if(t_end < 1 .or. t_end > nt0) then
       print*, 'ERROR: read_CCSM: time2 out of range'
       stop
    endif

    allocate(xc(nlonc,nlatc,nvar,4,ntc))

    do it=t_start,t_end
       do ivar=1,nvar
          start = [lon_start,lat_start,it]
          ncount = [nlonc,nlatc,1]
          if(lon_end >= lon_start) then
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(:,:,ivar,1,it-t_start+1),start,ncount))
          else
             start = [lon_start,lat_start,it]
             ncount = [nlon0-lon_start+1,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(1:nlon0-lon_start+1,:,ivar,1,it-t_start+1),start,ncount))
             start = [1,lat_start,it]
             ncount = [lon_end,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(nlonc-lon_end+1:nlonc,:,ivar,1,it-t_start+1),start,ncount))
          endif
       enddo
    enddo
    call handle_err(nf90_close(ncid(2)))
    call handle_err(nf90_close(ncid(1)))

    ! file numbers change for CLM:
    do ivar=1,nvar
       if(filen(ivar)==2) then
          varn(ivar) = varn(ivar) + 3
       endif
    enddo

    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.mam.Dave.nc',nf90_nowrite,ncid(1)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.clm2.ncrcat.mam.Dave.nc',nf90_nowrite,ncid(2)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.mam.PS.nc',nf90_nowrite,ncid(3)))
    do it=t_start,t_end
       do ivar=1,nvar
          start = [lon_start,lat_start,it]
          ncount = [nlonc,nlatc,1]
          if(lon_end >= lon_start) then
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(:,:,ivar,2,it-t_start+1),start,ncount))
          else
             start = [lon_start,lat_start,it]
             ncount = [nlon0-lon_start+1,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(1:nlon0-lon_start+1,:,ivar,2,it-t_start+1),start,ncount))
             start = [1,lat_start,it]
             ncount = [lon_end,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(nlonc-lon_end+1:nlonc,:,ivar,2,it-t_start+1),start,ncount))
          endif
       enddo
    enddo
    call handle_err(nf90_close(ncid(2)))
    call handle_err(nf90_close(ncid(1)))

    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.jja.Dave.nc',nf90_nowrite,ncid(1)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.clm2.ncrcat.jja.Dave.nc',nf90_nowrite,ncid(2)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.jja.PS.nc',nf90_nowrite,ncid(3)))
    do it=t_start,t_end
       do ivar=1,nvar
          start = [lon_start,lat_start,it]
          ncount = [nlonc,nlatc,1]
          if(lon_end >= lon_start) then
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(:,:,ivar,3,it-t_start+1),start,ncount))
          else
             start = [lon_start,lat_start,it]
             ncount = [nlon0-lon_start+1,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(1:nlon0-lon_start+1,:,ivar,3,it-t_start+1),start,ncount))
             start = [1,lat_start,it]
             ncount = [lon_end,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(nlonc-lon_end+1:nlonc,:,ivar,3,it-t_start+1),start,ncount))
          endif
       enddo
    enddo
    call handle_err(nf90_close(ncid(2)))
    call handle_err(nf90_close(ncid(1)))

    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.son.Dave.nc',nf90_nowrite,ncid(1)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.clm2.ncrcat.son.Dave.nc',nf90_nowrite,ncid(2)))
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/CCSM/special/'// &
         'b30.00_4kaDVTd.cam2.ncrcat.son.PS.nc',nf90_nowrite,ncid(3)))
    do it=t_start,t_end
       do ivar=1,nvar
          start = [lon_start,lat_start,it]
          ncount = [nlonc,nlatc,1]
          if(lon_end >= lon_start) then
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(:,:,ivar,4,it-t_start+1),start,ncount))
          else
             start = [lon_start,lat_start,it]
             ncount = [nlon0-lon_start+1,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(1:nlon0-lon_start+1,:,ivar,4,it-t_start+1),start,ncount))
             start = [1,lat_start,it]
             ncount = [lon_end,nlatc,1]
             call handle_err(nf90_get_var(ncid(filen(ivar)),varn(ivar), &
                  xc(nlonc-lon_end+1:nlonc,:,ivar,4,it-t_start+1),start,ncount))
          endif
       enddo
    enddo
    call handle_err(nf90_close(ncid(2)))
    call handle_err(nf90_close(ncid(1)))

    xc = merge(xc,missing,abs(xc)<=1e30)

    do ivar=1,nvar
       if(vars0(ivar)(1:1) == 'T') then
          xc(:,:,ivar,:,:) = merge(xc(:,:,ivar,:,:)-273.15,missing, &
               xc(:,:,ivar,:,:) /= missing)
       elseif(vars0(ivar) == 'PRECT') then
          xc(:,:,ivar,1,:) = merge(xc(:,:,ivar,1,:)*86400.*1000.*30.0833,missing, &
               xc(:,:,ivar,1,:) /= missing)
          xc(:,:,ivar,2,:) = merge(xc(:,:,ivar,2,:)*86400.*1000.*30.6667,missing, &
               xc(:,:,ivar,2,:) /= missing)
          xc(:,:,ivar,3,:) = merge(xc(:,:,ivar,3,:)*86400.*1000.*30.6667,missing, &
               xc(:,:,ivar,3,:) /= missing)
          xc(:,:,ivar,4,:) = merge(xc(:,:,ivar,4,:)*86400.*1000.*30.3333,missing, &
               xc(:,:,ivar,4,:) /= missing)
       endif
    enddo

  end subroutine get_CCSM
  
  subroutine handle_err(status)
    implicit none
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
    endif
  end subroutine handle_err

end module read_CCSM

