!gfortran -O3 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -c read_final.f90
module read_final

  use netcdf
  implicit none

  private
  public get_final
  
  integer :: nlon0, nlat0
  real, allocatable :: lons0(:),lats0(:)
  integer, public :: nlonc,nlatc
  real, allocatable, public :: lonsc(:),latsc(:)
  real, allocatable, public :: xc(:,:,:,:) ! lon,lat,var,month

contains

  ! time from = -2200 to 3 (unit == decades)

  subroutine get_final(lon11,lon22,lat1,lat2,decade,vars)

    use constants
    use strings

    real, intent(in) :: lon11, lon22
    real, intent(in) :: lat1, lat2
    integer, intent(in) :: decade
    character(len=*), intent(in) :: vars

    real:: lon1,lon2
    integer :: ncid,did,varid
    integer :: start(4),ncount(4)
    integer :: lon_start,lon_end,lat_start,lat_end,time
    integer :: dmin,dmax
    real :: test41,test42,test4
    real :: sf,ao
    real :: missing_value

    integer :: nvar
    integer, parameter :: maxvar = 30 ! maximum number of variables
    character (len=1000) :: vars2
    character (len=20) :: vars0(maxvar)
    integer :: ivar,ilon,ilat,jj,it

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
    
    call handle_err(nf90_open('/Users/djlorenz/Research/jack/final/CCSM/'// &
         'prcp.nc',nf90_nowrite,ncid))

    ! get longitudes:
    call handle_err(nf90_inq_dimid(ncid,'lon',did))
    call handle_err(nf90_inquire_dimension(ncid,did,len=nlon0))
    allocate(lons0(nlon0))
    call handle_err(nf90_inq_varid(ncid,'lon',varid))
    call handle_err(nf90_get_var(ncid,varid,lons0,[1],[nlon0]))

    ! get latitudes:
    call handle_err(nf90_inq_dimid(ncid,'lat',did))
    call handle_err(nf90_inquire_dimension(ncid,did,len=nlat0))
    allocate(lats0(nlat0))
    call handle_err(nf90_inq_varid(ncid,'lat',varid))
    call handle_err(nf90_get_var(ncid,varid,lats0,[1],[nlat0]))

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
          print*, 'ERROR: get_final: longitude range outside dataset'
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
          print*, 'ERROR: get_final: latitude range outside dataset'
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
    
    ! decade
    if(decade<-2200 .or. decade>3) then
       print*, 'ERROR: get_final: decade out of range'
       print*, 'decade = ',decade
       stop
    endif
    time = decade + 2201

    call handle_err(nf90_close(ncid))


    allocate(xc(nlonc,nlatc,nvar,12))

    do ivar=1,nvar
       if(vars0(ivar)/='tmax' .and. vars0(ivar)/='tmin') then
          call handle_err(nf90_open('/Users/djlorenz/Research/jack/final/CCSM/'// &
               trim(vars0(ivar))//'.nc',nf90_nowrite,ncid))
       else
          call handle_err(nf90_open('/Users/djlorenz/Research/jack/final/CCSM/'// &
               'temp.nc',nf90_nowrite,ncid))
       endif
       start = [lon_start,lat_start,1,time]
       ncount = [nlonc,nlatc,12,1]
       call handle_err(nf90_inq_varid(ncid,trim(vars0(ivar)),varid))
       if(lon_end >= lon_start) then
          call handle_err(nf90_get_var(ncid,varid, &
               xc(:,:,ivar,:),start,ncount))
       else
          start = [lon_start,lat_start,1,time]
          ncount = [nlon0-lon_start+1,nlatc,12,1]
          call handle_err(nf90_get_var(ncid,varid, &
               xc(1:nlon0-lon_start+1,:,ivar,:),start,ncount))
          start = [1,lat_start,1,time]
          ncount = [lon_end,nlatc,12,1]
          call handle_err(nf90_get_var(ncid,varid, &
               xc(nlonc-lon_end+1:nlonc,:,ivar,:),start,ncount))
       endif
       call handle_err(nf90_get_att(ncid,varid,'scale_factor',sf))
       call handle_err(nf90_get_att(ncid,varid,'add_offset',ao))
       call handle_err(nf90_get_att(ncid,varid,'missing_value',missing_value))
       xc(:,:,ivar,:) = merge(xc(:,:,ivar,:)*sf + ao,missing,xc(:,:,ivar,:)/=missing_value)
       call handle_err(nf90_close(ncid))
    enddo

    deallocate(lons0,lats0)

  end subroutine get_final
  
  subroutine handle_err(status)
    implicit none
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
    endif
  end subroutine handle_err

end module read_final

