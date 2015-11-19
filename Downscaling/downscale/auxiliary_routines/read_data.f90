module read_data
  
  use netcdf
  implicit none

  private
  public read_nc, read_no_time
  ! character string for files to find:
  integer, parameter :: lfile_expr = 1024
  character (len=lfile_expr), public :: file_expr
  ! character string for variables to get
  integer, parameter :: lvars = 100
  character (len=lvars), public :: vars
  ! dimensions of arrays for reading data:
  integer, public :: nlon, nlat, nlev, nt, nyear
  integer, public :: nvarx1,nvarxm
  ! nt is number of times within one year (so total times = nt*nyear)
  real, allocatable, public :: lons(:),lats(:),levs(:)
  ! order in arrays = lons, lats, levs, vars, ts, year
  real, allocatable, public :: x1(:,:,:,:,:)   ! one vertical level variables
  real, allocatable, public :: xm(:,:,:,:,:,:) ! multiple vertical levels vars
  ! for read_no_time:
  real, allocatable, public :: xnt1(:,:,:) ! no time vars 
  real, allocatable, public :: xntm(:,:,:,:) ! no time vars 

contains
  
  subroutine read_nc(lon11,lon22,ave_lon, lat1,lat2, lev1,lev2, &
       time1,time2,n_ave,unit_in,continuous)
    
    use constants

    character (len=lvars) :: vars1
    real, intent(in) :: lon11, lon22
    real, intent(in) :: lat1, lat2
    real, intent(in) :: lev1, lev2
    logical, intent(in) :: ave_lon
    integer, intent(in) :: time1(4),time2(4)
    integer, intent(in) :: n_ave
    logical, optional :: continuous ! continuous in time = .true.
                                    ! broken up into nt & nyear = .false. or missing 
    logical :: cont
    character (*) :: unit_in    
    character (len=8)  :: unit
    ! nt is number of times within one year (so total times = nt*nyear)
    real :: lon1, lon2
    ! put files & variables in these:
    integer :: nfile,nvar
    character (len=200), allocatable :: files(:)
    character (len=30), allocatable :: avars(:), varsx1(:), varsxm(:)
    !! for netcdf:
    integer, allocatable :: id(:),types(:),atts(:),vndims(:)
    real, allocatable :: scalef(:),add(:)
    integer :: ncid,ndims,nvars,natts,recdim,n,itype
    integer :: attype,attlen
    integer*2, allocatable :: miss3(:)
    real*4, allocatable :: miss5(:)
    integer :: dims(32) ! 32 hard wired into code
    integer, allocatable :: vdims(:,:)
    character (len=30) :: name,namelon,namelat,namelev,namet,attnam
    character (len=100) :: calendar,calendar0
    integer :: start(4),ncount(4)
    integer :: lon_start,lon_end,lat_start,lat_end,lev_start,lev_end
    integer :: dlonid,dlatid,dlevid,dtid
    integer :: lonid=-1,latid=-1,levid=-1,tid=-1,tbid=-1
    integer :: lontype,lattype,levtype,ttype,tbtype
    integer :: mlon=-1,mlat=-1,mlev=-1,mt=-1 
    integer :: tatts
    integer*2, allocatable :: data3(:,:,:,:)
    integer, allocatable :: levs4(:),times4(:),tb4(:,:)
    real,  allocatable :: data5(:,:,:,:),lons5(:),lats5(:),levs5(:),times5(:),tb5(:,:)
    real*8, allocatable :: lons6(:),lats6(:),levs6(:),times6(:),tb6(:,:)
    ! to test that all variables & files are on the same grid:
    real, allocatable :: lons_test(:),lats_test(:),levs_test(:)
    integer :: dmin,dmax,lev_flag
    real :: test41,test42,test4
    real*8 :: temp8
    !! for times:
    integer, allocatable :: time_out(:,:,:)
    character (len=100) :: time_units
    integer :: hour_day,year_cal,month_cal,day_cal,hour_cal,hour_diff
    integer :: n0,n1,nm,nm0,tmap,ntime
    integer, allocatable :: tt(:),n_ave_x1(:,:,:),n_ave_xm(:,:,:)
    integer :: ivar,dim
    !! for testing
    integer hour1,day1,month1,year1,hour2,day2,month2,year2
    
    integer :: ilon,ilat,ilev
    integer :: i,j,ii,jj,k,l,nn,it,jt,i1,i2,ifile
    
    character (len=8), save :: date=' '
    character (len=10), save :: time=' '

    if(date==' '.and.time==' ') then
       call date_and_time(date,time)
    endif

    if( present(continuous) ) then
       if(continuous) then
          cont = .true.
       else
          cont = .false.
       endif
    else
       cont = .false.
    endif

    unit=unit_in
    vars1=vars
    if(lat1/=missing.and.lat2/=missing.and.lat2<lat1) then
       print*, 'ERROR: lat2 < lat1'
       stop
    endif
    if(lev1/=missing.and.lev2/=missing.and.lev2<lev1) then
       print*, 'ERROR: lev2 < lev1'
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
    if(allocated(x1)) then
       deallocate(x1)
    endif
    if(allocated(xm)) then
       deallocate(xm)
    endif
    if(allocated(lons)) then
       deallocate(lons)
    endif
    if(allocated(lats)) then
       deallocate(lats)
    endif
    if(allocated(levs)) then
       deallocate(levs)
    endif
    
    ! determine # variables & put variables in character array:
    nvar = 0
    j = 0
    do i=1,lvars
       if(vars1(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          nvar = nvar + 1
          j = 0
       endif
    enddo
    if(nvar==0) then
       print*, 'ERROR: no variables entered'
       stop
    endif
    allocate(avars(nvar))
    allocate(id(nvar))
    allocate(types(nvar))
    allocate(atts(nvar))
    allocate(vndims(nvar))
    allocate(vdims(32,nvar)) ! 32 hard wired into code
    allocate(scalef(nvar))
    allocate(add(nvar))
    allocate(miss3(nvar))
    allocate(miss5(nvar))
    miss3=-32768
    miss5=missing
    vndims = -1
    do i=1,nvar
       vars1 = adjustl(vars1)
       l = index(vars1,' ') - 1
       avars(i) = vars1(1:l)
       vars1(1:l) = ' '
    enddo
    ! put list of files needed in "temp"
    call system('ls '//file_expr//' > temp'//date//time)
    ! find number of lines in "temp":
    call system('wc temp'//date//time//' > tempn'//date//time)
    open(98,file='tempn'//date//time)
    read(98,*) nfile
    close(98)
    call system('rm tempn'//date//time)
    
    ! get file names from "temp":
    allocate( files(nfile) )
    open(98,file='temp'//date//time)
    do i=1,nfile
       read(98,'(a)') files(i)
    enddo
    close(98)
    call system('rm temp'//date//time)
    
    ! find dimension of all variables first:
    ! (go through files in a way that might minimize the # files opened)
    ! also get the calendar
    id = -1
    if(mod(nfile,nvar)==0) then
       nn = nfile/nvar
    else
       nn = nfile/nvar + 1
    endif
    calendar = ' '
    do i=1,nn
       do j=0,nvar - 1
          k = nn*j + i
          if(k<=nfile) then
             !open file
             call handle_err(nf90_open(trim(files(k)),nf90_nowrite,ncid))
             ! get # variables:
             call handle_err(nf90_inquire(ncid, ndims, nvars, natts, recdim))
             ! get variable info:
             do ii=1,nvars
                ! get variable name
                call handle_err(nf90_inquire_variable(ncid,ii,name,itype,n,dims,natts))
                ! get calendar info (check exact match):
                if(k==1 .and. name=='time') then
                   do jj=1,natts
                      call handle_err(nf90_inq_attname(ncid,ii,jj,attnam))
                      if(attnam=='calendar') then
                         call handle_err(nf90_get_att(ncid,ii,attnam,calendar))
                      endif
                   enddo
                endif
                ! use keywords if calendar is not defined yet
                if(calendar==' ') then
                   if(k==1.and.index(name,'time')/=0) then
                      do jj=1,natts
                         call handle_err(nf90_inq_attname(ncid,ii,jj,attnam))
                         if(attnam=='calendar') then
                            call handle_err(nf90_get_att(ncid,ii,attnam,calendar))
                         endif
                      enddo
                   endif
                endif
                do jj=1,nvar
                   if(name==avars(jj)) then
                      id(jj) = ii
                      vndims(jj) = n
                   endif
                enddo
             enddo
             call handle_err(nf90_close(ncid))
             if (minval(id)>0) goto 10
          endif
       enddo
    enddo
10  continue
    ! default calendar is the standard gregorian:
    if(calendar.eq.' ') then
       calendar = 'gregorian'
    endif
    ! find the number of 3-D & 4-D variables 
    nvarx1 = sum(vndims,mask=vndims==3)/3
    nvarxm = sum(vndims,mask=vndims==4)/4
    if(minval(vndims)<0) then
       print*, 'ERROR: at least one variable does not exist in files'
       stop
    endif
    if((nvarx1 + nvarxm)/=nvar) then
       print*, 'ERROR: at least one variable does not have 3 or 4 dimensions'
       stop
    endif
    allocate(varsx1(nvarx1),varsxm(nvarxm))
    ! put 1-level variables in varsx1 & multi-level in varsxm
    i1 = 0
    i2 = 0
    do i=1,nvar
       if(vndims(i)==3) then
          i1 = i1 + 1
          varsx1(i1) = avars(i)
       elseif(vndims(i)==4) then
          i2 = i2 + 1
          varsxm(i2) = avars(i)
       endif
    enddo
    call maketimes(time1,time2,n_ave,unit,calendar,time_out,nt,nyear,cont)
    allocate(n_ave_x1(nvarx1,nt,nyear))
    allocate(n_ave_xm(nvarxm,nt,nyear))
    ! for avoiding missing times when searching for times 
    ! (new array "tt" has no missing values):
    allocate(tt(nt*nyear))
    ntime = 0
    do it=1,nt*nyear
       if(time_out(it,1,1)/=all) then
          ntime = ntime + 1
          tt(ntime) = it
       endif
    enddo
    ! for testing:
    !    do i = 1,nyear
    !       do it = 1,nt
    !          call cal(calendar,time_out(it,i,1),hour1,day1,month1,year1)
    !          call cal(calendar,time_out(it,i,2),hour2,day2,month2,year2)
    !          write(*,'(i3,i3,i5,i3,i7,i3,i5,i3)')  &
    !               month1,day1,year1,hour1,month2,day2,year2,hour2
    !       enddo
    !    enddo
    ! open netcdf file:
    lev_flag=-1
    do ifile=1,nfile
       call handle_err(nf90_open(trim(files(ifile)),nf90_nowrite,ncid))
       id = -1
       ! get # dimensions & variables:
       call handle_err(nf90_inquire(ncid, ndims, nvars, natts, recdim))
       ! get dimension sizes:
       ! first check for exact match:
       mlon = -1; mlat = -1; mlev = -1; mt = -1
       do i=1,ndims
          call handle_err(nf90_inquire_dimension(ncid,i,name,n))
          if(name=='lon') then
             mlon = n
             dlonid = i
          elseif(name=='lat') then
             mlat = n
             dlatid = i
          elseif(name=='lev' .or. name=='levsoi') then
             mlev = n
             dlevid = i
          elseif(name=='time') then
             mt = n
             dtid = i
          endif
       enddo
       ! then check for key words within dimension name
       do i=1,ndims
          call handle_err(nf90_inquire_dimension(ncid,i,name,n))
          if(index(name,'lon')/=0 .and. mlon==-1) then
             mlon = n
             dlonid = i
          elseif(index(name,'lat')/=0 .and. mlat==-1) then
             mlat = n
             dlatid = i
          elseif(index(name,'lev')/=0 .and. mlev==-1) then
             mlev = n
             dlevid = i
          elseif(index(name,'time')/=0 .and. mt==-1) then
             mt = n
             dtid = i
          endif
       enddo
       if(mt==0) goto 1 ! if no time steps in file go to next file
       ! get variable info:
       ! first check for exact match:
       lonid = -1; latid = -1; levid = -1; tid = -1; tbid = -1
       do i=1,nvars
          call handle_err(nf90_inquire_variable(ncid,i,name,itype,n,dims,natts))
          !          print*, i,trim(name),itype,n,dims(1:n),natts,rcode
          if(name=='lon') then
             lonid = i
             namelon = name
             lontype = itype
          elseif(name=='lat') then
             latid = i
             namelat = name
             lattype = itype
          elseif(name=='lev' .or. name=='levsoi') then
             levid = i
             namelev = name
             levtype = itype
          elseif(name=='time') then
             tid = i
             namet = name
             ttype = itype
             tatts = natts
             ! for ncar model output (because time is at end of step):
          elseif(name=='time_bnds'.or.name=='time_bounds') then
             tbid = i
             ! for ncar land model time_bnds type not the same as time type !!!!
             tbtype = itype 
          endif
          do j=1,nvar
             if(name==avars(j)) then
                id(j) = i
                types(j) = itype
                atts(j) = natts
                vdims(1:n,j) =  dims(1:n)
                ! get add_offset, scale_factor and missing values
                do jj=1,natts
                   call handle_err(nf90_inq_attname(ncid,id(j),jj,attnam))
                   if(attnam=='add_offset') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'add_offset',add(j)))
                      elseif(attype==nf90_double) then
                         call handle_err(nf90_get_att(ncid,id(j),'add_offset',temp8))
                         add(j)=real(temp8)
                      else
                         print*, 'ERROR: add_offset is not real*4 or real*8'
                         stop
                      endif
                   endif
                   if(attnam=='scale_factor') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'scale_factor',scalef(j)))
                      elseif(attype==nf90_double) then
                         call handle_err(nf90_get_att(ncid,id(j),'scale_factor',temp8))
                         scalef(j)=real(temp8)
                      else
                         print*, 'ERROR: scale_factor is not real*4 or real*8'
                         stop
                      endif
                   endif
                   if(attnam=='missing_value') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype/=types(j)) then
                         print*, 'ERROR: missing_value not same type as variable'
                         stop
                      endif
                      if(attype==nf90_short) then
                         call handle_err(nf90_get_att(ncid,id(j),'missing_value',miss3(j)))
                      elseif(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'missing_value',miss5(j)))
                      else
                         print*, 'ERROR: variable type is not real*4 or integer*2'
                         stop
                      endif
                   endif
                enddo
             endif
          enddo
       enddo
       ! then check for key words within variable name:
       do i=1,nvars
          call handle_err(nf90_inquire_variable(ncid,i,name,itype,n,dims,natts))
          if(index(name,'lon')/=0 .and. lonid==-1) then
             lonid = i
             namelon = name
             lontype = itype
          elseif(index(name,'lat')/=0 .and. latid==-1) then
             latid = i
             namelat = name
             lattype = itype
          elseif(index(name,'lev')/=0 .and. levid==-1) then
             levid = i
             namelev = name
             levtype = itype
          elseif(index(name,'time')/=0 .and. tid==-1) then
             tid = i
             namet = name
             ttype = itype
             tatts = natts
          endif
       enddo
       ! check that dimensions are in the default order
       do j=1,nvar
          if(id(j)>0) then
             if(vdims(1,j)/=dlonid) then
                print*, 'ERROR: longitude not first dimension in variable'
                stop
             endif
             if(vdims(2,j)/=dlatid) then
                print*, 'ERROR: latitude not second dimension in variable'
                stop
             endif
             if(vdims(3,j)/=dlevid .and. vndims(j) == 4) then
                print*, 'ERROR: level not third dimension in 4-D variable'
                stop
             endif
             if(vdims(4,j)/=dtid .and. vndims(j) == 4) then
                print*, 'ERROR: time not fourth dimension in 4-D variable'
                stop
             endif
             if(vdims(3,j)/=dtid .and. vndims(j) == 3) then
                print*, 'ERROR: time not third dimension in 3-D variable'
                stop
             endif
          endif
       enddo
       ! find lon, lat and level range
       ! first: longitudes
       if(lontype==nf90_float) then
          allocate(lons5(mlon))
          call handle_err(nf90_get_var(ncid,lonid,lons5,(/1/),(/mlon/)))
       elseif(lontype==nf90_double) then
          allocate(lons6(mlon))
          allocate(lons5(mlon))
          call handle_err(nf90_get_var(ncid,lonid,lons6,(/1/),(/mlon/)))
          do ilon=1,mlon
             lons5(ilon) = real(lons6(ilon))
          enddo
          deallocate(lons6)
       else
          print*, 'ERROR: lon type not real*4 or real*8'
          stop
       endif
       ! if lons == missing then do all longitudes
       if(lon1==missing .or. lon2==missing) then
          dmin = 1; dmax = mlon
       else
          ! make lons between 0 & 360
          do ilon=1,mlon
             if(lons5(ilon)<0.) then
                lons5(ilon) = lons5(ilon) + 360.
             elseif(lons5(ilon)>=360.) then
                lons5(ilon) = lons5(ilon) - 360.
             endif
          enddo
          test42 = -huge(0.0)
          test41 = huge(0.0)
          dmax = 0; dmin = 0
          do ilon=1,mlon
             if(lons5(ilon)<=lon2) then
                if(lons5(ilon)>test42) then
                   test42 = lons5(ilon)
                   dmax = ilon
                endif
             endif
             if(lons5(ilon)>=lon1) then
                if(lons5(ilon)<test41) then
                   test41 = lons5(ilon)
                   dmin = ilon
                endif
             endif
          enddo
          if(dmin==0 .or. dmax==0) then
             print*, 'ERROR: longitude range outside dataset'
             stop
          endif
       endif
       if(dmax>=dmin) then
          nlon = dmax - dmin + 1
       else
          nlon = (mlon - dmin + 1) + dmax
       endif
       allocate(lons_test(nlon))
       test4 = 0.
       if(lons5(dmin)>=180.) then
          test4 = -360.
       endif
       if(dmax>=dmin) then
          do ilon=dmin,dmax
             if(lons5(ilon)<180.) then
                lons_test(ilon - dmin + 1) = lons5(ilon)
             else
                lons_test(ilon - dmin + 1) = lons5(ilon) + test4
             endif
          enddo
       else
          do ilon=dmin,mlon
             lons_test(ilon - dmin + 1) = lons5(ilon) - 360.
          enddo
          do ilon=1,dmax
             lons_test(mlon - dmin + 1 + ilon) = lons5(ilon)
          enddo
       endif
       if(ifile>1) then
          do ilon=1,nlon
             if(lons_test(ilon)/=lons(ilon)) then
                print*, 'ERROR: longitude grids not the same among files'
                stop
             endif
          enddo
       else
          allocate(lons(nlon))
          lons = lons_test
       endif
       lon_start=dmin
       lon_end=dmax
       deallocate(lons5)          
       deallocate(lons_test)
       ! second: latitudes
       if(lattype==nf90_float) then
          allocate(lats5(mlat))
          call handle_err(nf90_get_var(ncid,latid,lats5,(/1/),(/mlat/)))
       elseif(lattype==nf90_double) then
          allocate(lats6(mlat))
          allocate(lats5(mlat))
          call handle_err(nf90_get_var(ncid,latid,lats6,(/1/),(/mlat/)))
          do ilat=1,mlat
             lats5(ilat) = real(lats6(ilat))
          enddo
          deallocate(lats6)
       else
          print*, 'ERROR: lat type not real*4 or real*8'
          stop
       endif
       if(lat1==missing .or. lat2==missing) then
          dmin = 1; dmax = mlat
       else
          test42 = -huge(0.0)
          test41 = huge(0.0)
          dmax = 0; dmin = 0
          do ilat=1,mlat
             if(lats5(ilat)<=lat2) then
                if(lats5(ilat)>test42) then
                   test42 = lats5(ilat)
                   dmax = ilat
                endif
             endif
             if(lats5(ilat)>=lat1) then
                if(lats5(ilat)<test41) then
                   test41 = lats5(ilat)
                   dmin = ilat
                endif
             endif
          enddo
          if(dmin==0 .or. dmax==0) then
             print*, 'ERROR: latitude range outside dataset'
             stop
          endif
       endif
       if(dmax<dmin) then
          jj=dmin
          dmin=dmax
          dmax=jj
       endif
       nlat = dmax - dmin + 1
       allocate(lats_test(nlat))
       lats_test = lats5(dmin:dmax)
       if(ifile>1) then
          do ilat=1,nlat
             if(lats_test(ilat)/=lats(ilat)) then
                print*, 'ERROR: latitude grids not the same among files'
                stop
             endif
          enddo
       else
          allocate(lats(nlat))
          lats = lats_test
       endif
       lat_start=dmin
       lat_end=dmax
       deallocate(lats5)          
       deallocate(lats_test)
       ! third: vertical levels
       if(mlev/=-1) then
          if(levtype==nf90_float) then
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs5,(/1/),(/mlev/)))
          elseif(levtype==nf90_double) then
             allocate(levs6(mlev))
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs6,(/1/),(/mlev/)))
             do ilev=1,mlev
                levs5(ilev) = real(levs6(ilev))
             enddo
             deallocate(levs6)
          elseif(levtype==nf90_int) then
             allocate(levs4(mlev))
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs4,(/1/),(/mlev/)))
             do ilev=1,mlev
                levs5(ilev) = real(levs4(ilev))
             enddo
             deallocate(levs4)
          else
             print*, 'ERROR: lev type not integer, real*4 or real*8'
             stop
          endif
          if(lev1==missing .or. lev2==missing) then
             dmin = 1; dmax = mlev
          else
             test42 = -huge(0.0)
             test41 = huge(0.0)
             dmax = 0; dmin = 0
             do ilev=1,mlev
                if(levs5(ilev)<=lev2) then
                   if(levs5(ilev)>test42) then
                      test42 = levs5(ilev)
                      dmax = ilev
                   endif
                endif
                if(levs5(ilev)>=lev1) then
                   if(levs5(ilev)<test41) then
                      test41 = levs5(ilev)
                      dmin = ilev
                   endif
                endif
             enddo
             if(dmin==0 .or. dmax==0) then
                print*, 'ERROR: vertical level range outside dataset'
                stop
             endif
          endif
          if(dmax<dmin) then
             jj=dmin
             dmin=dmax
             dmax=jj
          endif
          nlev = dmax - dmin + 1
          allocate(levs_test(nlev))
          levs_test = levs5(dmin:dmax)
          if(lev_flag>0) then
             do ilev=1,nlev
                if(levs_test(ilev)/=levs(ilev)) then
                   print*, 'ERROR: vertical levels not the same among files'
                   stop
                endif
             enddo
          else
             allocate(levs(nlev))
             levs = levs_test
             lev_flag=1
          endif
          lev_start=dmin
          lev_end=dmax
          deallocate(levs5)          
          deallocate(levs_test)
       endif
       ! get time units:
       time_units=' '
       do j=1,tatts
          call handle_err(nf90_inq_attname(ncid,tid,j,attnam))
          if(attnam=='units') then
             call handle_err(nf90_get_att(ncid,tid,attnam,time_units))
          endif
       enddo
       if(time_units==' ') then
          print*, 'ERROR: no time units in file'
          stop
       endif
       ! get time units
       ! are time units in hours or days?:
       if(index(time_units,'hours')/=0) then
          hour_day=1
       elseif(index(time_units,'days')/=0) then
          hour_day=24
       else
          print*, 'ERROR: time units not in hours or days'
          stop
       endif
       call get_time_units(time_units,year_cal,month_cal,day_cal,hour_cal)
       ! find time (in hours) between reference date in file with the 
       ! default reference date in this program:
       ! year_cal > year_ref  ==> positive time difference
       ! year_cal < year_ref  ==> negative time difference
       call get_time_diff(calendar,year_cal,month_cal,day_cal,hour_cal,hour_diff)
       ! fourth: get times
       if(ttype==nf90_int) then
          allocate(times4(mt))
          call handle_err(nf90_get_var(ncid,tid,times4,(/1/),(/mt/)))
          ! if times are at the end of time_bnds,
          ! then make times at the beginning of time_bnds instead:
          if(tbid/=-1) then
             allocate(tb4(2,mt))
             start(1:2)=(/1,1/)
             ncount(1:2)=(/2,mt/)
             if(tbtype==nf90_int) then
                call handle_err(nf90_get_var(ncid,tbid,tb4,start,ncount))
             elseif(tbtype==nf90_float) then
                allocate(tb5(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb5,start,ncount))
                tb4=int(tb5)
                deallocate(tb5)
             elseif(tbtype==nf90_double) then
                allocate(tb6(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb6,start,ncount))
                tb4=int(tb6)
                deallocate(tb6)
             endif
             if(tb4(2,1)==times4(1)) then
                do it=1,mt
                   times4(it)=tb4(1,it)
                enddo
             endif
             deallocate(tb4)
          endif
       elseif(ttype==nf90_double) then
          allocate(times6(mt))
          allocate(times4(mt))
          call handle_err(nf90_get_var(ncid,tid,times6,(/1/),(/mt/)))
          ! if times are at the end of time_bnds,
          ! then make times at the beginning of time_bnds instead:
          if(tbid/=-1) then
             allocate(tb6(2,mt))
             start(1:2)=(/1,1/)
             ncount(1:2)=(/2,mt/)
             if(tbtype==nf90_double) then
                call handle_err(nf90_get_var(ncid,tbid,tb6,start,ncount))
             elseif(tbtype==nf90_float) then
                allocate(tb5(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb5,start,ncount))
                tb6=dble(tb5)
                deallocate(tb5)
             elseif(tbtype==nf90_int) then
                allocate(tb4(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb4,start,ncount))
                tb6=dble(tb4)
                deallocate(tb4)
             endif
             if(tb6(2,1)==times6(1)) then
                do it=1,mt
                   times6(it)=tb6(1,it)
                enddo
             endif
             deallocate(tb6)
          endif
          ! special stuff for fractional days (convert everything to hours): 
          if(hour_day==1) then
             do it=1,mt
                times4(it) = nint(times6(it))
             enddo
          else
             do it=1,mt
                if(times6(it)/=dnint(times6(it))) then
                   do jt=1,mt
                      times4(jt) = nint(24*times6(jt))
                   enddo
                   hour_day=1
                   exit
                else
                   times4(it) = nint(times6(it))
                endif
             enddo
          endif
          deallocate(times6)
       elseif(ttype==nf90_float) then
          allocate(times5(mt))
          allocate(times4(mt))
          call handle_err(nf90_get_var(ncid,tid,times5,(/1/),(/mt/)))
          ! if times are at the end of time_bnds,
          ! then make times at the beginning of time_bnds instead:
          if(tbid/=-1) then
             allocate(tb5(2,mt))
             start(1:2)=(/1,1/)
             ncount(1:2)=(/2,mt/)
             if(tbtype==nf90_float) then
                call handle_err(nf90_get_var(ncid,tbid,tb5,start,ncount))
             elseif(tbtype==nf90_double) then
                allocate(tb6(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb6,start,ncount))
                tb5=real(tb6)
                deallocate(tb6)
             elseif(tbtype==nf90_int) then
                allocate(tb4(2,mt))
                call handle_err(nf90_get_var(ncid,tbid,tb4,start,ncount))
                tb5=real(tb4)
                deallocate(tb4)
             endif
             if(tb5(2,1)==times5(1)) then
                do it=1,mt
                   times5(it)=tb5(1,it)
                enddo
             endif
             deallocate(tb5)
          endif
          ! special stuff for fractional days (convert everything to hours): 
          if(hour_day==1) then
             do it=1,mt
                times4(it) = nint(times5(it))
             enddo
          else
             do it=1,mt
                if(times5(it)/=anint(times5(it))) then
                   do jt=1,mt
                      times4(jt) = nint(24*times5(jt)) 
                   enddo
                   hour_day=1
                   exit
                else
                   times4(it) = nint(times5(it)) 
                endif
             enddo
          endif
          deallocate(times5)
       else
          print*, 'ERROR: time type not integer, real*4 or real*8'
          stop
       endif
       ! scale if units are days & then shift to reference time in program (year_ref):
       times4 = times4*hour_day + hour_diff
       ! for testing:
       !       print *, hour_diff
       !       print *, year_cal,month_cal,day_cal,hour_cal
       !       do it=1,mt
       !          call cal(calendar,times4(it),hour1,day1,month1,year1)
       !          write(*,'(i3,i3,i5,i3)')  &
       !               month1,day1,year1,hour1
       !       enddo
       ! get output data arrays ready:
       if( (.not.allocated(x1)) ) then
          if(nvarx1>0) then
             if(ave_lon) then
                allocate(x1(1,nlat,nvarx1,nt,nyear))
             else
                allocate(x1(nlon,nlat,nvarx1,nt,nyear))
             endif
             x1 = 0.
          endif
          n_ave_x1 = 0
       endif
       if( (nlev>0) .and. (.not.allocated(xm))) then 
          if(nvarxm>0) then
             if(ave_lon) then
                allocate(xm(1,nlat,nlev,nvarxm,nt,nyear))
             else
                allocate(xm(nlon,nlat,nlev,nvarxm,nt,nyear))
             endif
             xm = 0.
          endif
          n_ave_xm = 0
       endif
       ! map times in file to correct times specified in input 
       ! (but first allocate data read arrays):
       if(nlev>0) then
          allocate(data3(nlon,nlat,nlev,1))
          allocate(data5(nlon,nlat,nlev,1))
       else
          allocate(data3(nlon,nlat,1,1))
          allocate(data5(nlon,nlat,1,1))
       endif
       do it=1,mt
          tmap=-1
          nm0=-1
          n0 = 1; n1 = ntime+1; nm = (n1 + n0)/2
          do i=1,ntime
             if(times4(it)<time_out(tt(nm),1,1)) then
                n1 = nm
                nm = (n1 + n0)/2
             elseif(times4(it)>time_out(tt(nm),1,2)) then
                n0 = nm
                nm = (n1 + n0)/2
             else
                tmap = tt(nm)
                ! for testing
                !                print*, tt(nm),time_out(tt(nm),1,1),times4(it), &
                !                     time_out(tt(nm),1,2),i,ifile
                exit
             endif
             if(nm==nm0) exit
             nm0=nm
          enddo
          ! get field variables (FINALLY!!!)
          if(tmap>0) then
             do j=1,nvar
                if(id(j)>0) then  ! check if variable is in this file
                   ivar=-1
                   do i=1,nvarx1  ! find # dimensions and order to put in output
                      if(avars(j)==varsx1(i)) then
                         ivar=i
                         dim=3
                         exit
                      endif
                   enddo
                   do i=1,nvarxm
                      if(avars(j)==varsxm(i)) then
                         ivar=i
                         dim=4
                         exit
                      endif
                   enddo
                   if(dim==3) then
                      if(types(j)==nf90_float) then
                         start(1:dim)=(/lon_start,lat_start,it/)
                         ncount(1:dim)=(/nlon,nlat,1/)
                         if(lon_end>=lon_start) then
                            call handle_err(nf90_get_var(ncid,id(j),data5,start,ncount))
                         else
                            start(1:dim)=(/lon_start,lat_start,it/)
                            ncount(1:dim)=(/mlon-lon_start+1,nlat,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data5(1:mlon-lon_start+1,:,:,:), &
                                 start,ncount))
                            start(1:dim)=(/1,lat_start,it/)
                            ncount(1:dim)=(/lon_end,nlat,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data5(nlon-lon_end+1:nlon,:,:,:), &
                                 start,ncount))
                         endif
                         if(lats(1)>lats(nlat)) then
                            data5=data5(:,nlat:1:-1,:,:)
                         endif
                         if(ave_lon) then
                            data5(1,:,:,:) = sum(data5,1,data5/=miss5(j))/ &
                                 count(data5/=miss5(j),1)
                         endif
                      elseif(types(j)==nf90_short) then
                         start(1:dim)=(/lon_start,lat_start,it/)
                         ncount(1:dim)=(/nlon,nlat,1/)
                         if(lon_end>=lon_start) then
                            call handle_err(nf90_get_var(ncid,id(j),data3,start,ncount))
                         else
                            start(1:dim)=(/lon_start,lat_start,it/)
                            ncount(1:dim)=(/mlon-lon_start+1,nlat,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data3(1:mlon-lon_start+1,:,:,:), &
                                 start,ncount))
                            start(1:dim)=(/1,lat_start,it/)
                            ncount(1:dim)=(/lon_end,nlat,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data3(nlon-lon_end+1:nlon,:,:,:), &
                                 start,ncount))
                         endif
                         if(lats(1)>lats(nlat)) then
                            data3=data3(:,nlat:1:-1,:,:)
                         endif
                         miss5(j)=scalef(j)*miss3(j) + add(j)
                         data5 = scalef(j)*data3 + add(j)
                         if(ave_lon) then
                            data5(1,:,:,:) = sum(data5,1,data5/=miss5(j))/ &
                                 count(data5/=miss5(j),1)
                         endif
                      else
                         print*, 'ERROR: variable not of integer*2 or real*4 type'
                         stop
                      endif
                      if(ave_lon) then
                         x1(1,:,ivar,tmap,1) = merge(x1(1,:,ivar,tmap,1) + &
                              merge(data5(1,:,1,1),missing,data5(1,:,1,1)/=miss5(j)) &
                              ,missing,x1(1,:,ivar,tmap,1)/=missing)
                      else
                         x1(:,:,ivar,tmap,1) = merge(x1(:,:,ivar,tmap,1) + &
                              merge(data5(:,:,1,1),missing,data5(:,:,1,1)/=miss5(j)) &
                              ,missing,x1(:,:,ivar,tmap,1)/=missing)
                      endif
                      n_ave_x1(ivar,tmap,1) = n_ave_x1(ivar,tmap,1) + 1
                   else ! 4 dimensional variables:
                      if(types(j)==nf90_float) then
                         start(1:dim)=(/lon_start,lat_start,lev_start,it/)
                         ncount(1:dim)=(/nlon,nlat,nlev,1/)
                         if(lon_end>=lon_start) then
                            call handle_err(nf90_get_var(ncid,id(j),data5,start,ncount))
                         else
                            start(1:dim)=(/lon_start,lat_start,lev_start,it/)
                            ncount(1:dim)=(/mlon-lon_start+1,nlat,nlev,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data5(1:mlon-lon_start+1,:,:,:), &
                                 start,ncount))
                            start(1:dim)=(/1,lat_start,lev_start,it/)
                            ncount(1:dim)=(/lon_end,nlat,nlev,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data5(nlon-lon_end+1:nlon,:,:,:),start,ncount))
                         endif
                         if(lats(1)>lats(nlat)) then
                            data5=data5(:,nlat:1:-1,:,:)
                         endif
                         if(levs(1)>levs(nlev)) then
                            data5=data5(:,:,nlev:1:-1,:)
                         endif
                         if(ave_lon) then
                            data5(1,:,:,:) = sum(data5,1,data5/=miss5(j))/ &
                                 count(data5/=miss5(j),1)
                         endif
                      elseif(types(j)==nf90_short) then
                         start(1:dim)=(/lon_start,lat_start,lev_start,it/)
                         ncount(1:dim)=(/nlon,nlat,nlev,1/)
                         if(lon_end>=lon_start) then
                            call handle_err(nf90_get_var(ncid,id(j),data3,start,ncount))
                         else
                            start(1:dim)=(/lon_start,lat_start,lev_start,it/)
                            ncount(1:dim)=(/mlon-lon_start+1,nlat,nlev,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data3(1:mlon-lon_start+1,:,:,:), &
                                 start,ncount))
                            start(1:dim)=(/1,lat_start,lev_start,it/)
                            ncount(1:dim)=(/lon_end,nlat,nlev,1/)
                            call handle_err(nf90_get_var(ncid,id(j),data3(nlon-lon_end+1:nlon,:,:,:), &
                                 start,ncount))
                         endif
                         if(lats(1)>lats(nlat)) then
                            data3=data3(:,nlat:1:-1,:,:)
                         endif
                         if(levs(1)>levs(nlev)) then
                            data3=data3(:,:,nlev:1:-1,:)
                         endif
                         miss5(j)=scalef(j)*miss3(j) + add(j)
                         data5 = scalef(j)*data3 + add(j)
                         if(ave_lon) then
                            data5(1,:,:,:) = sum(data5,1,data5/=miss5(j))/ &
                                 count(data5/=miss5(j),1)
                         endif
                      else
                         print*, 'ERROR: variable not of integer*2 or real*4 type'
                         stop
                      endif
                      if(ave_lon) then
                         xm(1,:,:,ivar,tmap,1) = merge(xm(1,:,:,ivar,tmap,1) + &
                              merge(data5(1,:,:,1),missing,data5(1,:,:,1)/=miss5(j)) &
                              ,missing,xm(1,:,:,ivar,tmap,1)/=missing)
                      else
                         xm(:,:,:,ivar,tmap,1) = merge(xm(:,:,:,ivar,tmap,1) + &
                              merge(data5(:,:,:,1),missing,data5(:,:,:,1)/=miss5(j)) &
                              ,missing,xm(:,:,:,ivar,tmap,1)/=missing)
                      endif
                      n_ave_xm(ivar,tmap,1) = n_ave_xm(ivar,tmap,1) + 1
                   endif ! end getting one variable at one file's time step
                endif
             enddo
          endif
       enddo
       
       deallocate(data3)
       deallocate(data5)
       deallocate(times4)
       
1      continue
       call handle_err(nf90_close(ncid))
              
    enddo
    
    ! divide by number of times averaged together (missing if no times):
    if(nvarx1>0) then
       do it=1,nt*nyear
          do j=1,nvarx1
             if(n_ave_x1(j,it,1)==0) then
                x1(:,:,j,it,1) = missing
             else
                if(ave_lon) then
                   x1(1,:,j,it,1) = merge(x1(1,:,j,it,1)/n_ave_x1(j,it,1),missing,x1(1,:,j,it,1)/=missing)
                else
                   x1(:,:,j,it,1) = merge(x1(:,:,j,it,1)/n_ave_x1(j,it,1),missing,x1(:,:,j,it,1)/=missing)
                endif
             endif
          enddo
       enddo
    endif
    if(nvarxm>0) then
       do it=1,nt*nyear
          do j=1,nvarxm
             if(n_ave_xm(j,it,1)==0) then
                xm(:,:,:,j,it,1) = missing
             else
                if(ave_lon) then
                   xm(1,:,:,j,it,1) = merge(xm(1,:,:,j,it,1)/n_ave_xm(j,it,1),missing,xm(1,:,:,j,it,1)/=missing)
                else
                   xm(:,:,:,j,it,1) = merge(xm(:,:,:,j,it,1)/n_ave_xm(j,it,1),missing,xm(:,:,:,j,it,1)/=missing)
                endif
             endif
          enddo
       enddo
    endif
    
    if(lats(1)>lats(nlat)) then
       lats=lats(nlat:1:-1)
    endif
    if(nvarxm>0) then
       if(levs(1)>levs(nlev)) then
          levs=levs(nlev:1:-1)
       endif
    endif
    
    deallocate(avars)
    deallocate(id)
    deallocate(types)
    deallocate(atts)
    deallocate(scalef)
    deallocate(add)
    deallocate(miss3)
    deallocate(miss5)
    deallocate(files)
    deallocate(vndims)
    deallocate(varsx1,varsxm)
    deallocate(time_out)
    deallocate(n_ave_x1)
    deallocate(n_ave_xm)
    deallocate(tt)
    
  end subroutine read_nc
  
  subroutine maketimes(time1,time2,n_ave,unit2,calendar,time_out,nt,nyear,cont)
    
    use constants
    
    integer :: time1(4),time2(4)
    integer :: n_ave
    character (len=8) :: unit2,unit
    character (len=100) :: calendar
    integer, allocatable :: time_out(:,:,:),time_out2(:,:,:)  
    logical :: cont
    integer :: nt,nyear
    integer :: nt2,nyear2
    integer :: m0,m1,d0,d1,h0,h1
    integer :: months(12),ileap
    ! hours of initial year (measured from start of year_ref)    
    integer :: hours0
    integer :: step0
    integer :: hour_start(0:1),hour_end(0:1),hour_incr(0:1),hour_acc
    integer :: rem,shift
    integer :: beg_year,end_year,hours_per_year
    integer :: ical
    integer :: imonth,iyear,it,i
    
    beg_year = time1(4)
    end_year = time2(4)
    nyear = end_year - beg_year + 1
    
    do i=1,100
       if(calendar(i:i)=='\0') calendar(i:i) = ' '
    enddo
    if ((calendar.eq.'gregorian').or.(calendar.eq.'standard')  &
         .or.(calendar.eq.'proleptic_gregorian')) then
       months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       ical = 1
       hours_per_year = 24*365
    elseif ((calendar.eq.'noleap').or.(calendar.eq.'365_day'))  then
       months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       ical = 2
       hours_per_year = 24*365
    elseif ((calendar.eq.'360_day'))  then
       months = (/30,30,30,30,30,30,30,30,30,30,30,30/)
       ical = 3
       hours_per_year = 24*360
    else
       print *, 'ERROR:',trim(calendar),' is not on calendar list'
       stop
    endif
    
    ! translate the "all" specifications to actual ranges:
    m0 = time1(3)
    m1 = time2(3)
    if(m0==all) m0 = 1
    if(m1==all) m1 = 12
    d0 = time1(2)
    d1 = time2(2)
    if(d0==all) d0 = 1
    if(d1==all) d1 = months(m1)
    h0 = time1(1)
    h1 = time2(1)
    if(h0==all) h0 = 0
    if(h1==all) h1 = 23
    if( h0<0 .or. h0>23 .or. h1<0 .or. h1>23) then
       print*, 'ERROR: input hours beyond valid range (valid range = 0 to 23)'
       stop
    endif
    if( m0<1 .or. m0>12 .or. m1<1 .or. m1>12) then
       print*, 'ERROR: input months beyond valid range (valid range = 1 to 12)'
       stop
    endif
    if( d0<1 .or. d0>months(m0) .or. d1<1 .or. d1>months(m1)) then
       print*, 'ERROR: input days beyond valid range for calendar = ', &
            trim(calendar)
       stop
    endif
    ! get hours from start of year_ref to start of beg_year
    hours0 = 0
    if(year_ref>beg_year) then
       print*, 'ERROR: change year_ref so that it is less than starting year'
       stop
    endif
    do iyear=year_ref,beg_year - 1
       if(ical==1) then
          if(leap(iyear)==0) then
             hours0 = hours0 + 24*365
          else
             hours0 = hours0 + 24*366
          endif
       elseif(ical==2) then
          hours0 = hours0 + 24*365
       else
          hours0 = hours0 + 24*360
       endif
    enddo
    ! get time unit & initialize some things related to this
    unit = adjustl(unit2)
    select case (unit(1:1))
    case("h")
       step0 = 1*n_ave
    case("d")
       step0 = 24*n_ave
    case("m")
       step0 = 0
       d0 = 1
       d1 = months(m1)
       h0 = 0
       h1 = 23
    case("y")
       step0 = 0
       d0 = 1
       d1 = months(m1)
       h0 = 0
       h1 = 23
    case default
       print*, 'ERROR: do not recognize time unit'
       stop
    end select
    ! find start and end hours (time1(1:3) & time2(1:3)) within standard year
    hour_incr = hours_per_year
    hour_start = 24*(d0 - 1) + h0
    do imonth=1,m0-1
       hour_start = hour_start + 24*months(imonth)
    enddo
    hour_end = 24*(d1 - 1) + h1
    do imonth=1,m1-1
       hour_end = hour_end + 24*months(imonth)
    enddo
    if(hour_end(0)<hour_start(0)) then
       hour_end = hour_end + hours_per_year
       m1 = m1 + 12
       nyear = nyear - 1
       end_year = end_year - 1
       if(nyear<1) then
          print*, 'ERROR: end time is before start time'
          stop
       endif
    endif
    ! adjust hours for leap years:
    if(ical==1) then
       if( (time1(3)==all.or.time2(3)==all) .and. &
            (time1(2)==all.or.time2(2)==all) .and. &
            (time1(1)==all.or.time2(1)==all) .and. &
            (step0==24.or.step0==12.or.step0==8.or.step0==6.or. &
            step0== 4.or.step0== 3.or.step0==2.or.step0==1) ) then
          hour_end(1) =  hour_end(0) + 24
          hour_incr(1) = hour_incr(0) + 24
       else
          if( m0>2 ) then
             hour_start(1) = hour_start(0) + 24
             hour_end(1) = hour_end(0) + 24
          endif
          hour_incr(1) = hour_incr(0) + 24
       endif
    endif
    ! do month & year time steps 
    ! (right now time step = 1 year or 1 month, not multiple years or months) 
    if(step0==0) then
       if(unit(1:1)=='m') then
          nt = m1 - m0 + 1
          allocate(time_out(nt,nyear,2))
          hour_acc = 0
          do iyear=beg_year,end_year
             if(ical==1 .and. m0==2 .and. leap(iyear)==1) then
                ileap = 24
             else
                ileap = 0
             endif
             time_out(1,iyear - beg_year + 1,1) = hours0 + hour_acc +  & 
                  hour_start(leap(iyear))
             time_out(1,iyear - beg_year + 1,2) =  &
                  time_out(1,iyear - beg_year + 1,1) + 24*months(m0) + ileap - 1
             do imonth=m0+1,m1
                if(ical==1 .and. mod(imonth,12)==2 .and. & 
                     leap(iyear + imonth/12)==1) then
                   ileap = 24
                else
                   ileap = 0
                endif
                time_out(imonth - m0 + 1,iyear - beg_year + 1,1) = &
                     time_out(imonth - m0,iyear - beg_year + 1,2) + 1 
                time_out(imonth - m0 + 1,iyear - beg_year + 1,2) = &
                     time_out(imonth - m0 + 1,iyear - beg_year + 1,1) + & 
                     24*months(modulo(imonth-1,12)+1) + ileap - 1
             enddo
             hour_acc =  hour_acc + hour_incr(leap(iyear))
          enddo
       else
          nt = 1
          allocate(time_out(nt,nyear,2))
          hour_acc = 0
          do iyear=beg_year,end_year
             time_out(1,iyear - beg_year + 1,1) = hours0 + hour_acc +  & 
                  hour_start(leap(iyear))
             if(m0<=2) then
                time_out(1,iyear - beg_year + 1,2) =  &
                     time_out(1,iyear - beg_year + 1,1) + &
                     hour_incr(leap(iyear)) - 1
             else
                time_out(1,iyear - beg_year + 1,2) =  &
                     time_out(1,iyear - beg_year + 1,1) + &
                     hour_incr(leap(iyear + 1)) - 1
             endif
             hour_acc =  hour_acc + hour_incr(leap(iyear))
          enddo
       endif
    else
       ! do day & hour time steps
       ! find # time steps
       nt = (hour_end(1) - hour_start(1) + 1)/step0
       rem = hour_end(1) - hour_start(1) + 1 - nt*step0
       shift = ((rem/24)/2)*24
       hour_start = hour_start + shift
       hour_end = hour_end + shift
       
       allocate(time_out(nt,nyear,2))
       if( (time1(3)==all.or.time2(3)==all) .and. &
            (time1(2)==all.or.time2(2)==all) .and. &
            (time1(1)==all.or.time2(1)==all) .and. &
            (step0==24.or.step0==12.or.step0==8.or.step0==6.or. &
            step0== 4.or.step0== 3.or.step0==2.or.step0==1) ) then
          hour_acc=0
          do iyear=beg_year,end_year
             time_out(1,iyear - beg_year + 1,1) = hours0 + hour_acc +  & 
                  hour_start(leap(iyear))
             time_out(1,iyear - beg_year + 1,2) =  &
                  time_out(1,iyear - beg_year + 1,1) + step0 - 1
             do it=2,nt
                time_out(it,iyear - beg_year + 1,1) =  &
                     time_out(it - 1,iyear - beg_year + 1,1) + step0 
                time_out(it,iyear - beg_year + 1,2) =  &
                     time_out(it,iyear - beg_year + 1,1) + step0 - 1
             enddo
             ! extra days in non-leap years get times that do not exist
             if(ical==1 .and. leap(iyear)==0) then
                shift = 24/step0 - 1
                time_out(nt - shift:nt,iyear - beg_year + 1,1) = all
                time_out(nt - shift:nt,iyear - beg_year + 1,2) = all
             endif
             hour_acc =  hour_acc + hour_incr(leap(iyear))
          enddo
       else
          hour_acc = 0
          do iyear=beg_year,end_year
             time_out(1,iyear - beg_year + 1,1) = hours0 + hour_acc +  & 
                  hour_start(leap(iyear))
             time_out(1,iyear - beg_year + 1,2) =  &
                  time_out(1,iyear - beg_year + 1,1) + step0 - 1
             do it=2,nt
                time_out(it,iyear - beg_year + 1,1) =  &
                     time_out(it - 1,iyear - beg_year + 1,1) + step0 
                time_out(it,iyear - beg_year + 1,2) =  &
                     time_out(it,iyear - beg_year + 1,1) + step0 - 1
             enddo
             hour_acc =  hour_acc + hour_incr(leap(iyear))
          enddo
       endif
       
    endif

    ! for all times in nt & 1 in nyear
    ! also for continuous data with no missing data on day 366 of non leap years
    if(cont) then
       nt2 = count(time_out(:,:,1)/=all)
       allocate(time_out2(nt2,1,2))
       i = 0
       do iyear=1,nyear
          do it=1,nt
             if(time_out(it,iyear,1)/=all) then
                i = i + 1
                time_out2(i,1,:) = time_out(it,iyear,:)
             endif
          enddo
       enddo
       deallocate(time_out)
       allocate(time_out(nt2,1,2))
       time_out = time_out2
       deallocate(time_out2)
       nt = nt2
       nyear = 1
    endif

  end subroutine maketimes
  
  subroutine get_time_units(time_units,year_cal,month_cal,day_cal,hour_cal)
    
    character (len=100) :: time_units
    character (len=30) :: time_format
    character (len=1) :: junk
    integer year_cal,month_cal,day_cal,hour_cal
    integer :: l,ll,lll,llll
    
    ! remove "days since" or "hours since"
    l=index(time_units,'since')+5
    time_units(1:l)=' '
    time_units=adjustl(time_units)
    ! get location of dashes, etc. so that you know the formatting of the date:
    l=index(time_units,'-')
    ll=index(time_units(1:12),'-',.true.)
    lll=index(time_units(1:12),' ')
    llll=index(time_units(1:12),'\0')
    if(lll/=0 .and. llll/=0) then ! for stupid people who terminate netcdf strings with "null"
       lll = min(lll,llll)
    endif
    if(lll==0) then
       lll=llll ! for stupid people who terminate netcdf strings with "null"
    endif
    write(time_format,'(a,i1,a,i1,a,i1,a)')  &
         '(i',l-1,',a1,i',ll-l-1,',a1,i',lll-ll-1,')'
    read(time_units,time_format) year_cal,junk,month_cal,junk,day_cal
    ! remove date:
    time_units(1:lll)=' '
    time_units=adjustl(time_units)
    ! get hours:
    l=index(time_units,':')
    if(l/=0) then
       write(time_format,'(a,i1,a)')  &
            '(i',l-1,')'
       read(time_units,time_format) hour_cal
    else
       hour_cal=0
    endif
    
    return
  end subroutine get_time_units
  
  subroutine get_time_diff(calendar,year_cal,month_cal,day_cal,hour_cal,hour_diff)
    
    use constants
    
    character (len=100) :: calendar
    integer :: year_cal,month_cal,day_cal,hour_cal,hour_diff
    integer :: nday,months(12)
    integer :: ical
    integer :: imonth,iyear
    
    if ((calendar.eq.'gregorian').or.(calendar.eq.'standard')  &
         .or.(calendar.eq.'proleptic_gregorian')) then
       months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       ical = 1
    elseif ((calendar.eq.'noleap').or.(calendar.eq.'365_day'))  then
       months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       ical = 2
    else
       months = (/30,30,30,30,30,30,30,30,30,30,30,30/)
       ical = 3
    endif
    
    hour_diff=0
    
    if(year_cal < year_ref) then
       
       do iyear=year_cal+1,year_ref - 1
          if(ical==1) then
             if(leap(iyear)==0) then
                hour_diff = hour_diff - 24*365
             else
                hour_diff = hour_diff - 24*366
             endif
          elseif(ical==2) then
             hour_diff = hour_diff - 24*365
          else
             hour_diff = hour_diff - 24*360
          endif
       enddo
       iyear=year_cal
       do imonth=month_cal+1,12
          if(ical==1) then
             if(leap(iyear)==1 .and. imonth==2) then
                hour_diff = hour_diff - 24*months(imonth) - 24
             else
                hour_diff = hour_diff - 24*months(imonth)
             endif
          else
             hour_diff = hour_diff - 24*months(imonth)
          endif
       enddo
       imonth=month_cal
       if(ical==1) then
          if(leap(iyear)==1 .and. imonth==2) then
             nday = months(imonth) + 1
          else
             nday = months(imonth)
          endif
       else
          nday = months(imonth)
       endif
       
       hour_diff = hour_diff - 24*(nday - day_cal) - (24-hour_cal)
       
    else
       
       do iyear=year_ref,year_cal - 1
          if(ical==1) then
             if(leap(iyear)==0) then
                hour_diff = hour_diff + 24*365
             else
                hour_diff = hour_diff + 24*366
             endif
          elseif(ical==2) then
             hour_diff = hour_diff + 24*365
          else
             hour_diff = hour_diff + 24*360
          endif
       enddo
       iyear=year_cal
       do imonth=1,month_cal-1
          if(ical==1) then
             if(leap(iyear)==1 .and. imonth==2) then
                hour_diff = hour_diff + 24*months(imonth) + 24
             else
                hour_diff = hour_diff + 24*months(imonth)
             endif
          else
             hour_diff = hour_diff + 24*months(imonth)
          endif
       enddo
       hour_diff = hour_diff + 24*(day_cal-1) + hour_cal
       
    endif

    ! for stupid people who use gregorian calendar with a reference date 
    ! before or during 1582
    ! WHAT ARE YOU THINKING!!!!
    ! for those who don't know gregorian means gregorian after Oct 15 1582 BUT  
    ! Julian before with the ADDED catch that 10 days were SKIPPED with the 
    ! implementation of the calendar. WHAT ARE YOU THINKING!!: holding on to 
    ! some DUMB historical accident!!!!!!!!! (Are you listening CDC!!)
    ! PLEASE use the proleptic_gregorian instead or a later reference date !!
    if ((calendar.eq.'gregorian').or.(calendar.eq.'standard') ) then
       if(year_cal<=1582) then
          hour_diff = hour_diff - 48 ! two days of difference after julian
                                     ! added 12 extra leap days between 
                                     ! 1 and 1582 but then
                                     ! 10 days were skipped in October 1582
       endif
    endif
    
  end subroutine get_time_diff
  
  function leap(iyear)
    
    integer :: iyear,leap
    
    if((mod(iyear,4)/=0)) then
       leap = 0
    elseif((mod(iyear,400)==0)) then
       leap = 1
    elseif((mod(iyear,100)==0)) then
       leap = 0
    else
       leap = 1
    endif
    
    return
  end function leap
  
  subroutine read_no_time(lon11,lon22,ave_lon, lat1,lat2, lev1,lev2)
    
    use constants

    character (len=lvars) :: vars1
    real, intent(in) :: lon11, lon22
    real, intent(in) :: lat1, lat2
    real, intent(in) :: lev1, lev2
    logical, intent(in) :: ave_lon
    real :: lon1, lon2
    ! put files & variables in these:
    integer :: nfile,nvar
    character (len=200), allocatable :: files(:)
    character (len=30), allocatable :: avars(:), varsx1(:), varsxm(:)
    !! for netcdf:
    integer, allocatable :: id(:),types(:),atts(:),vndims(:)
    real, allocatable :: scalef(:),add(:)
    integer :: ncid,ndims,nvars,natts,recdim,n,itype
    integer :: attype,attlen
    integer*2, allocatable :: miss3(:)
    real*4, allocatable :: miss5(:)
    integer :: dims(32) ! 32 hard wired into code
    integer, allocatable :: vdims(:,:)
    character (len=30) :: name,namelon,namelat,namelev,attnam
    integer :: start(3),ncount(3)
    integer :: lon_start,lon_end,lat_start,lat_end,lev_start,lev_end
    integer :: dlonid,dlatid,dlevid
    integer :: lonid=-1,latid=-1,levid=-1
    integer :: lontype,lattype,levtype
    integer :: mlon=-1,mlat=-1,mlev=-1 
    integer*2, allocatable :: data3(:,:,:)
    integer, allocatable :: levs4(:)
    real,  allocatable :: data5(:,:,:),lons5(:),lats5(:),levs5(:)
    real*8, allocatable :: lons6(:),lats6(:),levs6(:)
    ! to test that all variables & files are on the same grid:
    real, allocatable :: lons_test(:),lats_test(:),levs_test(:)
    integer :: dmin,dmax,lev_flag
    real :: test41,test42,test4
    real*8 :: temp8
!    integer :: n0,n1,nm,nm0,tmap,ntime
    integer :: ivar,dim
    
    integer :: ilon,ilat,ilev
    integer :: i,j,ii,jj,k,l,nn,i1,i2,ifile

    character (len=8), save :: date=' '
    character (len=10), save :: time=' '
    
    if(date==' '.and.time==' ') then
       call date_and_time(date,time)
    endif
    
    vars1=vars
    if(lat1/=missing.and.lat2/=missing.and.lat2<lat1) then
       print*, 'ERROR: lat2 < lat1'
       stop
    endif
    if(lev1/=missing.and.lev2/=missing.and.lev2<lev1) then
       print*, 'ERROR: lev2 < lev1'
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
    if(allocated(xnt1)) then
       deallocate(xnt1)
    endif
    if(allocated(xntm)) then
       deallocate(xntm)
    endif
    if(allocated(lons)) then
       deallocate(lons)
    endif
    if(allocated(lats)) then
       deallocate(lats)
    endif
    if(allocated(levs)) then
       deallocate(levs)
    endif
    
    ! determine # variables & put variables in character array:
    nvar = 0
    j = 0
    do i=1,lvars
       if(vars1(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          nvar = nvar + 1
          j = 0
       endif
    enddo
    if(nvar==0) then
       print*, 'ERROR: no variables entered'
       stop
    endif
    allocate(avars(nvar))
    allocate(id(nvar))
    allocate(types(nvar))
    allocate(atts(nvar))
    allocate(vndims(nvar))
    allocate(vdims(32,nvar)) ! 32 hard wired into code
    allocate(scalef(nvar))
    allocate(add(nvar))
    allocate(miss3(nvar))
    allocate(miss5(nvar))
    miss3=-32768
    miss5=missing
    vndims = -1
    do i=1,nvar
       vars1 = adjustl(vars1)
       l = index(vars1,' ') - 1
       avars(i) = vars1(1:l)
       vars1(1:l) = ' '
    enddo
    ! put list of files needed in "temp"
    call system('ls '//file_expr//' > temp'//date//time)
    ! find number of lines in "temp":
    call system('wc temp'//date//time//' > tempn'//date//time)
    open(98,file='tempn'//date//time)
    read(98,*) nfile
    close(98)
    call system('rm tempn'//date//time)
    
    ! get file names from "temp":
    allocate( files(nfile) )
    open(98,file='temp'//date//time)
    do i=1,nfile
       read(98,'(a)') files(i)
    enddo
    close(98)
    call system('rm temp'//date//time)
    
    ! find dimension of all variables first:
    ! (go through files in a way that might minimize the # files opened)
    ! also get the calendar
    id = -1
    if(mod(nfile,nvar)==0) then
       nn = nfile/nvar
    else
       nn = nfile/nvar + 1
    endif
    do i=1,nn
       do j=0,nvar - 1
          k = nn*j + i
          if(k<=nfile) then
             !open file
             call handle_err(nf90_open(trim(files(k)),nf90_nowrite,ncid))
             ! get # variables:
             call handle_err(nf90_inquire(ncid, ndims, nvars, natts, recdim))
             ! get variable info:
             do ii=1,nvars
                ! get variable name
                call handle_err(nf90_inquire_variable(ncid,ii,name,itype,n,dims,natts))
                do jj=1,nvar
                   if(name==avars(jj)) then
                      id(jj) = ii
                      vndims(jj) = n
                   endif
                enddo
             enddo
             call handle_err(nf90_close(ncid))
             if (minval(id)>0) goto 11
          endif
       enddo
    enddo
11  continue
    ! find the number of 2-D & 3-D variables 
    nvarx1 = sum(vndims,mask=vndims==2)/2
    nvarxm = sum(vndims,mask=vndims==3)/3
    if(minval(vndims)<0) then
       print*, 'ERROR: at least one variable does not exist in files'
       stop
    endif
    if((nvarx1 + nvarxm)/=nvar) then
       print*, 'ERROR: at least one variable does not have 2 or 3 dimensions'
       stop
    endif
    allocate(varsx1(nvarx1),varsxm(nvarxm))
    ! put 1-level variables in varsx1 & multi-level in varsxm
    i1 = 0
    i2 = 0
    do i=1,nvar
       if(vndims(i)==2) then
          i1 = i1 + 1
          varsx1(i1) = avars(i)
       elseif(vndims(i)==3) then
          i2 = i2 + 1
          varsxm(i2) = avars(i)
       endif
    enddo
    lev_flag=-1
    do ifile=1,nfile
       call handle_err(nf90_open(trim(files(ifile)),nf90_nowrite,ncid))
       id = -1
       ! get # dimensions & variables:
       call handle_err(nf90_inquire(ncid, ndims, nvars, natts, recdim))
       ! get dimension sizes:
       ! first check for exact match:
       mlon = -1; mlat = -1; mlev = -1
       do i=1,ndims
          call handle_err(nf90_inquire_dimension(ncid,i,name,n))
          if(name=='lon') then
             mlon = n
             dlonid = i
          elseif(name=='lat') then
             mlat = n
             dlatid = i
          elseif(name=='lev' .or. name=='levsoi') then
             mlev = n
             dlevid = i
          endif
       enddo
       ! then check for key words within dimension name
       do i=1,ndims
          call handle_err(nf90_inquire_dimension(ncid,i,name,n))
          if(index(name,'lon')/=0 .and. mlon==-1) then
             mlon = n
             dlonid = i
          elseif(index(name,'lat')/=0 .and. mlat==-1) then
             mlat = n
             dlatid = i
          elseif(index(name,'lev')/=0 .and. mlev==-1) then
             mlev = n
             dlevid = i
          endif
       enddo
       ! get variable info:
       ! first check for exact match:
       lonid = -1; latid = -1; levid = -1
       do i=1,nvars
          call handle_err(nf90_inquire_variable(ncid,i,name,itype,n,dims,natts))
          !          print*, i,trim(name),itype,n,dims(1:n),natts,rcode
          if(name=='lon') then
             lonid = i
             namelon = name
             lontype = itype
          elseif(name=='lat') then
             latid = i
             namelat = name
             lattype = itype
          elseif(name=='lev' .or. name=='levsoi') then
             levid = i
             namelev = name
             levtype = itype
          endif
          do j=1,nvar
             if(name==avars(j)) then
                id(j) = i
                types(j) = itype
                atts(j) = natts
                vdims(1:n,j) =  dims(1:n)
                ! get add_offset, scale_factor and missing values
                do jj=1,natts
                   call handle_err(nf90_inq_attname(ncid,id(j),jj,attnam))
                   if(attnam=='add_offset') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'add_offset',add(j)))
                      elseif(attype==nf90_double) then
                         call handle_err(nf90_get_att(ncid,id(j),'add_offset',temp8))
                         add(j)=real(temp8)
                      else
                         print*, 'ERROR: add_offset is not real*4 or real*8'
                         stop
                      endif
                   endif
                   if(attnam=='scale_factor') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'scale_factor',scalef(j)))
                      elseif(attype==nf90_double) then
                         call handle_err(nf90_get_att(ncid,id(j),'scale_factor',temp8))
                         scalef(j)=real(temp8)
                      else
                         print*, 'ERROR: scale_factor is not real*4 or real*8'
                         stop
                      endif
                   endif
                   if(attnam=='missing_value') then
                      call handle_err(nf90_inquire_attribute(ncid,id(j),attnam,attype,attlen))
                      if(attype/=types(j)) then
                         print*, 'ERROR: missing_value not same type as variable'
                         stop
                      endif
                      if(attype==nf90_short) then
                         call handle_err(nf90_get_att(ncid,id(j),'missing_value',miss3(j)))
                      elseif(attype==nf90_float) then
                         call handle_err(nf90_get_att(ncid,id(j),'missing_value',miss5(j)))
                      else
                         print*, 'ERROR: variable type is not real*4 or integer*2'
                         stop
                      endif
                   endif
                enddo
             endif
          enddo
       enddo
       ! then check for key words within variable name:
       do i=1,nvars
          call handle_err(nf90_inquire_variable(ncid,i,name,itype,n,dims,natts))
          if(index(name,'lon')/=0 .and. lonid==-1) then
             lonid = i
             namelon = name
             lontype = itype
          elseif(index(name,'lat')/=0 .and. latid==-1) then
             latid = i
             namelat = name
             lattype = itype
          elseif(index(name,'lev')/=0 .and. levid==-1) then
             levid = i
             namelev = name
             levtype = itype
          endif
       enddo
       ! check that dimensions are in the default order
       do j=1,nvar
          if(id(j)>0) then
             if(vdims(1,j)/=dlonid) then
                print*, 'ERROR: longitude not first dimension in variable'
                stop
             endif
             if(vdims(2,j)/=dlatid) then
                print*, 'ERROR: latitude not second dimension in variable'
                stop
             endif
             if(vdims(3,j)/=dlevid .and. vndims(j) == 3) then
                print*, 'ERROR: level not third dimension in 3-D variable'
                stop
             endif
          endif
       enddo
       ! find lon, lat and level range
       ! first: longitudes
       if(lontype==nf90_float) then
          allocate(lons5(mlon))
          call handle_err(nf90_get_var(ncid,lonid,lons5,(/1/),(/mlon/)))
       elseif(lontype==nf90_double) then
          allocate(lons6(mlon))
          allocate(lons5(mlon))
          call handle_err(nf90_get_var(ncid,lonid,lons6,(/1/),(/mlon/)))
          do ilon=1,mlon
             lons5(ilon) = real(lons6(ilon))
          enddo
          deallocate(lons6)
       else
          print*, 'ERROR: lon type not real*4 or real*8'
          stop
       endif
       ! if lons == missing then do all longitudes
       if(lon1==missing .or. lon2==missing) then
          dmin = 1; dmax = mlon
       else
          ! make lons between 0 & 360
          do ilon=1,mlon
             if(lons5(ilon)<0.) then
                lons5(ilon) = lons5(ilon) + 360.
             elseif(lons5(ilon)>=360.) then
                lons5(ilon) = lons5(ilon) - 360.
             endif
          enddo
          test42 = -huge(0.0)
          test41 = huge(0.0)
          dmax = 0; dmin = 0
          do ilon=1,mlon
             if(lons5(ilon)<=lon2) then
                if(lons5(ilon)>test42) then
                   test42 = lons5(ilon)
                   dmax = ilon
                endif
             endif
             if(lons5(ilon)>=lon1) then
                if(lons5(ilon)<test41) then
                   test41 = lons5(ilon)
                   dmin = ilon
                endif
             endif
          enddo
          if(dmin==0 .or. dmax==0) then
             print*, 'ERROR: longitude range outside dataset'
             stop
          endif
       endif
       if(dmax>=dmin) then
          nlon = dmax - dmin + 1
       else
          nlon = (mlon - dmin + 1) + dmax
       endif
       allocate(lons_test(nlon))
       test4 = 0.
       if(lons5(dmin)>=180.) then
          test4 = -360.
       endif
       if(dmax>=dmin) then
          do ilon=dmin,dmax
             if(lons5(ilon)<180.) then
                lons_test(ilon - dmin + 1) = lons5(ilon)
             else
                lons_test(ilon - dmin + 1) = lons5(ilon) + test4
             endif
          enddo
       else
          do ilon=dmin,mlon
             lons_test(ilon - dmin + 1) = lons5(ilon) - 360.
          enddo
          do ilon=1,dmax
             lons_test(mlon - dmin + 1 + ilon) = lons5(ilon)
          enddo
       endif
       if(ifile>1) then
          do ilon=1,nlon
             if(lons_test(ilon)/=lons(ilon)) then
                print*, 'ERROR: longitude grids not the same among files'
                stop
             endif
          enddo
       else
          allocate(lons(nlon))
          lons = lons_test
       endif
       lon_start=dmin
       lon_end=dmax
       deallocate(lons5)          
       deallocate(lons_test)
       ! second: latitudes
       if(lattype==nf90_float) then
          allocate(lats5(mlat))
          call handle_err(nf90_get_var(ncid,latid,lats5,(/1/),(/mlat/)))
       elseif(lattype==nf90_double) then
          allocate(lats6(mlat))
          allocate(lats5(mlat))
          call handle_err(nf90_get_var(ncid,latid,lats6,(/1/),(/mlat/)))
          do ilat=1,mlat
             lats5(ilat) = real(lats6(ilat))
          enddo
          deallocate(lats6)
       else
          print*, 'ERROR: lat type not real*4 or real*8'
          stop
       endif
       if(lat1==missing .or. lat2==missing) then
          dmin = 1; dmax = mlat
       else
          test42 = -huge(0.0)
          test41 = huge(0.0)
          dmax = 0; dmin = 0
          do ilat=1,mlat
             if(lats5(ilat)<=lat2) then
                if(lats5(ilat)>test42) then
                   test42 = lats5(ilat)
                   dmax = ilat
                endif
             endif
             if(lats5(ilat)>=lat1) then
                if(lats5(ilat)<test41) then
                   test41 = lats5(ilat)
                   dmin = ilat
                endif
             endif
          enddo
          if(dmin==0 .or. dmax==0) then
             print*, 'ERROR: latitude range outside dataset'
             stop
          endif
       endif
       if(dmax<dmin) then
          jj=dmin
          dmin=dmax
          dmax=jj
       endif
       nlat = dmax - dmin + 1
       allocate(lats_test(nlat))
       lats_test = lats5(dmin:dmax)
       if(ifile>1) then
          do ilat=1,nlat
             if(lats_test(ilat)/=lats(ilat)) then
                print*, 'ERROR: latitude grids not the same among files'
                stop
             endif
          enddo
       else
          allocate(lats(nlat))
          lats = lats_test
       endif
       lat_start=dmin
       lat_end=dmax
       deallocate(lats5)          
       deallocate(lats_test)
       ! third: vertical levels
       if(mlev/=-1) then
          if(levtype==nf90_float) then
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs5,(/1/),(/mlev/)))
          elseif(levtype==nf90_double) then
             allocate(levs6(mlev))
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs6,(/1/),(/mlev/)))
             do ilev=1,mlev
                levs5(ilev) = real(levs6(ilev))
             enddo
             deallocate(levs6)
          elseif(levtype==nf90_int) then
             allocate(levs4(mlev))
             allocate(levs5(mlev))
             call handle_err(nf90_get_var(ncid,levid,levs4,(/1/),(/mlev/)))
             do ilev=1,mlev
                levs5(ilev) = real(levs4(ilev))
             enddo
             deallocate(levs4)
          else
             print*, 'ERROR: lev type not integer, real*4 or real*8'
             stop
          endif
          if(lev1==missing .or. lev2==missing) then
             dmin = 1; dmax = mlev
          else
             test42 = -huge(0.0)
             test41 = huge(0.0)
             dmax = 0; dmin = 0
             do ilev=1,mlev
                if(levs5(ilev)<=lev2) then
                   if(levs5(ilev)>test42) then
                      test42 = levs5(ilev)
                      dmax = ilev
                   endif
                endif
                if(levs5(ilev)>=lev1) then
                   if(levs5(ilev)<test41) then
                      test41 = levs5(ilev)
                      dmin = ilev
                   endif
                endif
             enddo
             if(dmin==0 .or. dmax==0) then
                print*, 'ERROR: vertical level range outside dataset'
                stop
             endif
          endif
          if(dmax<dmin) then
             jj=dmin
             dmin=dmax
             dmax=jj
          endif
          nlev = dmax - dmin + 1
          allocate(levs_test(nlev))
          levs_test = levs5(dmin:dmax)
          if(lev_flag>0) then
             do ilev=1,nlev
                if(levs_test(ilev)/=levs(ilev)) then
                   print*, 'ERROR: vertical levels not the same among files'
                   stop
                endif
             enddo
          else
             allocate(levs(nlev))
             levs = levs_test
             lev_flag=1
          endif
          lev_start=dmin
          lev_end=dmax
          deallocate(levs5)          
          deallocate(levs_test)
       endif
       ! get output data arrays ready:
       if( (.not.allocated(xnt1)) ) then 
          if(ave_lon) then
             allocate(xnt1(1,nlat,nvarx1))
          else
             allocate(xnt1(nlon,nlat,nvarx1))
          endif
       endif
       if( (nlev>0) .and. (.not.allocated(xntm))) then 
          if(ave_lon) then
             allocate(xntm(1,nlat,nlev,nvarxm))
          else
             allocate(xntm(nlon,nlat,nlev,nvarxm))
          endif
       endif
       ! map times in file to correct times specified in input 
       ! (but first allocate data read arrays):
       if(nlev>0) then
          allocate(data3(nlon,nlat,nlev))
          allocate(data5(nlon,nlat,nlev))
       else
          allocate(data3(nlon,nlat,1))
          allocate(data5(nlon,nlat,1))
       endif
       ! get field variables (FINALLY!!!)
       do j=1,nvar
          if(id(j)>0) then  ! check if variable is in this file
             ivar=-1
             do i=1,nvarx1  ! find # dimensions and order to put in output
                if(avars(j)==varsx1(i)) then
                   ivar=i
                   dim=2
                   exit
                endif
             enddo
             do i=1,nvarxm
                if(avars(j)==varsxm(i)) then
                   ivar=i
                   dim=3
                   exit
                endif
             enddo
             if(dim==2) then
                if(types(j)==nf90_float) then
                   start(1:dim)=(/lon_start,lat_start/)
                   ncount(1:dim)=(/nlon,nlat/)
                   if(lon_end>=lon_start) then
                      call handle_err(nf90_get_var(ncid,id(j),data5,start,ncount))
                   else
                      start(1:dim)=(/lon_start,lat_start/)
                      ncount(1:dim)=(/mlon-lon_start+1,nlat/)
                      call handle_err(nf90_get_var(ncid,id(j),data5(1:mlon-lon_start+1,:,:), &
                           start,ncount))
                      start(1:dim)=(/1,lat_start/)
                      ncount(1:dim)=(/lon_end,nlat/)
                      call handle_err(nf90_get_var(ncid,id(j),data5(nlon-lon_end+1:nlon,:,:),start,ncount))
                   endif
                   if(lats(1)>lats(nlat)) then
                      data5=data5(:,nlat:1:-1,:)
                   endif
                   if(ave_lon) then
                      data5(1,:,:) = sum(data5,1,data5/=miss5(j))/ &
                           count(data5/=miss5(j),1)
                   endif
                elseif(types(j)==nf90_short) then
                   start(1:dim)=(/lon_start,lat_start/)
                   ncount(1:dim)=(/nlon,nlat/)
                   if(lon_end>=lon_start) then
                      call handle_err(nf90_get_var(ncid,id(j),data3,start,ncount))
                   else
                      start(1:dim)=(/lon_start,lat_start/)
                      ncount(1:dim)=(/mlon-lon_start+1,nlat/)
                      call handle_err(nf90_get_var(ncid,id(j),data3(1:mlon-lon_start+1,:,:), &
                           start,ncount))
                      start(1:dim)=(/1,lat_start/)
                      ncount(1:dim)=(/lon_end,nlat/)
                      call handle_err(nf90_get_var(ncid,id(j),data3(nlon-lon_end+1:nlon,:,:), &
                           start,ncount))
                   endif
                   if(lats(1)>lats(nlat)) then
                      data3=data3(:,nlat:1:-1,:)
                   endif
                   miss5(j)=scalef(j)*miss3(j) + add(j)
                   data5 = scalef(j)*data3 + add(j)
                   if(ave_lon) then
                      data5(1,:,:) = sum(data5,1,data5/=miss5(j))/ &
                           count(data5/=miss5(j),1)
                   endif
                else
                   print*, 'ERROR: variable not of integer*2 or real*4 type'
                   stop
                endif
                if(ave_lon) then
                   xnt1(1,:,ivar) = &
                        merge(data5(1,:,1),missing,data5(1,:,1)/=miss5(j))
                else
                   xnt1(:,:,ivar) = &
                        merge(data5(:,:,1),missing,data5(:,:,1)/=miss5(j))
                endif
             else ! 3 dimensional variables:
                if(types(j)==nf90_float) then
                   start(1:dim)=(/lon_start,lat_start,lev_start/)
                   ncount(1:dim)=(/nlon,nlat,nlev/)
                   if(lon_end>=lon_start) then
                      call handle_err(nf90_get_var(ncid,id(j),data5,start,ncount))
                   else
                      start(1:dim)=(/lon_start,lat_start,lev_start/)
                      ncount(1:dim)=(/mlon-lon_start+1,nlat,nlev/)
                      call handle_err(nf90_get_var(ncid,id(j),data5(1:mlon-lon_start+1,:,:), &
                           start,ncount))
                      start(1:dim)=(/1,lat_start,lev_start/)
                      ncount(1:dim)=(/lon_end,nlat,nlev/)
                      call handle_err(nf90_get_var(ncid,id(j),data5(nlon-lon_end+1:nlon,:,:), &
                           start,ncount))
                   endif
                   if(lats(1)>lats(nlat)) then
                      data5=data5(:,nlat:1:-1,:)
                   endif
                   if(levs(1)>levs(nlev)) then
                      data5=data5(:,:,nlev:1:-1)
                   endif
                   if(ave_lon) then
                      data5(1,:,:) = sum(data5,1,data5/=miss5(j))/ &
                           count(data5/=miss5(j),1)
                   endif
                elseif(types(j)==nf90_short) then
                   start(1:dim)=(/lon_start,lat_start,lev_start/)
                   ncount(1:dim)=(/nlon,nlat,nlev/)
                   if(lon_end>=lon_start) then
                      call handle_err(nf90_get_var(ncid,id(j),data3,start,ncount))
                   else
                      start(1:dim)=(/lon_start,lat_start,lev_start/)
                      ncount(1:dim)=(/mlon-lon_start+1,nlat,nlev/)
                      call handle_err(nf90_get_var(ncid,id(j),data3(1:mlon-lon_start+1,:,:), &
                           start,ncount))
                      start(1:dim)=(/1,lat_start,lev_start/)
                      ncount(1:dim)=(/lon_end,nlat,nlev/)
                      call handle_err(nf90_get_var(ncid,id(j),data3(nlon-lon_end+1:nlon,:,:), &
                           start,ncount))
                   endif
                   if(lats(1)>lats(nlat)) then
                      data3=data3(:,nlat:1:-1,:)
                   endif
                   if(levs(1)>levs(nlev)) then
                      data3=data3(:,:,nlev:1:-1)
                   endif
                   miss5(j)=scalef(j)*miss3(j) + add(j)
                   data5 = scalef(j)*data3 + add(j)
                   if(ave_lon) then
                      data5(1,:,:) = sum(data5,1,data5/=miss5(j))/ &
                           count(data5/=miss5(j),1)
                   endif
                else
                   print*, 'ERROR: variable not of integer*2 or real*4 type'
                   stop
                endif
                if(ave_lon) then
                   xntm(1,:,:,ivar) = &
                        merge(data5(1,:,:),missing,data5(1,:,:)/=miss5(j))
                else
                   xntm(:,:,:,ivar) =  &
                        merge(data5(:,:,:),missing,data5(:,:,:)/=miss5(j))
                endif
             endif
          endif
       enddo
       
       deallocate(data3)
       deallocate(data5)
       
       call handle_err(nf90_close(ncid))
       
    enddo
    
    if(lats(1)>lats(nlat)) then
       lats=lats(nlat:1:-1)
    endif
    if(nvarxm>0) then
       if(levs(1)>levs(nlev)) then
          levs=levs(nlev:1:-1)
       endif
    endif
    
    deallocate(avars)
    deallocate(id)
    deallocate(types)
    deallocate(atts)
    deallocate(scalef)
    deallocate(add)
    deallocate(miss3)
    deallocate(miss5)
    deallocate(files)
    deallocate(vndims)
    deallocate(varsx1,varsxm)
    
  end subroutine read_no_time

  subroutine handle_err(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
    endif
  end subroutine handle_err
  
end module read_data
