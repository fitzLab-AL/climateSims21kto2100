module write_data
  
  use netcdf
  implicit none
  private
  public scatter,map,line,make_nc,make_nc_short,make_ctl

contains
  
  subroutine make_nc(filename,title,lons,lats,levs,time0,ndt,dt1,ntime, &
       vars1,nvlevs,units1,vtitles1,x3,x4,x5,x6,calendar,deflate)
    
    use constants
    
    character (*) :: filename,title,vars1,units1,vtitles1,time0,dt1
    ! time0 = start time: format = year-month-day 00:00:0.0
    ! ndt = number of time increments per time step (i.e. for six hourly
    !            data ndt = 6 and dt1 = hours )
    ! dt = time increment (i.e. days, months etc.)
    character (len=10000) :: vars,units,vtitles
    character (len=40) :: dt
    character (len=1) :: ttest
    integer :: nvlevs(:)
    integer :: nlon,nlat,nlev,ndt,ntime
    real :: lons(:),lats(:),levs(:)
    real, optional :: x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:),x6(:,:,:,:,:,:)
    character (*), optional :: calendar 
    integer, optional :: deflate ! deflate level (0 to 9, 9 = most compression, 0 = none(=default))
    integer :: deflate_level
    logical :: pres(4)
    integer :: dim_size(6)

    integer nvar,nvlev
    character (len=80), allocatable :: avars(:),avtitles(:),aunits(:)
    integer, allocatable :: vlevs(:)
    integer, allocatable :: times(:)
    integer :: months(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    integer :: i,j,l,ll,iyear,imonth
    ! for netcdf:
    integer :: ncid
    integer :: start3(3),ncount3(3),start4(4),ncount4(4)
    integer :: vdims3(3),vdims4(4),did(4),id(4)
    integer, allocatable :: vid(:),vid1(:)
    integer :: nvar1

    if(present(deflate)) then
       if(deflate > 9 .or. deflate < 0) then
          print*, 'ERROR: read_nc: deflate must be >= 0 and <= 9'
          stop
       endif
       deflate_level = deflate
    else
       deflate_level = 0
    endif

    vars=vars1
    vtitles=vtitles1
    units=units1
    dt = dt1
    ! determine # variables & put variables in character array:
    nvar = 0
    j = 0
    do i=1,1000
       if(vars(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          nvar = nvar + 1
          j = 0
       endif
    enddo
    if(nvar==0) then
       print*, 'ERROR: make_nc: no variables entered'
       stop
    endif
    allocate(avars(nvar))
    allocate(avtitles(nvar))
    allocate(aunits(nvar))
    vars=adjustl(vars)
    vtitles=adjustl(vtitles)
    if(index(vtitles,';')==0) then
       ttest=' '
    else
       ttest=';'
    endif
    do i=1,nvar
       l=index(vars,' ')
       avars(i)=vars(:l)
       vars=adjustl(vars(l:))
       l=index(vtitles,ttest)
       avtitles(i)=vtitles(:l-1)
       vtitles=adjustl(vtitles(l+1:))
       if(i==1.and.avtitles(i)==' ') then
          print*, 'ERROR: make_nc: must enter long name for variable'
          stop
       elseif(avtitles(i)==' ') then
          avtitles(i)=avtitles(i-1)
       endif
       l=index(units,' ')
       aunits(i)=units(:l-1)
       units=adjustl(units(l+1:))
       if(i==1.and.aunits(i)==' ') then
          print*, 'ERROR: make_nc: must enter units for variable'
          stop
       elseif(aunits(i)==' ') then
          aunits(i)=aunits(i-1)
       endif
    enddo
    allocate(vlevs(nvar))
    allocate(vid(nvar))
    ! get size of lons, lats, levels
    nlon=size(lons)
    nlat=size(lats)
    nlev=size(levs)
    nvlev=size(nvlevs)
    if( count(nvlevs==0)==nvlev ) then
       nlev=-1
       vlevs = 0
    else
       if(nvlev==1) then
          if(nvlevs(1)==0) then
             vlevs=0
          else
             vlevs=nlev
          endif
       elseif(nvlev==nvar) then
          vlevs=merge(nvlevs,nlev,nvlevs==0)
       else
          print*, 'ERROR: make_nc: size of nvlevs /= # of variables'
          stop
       endif
    endif

    ! create times
    allocate(times(ntime))
    if(index(dt,'month')/=0) then
       l=index(time0,'-')-1
       read(time0(:l),*) iyear
       ll=index(time0(l+2:),'-') + l
       read(time0(l+2:ll),*) imonth
       times(1) = 0
       do i=2,ntime
          times(i) = times(i-1)
          do j=1,ndt
             if(imonth==2 .and. leap0(iyear)==1) then
                times(i) = times(i) + months(imonth) + 1
             else
                times(i) = times(i) + months(imonth)
             endif
             imonth =  imonth + 1
             if(imonth==13) then
                iyear =  iyear + 1
                imonth = 1
             endif
          enddo
       enddo
       dt='days'
    else
       do i=1,ntime
          times(i) = ndt*(i-1)
       enddo
    endif
    
    ! create nc file:
    if(deflate_level == 0) then
       call handle_err(nf90_create(trim(filename),nf90_clobber,ncid))
    else
       call handle_err(nf90_create(trim(filename),nf90_hdf5,ncid))
    endif
    ! define dimensions
    call handle_err(nf90_def_dim(ncid, 'lon', nlon, did(1)))
    call handle_err(nf90_def_dim(ncid, 'lat', nlat, did(2)))
    if(nlev/=-1) call handle_err(nf90_def_dim(ncid, 'lev', nlev, did(3)))
    call handle_err(nf90_def_dim(ncid, 'time', nf90_unlimited, did(4)))
    vdims3 = (/did(1),did(2),did(4)/)
    if(nlev/=-1) vdims4=did
    call handle_err(nf90_def_var(ncid,'lon',nf90_float,did(1),id(1)))
    call handle_err(nf90_def_var(ncid,'lat',nf90_float,did(2),id(2)))
    if(nlev/=-1) call handle_err(nf90_def_var(ncid,'lev',nf90_float,did(3),id(3)))
    call handle_err(nf90_def_var(ncid,'time',nf90_int,did(4),id(4)))
    do i=1,nvar
       if(vlevs(i)==0) then
          call handle_err(nf90_def_var(ncid,trim(avars(i)),nf90_float,vdims3,vid(i)))
       else
          call handle_err(nf90_def_var(ncid,trim(avars(i)),nf90_float,vdims4,vid(i)))
       endif
       if(deflate_level /= 0) then
          call handle_err(nf90_def_var_deflate(ncid,vid(i),shuffle=0,deflate=1, &
               deflate_level=deflate_level))
       endif
       ! attributes
       call handle_err(nf90_put_att(ncid,vid(i),'units',trim(aunits(i))))
       call handle_err(nf90_put_att(ncid,vid(i),'long_name',trim(avtitles(i))))
       call handle_err(nf90_put_att(ncid,vid(i),'missing_value',missing))
       call handle_err(nf90_put_att(ncid,vid(i),'_FillValue',fillvalue))
    enddo
    ! global attribute
    call handle_err(nf90_put_att(ncid,nf90_global,'title',trim(title)))
    ! attributes for lon, lat, lev and time
    call handle_err(nf90_put_att(ncid,id(1),'units','degrees_east'))
    call handle_err(nf90_put_att(ncid,id(2),'units','degrees_north'))
    if(nlev/=-1) call handle_err(nf90_put_att(ncid,id(3),'units','level'))
    write(units,'(a,a,a)') trim(dt),' since ',time0
    call handle_err(nf90_put_att(ncid,id(4),'units',trim(units)))
    if(present(calendar)) then
       call handle_err(nf90_put_att(ncid,id(4),'calendar',trim(calendar)))
    endif
    call handle_err(nf90_enddef(ncid)) ! end define mode
    ! data mode:
    call handle_err(nf90_put_var(ncid,id(1),lons))
    call handle_err(nf90_put_var(ncid,id(2),lats))
    if(nlev/=-1)  call handle_err(nf90_put_var(ncid,id(3),levs))
    call handle_err(nf90_put_var(ncid,id(4),times))
    pres=(/present(x3),present(x4),present(x5),present(x6)/)
    if( count(pres) >= 3 ) then
       print*, 'ERROR: make_nc: maximum of two out of four of '
       print*, '                x3, x4, x5 & x6 are allowed as input'
       stop
    endif
    ! two varables:
    if( any(vlevs==0) .and. any(vlevs==nlev) ) then
       if(count(pres)/=2) then
          print*, 'ERROR: make_nc: if there are both single & multi-level variables, '
          print*, '                then two of x3, x4, x5 & x6 must be input'
          stop
       endif
       nvar1=count(vlevs==0)
       allocate(vid1(nvar1))
       j=0
       do i=1,nvar
          if(vlevs(i)==0) then
             j=j+1
             vid1(j) = vid(i)
          endif
       enddo
       ! write single-level variables:
       if(present(x3)) then
          if(nvar1>1) then
             print*, 'ERROR: make_nc: if # of variables > 1, then variables must be '
             print*, '                in either arrays x4 or x5 (single level variables)'
             stop
          endif
          dim_size(1:3) = shape(x3)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x3 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x3 is not same as lats'
             stop
          endif
          start3=(/1,1,1/)
          ncount3=(/nlon,nlat,dim_size(3)/)
          if(start3(3)+ncount3(3)-1 > ntime) then
             print*, 'ERROR: make_nc: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid1(1),x3))
       elseif(present(x4)) then
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x4(:,:,i,:)))
             enddo
          elseif(nvar1==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x4))
          else
             print*, 'ERROR: make_nc: single level vars in x4 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)*dim_size(5)/)
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x5(:,:,i,:,:)))
             enddo
          elseif(nvar1==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x5))
          else
             print*, 'ERROR: make_nc: single level vars in x5 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       endif
       deallocate(vid1)
       nvar1=count(vlevs==nlev)
       allocate(vid1(nvar1))
       j=0
       do i=1,nvar
          if(vlevs(i)==nlev) then
             j=j+1
             vid1(j) = vid(i)
          endif
       enddo
       ! write multi-level variables:
       if(present(x6)) then
          dim_size(1:6) = shape(x6)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x6 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x6 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x6 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x6(:,:,:,i,:,:)))
             enddo
          elseif(nvar1==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x6))
          else
             print*, 'ERROR: make_nc: multi-level vars in x6 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x5 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x5(:,:,:,i,:)))
             enddo
          elseif(nvar1==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x5))
          else
             print*, 'ERROR: make_nc: multi-level vars in x5 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x4)) then
          if(nvar1>1) then
             print*, 'ERROR: make_nc: if # of variables > 1, then variables must be '
             print*, '                in either arrays x5 or x6 (for multi-level)'
             stop
          endif
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x4 is not same as levs'
             stop
          endif
          start4=(/1,1,1,1/)
          ncount4=(/nlon,nlat,nlev,dim_size(4)/)
          if(start4(4)+ncount4(4)-1 > ntime) then
             print*, 'ERROR: make_nc: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid1(1),x4))
       endif
       deallocate(vid1)
    ! single level variables:
    elseif(count(vlevs==0)==nvar) then
       if(count(pres(1:3))/=1) then
          print*, 'ERROR: make_nc: if there are only single level variables, '
          print*, '                then only one of x3, x4, x5 must be input'
          stop
       endif
       if(present(x6)) then
          print*, 'ERROR: make_nc: if there are only single level variables, '
          print*, '                then the variables cannot be in array x6'
          stop
       endif
       if(present(x3)) then
          if(nvar>1) then
             print*, 'ERROR: make_nc: if # of variables > 1, then variables must be '
             print*, '                in either arrays x4 or x5'
             stop
          endif
          dim_size(1:3) = shape(x3)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x3 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x3 is not same as lats'
             stop
          endif
          start3=(/1,1,1/)
          ncount3=(/nlon,nlat,dim_size(3)/)
          if(start3(3)+ncount3(3)-1 > ntime) then
             print*, 'ERROR: make_nc: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid(1),x3))
       elseif(present(x4)) then
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x4(:,:,i,:)))
             enddo
          elseif(nvar==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x4))
          else
             print*, 'ERROR: make_nc: only single level vars & x4 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x5(:,:,i,:,:)))
             enddo
          elseif(nvar==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x5))
          else
             print*, 'ERROR: make_nc: only single level vars & x5 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       endif
    ! multi-level variables:
    elseif(count(vlevs==nlev)==nvar) then
       if(count(pres(2:4))/=1) then
          print*, 'ERROR: make_nc: if there are only multi-level variables, '
          print*, '                then only one of x4, x5, x6 must be input'
          stop
       endif
       if(present(x3)) then
          print*, 'ERROR: make_nc: if there are only multi-level variables, '
          print*, '                then the variables cannot be in array x3'
          stop
       endif
       if(present(x4)) then
          if(nvar>1) then
             print*, 'ERROR: make_nc: if # of variables > 1, then variables must be '
             print*, '                in either arrays x5 or x6'
             stop
          endif
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x4 is not same as levs'
             stop
          endif
          start4=(/1,1,1,1/)
          ncount4=(/nlon,nlat,nlev,dim_size(4)/)
          if(start4(4)+ncount4(4)-1 > ntime) then
             print*, 'ERROR: make_nc: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid(1),x4))
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x5 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x5(:,:,:,i,:)))
             enddo
          elseif(nvar==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x5))
          else
             print*, 'ERROR: make_nc: only multi-level vars & x5 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x6)) then
          dim_size(1:6) = shape(x6)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc: size of dimension 1 of x6 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc: size of dimension 2 of x6 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc: size of dimension 3 of x6 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x6(:,:,:,i,:,:)))
             enddo
          elseif(nvar==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x6))
          else
             print*, 'ERROR: make_nc: only multi-level vars & x6 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       endif
    endif

    call handle_err(nf90_close(ncid))

    deallocate(avars,avtitles,aunits,vlevs,times,vid)

  contains
    
    function leap0(iyear)
      
      integer :: iyear,leap0
      
      if((mod(iyear,4)/=0)) then
         leap0 = 0
      elseif((mod(iyear,400)==0)) then
         leap0 = 1
      elseif((mod(iyear,100)==0)) then
         leap0 = 0
      else
         leap0 = 1
      endif
      
      return
    end function leap0

  end subroutine make_nc

  subroutine make_nc_short(filename,title,lons,lats,levs,time0,ndt,dt1,ntime, &
       vars1,nvlevs,units1,vtitles1,scale_factor,add_offset,x3,x4,x5,x6,calendar,deflate)
    
    use constants
    
    character (*) :: filename,title,vars1,units1,vtitles1,time0,dt1
    ! time0 = start time: format = year-month-day 00:00:0.0
    ! ndt = number of time increments per time step (i.e. for six hourly
    !            data ndt = 6 and dt1 = hours )
    ! dt = time increment (i.e. days, months etc.)
    character (len=10000) :: vars,units,vtitles
    character (len=40) :: dt
    character (len=1) :: ttest
    integer :: nvlevs(:)
    real :: scale_factor(:), add_offset(:)
    integer :: nlon,nlat,nlev,ndt,ntime
    real :: lons(:),lats(:),levs(:)
    integer*2, optional :: x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:),x6(:,:,:,:,:,:)
    character (*), optional :: calendar 
    integer, optional :: deflate ! deflate level (0 to 9, 9 = most compression, 0 = none(=default))
    integer :: deflate_level
    logical :: pres(4)
    integer :: dim_size(6)

    integer nvar,nvlev
    character (len=80), allocatable :: avars(:),avtitles(:),aunits(:)
    integer, allocatable :: vlevs(:)
    integer, allocatable :: times(:)
    integer :: months(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    integer :: i,j,l,ll,iyear,imonth
    ! for netcdf:
    integer :: ncid
    integer :: start3(3),ncount3(3),start4(4),ncount4(4)
    integer :: vdims3(3),vdims4(4),did(4),id(4)
    integer, allocatable :: vid(:),vid1(:)
    integer :: nvar1

    if(present(deflate)) then
       if(deflate > 9 .or. deflate < 0) then
          print*, 'ERROR: read_nc_short: deflate must be >= 0 and <= 9'
          stop
       endif
       deflate_level = deflate
    else
       deflate_level = 0
    endif

    vars=vars1
    vtitles=vtitles1
    units=units1
    dt = dt1
    ! determine # variables & put variables in character array:
    nvar = 0
    j = 0
    do i=1,1000
       if(vars(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          nvar = nvar + 1
          j = 0
       endif
    enddo
    if(nvar==0) then
       print*, 'ERROR: make_nc_short: no variables entered'
       stop
    endif
    i = size(scale_factor)
    if(i/=nvar) then
       print*, 'ERROR: make_nc_short: # vars /= # scale factors'
       stop
    endif
    i = size(add_offset)
    if(i/=nvar) then
       print*, 'ERROR: make_nc_short: # vars /= # add offsets'
       stop
    endif
    allocate(avars(nvar))
    allocate(avtitles(nvar))
    allocate(aunits(nvar))
    vars=adjustl(vars)
    vtitles=adjustl(vtitles)
    if(index(vtitles,';')==0) then
       ttest=' '
    else
       ttest=';'
    endif
    do i=1,nvar
       l=index(vars,' ')
       avars(i)=vars(:l)
       vars=adjustl(vars(l:))
       l=index(vtitles,ttest)
       avtitles(i)=vtitles(:l-1)
       vtitles=adjustl(vtitles(l+1:))
       if(i==1.and.avtitles(i)==' ') then
          print*, 'ERROR: make_nc_short: must enter long name for variable'
          stop
       elseif(avtitles(i)==' ') then
          avtitles(i)=avtitles(i-1)
       endif
       l=index(units,' ')
       aunits(i)=units(:l-1)
       units=adjustl(units(l+1:))
       if(i==1.and.aunits(i)==' ') then
          print*, 'ERROR: make_nc_short: must enter units for variable'
          stop
       elseif(aunits(i)==' ') then
          aunits(i)=aunits(i-1)
       endif
    enddo
    allocate(vlevs(nvar))
    allocate(vid(nvar))
    ! get size of lons, lats, levels
    nlon=size(lons)
    nlat=size(lats)
    nlev=size(levs)
    nvlev=size(nvlevs)
    if( count(nvlevs==0)==nvlev ) then
       nlev=-1
       vlevs = 0
    else
       if(nvlev==1) then
          if(nvlevs(1)==0) then
             vlevs=0
          else
             vlevs=nlev
          endif
       elseif(nvlev==nvar) then
          vlevs=merge(nvlevs,nlev,nvlevs==0)
       else
          print*, 'ERROR: make_nc_short: size of nvlevs /= # of variables'
          stop
       endif
    endif

    ! create times
    allocate(times(ntime))
    if(index(dt,'month')/=0) then
       l=index(time0,'-')-1
       read(time0(:l),*) iyear
       ll=index(time0(l+2:),'-') + l
       read(time0(l+2:ll),*) imonth
       times(1) = 0
       do i=2,ntime
          times(i) = times(i-1)
          do j=1,ndt
             if(imonth==2 .and. leap0(iyear)==1) then
                times(i) = times(i) + months(imonth) + 1
             else
                times(i) = times(i) + months(imonth)
             endif
             imonth =  imonth + 1
             if(imonth==13) then
                iyear =  iyear + 1
                imonth = 1
             endif
          enddo
       enddo
       dt='days'
    else
       do i=1,ntime
          times(i) = ndt*(i-1)
       enddo
    endif
    
    ! create nc file:
    if(deflate_level == 0) then
       call handle_err(nf90_create(trim(filename),nf90_clobber,ncid))
    else
       call handle_err(nf90_create(trim(filename),nf90_hdf5,ncid))
    endif
    ! define dimensions
    call handle_err(nf90_def_dim(ncid, 'lon', nlon, did(1)))
    call handle_err(nf90_def_dim(ncid, 'lat', nlat, did(2)))
    if(nlev/=-1) call handle_err(nf90_def_dim(ncid, 'lev', nlev, did(3)))
    call handle_err(nf90_def_dim(ncid, 'time', nf90_unlimited, did(4)))
    vdims3 = (/did(1),did(2),did(4)/)
    if(nlev/=-1) vdims4=did
    call handle_err(nf90_def_var(ncid,'lon',nf90_float,did(1),id(1)))
    call handle_err(nf90_def_var(ncid,'lat',nf90_float,did(2),id(2)))
    if(nlev/=-1) call handle_err(nf90_def_var(ncid,'lev',nf90_float,did(3),id(3)))
    call handle_err(nf90_def_var(ncid,'time',nf90_int,did(4),id(4)))
    do i=1,nvar
       if(vlevs(i)==0) then
          call handle_err(nf90_def_var(ncid,trim(avars(i)),nf90_short,vdims3,vid(i)))
       else
          call handle_err(nf90_def_var(ncid,trim(avars(i)),nf90_short,vdims4,vid(i)))
       endif
       if(deflate_level /= 0) then
          call handle_err(nf90_def_var_deflate(ncid,vid(i),shuffle=0,deflate=1, &
               deflate_level=deflate_level))
       endif
       ! attributes
       call handle_err(nf90_put_att(ncid,vid(i),'units',trim(aunits(i))))
       call handle_err(nf90_put_att(ncid,vid(i),'long_name',trim(avtitles(i))))
       call handle_err(nf90_put_att(ncid,vid(i),'scale_factor',scale_factor(i)))
       call handle_err(nf90_put_att(ncid,vid(i),'add_offset',add_offset(i)))
       call handle_err(nf90_put_att(ncid,vid(i),'missing_value',missing_short))
       call handle_err(nf90_put_att(ncid,vid(i),'_FillValue',fillvalue_short))
    enddo
    ! global attribute
    call handle_err(nf90_put_att(ncid,nf90_global,'title',trim(title)))
    ! attributes for lon, lat, lev and time
    call handle_err(nf90_put_att(ncid,id(1),'units','degrees_east'))
    call handle_err(nf90_put_att(ncid,id(2),'units','degrees_north'))
    if(nlev/=-1) call handle_err(nf90_put_att(ncid,id(3),'units','level'))
    write(units,'(a,a,a)') trim(dt),' since ',time0
    call handle_err(nf90_put_att(ncid,id(4),'units',trim(units)))
    if(present(calendar)) then
       call handle_err(nf90_put_att(ncid,id(4),'calendar',trim(calendar)))
    endif
    call handle_err(nf90_enddef(ncid)) ! end define mode
    ! data mode:
    call handle_err(nf90_put_var(ncid,id(1),lons))
    call handle_err(nf90_put_var(ncid,id(2),lats))
    if(nlev/=-1)  call handle_err(nf90_put_var(ncid,id(3),levs))
    call handle_err(nf90_put_var(ncid,id(4),times))
    pres=(/present(x3),present(x4),present(x5),present(x6)/)
    if( count(pres) >= 3 ) then
       print*, 'ERROR: make_nc_short: maximum of two out of four of '
       print*, '                x3, x4, x5 & x6 are allowed as input'
       stop
    endif
    ! two varables:
    if( any(vlevs==0) .and. any(vlevs==nlev) ) then
       if(count(pres)/=2) then
          print*, 'ERROR: make_nc_short: if there are both single & multi-level variables, '
          print*, '                then two of x3, x4, x5 & x6 must be input'
          stop
       endif
       nvar1=count(vlevs==0)
       allocate(vid1(nvar1))
       j=0
       do i=1,nvar
          if(vlevs(i)==0) then
             j=j+1
             vid1(j) = vid(i)
          endif
       enddo
       ! write single-level variables:
       if(present(x3)) then
          if(nvar1>1) then
             print*, 'ERROR: make_nc_short: if # of variables > 1, then variables must be '
             print*, '                in either arrays x4 or x5 (single level variables)'
             stop
          endif
          dim_size(1:3) = shape(x3)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x3 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x3 is not same as lats'
             stop
          endif
          start3=(/1,1,1/)
          ncount3=(/nlon,nlat,dim_size(3)/)
          if(start3(3)+ncount3(3)-1 > ntime) then
             print*, 'ERROR: make_nc_short: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid1(1),x3))
       elseif(present(x4)) then
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x4(:,:,i,:)))
             enddo
          elseif(nvar1==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x4))
          else
             print*, 'ERROR: make_nc_short: single level vars in x4 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)*dim_size(5)/)
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x5(:,:,i,:,:)))
             enddo
          elseif(nvar1==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x5))
          else
             print*, 'ERROR: make_nc_short: single level vars in x5 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       endif
       deallocate(vid1)
       nvar1=count(vlevs==nlev)
       allocate(vid1(nvar1))
       j=0
       do i=1,nvar
          if(vlevs(i)==nlev) then
             j=j+1
             vid1(j) = vid(i)
          endif
       enddo
       ! write multi-level variables:
       if(present(x6)) then
          dim_size(1:6) = shape(x6)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x6 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x6 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x6 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x6(:,:,:,i,:,:)))
             enddo
          elseif(nvar1==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x6))
          else
             print*, 'ERROR: make_nc_short: multi-level vars in x6 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x5 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar1
                call handle_err(nf90_put_var(ncid,vid1(i),x5(:,:,:,i,:)))
             enddo
          elseif(nvar1==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid1(1),x5))
          else
             print*, 'ERROR: make_nc_short: multi-level vars in x5 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x4)) then
          if(nvar1>1) then
             print*, 'ERROR: make_nc_short: if # of variables > 1, then variables must be '
             print*, '                in either arrays x5 or x6 (for multi-level)'
             stop
          endif
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x4 is not same as levs'
             stop
          endif
          start4=(/1,1,1,1/)
          ncount4=(/nlon,nlat,nlev,dim_size(4)/)
          if(start4(4)+ncount4(4)-1 > ntime) then
             print*, 'ERROR: make_nc_short: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid1(1),x4))
       endif
       deallocate(vid1)
    ! single level variables:
    elseif(count(vlevs==0)==nvar) then
       if(count(pres(1:3))/=1) then
          print*, 'ERROR: make_nc_short: if there are only single level variables, '
          print*, '                then only one of x3, x4, x5 must be input'
          stop
       endif
       if(present(x6)) then
          print*, 'ERROR: make_nc_short: if there are only single level variables, '
          print*, '                then the variables cannot be in array x6'
          stop
       endif
       if(present(x3)) then
          if(nvar>1) then
             print*, 'ERROR: make_nc_short: if # of variables > 1, then variables must be '
             print*, '                in either arrays x4 or x5'
             stop
          endif
          dim_size(1:3) = shape(x3)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x3 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x3 is not same as lats'
             stop
          endif
          start3=(/1,1,1/)
          ncount3=(/nlon,nlat,dim_size(3)/)
          if(start3(3)+ncount3(3)-1 > ntime) then
             print*, 'ERROR: make_nc_short: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid(1),x3))
       elseif(present(x4)) then
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x4(:,:,i,:)))
             enddo
          elseif(nvar==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x4))
          else
             print*, 'ERROR: make_nc_short: only single level vars & x4 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)==nvar) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x5(:,:,i,:,:)))
             enddo
          elseif(nvar==1) then
             start3=(/1,1,1/)
             ncount3=(/nlon,nlat,dim_size(3)*dim_size(4)*dim_size(5)/)
             if(start3(3)+ncount3(3)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x5))
          else
             print*, 'ERROR: make_nc_short: only single level vars & x5 =>'
             print*, '                dimension 3 = nvar or nvar==1'
             stop
          endif
       endif
    ! multi-level variables:
    elseif(count(vlevs==nlev)==nvar) then
       if(count(pres(2:4))/=1) then
          print*, 'ERROR: make_nc_short: if there are only multi-level variables, '
          print*, '                then only one of x4, x5, x6 must be input'
          stop
       endif
       if(present(x3)) then
          print*, 'ERROR: make_nc_short: if there are only multi-level variables, '
          print*, '                then the variables cannot be in array x3'
          stop
       endif
       if(present(x4)) then
          if(nvar>1) then
             print*, 'ERROR: make_nc_short: if # of variables > 1, then variables must be '
             print*, '                in either arrays x5 or x6'
             stop
          endif
          dim_size(1:4) = shape(x4)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x4 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x4 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x4 is not same as levs'
             stop
          endif
          start4=(/1,1,1,1/)
          ncount4=(/nlon,nlat,nlev,dim_size(4)/)
          if(start4(4)+ncount4(4)-1 > ntime) then
             print*, 'ERROR: make_nc_short: times in variables > ntime' 
             stop
          endif
          call handle_err(nf90_put_var(ncid,vid(1),x4))
       elseif(present(x5)) then
          dim_size(1:5) = shape(x5)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x5 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x5 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x5 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x5(:,:,:,i,:)))
             enddo
          elseif(nvar==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x5))
          else
             print*, 'ERROR: make_nc_short: only multi-level vars & x5 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       elseif(present(x6)) then
          dim_size(1:6) = shape(x6)
          if(dim_size(1)/=nlon) then
             print*, 'ERROR: make_nc_short: size of dimension 1 of x6 is not same as lons'
             stop
          endif
          if(dim_size(2)/=nlat) then
             print*, 'ERROR: make_nc_short: size of dimension 2 of x6 is not same as lats'
             stop
          endif
          if(dim_size(3)/=nlev) then
             print*, 'ERROR: make_nc_short: size of dimension 3 of x6 is not same as levs'
             stop
          endif
          if(dim_size(4)==nvar) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             do i=1,nvar
                call handle_err(nf90_put_var(ncid,vid(i),x6(:,:,:,i,:,:)))
             enddo
          elseif(nvar==1) then
             start4=(/1,1,1,1/)
             ncount4=(/nlon,nlat,nlev,dim_size(4)*dim_size(5)*dim_size(6)/)
             if(start4(4)+ncount4(4)-1 > ntime) then
                print*, 'ERROR: make_nc_short: times in variables > ntime' 
                stop
             endif
             call handle_err(nf90_put_var(ncid,vid(1),x6))
          else
             print*, 'ERROR: make_nc_short: only multi-level vars & x6 =>'
             print*, '                dimension 4 = nvar or nvar==1'
             stop
          endif
       endif
    endif

    call handle_err(nf90_close(ncid))

    deallocate(avars,avtitles,aunits,vlevs,times,vid)

  contains
    
    function leap0(iyear)
      
      integer :: iyear,leap0
      
      if((mod(iyear,4)/=0)) then
         leap0 = 0
      elseif((mod(iyear,400)==0)) then
         leap0 = 1
      elseif((mod(iyear,100)==0)) then
         leap0 = 0
      else
         leap0 = 1
      endif
      
      return
    end function leap0

  end subroutine make_nc_short

  ! make control file for GrADS:
  subroutine make_ctl(ctlfile,datafile,title,options1,lons,lats,levs,times, &
       vars1,vlevs1,vtitles1)
    
    use constants
    
    character (*) :: ctlfile,datafile,title,options1,times,vars1,vtitles1
    character (2000) :: options,vars,vlevs,vtitles
    character (len=1) :: ttest
    real :: lons(:),lats(:),levs(:)
    integer :: vlevs1(:)
    integer :: nlon,nlat,nlev,nvlev
    real :: dlon,dlat,dlev,dl,dl1,dl2
    real :: slon,slat,slev
    
    integer nvar,noption
    character (len=80), allocatable :: avars(:),avlevs(:),avtitles(:),aoptions(:)
    
    integer :: i,j,l,ll,il,jl,maxl1,maxl2,maxl3
    character (len=2) :: iformat,if1,if2,if3
    character (len=500) :: work
    character (len=12) :: temp

    nvlev=size(vlevs1)
    write(temp,'(a,i8,a)') '(',nvlev,'i8)'
    write(vlevs,temp) vlevs1
    vars=vars1
    vtitles=vtitles1
    options=options1
    ! determine # variables & put variables in character array:
    nvar = 0
    j = 0
    do i=1,1000
       if(vars(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          nvar = nvar + 1
          j = 0
       endif
    enddo
    if(nvar==0) then
       print*, 'ERROR: make_ctl: no variables entered'
       stop
    endif
    allocate(avars(nvar))
    allocate(avlevs(nvar))
    allocate(avtitles(nvar))
    vars=adjustl(vars)
    vlevs=adjustl(vlevs)
    vtitles=adjustl(vtitles)
    if(index(vtitles,';')==0) then
       ttest=' '
    else
       ttest=';'
    endif
    maxl1=0
    maxl2=0
    maxl3=0
    do i=1,nvar
       l=index(vars,' ')
       avars(i)=vars(:l)
       vars=adjustl(vars(l:))
       l=index(vlevs,' ')
       avlevs(i)=vlevs(:l)
       vlevs=adjustl(vlevs(l:))
       if(i==1.and.avlevs(i)==' ') then
          print*, 'ERROR: make_ctl: number of levels for variables must be entered'
          stop
       elseif(avlevs(i)==' ') then
          avlevs(i)=avlevs(i-1)
       endif
       l=index(vtitles,ttest)
       avtitles(i)=vtitles(:l-1)
       vtitles=adjustl(vtitles(l+1:))
       if(i==1.and.avtitles(i)==' ') then
          print*, 'ERROR: make_ctl: must enter variable comments'
          stop
       elseif(avtitles(i)==' ') then
          avtitles(i)=avtitles(i-1)
       endif
       ! find max length:
       ll=len_trim(avars(i))
       if(ll>maxl1) maxl1=ll
       ll=len_trim(avlevs(i))
       if(ll>maxl2) maxl2=ll
       ll=len_trim(avtitles(i))
       if(ll>maxl3) maxl3=ll
    enddo
    ! determine # options & put options in character array:
    noption = 0
    j = 0
    do i=1,1000
       if(options(i:i)/=' ') then
          j = 1
       elseif(j==1) then
          noption = noption + 1
          j = 0
       endif
    enddo
    if(noption==0) then
       print*, 'ERROR: make_ctl: no options entered'
       stop
    endif
    allocate(aoptions(noption))
    options=adjustl(options)
    do i=1,noption
       l=index(options,' ')
       aoptions(i)=options(:l)
       options=adjustl(options(l:))
    enddo
    ! determine nlon and whether lons are 'LINEAR' or 'LEVELS'
    nlon=size(lons)
    slon=lons(1)
    if(nlon>1) then
       dlon=lons(2)-lons(1)
    else
       dlon=1.
    endif
    do i=3,nlon
       dl=lons(i)-lons(i-1)
       if(dl/=dlon) then
          dlon=missing
          exit
       endif
    enddo
    ! determine nlat and whether lats are 'LINEAR' or 'LEVELS'
    nlat=size(lats)
    slat=lats(1)
    if(nlat>1) then
       dlat=lats(2)-lats(1)
    else
       dlat=1.
    endif
    do i=3,nlat
       dl=lats(i)-lats(i-1)
       if(dl/=dlat) then
          dlat=missing
          exit
       endif
    enddo
    ! determine nlev and whether levs are 'LINEAR' or 'LEVELS'
    nlev=size(levs)
    slev=levs(1)
    if(nlev>1) then
       dlev=levs(2)-levs(1)
    else
       dlev=1.
    endif
    if(dlev>0.) then
       do i=3,nlev
          dl=levs(i)-levs(i-1)
          if(dl/=dlev) then
             dlev=missing
             exit
          endif
       enddo
    else
       dlev = missing
    endif

    open(99,file=ctlfile)
    
    write(99,'(a,t10,a)') 'DSET',trim(adjustl(datafile))
    write(99,'(a,t10,a)') 'TITLE',trim(adjustl(title))
    ! missing values
    write(work,'(a,t10,1PE13.6)') 'UNDEF',missing
    l=index(work,'E',.true.)
    work(l:l)='e'
    l=l-1
    do while(work(l:l)=='0')
       work=work(:l-1)//work(l+1:)
       l=l-1
    enddo
    ! options
    write(99,'(a)') trim(work)    
    do i=1,noption
       write(99,'(a,t10,a)') 'OPTIONS',trim(adjustl(aoptions(i)))
    enddo
    ! longitudes
    if(dlon/=missing) then
       write(iformat,'(a,i1)') 'i',int(log10(nlon+0.5)+1)
       write(work,'(a,t10,'//iformat//'a,F18.6,F18.6)') 'XDEF',nlon,' LINEAR ',slon,dlon
       l=len_trim(work)
       ll=index(work,'LINEAR')
       il=l
       do while(il>ll)
          jl=il
          do while(work(jl:jl)=='0'.or.work(jl:jl)==' ')
             work=work(:jl-1)//work(jl+1:)
             jl=jl-1
          enddo
          il=index(work(:jl),' ',.true.)-1
       enddo
       write(99,'(a)') trim(work)
    else
       write(iformat,'(a,i1)') 'i',int(log10(nlon+0.5)+1)
       write(99,'(a,t10,'//iformat//'a)') 'XDEF',nlon,' LEVELS'
       i=1
       do while(i<=nlon)
          j=min(i+4,nlon)
          write(iformat,'(a,i1)') ' ',j-i+1
          write(99,'('//iformat//'(1PG15.7))') lons(i:j)
          i=j+1
       enddo
    endif
    ! latitudes
    if(dlat/=missing) then
       write(iformat,'(a,i1)') 'i',int(log10(nlat+0.5)+1)
       write(work,'(a,t10,'//iformat//'a,F18.6,F18.6)') 'YDEF',nlat,' LINEAR ',slat,dlat
       l=len_trim(work)
       ll=index(work,'LINEAR')
       il=l
       do while(il>ll)
          jl=il
          do while(work(jl:jl)=='0'.or.work(jl:jl)==' ')
             work=work(:jl-1)//work(jl+1:)
             jl=jl-1
          enddo
          il=index(work(:jl),' ',.true.)-1
       enddo
       write(99,'(a)') trim(work)
    else
       write(iformat,'(a,i1)') 'i',int(log10(nlat+0.5)+1)
       write(99,'(a,t10,'//iformat//'a)') 'YDEF',nlat,' LEVELS'
       i=1
       do while(i<=nlat)
          j=min(i+4,nlat)
          write(iformat,'(a,i1)') ' ',j-i+1
          write(99,'('//iformat//'(1PG15.7))') lats(i:j)
          i=j+1
       enddo
    endif
    ! levels
    if(dlev/=missing) then
       write(iformat,'(a,i1)') 'i',int(log10(nlev+0.5)+1)
       write(work,'(a,t10,'//iformat//'a,F18.6,F18.6)') 'ZDEF',nlev,' LINEAR ',slev,dlev
       l=len_trim(work)
       ll=index(work,'LINEAR')
       il=l
       do while(il>ll)
          jl=il
          do while(work(jl:jl)=='0'.or.work(jl:jl)==' ')
             work=work(:jl-1)//work(jl+1:)
             jl=jl-1
          enddo
          il=index(work(:jl),' ',.true.)-1
       enddo
       write(99,'(a)') trim(work)
    else
       write(iformat,'(a,i1)') 'i',int(log10(nlev+0.5)+1)
       write(99,'(a,t10,'//iformat//'a)') 'ZDEF',nlev,' LEVELS'
       i=1
       do while(i<=nlev)
          j=min(i+4,nlev)
          write(iformat,'(a,i1)') ' ',j-i+1
          write(99,'('//iformat//'(1PG15.7))') levs(i:j)
          i=j+1
       enddo
    endif
    ! times
    write(99,'(a,t10,a)') 'TDEF',trim(adjustl(times))
    ! variables
    write(iformat,'(a,i1)') 'i',int(log10(nvar+0.5)+1)
    write(99,'(a,t10,'//iformat//')') 'VARS',nvar
    ! use max lengths
    write(if1,'(i2)') maxl1+2
    write(if2,'(i2)') maxl1+maxl2+3
    write(if3,'(i2)') maxl1+maxl2+5
    do i=1,nvar
       write(99,'(a,t'//if1//',a,t'//if2//',a,t'//if3//',a)') & 
            trim(avars(i)),trim(avlevs(i)),'0',trim(avtitles(i))
    enddo
    write(99,'(a)') 'ENDVARS'
    
    close(99)    
    
    deallocate(avars,avlevs,avtitles,aoptions)

  end subroutine make_ctl

  subroutine handle_err(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
    endif
  end subroutine handle_err
  
end module write_data
