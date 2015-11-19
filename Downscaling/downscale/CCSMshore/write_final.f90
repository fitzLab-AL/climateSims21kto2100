!gfortran -O3 -I/Users/djlorenz/Research/netcdf/include -I/Users/djlorenz/lib -c write_final.f90
module write_final

  use netcdf
  use constants, only: missing, fillvalue, missing_short, fillvalue_short
  use strings
  implicit none

  private
  public ncopen,ncwrite,ncclose

  integer, parameter :: maxopen = 30 ! maximum number of open files
  integer, parameter :: maxvar = 60 ! maximum number of variables
  integer, save :: ncids(maxopen) = -1
  integer, save :: tsteps(maxopen), times(maxopen), years(maxopen), months(maxopen), starts(maxopen)
  integer, save :: nlons(maxopen),nlats(maxopen)
  logical, save :: shorts(maxopen)
  integer, save :: vids(maxvar,maxopen),nvars(maxopen),tids(maxopen)

contains

  ! opens a netcdf file and defines all dimensions etc to get ready for writing
  ! files are identified by integer "nfile"  (nfile >= 1 & nfile <= maxopen)
  !
  ! bare minimum: nfile,filename,vars,   all the rest optional
  ! variable names in vars separated by spaces
  ! variable units in vunits separated by spaces
  ! variable long-names in vnames separated by SEMI-COLONS
  ! presence of both sf (scale_factor) and ao (add_offset) also like a flag for integer*2
  !
  ! defaults: vnames = vars, vunits = '_'
  !
  !           if # vunits < # vars then fill in rest with last available vunits
  !              same for nvlevs, sf, ao  (sf = scale factor, ao = add offset)
  !           if # vnames < # vars then fill in rest with vars
  !
  !           nlon = 1, minlon = 1., dlon = 1.
  !           same for lat
  !
  !           default: prepare to write real*4 output; if "sf", "ao" present prepare for integer*2
  !
  subroutine ncopen(nfile,filename,vars,vunits,vnames,lons,lats, &
       title,deflate,sf,ao)

    integer :: nfile ! file number (identifies file for future calls to this and other routines)
    character (*) :: filename,vars
    character (*), optional :: vunits,vnames
    real, optional :: lons(:),lats(:)
    character (*), optional :: title
    integer, optional :: deflate ! deflate level (0 to 9, 9 = most compression, 0 = none(=default))
    ! if both of the following are present then this is also a flag for integer*2:
    real, optional :: sf(:), ao(:) ! scale factor and add offset (for integer*2 storage)

    ! work
    integer :: nvar,nvar2
    character (len=80) :: vars0(maxvar),vunits0(maxvar),vnames0(maxvar)
    character (len=1000) :: vars2,vnames2,vunits2
    real :: sf0(maxvar),ao0(maxvar)

    real, allocatable :: lons0(:),lats0(:)
    integer :: nlon0,nlat0

    integer :: ncid,deflate_level
    logical :: short0

    integer :: vdims4(4),did(4),id(4)

    integer :: j

    if(nfile > maxopen .or. nfile < 1) then
       print*, 'ERROR: ncopen: file number must be >=1 and <= ',maxopen
       stop
    endif

    if(ncids(nfile) /= -1) then
       print*, 'ERROR: ncopen: file number already taken: nfile = ',nfile
       stop
    endif

    ! find longitudes & put in allocatable variable lons0:
    if(present(lons)) then
       nlon0 = size(lons)
       allocate(lons0(nlon0))
       lons0 = lons
    else
       print*, 'ERROR: ncopen: no lons'
       stop
    endif
    nlons(nfile) = nlon0
    ! find latitudes & put in allocatable variable lats0:
    if(present(lats)) then
       nlat0 = size(lats)
       allocate(lats0(nlat0))
       lats0 = lats
    else
       print*, 'ERROR: ncopen: no lats'
       stop
    endif
    nlats(nfile) = nlat0

    vars2 = vars
    call parse(vars2,' ',vars0,nvar)
    if(nvar == 0) then
       print*, 'ERROR: ncopen: the string "vars" is empty'
       stop
    endif
    nvars(nfile) = nvar
    if(present(vunits)) then
       vunits2 = vunits
       call parse(vunits2,' ',vunits0,nvar2)
       if(nvar2 < nvar) then
          vunits0(nvar2+1:nvar) = vunits0(nvar2)
       endif
    else
       vunits0(1:nvar) = '_'
    endif
    if(present(vnames)) then
       vnames2 = vnames
       call parse(vnames2,';',vnames0,nvar2)
       if(nvar2 < nvar) then
          vnames0(nvar2+1:nvar) = vars0(nvar2+1:nvar)
       endif
    else
       vnames0 = vars0
    endif

    if(present(sf) .and. present(ao)) then
       short0 = .true.
       nvar2 = size(sf)
       sf0(1:nvar2) = sf(1:nvar2)
       if(nvar2 < nvar) then
          sf0(nvar2+1:nvar) = sf(nvar2)
       endif
       nvar2 = size(ao)
       ao0(1:nvar2) = ao(1:nvar2)
       if(nvar2 < nvar) then
          ao0(nvar2+1:nvar) = ao(nvar2)
       endif
    else
       short0 = .false. ! default
    endif
    shorts(nfile) = short0

    ! current time:
    times(nfile) = -2200
    starts(nfile) = 1

    if(present(deflate)) then
       if(deflate > 9 .or. deflate < 0) then
          print*, 'ERROR: ncopen: deflate must be >= 0 and <= 9'
          stop
       endif
       deflate_level = deflate
    else
       deflate_level = 0
    endif

    ! create nc file:
    if(deflate_level == 0) then
       call handle_err(nf90_create(trim(filename),nf90_clobber,ncid))
    else
       call handle_err(nf90_create(trim(filename),nf90_hdf5,ncid))
    endif
    ncids(nfile) = ncid
    ! define dimensions
    call handle_err(nf90_def_dim(ncid,'lon',nlon0,did(1)))
    call handle_err(nf90_def_dim(ncid,'lat',nlat0,did(2)))
    call handle_err(nf90_def_dim(ncid,'month',12,did(3)))
    call handle_err(nf90_def_dim(ncid,'time',nf90_unlimited,did(4)))
    vdims4=did
    call handle_err(nf90_def_var(ncid,'lon',nf90_float,did(1),id(1)))
    call handle_err(nf90_def_var(ncid,'lat',nf90_float,did(2),id(2)))
    call handle_err(nf90_def_var(ncid,'month',nf90_int,did(3),id(3)))
    call handle_err(nf90_def_var(ncid,'time',nf90_int,did(4),id(4)))
    tids(nfile) = id(4)
    do j=1,nvar
       if(short0) then
          call handle_err(nf90_def_var(ncid,trim(vars0(j)),nf90_short,vdims4,vids(j,nfile)))
       else
          call handle_err(nf90_def_var(ncid,trim(vars0(j)),nf90_float,vdims4,vids(j,nfile)))
       endif
       if(deflate_level /= 0) then
          call handle_err(nf90_def_var_deflate(ncid,vids(j,nfile),shuffle=1,deflate=1, &
               deflate_level=deflate_level))
       endif
       ! attributes
       call handle_err(nf90_put_att(ncid,vids(j,nfile),'units',trim(vunits0(j))))
       call handle_err(nf90_put_att(ncid,vids(j,nfile),'long_name',trim(vnames0(j))))
       if(short0) then
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'scale_factor',sf0(j)))
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'add_offset',ao0(j)))
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'missing_value',missing_short))
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'_FillValue',fillvalue_short))
       else
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'missing_value',missing))
          call handle_err(nf90_put_att(ncid,vids(j,nfile),'_FillValue',fillvalue))
       endif
    enddo
    ! global attribute
    if(present(title)) then
       call handle_err(nf90_put_att(ncid,nf90_global,'title',trim(title)))
    endif
    ! attributes for lon, lat, lev and time
    call handle_err(nf90_put_att(ncid,id(1),'units','degrees_east'))
    call handle_err(nf90_put_att(ncid,id(2),'units','degrees_north'))
    call handle_err(nf90_put_att(ncid,id(3),'units','month'))
    call handle_err(nf90_put_att(ncid,id(4),'units','decades before present'))
     call handle_err(nf90_enddef(ncid)) ! end define mode
    ! data mode:
    call handle_err(nf90_put_var(ncid,id(1),lons0))
    call handle_err(nf90_put_var(ncid,id(2),lats0))
    call handle_err(nf90_put_var(ncid,id(3),[1,2,3,4,5,6,7,8,9,10,11,12]))

    deallocate(lons0,lats0)
    
  end subroutine ncopen
  
  ! assumes writing consecutive times to netcdf file
  ! no checking of correct size!!!
  ! if you get seg fault or bus error check size of r1,rm,i1,im compared to nlon,nlat,nvar, etc
  ! Must enter all variables at once
  ! ATTENTION:
  ! this routine does not "move" to next time step until all variables have been written 
  ! for current time step. 
  subroutine ncwrite(nfile,rm,im)

    integer :: nfile
    real, optional :: rm(nlons(nfile)*nlats(nfile)*12*nvars(nfile))
    integer*2, optional :: im(nlons(nfile)*nlats(nfile)*12*nvars(nfile))
    integer :: start4(4),ncount4(4)
    integer :: j,jm,nspm

    if(ncids(nfile) == -1) then
       print*, 'ERROR: ncwrite: file with nfile = ',nfile,' not open'
       stop
    endif
    
    if(shorts(nfile)) then
       if(.not.present(im)) then
          print*, 'ERROR: ncwrite: file opened as integer*2 but im not input'
          print*, 'file number = ',nfile
          stop
       endif
    else
       if((.not.present(rm))) then
          print*, 'ERROR: ncwrite: file opened as real*4 but rm not input'
          print*, 'file number = ',nfile
          stop
       endif
    endif
    
    start4 = (/1,1,1,starts(nfile)/)
    ncount4 = (/nlons(nfile),nlats(nfile),12,1/)

    nspm = nlons(nfile)*nlats(nfile)*12
    jm = 0
    
    if(shorts(nfile)) then
       do j=1,nvars(nfile)
          call handle_err(nf90_put_var(ncids(nfile),vids(j,nfile), &
               im(nspm*jm+1:nspm*(jm+1)), &
               start=start4,count=ncount4))
          jm = jm + 1
       enddo
    else
       do j=1,nvars(nfile)
          call handle_err(nf90_put_var(ncids(nfile),vids(j,nfile), &
               rm(nspm*jm+1:nspm*(jm+1)), &
               start=start4,count=ncount4))
          jm = jm + 1
       enddo
    endif
    
    ! write time
    call handle_err(nf90_put_var(ncids(nfile),tids(nfile),(/times(nfile)/), &
         start=(/starts(nfile)/),count=(/1/)))
    ! change times to next time:
    times(nfile) = times(nfile) + 1
    starts(nfile) = starts(nfile) + 1

  end subroutine ncwrite

  subroutine ncclose(nfile)

    integer :: nfile

    if(ncids(nfile) == -1) then
       print*, 'ERROR: ncclose: file with nfile = ',nfile,' not open'
       stop
    endif

    call handle_err(nf90_close(ncids(nfile)))
    ncids(nfile) = -1

  end subroutine ncclose
  

  subroutine handle_err(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
    endif
  end subroutine handle_err

end module write_final
