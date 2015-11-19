
subroutine bilinear(x,y,nx,ny,xf,yf,nxf,nyf,z,zf,nsam,periodic)

  use constants
  implicit none

  integer :: nx,ny,nxf,nyf,jx1,jx2,jy1,jy2
  integer :: kx(nxf,2),ky(nyf,2)
  real :: x(nx),y(ny),xf(nxf),yf(nyf) ! dimensions: x & y, f == final
  real :: tt(nxf),uu(nyf)
  real :: z(nx,ny,nsam),zf(nxf,nyf,nsam),u,t
  integer iy,jy,ix,jx,isam,nsam
  real sum0,z1,z2,z3,z4
  integer periodic

  ! initialize arrays:
  do iy=1,nyf
     jy=1
     do while((y(jy)<=yf(iy)).and.(jy<=ny))
        jy=jy+1
     enddo
     if(jy==1) then
        jy1=1
        jy2=1
     elseif(jy>ny) then
        jy1=ny
        jy2=ny
     else
        jy1=jy-1
        jy2=jy
     endif
     if(jy2==jy1) then
        uu(iy)=0.
     else
        uu(iy) = (yf(iy) - y(jy1))/(y(jy2) - y(jy1))
     endif
     ky(iy,1)=jy1
     ky(iy,2)=jy2
  enddo
  do ix=1,nxf
     jx=1
     do while((x(jx)<=xf(ix)).and.(jx<=nx))
        jx=jx+1
     enddo
     if(jx==1) then
        if(periodic==1) then
           jx1=nx
           jx2=1
        else
           jx1=1
           jx2=1
        endif
     elseif(jx>nx) then
        if(periodic==1) then
           jx1=nx
           jx2=1
        else
           jx1=nx
           jx2=nx
        endif
     else
        jx1=jx-1
        jx2=jx
     endif
     if(jx2==jx1) then
        tt(ix)=0.
     else
        tt(ix) = (xf(ix) - x(jx1))/(x(jx2) - x(jx1))
     endif
     kx(ix,1)=jx1
     kx(ix,2)=jx2
  enddo

  do isam=1,nsam
     do iy=1,nyf
        do ix=1,nxf
           sum0=0.
           jx1=kx(ix,1)
           jx2=kx(ix,2)
           jy1=ky(iy,1)
           jy2=ky(iy,2)
           t=tt(ix)
           u=uu(iy)
           z1=0.
           z2=0.
           z3=0.
           z4=0.
           if(z(jx1,jy1,isam).ne.missing) then
              z1=z(jx1,jy1,isam)
              sum0 = sum0 + (1.-t)*(1.-u)
           endif
           if(z(jx2,jy1,isam).ne.missing) then
              z2=z(jx2,jy1,isam)
              sum0 = sum0 + t*(1.-u)
           endif
           if(z(jx2,jy2,isam).ne.missing) then
              z3=z(jx2,jy2,isam)
              sum0 = sum0 + t*u
           endif
           if(z(jx1,jy2,isam).ne.missing) then
              z4=z(jx1,jy2,isam)
              sum0 = sum0 + (1.-t)*u
           endif
           if(sum0.gt.0.) then
              zf(ix,iy,isam) = ((1.-t)*(1.-u)*z1 + t*(1.-u)*z2 &
              + t*u*z3 + (1.-t)*u*z4)/sum0
           else
              zf(ix,iy,isam) = missing
           endif
        enddo
     enddo
  enddo

  return

end subroutine bilinear

subroutine box(x,y,nx,ny,xf,yf,nxf,nyf,z,zf,nsam,periodic)

  use constants
  implicit none

  integer :: nx,ny,nxf,nyf
  real :: x(nx),y(ny),xf(nxf),yf(nyf) ! dimensions: x & y, f == final
  real :: x2(nx)
  real :: xfm(nxf+1),yfm(nyf+1)
  integer :: xmin(nxf),xmax(nxf),xminb(nxf),xmaxb(nxf)
  integer :: ymin(nyf),ymax(nyf),yminb(nyf),ymaxb(nyf)
  real :: z(nx,ny,nsam),zf(nxf,nyf,nsam)
  real :: testmax,testmin
  integer iy,jy,ix,jx,isam,nsam
  integer periodic
  real :: sum0,count0
  real :: weight(ny)

  ! initialize arrays (latitude):
  do iy=1,ny
     weight(iy)=cos(pi*y(iy)/180.)
  enddo
  yminb=-1
  ymaxb=-1
  ymin=-1
  ymax=-1
!  yfm(1)=-huge(0.0)
!  yfm(nyf+1)=huge(0.0)
  do iy=2,nyf
     yfm(iy) = (yf(iy) + yf(iy-1))/2.
  enddo
  yfm(1) = yfm(2) - (yfm(3) - yfm(2))
  yfm(nyf+1) = yfm(nyf) + (yfm(nyf) - yfm(nyf-1))
  do iy=1,nyf
     testmax=-huge(0.0)
     testmin=huge(0.0)
     do jy=1,ny
        if( yfm(iy)<y(jy) .and. y(jy)<yfm(iy+1)) then
           if(y(jy)<testmin) then
              ymin(iy)=jy
              testmin=y(jy)
           endif
           if(y(jy)>testmax) then
              ymax(iy)=jy
              testmax=y(jy)
           endif
        elseif(y(jy)==yfm(iy)) then
           yminb(iy)=jy
        elseif(y(jy)==yfm(iy+1)) then
           ymaxb(iy)=jy
        endif
     enddo
  enddo
  ! initialize arrays (longitude):
  if(periodic==1) then
     xminb=-1
     xmaxb=-1
     xmin=-1
     xmax=-1
     xfm(1)=(xf(1) + xf(nxf)-360.)/2.
     xfm(nxf+1)=(xf(1) + 360. + xf(nxf))/2.
     do ix=2,nxf
        xfm(ix) = (xf(ix) + xf(ix-1))/2.
     enddo
     do ix=1,nx
        if(x(ix)>xfm(nxf+1)) then
           x2(ix) = x(ix) - 360.
        elseif(x(ix)<xfm(1)) then
           x2(ix) = x(ix) + 360.
        else
           x2(ix) = x(ix)
        endif
     enddo
     do ix=1,nxf
        testmax=-huge(0.0)
        testmin=huge(0.0)
        do jx=1,nx
           if( xfm(ix)<x2(jx) .and. x2(jx)<xfm(ix+1)) then
              if(x2(jx)<testmin) then
                 xmin(ix)=jx
                 testmin=x2(jx)
              endif
              if(x2(jx)>testmax) then
                 xmax(ix)=jx
                 testmax=x2(jx)
              endif
           elseif(x2(jx)==xfm(ix)) then
              xminb(ix)=jx
              if(ix==1) then ! for periodic
                 xmaxb(nxf)=jx
              endif
           elseif(x2(jx)==xfm(ix+1)) then
              xmaxb(ix)=jx
              if(ix==nxf) then ! for periodic
                 xminb(1)=jx
              endif
           endif
        enddo
     enddo
  else
     xminb=-1
     xmaxb=-1
     xmin=-1
     xmax=-1
!     xfm(1)=-huge(0.0)
!     xfm(nxf+1)=huge(0.0)
     do ix=2,nxf
        xfm(ix) = (xf(ix) + xf(ix-1))/2.
     enddo
     xfm(1) = xfm(2) - (xfm(3) - xfm(2))
     xfm(nxf+1) = xfm(nxf) + (xfm(nxf) - xfm(nxf-1))
     do ix=1,nxf
        testmax=-huge(0.0)
        testmin=huge(0.0)
        do jx=1,nx
           if( xfm(ix)<x(jx) .and. x(jx)<xfm(ix+1)) then
              if(x(jx)<testmin) then
                 xmin(ix)=jx
                 testmin=x(jx)
              endif
              if(x(jx)>testmax) then
                 xmax(ix)=jx
                 testmax=x(jx)
              endif
           elseif(x(jx)==xfm(ix)) then
              xminb(ix)=jx
           elseif(x(jx)==xfm(ix+1)) then
              xmaxb(ix)=jx
           endif
        enddo
     enddo
  endif
  
  do isam=1,nsam
     do ix=1,nxf
        if(xmin(ix)<=xmax(ix)) then
           do iy=1,nyf
              sum0=0.
              count0=0.
              do jy=ymin(iy),ymax(iy)
                 do jx=xmin(ix),xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+weight(jy)*z(jx,jy,isam)
                       count0=count0+weight(jy)
                    endif
                 enddo
              enddo
              ! y boundaries
              if(yminb(iy)/=-1) then
                 jy=yminb(iy)
                 do jx=xmin(ix),xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              if(ymaxb(iy)/=-1) then
                 jy=ymaxb(iy)
                 do jx=xmin(ix),xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              ! x boundaries
              if(xminb(ix)/=-1) then
                 jx=xminb(ix)
                 do jy=ymin(iy),ymax(iy)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              if(xmaxb(ix)/=-1) then
                 jx=xmaxb(ix)
                 do jy=ymin(iy),ymax(iy)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              ! corners:
              if(yminb(iy)/=-1 .and. xminb(ix)/=-1) then
                 if(z(xminb(ix),yminb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xminb(ix),yminb(iy),isam)*weight(yminb(iy))
                    count0=count0+0.25*weight(yminb(iy))
                 endif
              endif
              if(yminb(iy)/=-1 .and. xmaxb(ix)/=-1) then
                 if(z(xmaxb(ix),yminb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xmaxb(ix),yminb(iy),isam)*weight(yminb(iy))
                    count0=count0+0.25*weight(yminb(iy))
                 endif
              endif
              if(ymaxb(iy)/=-1 .and. xminb(ix)/=-1) then
                 if(z(xminb(ix),ymaxb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xminb(ix),ymaxb(iy),isam)*weight(ymaxb(iy))
                    count0=count0+0.25*weight(ymaxb(iy))
                 endif
              endif
              if(ymaxb(iy)/=-1 .and. xmaxb(ix)/=-1) then
                 if(z(xmaxb(ix),ymaxb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xmaxb(ix),ymaxb(iy),isam)*weight(ymaxb(iy))
                    count0=count0+0.25*weight(ymaxb(iy))
                 endif
              endif
              if(count0/=0.) then
                 zf(ix,iy,isam) = sum0/count0
              else
                 zf(ix,iy,isam) = missing
              endif
           enddo
        else
           do iy=1,nyf
              sum0=0.
              count0=0.
              do jy=ymin(iy),ymax(iy)
                 do jx=1,xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+weight(jy)*z(jx,jy,isam)
                       count0=count0+weight(jy)
                    endif
                 enddo
                 do jx=xmin(ix),nx
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+weight(jy)*z(jx,jy,isam)
                       count0=count0+weight(jy)
                    endif
                 enddo
              enddo
              ! y boundaries
              if(yminb(iy)/=-1) then
                 jy=yminb(iy)
                 do jx=1,xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
                 do jx=xmin(ix),nx
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              if(ymaxb(iy)/=-1) then
                 jy=ymaxb(iy)
                 do jx=1,xmax(ix)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
                 do jx=xmin(ix),nx
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              ! x boundaries
              if(xminb(ix)/=-1) then
                 jx=xminb(ix)
                 do jy=ymin(iy),ymax(iy)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              if(xmaxb(ix)/=-1) then
                 jx=xmaxb(ix)
                 do jy=ymin(iy),ymax(iy)
                    if(z(jx,jy,isam)/=missing) then
                       sum0=sum0+0.5*weight(jy)*z(jx,jy,isam)
                       count0=count0+0.5*weight(jy)
                    endif
                 enddo
              endif
              ! corners:
              if(yminb(iy)/=-1 .and. xminb(ix)/=-1) then
                 if(z(xminb(ix),yminb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xminb(ix),yminb(iy),isam)*weight(yminb(iy))
                    count0=count0+0.25*weight(yminb(iy))
                 endif
              endif
              if(yminb(iy)/=-1 .and. xmaxb(ix)/=-1) then
                 if(z(xmaxb(ix),yminb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xmaxb(ix),yminb(iy),isam)*weight(yminb(iy))
                    count0=count0+0.25*weight(yminb(iy))
                 endif
              endif
              if(ymaxb(iy)/=-1 .and. xminb(ix)/=-1) then
                 if(z(xminb(ix),ymaxb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xminb(ix),ymaxb(iy),isam)*weight(ymaxb(iy))
                    count0=count0+0.25*weight(ymaxb(iy))
                 endif
              endif
              if(ymaxb(iy)/=-1 .and. xmaxb(ix)/=-1) then
                 if(z(xmaxb(ix),ymaxb(iy),isam)/=missing) then
                    sum0=sum0+0.25*z(xmaxb(ix),ymaxb(iy),isam)*weight(ymaxb(iy))
                    count0=count0+0.25*weight(ymaxb(iy))
                 endif
              endif
              if(count0/=0.) then
                 zf(ix,iy,isam) = sum0/count0
              else
                 zf(ix,iy,isam) = missing
              endif
           enddo
        endif
     enddo
  enddo

end subroutine box
