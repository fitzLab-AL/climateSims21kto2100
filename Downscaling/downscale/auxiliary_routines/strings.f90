module strings
  
contains
  
  subroutine parse(str,delims,args,nargs)
    
    implicit none
    
    ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    ! the delimiters contained in the string 'delims'. 
    ! The integer output variable nargs contains the number of arguments found.
    
    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:) :: args
    integer :: nargs,na,i,lenstr
    
    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
       args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    
    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       if(nargs > na) then
          print*, 'ERROR: parse: nargs > size of args'
          stop
       endif
       call split(str,delims,args(nargs))
    end do
    str=strsav
    
  end subroutine parse
  
  !**********************************************************************
  
  subroutine compact(str)
    
    implicit none
    
    ! Converts multiple spaces and tabs to single spaces; deletes control characters;
    ! removes initial spaces.
    
    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str)):: outstr
    integer :: lenstr,isp,k,i,ich
    
    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    isp=0
    k=0
    
    do i=1,lenstr
       ch=str(i:i)
       ich=iachar(ch)
       
       select case(ich)
          
       case(9,32)     ! space or tab character
          if(isp==0) then
             k=k+1
             outstr(k:k)=' '
          end if
          isp=1
          
       case(33:)      ! not a space, quote, or control character
          k=k+1
          outstr(k:k)=ch
          isp=0
          
       end select
       
    end do
    
    str=adjustl(outstr)
    
  end subroutine compact
  
  !**********************************************************************
  subroutine split(str,delims,before,sep)
    
    implicit none
    
    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the 
    ! found delimiter. 
    
    character(len=*) :: str,delims,before
    character,optional :: sep
    logical :: pres
    character :: ch,cha
    integer :: lenstr,k,ibsl,i,ipos,iposa
    
    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if(lenstr == 0) return        ! string str is empty
    k=0
    ibsl=0                        ! backslash initially inactive
    before=' '
    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          before(k:k)=ch
          ibsl=0
          cycle
       end if
       ipos=index(delims,ch)         
       if(ipos == 0) then          ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
       end if
       if(ch /= ' ') then          ! character is a delimiter that is not a space
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
       cha=str(i+1:i+1)            ! character is a space delimiter
       iposa=index(delims,cha)
       if(iposa > 0) then          ! next character is a delimiter
          str=str(i+2:)
          if(pres) sep=cha
          exit
       else
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
    end do
    if(i >= lenstr) str=''
    str=adjustl(str)              ! remove initial spaces
    return
    
  end subroutine split
  
end module strings
