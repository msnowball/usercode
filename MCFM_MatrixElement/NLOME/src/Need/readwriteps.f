c--- subroutines to write out PS points or read them in

c--- p(mxpart,4): standard MCFM momentum array
c--- nmom: number of momentum entries to write out
c--- pswt: phase space weight of the point
      subroutine writeps(p,nmom,pswt)
      implicit none
      include 'constants.f'
      include 'outputstring.f' 
 !     integer mxpart
 !     parameter(mxpart=12)
      integer nmom,zero
      double precision p(mxpart,4),pswt
      character*200 strformat
      character*14 file_name
      logical first
      data first/.true./
      save first
      
    
      
      file_name=outputstring//'.DAT'
       
      if (first) then
        open(unit=51,file=file_name,status='unknown')
       
        first=.false.
      endif
      
   

      if ((nmom .ge. 1) .and. (nmom .le. 9)) then
	 zero=ichar('0')
	 strformat='('//
     &    char(nmom+zero)//
     &    '(E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X),E20.12)'
      else
        write(6,*) 'Unexpected number of momenta to write out: ',nmom
	write(6,*) '(must be greater than zero and <10 at present'
	stop
      endif
      
      write(6,strformat) p(1:nmom,1),p(1:nmom,2),
     &                    p(1:nmom,3),p(1:nmom,4),pswt
      
      
      write(51,strformat) p(1:nmom,1),p(1:nmom,2),
     &                    p(1:nmom,3),p(1:nmom,4),pswt
      
      return
      end
      
c--- this routine returns the array p and pswt, given input nmom
c--- returns flag 'nomore' TRUE if the end of the file is reached, else false
      subroutine readps(p,nmom,pswt,nomore)
      implicit none
      include 'constants.f'
      include 'outputstring.f' 
!      integer mxpart
!      parameter(mxpart=12)
      integer nmom,zero,errno
      double precision p(mxpart,4),pswt
      character*200 strformat
      character*14 file_name
      logical first,nomore
      data first/.true./
      save first
      
      file_name=outputstring//'.DAT'
      nomore=.false.
      if (first) then
        open(unit=51,file=file_name,status='unknown')
      endif
      
      if ((nmom .ge. 1) .and. (nmom .le. 9)) then
	 zero=ichar('0')
	 strformat='('//
     &    char(nmom+zero)//
     &    '(E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X),E20.12)'
      else
        write(6,*) 'Unexpected number of momenta to read in: ',nmom
	write(6,*) '(must be greater than zero and <10 at present'
	stop
      endif
      
      read(51,strformat,iostat=errno) p(1:nmom,1),p(1:nmom,2),
     &                          p(1:nmom,3),p(1:nmom,4),pswt
      
      if (errno .eq. -1) then 
         nomore=.true.
         close(51) 
      endif
      return
      
      end
      
      subroutine closeps(unit) 
      implicit none 
      integer unit

      close(unit) 
      return 
      end 


