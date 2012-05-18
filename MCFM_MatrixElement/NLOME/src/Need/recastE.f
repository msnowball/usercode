
      subroutine recastE(pin) 
      implicit none 
      include 'interface_settings.f'
      double precision pin(n_final,4) 
      logical Erec 
      double precision p_sq
      integer i,j

!----- Erec controls whether we recast E or |p| 
      Erec=.true.
!
!      write(6,*) '***** Before Recasting ********'
!      
!      do i=1,4
!         write(6,*) 'pin (',i,') = ',pin(i,1),pin(i,2)
!     &        ,pin(i,3),pin(i,4)
!         write(6,*) 'p**2 = ',pin(i,4)**2-(pin(i,1)**2
!!     &        +pin(i,2)**2+pin(i,3)**2)
!      enddo
 
    

      if(Erec) then 
!----- recast E
         do i=1,n_final
            p_sq=0d0
            do j=1,3
               p_sq=p_sq+pin(i,j)**2
!-------- make E such that mass=0d0 
               pin(i,4)=dsqrt(p_sq) 
            enddo
         enddo
      else
         write(6,*) 'Not yet implemented' 
         stop
      endif

      
!      write(6,*) '***** After Recasting ********'!
!
!      do i=1,4
!         write(6,*) 'pin (',i,') = ',pin(i,1),pin(i,2)
!     &        ,pin(i,3),pin(i,4)
!         write(6,*) 'p**2 = ',pin(i,4)**2-(pin(i,1)**2
!     &        +pin(i,2)**2+pin(i,3)**2)
!      enddo
 
    

      return 
      end 
