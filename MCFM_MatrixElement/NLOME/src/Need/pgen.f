      
      subroutine pgen(pin,pout) 
      implicit none 
      include 'interface_settings.f' 
      include 'constants.f' 
      double precision pin(n_particles,4),pout(mxpart,4)
      integer i,j 
      
      !-- Initailse MCFM matrix 
      do i=1,mxpart 
         do j=1,4
            pout(i,j)=0d0
         enddo
      enddo

      !-- Copy accross elements 
!--- Input (E,x,y,z), MCFM: (x,y,z,E) 
      do i=1,n_particles  
         pout(i,1)=pin(i,2) 
         pout(i,2)=pin(i,3) 
         pout(i,3)=pin(i,4) 
         pout(i,4)=pin(i,1) 
      enddo
         
      
      return 
      end subroutine
