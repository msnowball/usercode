

!--------- pT and eta in lab frame

    


      double precision function my_pt(j,p) 
      implicit none 
      include 'constants.f' 
      double precision p(mxpart,4)
      double precision taj,tjb,sab
      double precision my_dot
      double precision pa(4),pb(4) 
      common/papb/pa,pb 
      integer i,j

      my_pt=0d0 
      sab=0d0 
      taj=0d0 
      tjb=0d0 

      do i=1,3 
         sab=sab-(p(1,i)+p(2,i))**2 
         taj=taj-(p(1,i)+p(j,i))**2
         tjb=tjb-(p(2,i)+p(j,i))**2
      enddo

      sab=sab+(p(1,4)+p(2,4))**2 
      taj=taj+(p(1,4)+p(j,4))**2
      tjb=tjb+(p(2,4)+p(j,4))**2
  
     
      my_pt=dsqrt(taj*tjb/sab)

      return 
      end

      
      double precision function my_etarap(j,p) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      double precision p(mxpart,4)
      integer j 
      integer i,nu
      double precision taj,tjb
      double precision tiny,my_dot
      double precision x1,x2
      integer tag 
      double precision pa(4),pb(4)
      parameter(tiny=1d-8)
      common/tag/tag
      common/x1x2/x1,x2
!      common/papb/pa,pb

      my_etarap=1d8
      
      taj=0d0 
      tjb=0d0 

     
      do i=1,3 
         taj=taj-(p(1,i)+p(j,i))**2
         tjb=tjb-(p(2,i)+p(j,i))**2 
      enddo
      taj=taj+(p(1,4)+p(j,4))**2
      tjb=tjb+(p(2,4)+p(j,4))**2

      my_etarap=-0.5d0*dlog(x2*taj/x1/tjb)        
        
      return 
      end function 

      double precision function my_dot(a,b,p) 
      implicit none 
      include 'constants.f' 
      double precision p(mxpart,4) 
      integer a,b,i
      
      my_dot=0d0 
      
      do i=1,3 
         my_dot=my_dot-(p(a,i)+p(b,i))**2          
      enddo

      my_dot=my_dot+(p(a,4)+p(b,4))**2

      return 
      end 

  
