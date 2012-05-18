!---- So simple it should need no explanation! 
!---- Sets epinv1 to e1 and epinv2 to e2 
!-----Ciaran Williams 2011

      subroutine set_epinv(e1) 
      implicit none 
      include 'epinv.f' 
      include 'epinv2.f' 
      double precision e1 
      
      
      epinv=e1 
      epinv2=epinv
      
      return 
      end subroutine 
      
