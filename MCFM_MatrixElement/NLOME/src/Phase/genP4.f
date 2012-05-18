      subroutine genP4(p,smax,wt,r1,r2,r3)
      implicit none
      double precision p(4)
      double precision smax,wt,r1,r2,r3

      call genP3(p,smax,wt,r1,r2,r3)
      p(4)=sqrt(p(1)**2+p(2)**2+p(3)**2)

      return
      end subroutine
