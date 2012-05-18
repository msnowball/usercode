      subroutine genP3(p,smax,wt,r1,r2,r3)
      implicit none
      integer i1,i2,i3
      double precision p(4)
      double precision smax,wt,r1,r2,r3

      wt=1d0
      p(3)=smax*(2d0*r1-1d0)
      wt=wt*2d0*smax
      smax=sqrt(smax*smax-p(3)**2)
      p(2)=smax*(2d0*r2-1d0)
      wt=wt*2d0*smax
      smax=sqrt(smax*smax-p(2)**2)
      p(1)=smax*(2d0*r3-1d0)
      wt=wt*2d0*smax
      smax=sqrt(smax*smax-p(1)**2)

      return
      end subroutine
