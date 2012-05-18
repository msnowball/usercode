      subroutine gen2flat(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common_cr.f' 
      include 'limits.f' 
      include 'process.f' 
      double precision p(n_final,4),r(mxdim_cr),p1(4) 
      double precision sqrts,wt,pswt,smax,Qsq,mxmass
      integer i
      common/energy/sqrts
      common/Qsq/Qsq

      wt=1d0
      mxmass=wsqmax
      smax=mxmass/2d0
      call genP4(p1,smax,pswt,r(1),r(2),r(3))
      wt=wt*pswt
      p(1,4)=p1(4)
      do i=1,3
         p(1,i)=p1(i)
         p(2,i)=-p1(i)
      enddo
      p(2,4)=sqrt(p(2,1)**2+p(2,2)**2+p(2,3)**2)
      Qsq=(p(1,4)+p(2,4))**2
      do i=1,3
         Qsq=Qsq-(p(1,i)+p(2,i))**2
      enddo
      return
      end subroutine
