

      subroutine gen2l_gam(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common_cr.f' 
      include 'process.f' 
      include 'limits.f' 
      include 'cutparams.f' 
      include 'user_cuts.f'
      double precision p(n_final,4),r(mxdim_cr) 
      double precision sqrts,taumin,wt 
      double precision lntaum,tau,xjac
      double precision Qsq,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision pswt,rdk1,rdk2
      double precision mllga,rQsq
      integer i,nu
      common/energy/sqrts

      mllga=max(mll_s,mllgam_min)
      mllga=max(mllga,1d-2)
      mllga=mllga**2
      

      wt=0d0
      taumin=mllga/sqrts**2
!     taumin=1d-5
      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(4)))         
      xjac=-lntaum*tau      
      Qsq=tau*sqrts**2
      xjac=xjac*sqrts**2
   
      rQsq=dsqrt(max(Qsq,0d0))
      p1(3)=-rQsq*half 
      p1(4)=-rQsq*half 
      p2(3)=rQsq*half 
      p2(4)=-rQsq*half
      
      do i=1,2
         p1(i)=zip
         p2(i)=zip
      enddo

      rdk1=r(5)
      rdk2=r(6)
     

      call phase3(r(1),r(2),r(3),rdk1,rdk2,p1,p2,p3,p4,p5,p6,p7,pswt)

 

      wt=xjac*pswt
      if(wt.eq.0d0) return 1 
      
      do nu=1,4 
         p(1,nu)=p3(nu) 
         p(2,nu)=p4(nu) 
         p(3,nu)=p5(nu) 
      enddo

    
      return 
      end 
