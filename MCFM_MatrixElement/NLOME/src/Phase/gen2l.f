      subroutine gen2l(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common_cr.f' 
      include 'limits.f'
      include 'masses.f'
      include 'process.f' 
      double precision p(n_final,4),r(mxdim_cr) 
      double precision p34(4),tau,sqrts,Qsq,xjac,p3(4),p4(4),wt34,wt
      integer nu
      double precision wt0,taumin,lntaum
      common/energy/sqrts
      common/Qsq/Qsq

      parameter(wt0=1/(16d0*pi))
      
      wt=0d0
!------Old phase space
!      taumin=wsqmin/sqrts**2
!      lntaum=dlog(taumin)
!      tau=dexp(lntaum*(one-r(3)))         
!      xjac=-lntaum*tau      
!      Qsq=tau*sqrts**2
      
      call breitw(r(3),wsqmin,wsqmax,zmass,zwidth,Qsq,xjac)
!      xjac=xjac*sqrts**2
       do nu=1,3 
          p34(nu)=0d0 
      enddo
      p34(4)=dsqrt(max(Qsq,0d0)) 


      call phi3m0(r(1),r(2),p34,p3,p4,wt34,*99) 

      do nu=1,4 
         p(1,nu)=p3(nu)
         p(2,nu)=p4(nu) 
      enddo
      wt=xjac*wt34
!     wt=xjac

 99   continue
      return 
      end subroutine
