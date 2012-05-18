      subroutine gen4l(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common.f'
      include 'limits.f' 
      include 'masses.f'
      include 'process.f'
      double precision p(n_final,4),r(mxdim) 
      double precision p3(4),p4(4),p5(4),p6(4),p1(4),p2(4)
      double precision p34(4),p56(4)
      double precision tau,sqrts,Qsq,xjac,taumin,lntaum
      double precision wt,wt0,wt3456,wt34,wt56,ymax,rtshat
      double precision s3456
      integer nu,i
      common/energy/sqrts
      common/Qsq/Qsq
      parameter(wt0=1d0/twopi**2)
      double precision minZsq 


      wt=0d0 
      
       minZsq=min(wsqmin,bbsqmin) 
 

 !     if((case.eq.'ZZlept').or.(case.eq.'WZbbar')) then 
         taumin=(minZsq/sqrts**2)
         lntaum=dlog(taumin)
         tau=dexp(lntaum*(one-r(9)))         
         xjac=-lntaum*tau      
         Qsq=tau*sqrts**2
         xjac=xjac*sqrts**2
        
        
         
!      elseif(case.eq.'HZZ_4l') then 
!------ Generate around s = m_h 
!         call breitw(r(9),minZsq,sqrts**2,hmass,hwidth,s3456,wt3456)
!         rtshat=dsqrt(s3456)
!         Qsq=s3456
!         ymax=dlog(sqrts/rtshat)
!         xjac=wt3456
!      endif
         
     
        
         do i=1,4
            p1(i)=0d0 
            p2(i)=0d0 
         enddo

         p1(4)=-dsqrt(Qsq)/2d0 
         p2(4)=-dsqrt(Qsq)/2d0

!      call phi1_2(r(1),r(2),r(3),r(4),Qsq,p34,p56,wt3456,*99) 
!      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99) 
!      call phi3m0(r(7),r(8),p56,p5,p6,wt56,*99)

         call phase4(r,p1,p2,p3,p4,p5,p6,wt3456,*99)
      


      do nu=1,4 
         p(1,nu)=p3(nu)
         p(2,nu)=p4(nu)
         p(3,nu)=p5(nu)
         p(4,nu)=p6(nu)
      enddo
      

!      write(6,*) p 
!      pause
      
      wt=wt0*wt3456*xjac
      return 
!      write(6,*) wt  
!      pause
 99   continue
      return  1
      end 
