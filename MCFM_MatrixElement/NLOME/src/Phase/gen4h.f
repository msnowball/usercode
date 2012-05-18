      subroutine gen4h(p,r,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'interface_settings.f' 
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(n_final,4),sqrts,rtshat
      double precision pswt,xjac
      double precision xx(2),s3456,wt3456,ymax,yave
      double precision Qsq 
      common/Qsq/Qsq
      common/energy/sqrts
     

      wt4=0d0
      
      
      call breitw(r(9),0d0,sqrts**2,hmass,hwidth,Qsq,wt3456)
    
   
      xjac=wt3456
           
    

      p1(4)=-dsqrt(Qsq)/2d0 
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=0d0
      
      p2(4)=-dsqrt(Qsq)/2d0
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=0d0 

      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999) 

      do nu=1,4
      p(1,nu)=p3(nu)
      p(2,nu)=p4(nu)
      p(3,nu)=p5(nu)
      p(4,nu)=p6(nu)     
      enddo 

      wt4=xjac*pswt/pi**3
      
      if (debug) write(6,*) 'wt4 in gen4h',wt4
      return

 999  return 1
      end
