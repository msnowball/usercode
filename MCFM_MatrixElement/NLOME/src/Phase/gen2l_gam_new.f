


      subroutine gen2l_gam_new(p,r,wt,*) 
      implicit none 
      include 'constants.f'
      include 'vegas_common_cr.f' 
      include 'interface_settings.f' 
      include 'user_cuts.f' 
      include 'limits.f'
      double precision r(mxdim_cr),p(n_final,4),wt
      integer nu,i,j 
      double precision ptmin_part,etamax_part
      double precision sqrts 
      common/energy/sqrts
      double precision hmin,hmax,pbreak_part
      double precision delh,pt,etamax,y
      double precision phi,sinhy,coshy
      double precision h
      double precision wtbw,mv2,Q(4),p3(4),p34(4),pcm(4)
      double precision rdk1,rdk2,costh,sinth
      integer n2,n3
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      double precision Qsq
      double precision xx(2) 
      double precision my_w,rshat,q0sqt
      double precision xjac,lnhmin

      wt=1d0

!====== generate photon using min Pt and max eta 

      wt=wt/16d0/pi**3 

!      write(6,*) gammpt,eta_phot,mass3,width3 
!      write(6,*) wsqmin,wsqmax 
!      pause
      ptmin_part=gammpt 
      etamax_part=eta_phot 
      pbreak_part=0d0

C--- favour small pt region 
      hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak_part**2)
      hmax=1d0/dsqrt(ptmin_part**2+pbreak_part**2)
      delh=hmax-hmin
      h=hmin+r(1)*delh        
      pt=dsqrt(1d0/h**2-pbreak_part**2)
      wt=wt*delh/h**3
 
      
!       write(6,*) pt  
       
     
       etamax=sqrts/2d0/pt


      if (etamax**2 .le. 1d0) then
         write(6,*) 'etamax**2 .le. 1d0 in gen_phots_jets.f',etamax**2 
         wt=0d0
         return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))
      
      etamax=min(etamax,etamax_part)
      y=etamax*(2d0*r(2)-1d0)
   
      wt=wt*2d0*etamax
      
     
      sinhy=dsinh(y)
      coshy=dsqrt(1d0+sinhy**2)
        
      p(3,4)=pt*coshy
      phi=2d0*pi*r(3)
 
      wt=wt*2d0*pi
 
     
     
      p(3,1)=pt*dcos(phi)
      p(3,2)=pt*dsin(phi)
      p(3,3)=pt*sinhy
      

     
!======Generate BW,  
      call breitw(r(4),wsqmin,wsqmax,mass3,width3,mv2,wtbw)
 !     write(6,*) 'mv2 = ',dsqrt(mv2)
      wt=wt*wtbw/2d0/pi
      rdk1=r(5)
      rdk2=r(6)
! --- decay boson into leptons, in boson rest frame
      costh=2d0*rdk1-1d0
      sinth=dsqrt(1d0-costh**2)
      phi=2d0*pi*rdk2
      p34(4)=dsqrt(mv2)/2d0
      p34(1)=p34(4)*sinth*dcos(phi)
      p34(2)=p34(4)*sinth*dsin(phi)
      p34(3)=p34(4)*costh
     
      xx(1)=p(3,4)+dsqrt(p(3,4)**2+mv2)
      xx(1)=xx(1)/sqrts
      if(xx(1).gt.1d0) return 1 

     

      do j=1,3
         Q(j)=-p(3,j) 
      enddo
      Q(4)=-p(3,4)+xx(1)*sqrts
      
 !     write(6,*) Q(4)
     
 !     write(6,*) 'dsqrt(Qsq) = ',dsqrt(Q(4)**2-Q(3)**2-Q(2)**2-Q(1)**2)
      
c---  boost into lab frame    
      call boost(dsqrt(mv2),Q,p34,p3)

   
   
      do j=1,4
         p(1,j)=p3(j)
         p(2,j)=Q(j)-p3(j)
      enddo


      wt=wt/8d0/pi

    
     

      Qsq=0d0 
      do i=1,3 
         Qsq=Qsq-(p(1,i)+p(2,i)+p(3,i))**2
      enddo
      Qsq=Qsq+(p(1,4)+p(2,4)+p(3,4))**2 
      
     

      if(dsqrt(Qsq).gt.sqrts) return 1


!      write(6,*) p(1,1),p(1,2),p(1,3),p(1,4) 
!      write(6,*) p(2,1),p(2,2),p(2,3),p(2,4)
!      write(6,*) p(3,1),p(3,2),p(3,3),p(3,4) 
!      write(6,*) p(1,1)+p(2,1)+p(3,1),p(1,2)+p(2,2)+p(3,2),
!     & p(1,3)+p(2,3)+p(3,3),p(1,4)+p(2,4)+p(3,4)
!      pause

      return 
      end 
      
