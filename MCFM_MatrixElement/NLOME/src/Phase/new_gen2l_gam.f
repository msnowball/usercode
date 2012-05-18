



      subroutine new_gen2l_gam(p,r,wt,*) 
      implicit none
      include 'constants.f'
      include 'interface_settings.f'
      include 'vegas_common_cr.f'
      include 'user_cuts.f' 
      include 'limits.f'
      integer i,j,k,n2,n3
      double precision mv2
      double precision r(mxdim_cr),p(n_final,4),wt
      double precision pl1(4),pl2(4),pgam(4),pt
      double precision pz(4),Q(4)
      double precision rdk1,rdk2,costh,sinth
      double precision hmax,hmin,delh 
      double precision rQsq_min,Qsq_min
      double precision taumin,lntaum,tau
      double precision Qsq,xjac,sqrts
      common/energy/sqrts
      double precision ptmin_part,pbreak_part,etamax
      double precision coshy,sinhy,phi
      double precision h,y,etamax_part
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      double precision wtbw,rQsq
      double precision x_th,y_th,z_th,theta,rxyz
      double precision a,b,p_in(4)

      do i=1,3
         p_in(i)=0d0
      enddo

      wt=1d0

!====== this routine will generate  

!------ int d tau d^4 p_gam d^4 p_l1 d^4 p_l2 delta(p_l1^2)delta(p_l2^2)delta(p_gam^2)delta^4(Q+p_gam+p_l1+p_l2)
!------ 
!------ 4+4+4+1 = 13 integration variables 
!------ 3+4=7 constraints => 6 variables of integration 
!======= use Q_sq, (p_l1+p_l2)^2 (with breit wigner)  

!====== Fix p_l2 using delta function removes 4 variables, uses 4 deltas

!-------------------------------------------------------------
!- Qsq
!-------------------------------------------------------------

!======= Q_sq min is either mll or mllgam if user has defined this cut
      rQsq_min=max(mll_s,mllgam_min)
      Qsq_min=rQsq_min**2
      
      taumin=Qsq_min/sqrts**2
!     taumin=1d-5
      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(1)))    
      xjac=-lntaum*tau      
      Qsq=tau*sqrts**2
      rQsq=dsqrt(Qsq)
      xjac=xjac*sqrts**2  
      wt=wt*xjac

      p_in(4)=rQsq

!--------------------------------------------------------------
!-  pt (gamma) 
!-------------------------------------------------------------- 

      ptmin_part=gammpt 
      etamax_part=eta_phot 
      pbreak_part=0d0

C--- favour small pt region 
      hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak_part**2)
      hmax=1d0/dsqrt(ptmin_part**2+pbreak_part**2)
      delh=hmax-hmin
      h=hmin+r(2)*delh        
      pt=dsqrt(1d0/h**2-pbreak_part**2)
      wt=wt*delh/h**3

      etamax=sqrts/2d0/pt

      
      if (etamax**2 .le. 1d0) then
         write(6,*) 'etamax**2 .le. 1d0 in gen_phots_jets.f',etamax**2 
         wt=0d0
         return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))
      
      etamax=min(etamax,etamax_part)
      y=etamax*(2d0*r(3)-1d0)   
      wt=wt*2d0*etamax
      
     
      sinhy=dsinh(y)
      coshy=dsqrt(1d0+sinhy**2)
        
      p(3,4)=pt*coshy
      phi=2d0*pi*r(4)
 
      wt=wt*2d0*pi
 
      p(3,1)=pt*dcos(phi)
      p(3,2)=pt*dsin(phi)
      p(3,3)=pt*sinhy
      
      do i=1,4
         pgam(i)=p(3,i) 
      enddo
      
!------------------------------------------------------------------
!- Now have generated a photon and a total invariant mass 
!- 
!- Have 2 variables left, generate a invariant mass from a breit
!- wigner, 
!------------------------------------------------------------------
 
      call breitw(r(5),wsqmin,wsqmax,mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2d0/pi

!----- Now need to work out jacobian and constraint from delta 

!----- Decay boson to leptons in boson rest frame.       
      phi=2d0*pi*r(6) 
      
!----- Now costh is fixed by delta function
!      costh=2d0*rdk1-1d0
!------ delta function is delta(p_l1^4) = delta((Q-p_l1-p_gam)^2)
     
!---- parametrize as sinth*x+costh*y=z 
!
!---  x = sqrt(mv_2)*(cos(phi)*pgam_x+sin(phi)*pgam_y)
!---  y = sqrt(mv_2)*(pgam_z) 
!---  z = 1/2*(2d0*E_gam*Sqrt(Q)-Qsq)+sqrt(mv_2)*(E_gam+Sqrt(Q))
      

      x_th=dsqrt(mv2)*(dcos(phi)*pgam(1)+dsin(phi)*pgam(2))
      y_th=dsqrt(mv2)*(pgam(3))
      z_th=1d0/2d0*(2d0*pgam(4)*rQsq-Qsq)+dsqrt(mv2)*(pgam(4)+rQsq)

!------ check root 
      rxyz=x_th**4+x_th**2*y_th**2-x_th**2*z_th**2 
      
      if(rxyz.lt.0d0) return 1 

      costh=y_th*x_th-dsqrt(rxyz) 
      costh=costh/(x_th**2+y_th**2) 
      
!      if(dabs(theta).gt.1d0) then          
!         costh=y_th*x_th+dsqrt(rxyz) 
!         costh=costh/(x_th**2+y_th**2) 
!         if(dabs(theta).gt.0.999999999d0) return 1     
!      endif

      sinth=dsqrt(1d0-costh**2)
      
     
      
      pz(4)=dsqrt(mv2)/2d0

      
      pz(1)=pz(4)*sinth*dcos(phi)
      pz(2)=pz(4)*sinth*dsin(phi)
      pz(3)=pz(4)*costh

     
!------ define pz in lab frame 

      do i=1,2 
         Q(i)=-pgam(i)
      enddo
      

      a=  -(-(mv2*pgam(3)) + Qsq*pgam(3) + pgam(3)**3 - 
     &        pgam(3)*pgam(4)**2 + 
     &        Sqrt(pgam(4)**2*
     &          (mv2**2 + (Qsq + pgam(3)**2 - pgam(4)**2)**2 - 
     &            2*mv2*(Qsq - pgam(3)**2 + pgam(4)**2))))/
     &     (2.*(pgam(3)**2 - pgam(4)**2))
   
     
      
      b=dsqrt(mv2+a**2+pgam(1)**2+pgam(2)**2)
      Q(3)=a 
      Q(4)=b


      write(6,*) Q(1),Q(2),Q(3),Q(4)
      write(6,*) pgam(1),pgam(2),pgam(3),pgam(4)
      write(6,*) 0d0,0d0,0d0,rQsq

      write(6,*) dsqrt(Q(4)**2-Q(1)**2-Q(2)**2-Q(3)**2)
      write(6,*) dsqrt(mv2)

      call boost(dsqrt(mv2),Q,pz,pl1)

      
   
      do j=1,4
         p(1,j)=pl1(j)
         p(2,j)=Q(j)-pl1(j)
      enddo

       

      wt=wt/8d0/pi

      write(6,*) p_in(1),p_in(2),p_in(3),p_in(4)

      write(6,*) p(1,1),p(1,2),p(1,3),p(1,4)
      write(6,*) p(2,1),p(2,2),p(2,3),p(2,4) 
      write(6,*) p(3,1),p(3,2),p(3,3),p(3,4)
      write(6,*) p(1,1)+p(2,1)+p(3,1),p(1,2)+p(2,2)+p(3,2)
     &,p(1,3)+p(2,3)+p(3,3),p(1,4)+p(2,4)+p(3,4),rQsq
      write(6,*)dsqrt(-2d0*(p(1,1)*p(2,1)+p(1,2)*p(2,2)+p(1,3)*p(2,3))
     &+2d0*p(1,4)*p(2,4))

      pause

      return 
      end 
