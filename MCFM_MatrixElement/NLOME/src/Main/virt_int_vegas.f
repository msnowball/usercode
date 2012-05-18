
      double precision function virt_int_vegas(vector,wgt)
      implicit none 
      include 'constants.f' 
      include 'qcdcouple.f'
      include 'facscale.f' 
      include 'scale.f'
      include 'epinv.f' 
      include 'agq.f'
      include 'PR_new.f'
      include 'Gen_Q.f' 
      include 'b0.f'
      include 'n_born.f' 
      include 'vegas_common.f'
      include 'interface_settings.f' 
  
      double precision vector(mxdim),wgt 
      integer j,nu,ih1,ih2,k
      double precision pborn(mxpart,4)
      double precision virt_me(-nf:nf,-nf:nf),xmsq
      double precision xx(2),Qsq,Q_f(4),sqrts
      double precision fx1(-nf:nf),fx2(-nf:nf),W,flux
      double precision fx1z(-nf:nf),fx2z(-nf:nf)
      double precision lo_me(-nf:nf,-nf:nf),xzz,z
      double precision msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      double precision AP(-1:1,-1:1,3),x1onz,x2onz,omz
      double precision xmin,lnxmin,xjac,xmax
      double precision zjac
      double precision pa(4),pb(4),x1,x2
      common/papb/pa,pb
      integer ia,ib,ic,is
      double precision xa,xb
      common/x1x2/x1,x2
      logical cuts,passed_cuts 
      common/cuts/cuts
      common/virt_me/virt_me
      common/lo_me/lo_me 
      common/pborn/pborn 
      common/energy/sqrts
      common/density/ih1,ih2
      common/xminmax/xmin,xmax
      logical logxsam
      common/logxsam/logxsam

      virt_int_vegas=0d0 

      
!----- Then calculate Q (final state partons 
      do nu=1,4
         Q_f(nu)=0d0 
         pa(nu)=pborn(1,nu) 
         pb(nu)=pborn(2,nu)
         do j=3,n_particles            
            Q_f(nu)=Q_f(nu)-pborn(j,nu) 
         enddo 
      enddo
      Qsq=Q_f(4)**2-Q_f(3)**2-Q_f(2)**2-Q_f(1)**2
      
     
!---- Now use that (x1x2*sqrts**2)=Q**2 to set x2 given that x1 = r1
      if(logxsam.eqv..false.) then 

!---- Generate xx(1) linearly 
      xx(1)=xmax*(one-vector(1))+vector(1)*xmin 
      xjac=dabs(xmin-xmax)
      xx(2)=Qsq/(xx(1)*sqrts**2) 
      if(xx(2).gt.1d0) goto 999 
      x1=xx(1)
      x2=xx(2)
      else
    
   
!---- Generate xx(1) log     
      lnxmin=dlog(xmin/xmax)
      xx(1)=xmax*dexp(lnxmin*(one-vector(1)))
      xjac=-lnxmin*xx(1)
      xx(2)=Qsq/(xx(1)*sqrts**2) 
      x1=xx(1)
      x2=xx(2)
     
      if(xx(2).gt.0.99d0) goto 999 

      endif

 
!----- Calculate PDFS
      W=sqrts**2

 
      call docuts(0,pborn,passed_cuts) 
      if(passed_cuts.eqv..false.) goto 999 

      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)*xjac
    
      do ia=-1,+1
      do ib=-1,+1
      do ic=-1,+1
      do is=1,3
        Q1(ia,ib,ic,is)=0d0
        Q2(ia,ib,ic,is)=0d0
      enddo
      enddo
      enddo
      enddo
      
!---- generate z quadratically 
      z=vector(2)**2
      zjac=two*dsqrt(z)
      xmsq=0d0
      omz=one-z

      
      flux=flux*zjac
      call setzdips(z)
  
      epcorr=epinv+2d0*dlog(scale/facscale)
      
      AP(q,q,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(q,q,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(q,q,3)=+ason2pi*Cf*2d0/omz*epcorr
      AP(a,a,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(a,a,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(a,a,3)=+ason2pi*Cf*2d0/omz*epcorr
      
      AP(q,g,1)=0d0
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(q,g,3)=0d0
      AP(a,g,1)=0d0
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(a,g,3)=0d0
      
      AP(g,q,1)=0d0
      AP(g,q,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,q,3)=0d0
      AP(g,a,1)=0d0
      AP(g,a,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,a,3)=0d0

      AP(g,g,1)=+ason2pi*b0*epcorr
      AP(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epcorr
      AP(g,g,3)=+ason2pi*xn*2d0/omz*epcorr

      

      do j=-nf,nf
         fx1(j)=0d0 
         fx2(j)=0d0 
         fx1z(j)=0d0 
         fx2z(j)=0d0
      enddo
      
      if (z .gt. xx(1)) then
         x1onz=xx(1)/z
         call fdist(ih1,x1onz,facscale,fx1z)
      endif
      if (z .gt. xx(2)) then
         x2onz=xx(2)/z         
         call fdist(ih2,x2onz,facscale,fx2z)
      endif         
      
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih1,xx(2),facscale,fx2)
      
      xmsq=0d0 
      xzz=0d0
      
      
      do j=-nf,nf
         do k=-nf,nf
!     write(6,*) j,k,xmsq
C--QQ 
            if     ((j .gt. 0) .and. (k.gt.0)) then
   
               xmsq=xmsq+(virt_me(j,k)
     & + lo_me(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + lo_me(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + lo_me(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
               
C--QbarQbar
      elseif ((j .lt. 0) .and. (k.lt.0)) then
   
      xmsq=xmsq+(virt_me(j,k)
     & + lo_me(j,k)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + lo_me(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + lo_me(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
  

C--QQbar
      elseif ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(virt_me(j,k)
     & + lo_me(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + lo_me(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + lo_me(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z
      elseif ((j .lt. 0) .and. (k.gt.0)) then

C--QbarQ
      xmsq=xmsq+(virt_me(j,k)
     & +lo_me(j,k)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &               *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     & + lo_me(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     & + lo_me(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .eq. g) .and. (k.eq.g)) then
C--gg
 
     
       msq_qg=
     &       lo_me(+5,g)+lo_me(+4,g)+lo_me(+3,g)+lo_me(+2,g)+lo_me(+1,g)
     &      +lo_me(-5,g)+lo_me(-4,g)+lo_me(-3,g)+lo_me(-2,g)+lo_me(-1,g)
       msq_gq=
     &    lo_me(g,+5)+lo_me(g,+4)+lo_me(g,+3)+lo_me(g,+2)+lo_me(g,+1)
     &     +lo_me(g,-5)+lo_me(g,-4)+lo_me(g,-3)+lo_me(g,-2)+lo_me(g,-1)
       xmsq=xmsq+(virt_me(g,g)
     &  +lo_me(g,g)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     &                +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &                *fx1(g)*fx2(g)
     &  +(lo_me(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     &  +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     &  +(lo_me(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     &  +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z
      

      elseif (j .eq. g) then
C--gQ
       if    (k .gt. 0) then
       msq_aq=
     &      lo_me(-1,k)+lo_me(-2,k)+lo_me(-3,k)+lo_me(-4,k)+lo_me(-5,k)
       msq_qq=
     &      lo_me(+1,k)+lo_me(+2,k)+lo_me(+3,k)+lo_me(+4,k)+lo_me(+5,k)
       xmsq=xmsq+(virt_me(g,k)
     & +lo_me(g,k)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1(g)*fx2(k)
     & +(lo_me(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z(g)/z*fx2(k)
     & +(lo_me(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + lo_me(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z
 
C--gQbar
       elseif (k.lt.0) then
       msq_qa=
     &       lo_me(+1,k)+lo_me(+2,k)+lo_me(+3,k)+lo_me(+4,k)+lo_me(+5,k)
       msq_aa=
     &      lo_me(-1,k)+lo_me(-2,k)+lo_me(-3,k)+lo_me(-4,k)+lo_me(-5,k)
       xmsq=xmsq+(virt_me(g,k)
     & +lo_me(g,k)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1(g)*fx2(k)
     & +(lo_me(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z(g)/z*fx2(k)
     & +(lo_me(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + lo_me(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z
    
       endif
C--Qg
      elseif (k .eq. g) then
       if     (j.gt.0) then
       msq_qa=
     &      lo_me(j,-1)+lo_me(j,-2)+lo_me(j,-3)+lo_me(j,-4)+lo_me(j,-5)
       msq_qq=
     &      lo_me(j,+1)+lo_me(j,+2)+lo_me(j,+3)+lo_me(j,+4)+lo_me(j,+5)
       xmsq=xmsq+(virt_me(j,g)
     & +lo_me(j,g)*(
     &               +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1(j)*fx2(g)
     & +(lo_me(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + lo_me(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     & +(lo_me(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1(j)*fx2z(g)/z
  
C--Qbarg
       elseif (j.lt.0) then
       msq_aq=
     &       lo_me(j,+1)+lo_me(j,+2)+lo_me(j,+3)+lo_me(j,+4)+lo_me(j,+5)
       msq_aa=
     &      lo_me(j,-1)+lo_me(j,-2)+lo_me(j,-3)+lo_me(j,-4)+lo_me(j,-5)
       xmsq=xmsq+(virt_me(j,g)
     & +lo_me(j,g)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                *fx1(j)*fx2(g)
     & +(lo_me(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + lo_me(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     & +(lo_me(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     & + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1(j)*fx2z(g)/z
   
      endif
      endif

      
      enddo
      enddo

      
      virt_int_vegas=xmsq*flux
 999  continue 
      return 
      end 
