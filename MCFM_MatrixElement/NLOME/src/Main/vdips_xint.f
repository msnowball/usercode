      double precision function vdips_xint(x) 
      implicit none 
      include 'constants.f'
      include 'facscale.f' 
      include 'mxdim.f'
      include 'ptilde.f' 
      include 'epinv.f' 
      include 'agq.f'
      include 'PR_new.f'
      double precision p(mxpart,4),Qsq,z
      integer nd,ih1,ih2
      integer i,j,k
      double precision sqrts,xx(2),fx1(-nf:nf),fx2(-nf:nf),flux
      double precision fx1z(-nf:nf),fx2z(-nf:nf) 
      double precision AP(-1:1,-1:1,3)
      common/myAP/AP
      common/energy/sqrts
      common/density/ih1,ih2  
      double precision xmsq,x,x1onz,x2onz
      double precision lo_me(-nf:nf,-nf:nf)
      common/lo_me/lo_me
      double precision msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      common/Qsq/Qsq
      common/intdip_z/z 

      vdips_xint=0d0

     

      xx(1)=x 
      xx(2)=Qsq/(xx(1)*sqrts**2) 
      if(xx(2).gt.0.99d0) return 
      if(xx(1).gt.0.99d0) return 

      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)

   
      do j=-nf,nf
         fx1(j)=0d0 
         fx2(j)=0d0 
         fx1z(j)=0d0 
         fx2z(j)=0d0
      enddo
      
      if ((z .gt. xx(1))) then
         x1onz=xx(1)/z
         call fdist(ih1,x1onz,facscale,fx1z)
      endif
      if ((z .gt. xx(2))) then
         x2onz=xx(2)/z         
         call fdist(ih2,x2onz,facscale,fx2z)
      endif         
      
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih1,xx(2),facscale,fx2)
      
      xmsq=0d0 
   
      
      do j=-nf,nf
         do k=-nf,nf
!     write(6,*) j,k,xmsq
C--QQ 
            if     ((j .gt. 0) .and. (k.gt.0)) then
   
               xmsq=xmsq+(
     & + lo_me(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + lo_me(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + lo_me(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
               
C--QbarQbar
      elseif ((j .lt. 0) .and. (k.lt.0)) then
   
      xmsq=xmsq+(
     & + lo_me(j,k)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + lo_me(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + lo_me(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
  

C--QQbar
      elseif ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(
     & + lo_me(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(lo_me(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + lo_me(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(lo_me(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + lo_me(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z
      elseif ((j .lt. 0) .and. (k.gt.0)) then

C--QbarQ
      xmsq=xmsq+(
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
       xmsq=xmsq+(
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
       xmsq=xmsq+(
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
       xmsq=xmsq+(
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
       xmsq=xmsq+(
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
       xmsq=xmsq+(
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

      vdips_xint=xmsq*flux
     
      return 
      end
