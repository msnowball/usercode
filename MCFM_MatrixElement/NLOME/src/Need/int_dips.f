

      double precision function int_dips(p,z,flag) 
      implicit none 
      include 'constants.f' 
      include 'qcdcouple.f'
      include 'facscale.f' 
      include 'scale.f'
      include 'epinv.f' 
      include 'agq.f'
      include 'PR_new.f'
      double precision p(mxpart,4),z,msq(-nf:nf,-nf:nf)
      double precision msqv(-nf:nf,-nf:nf)
      double precision msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      double precision fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf)
      double precision AP(-1:1,-1:1,3),sqrts,flux,W,xmsq
      double precision xx(2),x1onz,x2onz,omz

      integer flag !--- chosses zz,hzz
      integer j,k,ih1,ih2
      common/energy/sqrts
      common/density/ih1,ih2

      int_dips=0d0 
      
      do j=-nf,nf
         do k=-nf,nf 
            msq(j,k)=0d0 
            msqv(j,k)=0d0 
         enddo
      enddo
      
      xmsq=0d0
      omz=one-z

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


!------ Fill appropriate common blocks based on flag 
      if(flag.eq.1) then 
         call qqb_zz_z(p,z) 
         call qqb_zz(p,msq) 
      elseif(flag.eq.2) then 
         call qqb_hzz_z(p,z) 
         call qqb_hzz(p,msq) 
      else
         write(6,*) 'Unknown option in int_dips = ',flag
         stop
      endif

      do j=-nf,nf
         fx1z(j)=0d0
         fx1(j)=0d0
         fx2z(j)=0d0 
         fx2(j)=0d0
      enddo
           
       xx(1)=-2d0*p(1,4)/sqrts
       xx(2)=-2d0*p(2,4)/sqrts
 
       
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
    
      do j=-nf,nf
         do k=-nf,nf
!           write(6,*) j,k,xmsq
C--QQ 
            if     ((j .gt. 0) .and. (k.gt.0)) then
   
               xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
               
C--QbarQbar
      elseif ((j .lt. 0) .and. (k.lt.0)) then
   
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
  

C--QQbar
      elseif ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z
      elseif ((j .lt. 0) .and. (k.gt.0)) then

C--QbarQ
      xmsq=xmsq+(msqv(j,k)
     & +msq(j,k)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &               *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .eq. g) .and. (k.eq.g)) then
C--gg
 
     
       msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     &       +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
       msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     &       +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)
       xmsq=xmsq+(msqv(g,g)
     &  +msq(g,g)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     &                +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &                *fx1(g)*fx2(g)
     &  +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     &  +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     &  +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     &  +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z
      

      elseif (j .eq. g) then
C--gQ
       if    (k .gt. 0) then
       msq_aq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       msq_qq=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z
 
C--gQbar
       elseif (k.lt.0) then
       msq_qa=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       msq_aa=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z
    
       endif
C--Qg
      elseif (k .eq. g) then
       if     (j.gt.0) then
       msq_qa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       msq_qq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(
     &               +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1(j)*fx2z(g)/z
  
C--Qbarg
       elseif (j.lt.0) then
       msq_aq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       msq_aa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     & + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1(j)*fx2z(g)/z
   
      endif
      endif

      enddo
      enddo
      
      int_dips=xmsq

      return 
      end 
