

      double precision function vdips_int(z)
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
      include 'interface_settings.f' 
  
     
      integer j,nu,ih1,ih2,k
      double precision xmin,xmax,my_simp_dips
      double precision pborn(mxpart,4)
     
      double precision xx(2),Qsq,Q_f(4),sqrts    
      double precision lo_me(-nf:nf,-nf:nf),xzz,z
      double precision msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      double precision AP(-1:1,-1:1,3),x1onz,x2onz,omz
      double precision zjac
      double precision pa(4),pb(4),x1,x2
      common/papb/pa,pb
      integer ia,ib,ic,is
      double precision xa,xb,vdips_xint,zdip
      double precision vout
      external vdips_xint
      common/x1x2/x1,x2
      logical cuts,passed_cuts 
      common/cuts/cuts
    
      common/lo_me/lo_me 
      common/pborn/pborn 
      common/energy/sqrts
      common/density/ih1,ih2
      common/myAP/AP
      common/xminmax/xmin,xmax
      common/intdip_z/zdip

      vdips_int=0d0 

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
      
      
      omz=one-z
     
      call setzdips(z)
      zdip=z
  
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

   
      call dromb(vdips_xint,xmin,xmax,vout)
      vdips_int=vout
     
      return 
 999  continue 
      return 
      end 

      
