

      double precision function LOME(p) 
!------ takes in p returns integral of ME over boosts at this point 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f'
      include 'vegas_common.f' 
      include 'scale.f' 
      include 'facscale.f' 
      include 'masses.f' 
      double precision p(n_final,4)
      double precision pborn(mxpart,4),ptemp(n_final,4) 
      common/pborn/pborn  
      double precision lo_int,midpnt 
      external lo_int,midpnt
      double precision xmin,xmax
      integer nu,j 
      double precision lo_me(-nf:nf,-nf:nf) 
      common/lo_me/lo_me 
      double precision Q_f(4),Qsq,sqrts,pa(4),pb(4)
      common/Qsq/Qsq
      common/energy/sqrts 
      logical passed
      common/papb/pa,pb
      logical cuts 
      common/cuts/cuts
      integer itmx1,ncall1,int_flag 
      common/xminmax/xmin,xmax
      double precision err_l1,my_pt,x1
      logical first
      data first /.true./ 
      save first
      
      int_flag=2

      LOME=0d0 
      
      do j=1,mxpart
         do nu=1,4 
            pborn(j,nu)=0d0 
         enddo
      enddo

      Qsq=0d0
!------ make a full phase space point
      call boost_Q(p,ptemp) 
      call recastE(ptemp) 
      call fixx1x2(ptemp,pborn)
      do nu=1,4
         Q_f(nu)=0d0 
         pa(nu)=pborn(1,nu) 
         pb(nu)=pborn(2,nu)
         do j=3,n_final+2            
            Q_f(nu)=Q_f(nu)-pborn(j,nu) 
         enddo 
      enddo
      Qsq=Q_f(4)**2-Q_f(3)**2-Q_f(2)**2-Q_f(1)**2
      if(cuts) call docuts(0,pborn,passed) 
     
      if(passed.eqv..false.) goto 999 

        
!     ======== setup dynamic scale, print warning 
!      if(first) then 
!         write(6,*) '****************************************'
!     write(6,*) '* Using dynamic scale : pt_ph**2+mz**2 *' 
!         write(6,*) '****************************************'
!         first=.false.
!         scale=my_pt(3,pborn)**2+zmass**2
!         scale=dsqrt(scale) 
!         facscale=scale
!      endif

!------ determine xmin and xmax 
      call determine_xminmax(xmin,xmax,Qsq/sqrts**2,pborn,*999)

!----- integrate over this region 
!----- calculate ME 
      call lord_ME
   
!----- OPTION 1) USE ROMBERG INTEGRATION,  
      if(int_flag.eq.1) then 
         call dromb(lo_int,xmin,xmax,LOME) 
!------OPTION 2) USE VEGAS INTEGRATION 
       elseif(int_flag.eq.2) then 
          ndim=1
          itmx1=1
          ncall1=400 
          call Lord_vegas(0,itmx1,ncall1,.false.,LOME,err_l1)
       endif



      return
 999  continue 
      LOME=0d0 
      return 
      end 

      double precision function lo_int(x)
      implicit none 
      include 'constants.f' 
      include 'facscale.f'
      double precision x,msq(-nf:nf,-nf:nf),fx1(-nf:nf),fx2(-nf:nf) 
      double precision Qsq,sqrts,xjac,xmsq 
      logical log_mode_x
      common/log_mode_x/log_mode_x
      common/Qsq/Qsq
      common/energy/sqrts 
      common/lo_me/msq
      integer i,j,k,ih1,ih2
      common/density/ih1,ih2
      double precision xx(2),flux,t


     
      lo_int=0d0 

      
      xx(1)=x
      
     
      do i=-nf,nf 
         fx1(i)=0d0 
         fx2(i)=0d0 
      enddo 

   
      xx(2)=Qsq/(xx(1)*sqrts**2) 

  
      if(xx(2).gt.1d0) return 
      if(xx(1).gt.1d0) return 
      
      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)
      
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
      
      xmsq=0d0 

      do j=-nf,nf
         do k=-nf,nf
            xmsq=fx1(j)*fx2(k)*msq(j,k)
            lo_int=lo_int+xmsq*flux
         enddo
      enddo
      
      

      return 
      end
