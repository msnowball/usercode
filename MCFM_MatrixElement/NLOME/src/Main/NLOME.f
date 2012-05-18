

      double precision function NLOME(p) 
!------ takes in p returns integral of ME over boosts at this point 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f'
      include 'vegas_common.f' 
      include 'process.f'
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
      logical log_mode_x 
      common/papb/pa,pb
      common/log_mode_x/log_mode_x 
      integer itmx1,ncall1,int_flag,int_fr
      common/xminmax/xmin,xmax
      double precision err_l1,err_v,v_1,lo_1,lord_virt
      double precision real_piece,err_v1,lv,vdips,err_r1
      double precision lv_int,vdips_int 
      external lv_int,vdips_int
      logical passed,cuts
      common/cuts/cuts
      integer itmxb,itmxv,itmxr
      integer ncallb,ncallv,ncallr
      character*4 part
      common/part/part
      common/my_itsNLOME/itmxb,itmxv,itmxr
      common/my_ncsNLOME/ncallb,ncallv,ncallr 
      common/int_setup/int_flag,int_fr
      integer ncallr_reset,n_fail
      common/n_fail/n_fail
      logical calc_real
      common/calc_real/calc_real
      logical first 
      data first /.true./
      save first 
      double precision my_pt

   

      calc_real=.true.
      

      NLOME=0d0 
      lord_virt=0d0 
      lv=0d0 
      vdips=0d0 
      real_piece=0d0 
      lo_1=0d0 
      v_1=0d0 

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
      
!      call writeout(pborn) 
      
!     
!------ CHECK that we will have a contribution from born and virt 
      passed=.true.
      if(cuts) call docuts(0,pborn,passed)
      if(calc_real.eqv..false.) goto 999 
      if(passed.eqv..false.) goto 101
     
      
!------ determine xmin and xmax 
      call determine_xminmax(xmin,xmax,Qsq/sqrts**2,pborn,*101)
!      write(6,*) xmin,xmax
        
!----- integrate over this region 
!----- calculate lord and virt ME 
      call lord_ME
      call virtME
     
!----- OPTION 1) USE ROMBERG INTEGRATION,  
      if(int_flag.eq.1) then 
!--------- virt and lord no z dependence
         call dromb(lv_int,xmin,xmax,lv)      
!---------- open integral for dipoles 
         call dromo(vdips_int,0d0,1d0,vdips,midpnt) 
         lord_virt=lv+vdips
!------OPTION 2) USE VEGAS INTEGRATION 
       elseif(int_flag.eq.2) then 
          ndim=1
          call Lord_vegas(0,itmxb,ncallb,.false.,lo_1,err_l1)
          ndim=2        
          call virt_vegas(0,itmxv,ncallv,.false.,v_1,err_v1)
          lord_virt=lo_1+v_1
       endif
      
 101   continue 
       
!--------- REAL PIECES  
       real_piece=0d0 
     
       part='real'   
       n_fail=0 
     
       call real_vegas(int_fr,0,itmxr,ncallr,.false.,real_piece,err_l1)
    
!-------ensure correct dimension 
       if((case.eq.'HZZ_4l').or.(case.eq.'ZZlept')) then
          real_piece=real_piece/sqrts**6
       endif

 102   continue

      NLOME=lord_virt+real_piece
!      write(6,*) lord_virt
!      write(6,*) real_piece 
!      pause
      return
 999  continue 
      NLOME=0d0 
      return 
      end 
