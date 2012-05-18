

      subroutine cuts_Vgamma(nd,p,passed) 
      implicit none 
      include 'interface_settings.f' 
      include 'constants.f' 
      include 'frag.f'
      include 'user_cuts.f' 
      double precision p(mxpart,4),my_pt,my_etarap
      double precision pt_gam_min,mll_min
      double precision Qsq,my_R,mll_max
      double precision pt_l_s,pt_l_h
      logical passed 
      integer i,nd
      logical first 
      data first/.true. /
      save first 
      double precision x1,x2 
      common/x1x2/x1,x2 
      double precision dot
      
    
     
      passed=.true. 

      
      if(first) then 
         write(6,*) '**************** V GAMMA CUTS **************' 
         write(6,*) '*    pt (gamma) > ',gammpt, '      *' 
         write(6,*) '*    R(l,gamma) > ',R_lgam, '      *' 
         write(6,*) '*    pt_l (hard ) > ',lept_pt_1,'       *' 
         write(6,*) '*    pt_l (soft) >  ',lept_pt_2,'        *'
!         write(6,*) '*   ',mll_s, ' < m_ll  <  ',mll_h,' *' 
!         write(6,*) '*   ',mll_s, ' < m_llgam  <  ',mll_h,' *' 
!         write(6,*) '*
         write(6,*) '**************************************************'
         first=.false. 
      endif

     
 !     Qsq=0d0 
 !     do i=1,3 
 !        Qsq=Qsq-(p(3,i)+p(4,i))**2
 !     enddo
 !     Qsq=Qsq+(p(3,4)+p(4,4))**2
 !     Qsq=dsqrt(Qsq)
      
!-   Make sure both are hard enough to pass hard 
      if((my_pt(3,p).lt.lept_pt_1).and.(my_pt(4,p).lt.lept_pt_1)) then 
         passed = .false. 
         return 
      endif

!--   Make sure at least one is hard enough to pass soft 
      if((my_pt(3,p).lt.lept_pt_2).or.(my_pt(4,p).lt.lept_pt_2)) then 
         passed = .false. 
         return 
      endif

      if((my_R(3,5,p).lt.R_lgam).or.(my_R(4,5,p).lt.R_lgam)) then 
         passed=.false.
         return 
      endif

!      if((Qsq.lt.mll_s).or.(Qsq.gt.mll_h)) then 
!         passed=.false.
        
!         return 
!      endif
      

      if(my_pt(5,p).lt.gammpt) then 
         passed=.false.       
         return 
      endif
      
      return 
      end 


      double precision function dphi(i,j,p) 
      implicit none 
      include 'interface_settings.f'
      include 'constants.f' 
      integer i,j
      double precision p(mxpart,4) 
      double precision my_pt,my_etarap
      double precision yi,yj,pti,ptj,sij 
      integer nu

      sij=0d0 
      do nu=1,3 
         sij=sij-(p(i,nu)+p(j,nu))**2 
      enddo
      sij=sij+(p(i,4)+p(j,4))**2 

      yi=my_etarap(i,p) 
      yj=my_etarap(j,p) 
      yi=0d0 
      yj=0d0
      pti=my_pt(i,p) 
      ptj=my_pt(j,p) 

      dphi=dcosh(yi-yj)-sij/(2d0*pti*ptj)     
      if(dphi.lt.-0.99999999999d0) then 
         dphi=-1d0 
      elseif(dphi.gt.0.999999999999d0) then 
         dphi=1d0 
      endif
 
      dphi=dacos(dphi) 

      return 
      end 


      
      double precision function my_R(i,j,p) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      integer i,j 
      double precision p(mxpart,4) 
      double precision my_etarap,dphi 
      double precision dely

      my_R=0d0
  !    call writeout(p) 
      dely=my_etarap(i,p)-my_etarap(j,p) 
      dely=0d0 
!     write(6,*) 'dely',dely 
      my_R=dely**2+dphi(i,j,p)**2
 !     write(6,*) 'dphi ',dphi(i,j,p) 
 !     pause
      my_R=dsqrt(my_R) 
      
      return 
      end

      
