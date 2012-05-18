
      subroutine RealME_newdips(p,msq_r,msq_c,includereal) 
      implicit none 
      include 'constants.f' 
      include 'ptilde.f' 
      include 'process.f' 
!-----fills lo_me in common block 
      double precision p(mxpart,4),msq_r(-nf:nf,-nf:nf) 
      double precision msq_c(-nf:nf,-nf:nf)
      logical first,includereal
      integer j,k 
      data first/.true./
      save first 
      

      do j=-nf,nf
         do k=-nf,nf
            msq_r(j,k)=0d0 
         enddo
      enddo
      
      if(case.eq.'W_only') then 
!     if(includereal) call qqb_w_g(p,msq_r)
!     call qqb_w(p,msq_c)
      elseif (case .eq. 'Wgamma') then
!     if(includereal) call qqb_wgam_g(p,msq_r)
!         call qqb_wgam(p,msq_c)
      elseif (case .eq. 'Z_only') then
         if(includereal)   call qqb_z1jet(p,msq_r)
         call qqb_z(p,msq_c)
      elseif (case .eq. 'Zgamma') then
         if(includereal)   call qqb_zgam_g(p,msq_r)
         call qqb_zgam(p,msq_c)
      elseif (case .eq. 'WWqqbr') then
!        if(includereal)   call qqb_ww_g(p,msq_r)
!         call qqb_ww(p,msq_c)
      elseif (case .eq. 'ZZlept') then
         if(includereal) call qqb_zz_g(p,msq_r)
        call qqb_zz(p,msq_c)
      
      elseif (case .eq. 'gamgam') then
!     if(includereal)    call qqb_gamgam_g(p,msq_r)
!         call qqb_gamgam(p,msq_c)
      elseif (case .eq. 'WH__WW') then
!       if(includereal)    call qqb_wh_ww_g(p,msq_r)
!         call qqb_wh_ww(p,msq_c)
      elseif (case .eq. 'WH__ZZ') then
!        if(includereal)   call qqb_wh_zz_g(p,msq_r)
!         call qqb_wh_zz(p,msq_c)
      elseif (case .eq. 'WHgaga') then
!         if(includereal)  call qqb_wh_gaga_g(p,msq_r)
!         call qqb_wh_gaga(p,msq_c)
      elseif (case .eq. 'ZH__WW') then
!         if(includereal) call qqb_zh_ww_g(p,msq_r)
!        call qqb_zh_ww(p,msq_c)
      elseif (case .eq. 'ZH__ZZ') then
!         if(includereal)  call qqb_zh_zz_g(p,msq_r)
!        call qqb_zh_zz(p,msq_c)
      elseif (case .eq. 'ZHgaga') then
!         if(includereal) call qqb_zh_gaga_g(p,msq_r)
!        call qqb_zh_gaga(p,msq_c)
      elseif (case .eq. 'Higaga') then
!         if(includereal) call gg_hgamgam_g(p,msq_r)
!        call gg_hgamgam(p,msq_c)
      elseif (case .eq. 'HWW_4l') then
!         if(includereal) call qqb_hww_g(p,msq_r)
!        call qqb_hww(p,msq_c)
      elseif (case .eq. 'HWW_tb') then
!         if(includereal) call qqb_hww_tb_g(p,msq_r)
!        call qqb_hww_tb(p,msq_c)
      elseif (case .eq. 'HWWint') then
!         if(includereal) call gg_ww_int_g(p,msq_r)
!        call qqb_hww_tb(p,msq_c)
      elseif (case .eq. 'HZZ_4l') then
         if(includereal) call qqb_hzz_g(p,msq_r)
         call qqb_hzz(p,msq_c)
      else
         if(first) then
         write(6,*) 'Unimplmented process = ',case  
         write(6,*) 'Setting REAL Matrix Element to 0' 
         do j=-nf,nf
            do k=-nf,nf
               msq_r(j,k)=0d0 
            enddo
         enddo
         first=.false. 
         endif
      endif
      

      return 
      end subroutine 

