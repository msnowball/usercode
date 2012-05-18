

      subroutine VirtME 
      implicit none 
      include 'constants.f' 
      include 'process.f' 
!-----fills lo_me in common block 
      double precision pborn(mxpart,4),virt_me(-nf:nf,-nf:nf) 
      logical first
      integer j,k 
      common/virt_me/virt_me
      common/pborn/pborn
      data first/.true./ 
      save first 
      
      if(case.eq.'W_only') then 
!          call qqb_w_v(pborn,virt_me)
      elseif (case .eq. 'Wgamma') then
!         call qqb_wgam_v(pborn,virt_me)
      elseif (case .eq. 'Z_only') then
         call qqb_z_v(pborn,virt_me)
      elseif (case .eq. 'Zgamma') then
         call qqb_zgam_v(pborn,virt_me)
      elseif (case .eq. 'WWqqbr') then
!         call qqb_ww_v(pborn,virt_me)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz_v(pborn,virt_me)
        call gg_ZZ(pborn,virt_me(0,0))
      elseif (case .eq. 'gamgam') then
!         call qqb_gamgam_v(pborn,virt_me)
      elseif (case .eq. 'WH__WW') then
!         call qqb_wh_ww_v(pborn,virt_me)
      elseif (case .eq. 'WH__ZZ') then
!         call qqb_wh_zz_v(pborn,virt_me)
      elseif (case .eq. 'WHgaga') then
!         call qqb_wh_gaga_v(pborn,virt_me)
      elseif (case .eq. 'ZH__WW') then
!        call qqb_zh_ww_v(pborn,virt_me)
      elseif (case .eq. 'ZH__ZZ') then
!        call qqb_zh_zz_v(pborn,virt_me)
      elseif (case .eq. 'ZHgaga') then
!        call qqb_zh_gaga_v(pborn,virt_me)
      elseif (case .eq. 'Higaga') then
!        call gg_hgamgam_v(pborn,virt_me)
      elseif (case .eq. 'HWW_4l') then
!        call qqb_hww_v(pborn,virt_me)
      elseif (case .eq. 'HWW_tb') then
!        call qqb_hww_tb_v(pborn,virt_me)
      elseif (case .eq. 'HWWint') then
!        call gg_ww_int_v(pborn,virt_me)
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz_v(pborn,virt_me)
 
      else
         if(first) then
         write(6,*) 'Unimplmented process = ',case  
         write(6,*) 'Setting VIRT Matrix Element to 0' 
         do j=-nf,nf
            do k=-nf,nf
               virt_me(j,k)=0d0 
            enddo
         enddo
         first=.false. 
         endif
      endif
      

      return 
      end subroutine 
