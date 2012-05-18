

      subroutine setzdips(z) 
      implicit none 
      include 'constants.f'
      include 'process.f'  
      double precision pborn(mxpart,4),z 
      logical first
      integer j,k 
      common/pborn/pborn
      data first/.true./ 
      save first 

      
      if(case.eq.'W_only') then 
!          call qqb_w_z(pborn,z)
      elseif (case .eq. 'Wgamma') then
!         call qqb_wgam_z(pborn,z)
      elseif (case .eq. 'Z_only') then
         call qqb_z_z(pborn,z)
      elseif (case .eq. 'Zgamma') then
          call qqb_zgam_z(pborn,z)
      elseif (case .eq. 'WWqqbr') then
!         call qqb_ww_z(pborn,z)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz_z(pborn,z)
      elseif (case .eq. 'gamgam') then
!         call qqb_gamgam_z(pborn,z)
      elseif (case .eq. 'WH__WW') then
!         call qqb_wh_ww_z(pborn,z)
      elseif (case .eq. 'WH__ZZ') then
!         call qqb_wh_zz_z(pborn,z)
      elseif (case .eq. 'WHgaga') then
!         call qqb_wh_gaga_z(pborn,z)
      elseif (case .eq. 'ZH__WW') then
!        call qqb_zh_ww_z(pborn,z)
      elseif (case .eq. 'ZH__ZZ') then
!        call qqb_zh_zz_z(pborn,z)
      elseif (case .eq. 'ZHgaga') then
!        call qqb_zh_gaga_z(pborn,z)
      elseif (case .eq. 'Higaga') then
!        call gg_hgamgam_z(pborn,z)
      elseif (case .eq. 'HWW_4l') then
!        call qqb_hww_z(pborn,z)
      elseif (case .eq. 'HWW_tb') then
!        call qqb_hww_tb_z(pborn,z)
      elseif (case .eq. 'HWWint') then
!        call gg_ww_int_z(pborn,z)
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz_z(pborn,z)
  
      else
         if(first) then
         write(6,*) 'Unimplmented process = ',case  
      
         first=.false. 
         endif
      endif
      

      return 
      end subroutine 
