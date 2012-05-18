

      subroutine lord_ME
      implicit none 
      include 'constants.f' 
      include 'process.f' 
!-----fills lo_me in common block 
      double precision pborn(mxpart,4),lo_me(-nf:nf,-nf:nf) 
      logical first
      integer j,k 
      common/lo_me/lo_me 
      common/pborn/pborn
      data first/.true./ 
      save first 
      
      if(case.eq.'W_only') then 
          call qqb_w(pborn,lo_me)
      elseif (case .eq. 'Wgamma') then
!         call qqb_wgam(pborn,lo_me)
      elseif (case .eq. 'Z_only') then
         call qqb_z(pborn,lo_me)
      elseif (case .eq. 'Zgamma') then
         call qqb_zgam(pborn,lo_me)
      elseif (case .eq. 'WWqqbr') then
!         call qqb_ww(pborn,lo_me)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz(pborn,lo_me)
!        pause
!        call writeout(pborn) 
      elseif (case .eq. 'gamgam') then
!         call qqb_gamgam(pborn,lo_me)
      elseif (case .eq. 'WH__WW') then
!         call qqb_wh_ww(pborn,lo_me)
      elseif (case .eq. 'WH__ZZ') then
!         call qqb_wh_zz(pborn,lo_me)
      elseif (case .eq. 'WHgaga') then
!         call qqb_wh_gaga(pborn,lo_me)
      elseif (case .eq. 'ZH__WW') then
!        call qqb_zh_ww(pborn,lo_me)
      elseif (case .eq. 'ZH__ZZ') then
!        call qqb_zh_zz(pborn,lo_me)
      elseif (case .eq. 'ZHgaga') then
!        call qqb_zh_gaga(pborn,lo_me)
      elseif (case .eq. 'Higaga') then
!        call gg_hgamgam(pborn,lo_me)
      elseif (case .eq. 'HWW_4l') then
!        call qqb_hww(pborn,lo_me)
      elseif (case .eq. 'HWW_tb') then
!        call qqb_hww_tb(pborn,lo_me)
      elseif (case .eq. 'HWWint') then
!        call gg_ww_int(pborn,lo_me)
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz(pborn,lo_me)
      
      else
         if(first) then
         write(6,*) 'Unimplmented process = ',case  
         write(6,*) 'Setting LO Matrix Element to 0' 
         do j=-nf,nf
            do k=-nf,nf
               lo_me(j,k)=0d0 
            enddo
         enddo
         first=.false. 
         endif
      endif
      

      return 
      end subroutine 
