
      subroutine docuts(nd,p,passed)
      implicit none 
      include 'constants.f' 
      include 'process.f' 
      include 'frag.f' 
      character*4 part
      common/part/part
      integer nd,isub
      logical cuts
      common/cuts/cuts
      double precision p(mxpart,4)
      logical phot_dip(mxpart)
      logical passed,iso
      common/phot_dip/phot_dip
      if(nd.gt.0) then 
         isub=1 
      else
         isub=0 
      endif

      passed=.true. 
      if(cuts.eqv..false.) return
     
      if(case.eq.'Z_only') then 
         call cuts_drell_yan(p,passed) 
         return 
      elseif(case.eq.'Zgamma') then     
!         if(rescale) call rescale_pjet(p) 
        
         call cuts_Vgamma(nd,p,passed)
!         if(passed.and.(nd.eq.0)) stop
         if(passed.eqv..false.) return
!-------- Check isolation of photon
         if(frag) then 
            passed=iso(p,isub,phot_dip(nd),5,nd)         
         else
          
!           pause
            if((part.eq.'real').and.(nd.eq.0)) then 
     
               call frix(p,passed,5,isub)
            endif
         endif
!        stop
        
        return 
      elseif((case.eq.'ZZlept').or.(case.eq.'HZZ_4l')) then 
         call cuts_ZZ_4l(p,passed)
         
         return 
      endif

      return 
      end 
