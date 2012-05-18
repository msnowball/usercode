      subroutine scaleset(rscalestart,fscalestart,p)
c--- wrapper routine to set a dynamic scale; please refer to individual
c--- routines for exact definitions of the scales.
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'frag.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
!      include 'inp_scales.f'
      double precision rscalestart,fscalestart,p(mxpart,4),mu0,
     & alphas,amz
      logical first
      common/couple/amz
      data first/.true./  
      save first
      
!      inp_scale=rscalestart
!      inp_fscale=fscalestart

      return
c$$$
c$$$      if     (dynstring .eq. 'm(34)') then
c$$$        call scaleset_m34(p,mu0)
c$$$      elseif (dynstring .eq. 'm(345)') then
c$$$        call scaleset_m345(p,mu0)
c$$$      elseif (dynstring .eq. 'm(3456)') then
c$$$        call scaleset_m3456(p,mu0)
c$$$      elseif (dynstring .eq. 'sqrt(M^2+pt34^2)') then
c$$$        call scaleset_Msqpt34sq(p,mu0)
c$$$      elseif (dynstring .eq. 'sqrt(M^2+pt5^2)') then
c$$$        call scaleset_Msqpt5sq(p,mu0)
c$$$      elseif (dynstring .eq. 'pt(photon)') then
c$$$        call scaleset_ptphoton(p,mu0)
c$$$      elseif (dynstring .eq. 'HT') then
c$$$        call scaleset_HT(p,mu0)
c$$$      elseif (dynstring .eq. 'DDIS') then
c$$$        call scaleset_ddis(p,mu0)
c$$$      else
c$$$        write(6,*) 'Dynamic scale choice not recognized'
c$$$	write(6,*) '   dynamicscale = ',dynstring
c$$$	stop
c$$$      endif
c$$$      
c$$$      scale=rscalestart*mu0
c$$$      facscale=fscalestart*mu0
c$$$c--- piggy-back renomalization scale for fragmentation scale 
c$$$      frag_scale=scale
c$$$	    
c$$$      if (first) then
c$$$        write(6,*)
c$$$        write(6,*)'************** Dynamic scale choice ****************'
c$$$        write(6,*)'*                                                  *'
c$$$        write(6,*)'*                 RENORMALIZATION                  *'
c$$$        write(6,45) ' mu_ren  =',rscalestart,dynstring
c$$$        write(6,*)'*                                                  *'
c$$$        write(6,*)'*                  FACTORIZATION                   *'
c$$$        write(6,45) ' mu_fac  =',fscalestart,dynstring
c$$$	if (frag) then
c$$$        write(6,*)'*                                                  *'
c$$$        write(6,*)'*                  FRAGMENTATION                   *'
c$$$        write(6,45) ' mu_frag =',rscalestart,dynstring
c$$$	endif
c$$$        write(6,*)'*                                                  *'
c$$$        write(6,*)'****************************************************'
c$$$        first=.false.      
c$$$      endif
c$$$  
c$$$c--- catch absurdly large and small scales      
c$$$      if  (scale .gt. 3000d0) scale=3000d0
c$$$      if  (facscale .gt. 3000d0) facscale=3000d0
c$$$      if  (frag_scale .gt. 3000d0) frag_scale=3000d0
c$$$      if  (scale .lt. 1d0) scale=1d0
c$$$      if  (facscale .lt. 1d0) facscale=1d0
c$$$      if  (frag_scale .lt. 1d0) frag_scale=1d0
c$$$
c$$$c--- run alpha_s
c$$$      as=alphas(scale,amz,nlooprun)
c$$$	
c$$$      ason2pi=as/twopi
c$$$      ason4pi=as/fourpi
c$$$      gsq=fourpi*as
c$$$      musq=scale**2
c$$$      
c$$$      return
c$$$
c$$$ 45   format(1x,'* ',a15,f6.2,' x ',a24,' *')

      end
      
