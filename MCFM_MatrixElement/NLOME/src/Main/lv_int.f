      double precision function lv_int(x)
      implicit none 
      include 'constants.f' 
      include 'facscale.f'
      double precision x,msq(-nf:nf,-nf:nf),fx1(-nf:nf),fx2(-nf:nf) 
      double precision msqv(-nf:nf,-nf:nf)
      double precision Qsq,sqrts,xjac,xmsq 
      logical log_mode_x
      common/log_mode_x/log_mode_x
      common/Qsq/Qsq
      common/energy/sqrts 
      common/lo_me/msq
      common/virt_me/msqv
      integer i,j,k,ih1,ih2
      common/density/ih1,ih2
      double precision xx(2),flux,t

      lv_int=0d0      
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
            xmsq=fx1(j)*fx2(k)*(msq(j,k)+msqv(j,k))
            lv_int=lv_int+xmsq*flux
         
         enddo
      enddo
      
      

      return 
      end
