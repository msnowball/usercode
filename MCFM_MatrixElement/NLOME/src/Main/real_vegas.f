      subroutine real_vegas(int,myinit,myitmx,myncall,mybin,xinteg,xerr)!
!----- NLOME Branch 
************************************************************************
*                                                                      *
*  This routine should perform the sweeps of vegasnr                   *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                     *
*       myitmx  :  the number of vegasnr sweeps                        *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error                                   *
*                                                                      *
************************************************************************

      implicit none 
      include 'gridinfo.f' 
      include 'realwt.f' 
      include 'vegas_common.f' 
      include 'constants.f' 
      integer int 
      integer myitmx,myncall,myinit,i,j,k,nproc,mynproc
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsd,sumsdr
      double precision pborn(mxpart,4)
      character*4 part,mypart
       
    
      double precision region(2*mxdim)
      logical first,myreadin
      integer it_c 
!      logical theory_mode,total_cr
      double precision theory_mode_lord 
!      common/techP/theory_mode,total_cr
      common/pborn/pborn 
      common/nproc/nproc
      common/bin/bin
      common/xreal/xreal,xreal2
      common/it_c/it_c
      external real_int
      data first/.true./
      save first
      integer ndim_in 
      external real_int_2dv,real_int_newdips
     
      ndim_in=ndim 

      xerr=0d0 
      xinteg=0d0 
      
      itmx=myitmx
      bin=mybin 
     
      if(int.eq.1) then 
!----- ndim = 2 remaining int is done romb 
         ndim=2
         call boundregion(ndim,region)
         call vegasnr(region,ndim,real_int_2dv,myinit,myncall,myitmx,
     &        0,xinteg,xerr,chi)
      else
!-----all vegas
         ndim=3
         call boundregion(ndim,region)
         call vegasnr(region,ndim,real_int,myinit,myncall,myitmx
     &        ,0,xinteg,xerr,chi)
      endif

      
      ndim=ndim_in 
      return 
      end 
          
    
