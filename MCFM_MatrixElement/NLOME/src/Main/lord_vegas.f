      subroutine lord_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)!
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
      integer myitmx,myncall,myinit,i,j,k,nproc,mynproc
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsd,sumsdr
      double precision pborn(mxpart,4)
      character*4 part,mypart
       
    
      double precision region(2*mxdim)
      logical first,myreadin
!      logical theory_mode,total_cr
      double precision theory_mode_lord 
!      common/techP/theory_mode,total_cr
      common/pborn/pborn 
      common/nproc/nproc
      common/bin/bin
      common/xreal/xreal,xreal2

      external LO_int_vegas
      data first/.true./
      save first
      integer ndim_in 

      ndim_in=ndim 

      xerr=0d0 
      xinteg=0d0 
      
      itmx=myitmx
      bin=mybin 
     
!-----one dim int for PDFS
      ndim=1
      call boundregion(ndim,region)
      call vegasnr(region,ndim,LO_int_vegas,myinit,myncall,myitmx,
     .     0,xinteg,xerr,chi)
      
      
      ndim=ndim_in 
      return 
      end 
          
    
