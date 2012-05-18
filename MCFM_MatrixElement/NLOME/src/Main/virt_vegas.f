      subroutine virt_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
!----- NLO ME  
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
      integer myitmx,myncall,myinit,i,j,k
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsd,sumsdr
     
      common/bin/bin
      double precision region(2*mxdim)
      logical first,myreadin
      external virt_int_vegas
     
      integer ndim_in 

      ndim_in=ndim 
      xerr=0d0 
      xinteg=0d0 

      
      itmx=myitmx
      bin=mybin 
      
!      for PDFS + integrated dipoles, increase ndim=1 by 1 
      ndim=2
      call boundregion(ndim,region)
      call vegasnr(region,ndim,Virt_int_vegas,myinit,myncall,myitmx,
     .     0,xinteg,xerr,chi)
      
      
      ndim=ndim_in
      return 
      end 
          
