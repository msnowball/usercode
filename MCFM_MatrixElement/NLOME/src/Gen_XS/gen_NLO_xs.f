!----- C. Williams Jan 2012

      subroutine gen_NLO_XS(mu_rin,mu_fin,LO_xs,n_w) 
      implicit none 
!----- subroutine will generate LO cross section 
!----- and write weighted events to file 
!----- takes in scales and number of points for writoute iteration as inputs 
!----- and returns LO cross section as an output  
      include 'vegas_common_cr.f' 
      include 'outputstring.f'
      include 'interface_settings.f'
      include 'scale.f' 
      include 'facscale.f'
      double precision mu_rin,mu_fin,LO_xs
      double precision lo_1,err_l1
      integer n_i,n_w,itmx,ncall
      logical write_out
      common/write_out/write_out
    
!--------- INITIALIZE CODE 

      call init_NLOME(mu_rin,mu_fin) 
      scale=mu_rin
      facscale=mu_fin
      itmx=10
      ncall=5000
!-------- VEGAS 
      call gen_NLO_vegas(0,itmx,ncall,.false.,lo_1,err_l1)

!---------WRITE TO FILE 
      if(write_out) then 
         call gen_NLO_vegas(1,10,n_w,.true.,LO_xs,err_l1)
      
         write(6,*) '*****************************************' 
         write(6,*) '* Used ,',n_w,' events                  *' 
         write(6,*) '* Calculated cross section, ',LO_xs,'   *' 
         write(6,*) '* Wrote events to file ',outputstring,'.DAT   *'
         write(6,*) '*****************************************' 
         
      endif
!-------------- EXIT 

      return 
      end 


      subroutine gen_NLO_vegas(myinit,myitmx,
     &     myncall,mybin,xinteg,xerr)
  
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
      include 'gridinfo_cr.f' 
      include 'realwt.f' 
      include 'vegas_common_cr.f' 
      include 'constants.f' 
      integer myitmx,myncall,myinit,i,j,k,nproc,mynproc
      logical mybin,bin_cr
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsd,sumsdr
      double precision pborn(mxpart,4)
      character*4 part,mypart
       
    
      double precision region(2*mxdim_cr)
      logical first,myreadin
      logical theory_mode,total_cr
      double precision theory_mode_lord 
      common/techP/theory_mode,total_cr
      common/pborn/pborn 
      common/nproc/nproc
      common/bin_cr/bin_cr
      common/xreal/xreal,xreal2

      external genNLO_int
      data first/.true./
      save first
      integer ndim_in 

      ndim_in=ndim_cr

      xerr=0d0 
      xinteg=0d0 
      
      itmx_cr=myitmx
      bin_cr=mybin 
      
 !     write(6,*) 'Ndim in vegas =' ,ndim_Cr

      call boundregion(ndim_cr,region) 
      
      call vegasnr_cr(region,ndim_cr,genNLO_int,myinit,myncall,myitmx,
     .        0,xinteg,xerr,chi)
     
      ndim_cr=ndim_in 
      return 
      end 

