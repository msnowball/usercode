!------ This routine will initiate NLOME by reading the interface_settings and 
!------ Reading the file NLOME_params.DAT


      subroutine init_NLOME(mu_R,mu_F)
      implicit none 
      include 'interface_settings.f' 
      include 'alfacut.f' 
      include 'betacut.f'
      include 'zerowidth.f' 
      include 'removebr.f'
      include 'scale.f' 
      include 'facscale.f' 
      double precision mu_R,mu_F
      
      
!------ setup some non-zero scale for correct initialisation 
      scale=mu_R
      facscale=mu_F

!----- inital setup 
      call read_interface_settings
!----- Read user params
      call read_params 
    
!------Setup arrays 
      call basic_setup(nproc_1) 
      call read_cuts
      call readcoup
!------- MCFM parametrs (shouldnt normally be changed unless you know what your doing!) 

!----- REMOVE BR, we dont want to do this
      removebr=.false. 
!------ ZEROWIDTH We dont want to do this either 
      zerowidth=.false. 

!-----alpha parameters 
      aii=0.5d0
      aif=1d0 
      afi=1d0 
      aff=1d0 
      bfi=1d0 
      bff=1d0



      return 
      end 

