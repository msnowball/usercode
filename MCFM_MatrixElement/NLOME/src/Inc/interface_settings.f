
      
      
      logical efficient,params 
      integer n_particles,n_final
      integer n_out
      integer nproc_1,nproc_2
      logical inc_sig
      logical read_ew,read_mv,read_mq,read_ml,read_ckm,read_as
      logical read_scheme
      logical inherit_as,ret_poles
      double precision inh_as 
      integer pid
      logical alphas_eq1,alphaew_eq1
      character*4 inscheme 
      logical inc_MET 
      integer no_neutrinos 
      common/MET_params/inc_MET,no_neutrinos

      common/inc_sig/inc_sig
      common/nprocs/nproc_1,nproc_2
      common/n_out/n_out 
      common/inscheme/inscheme
      common/efficient/efficient
      common/pid/pid
     
      common/n_particles/n_particles,n_final
      common/params/params
      common/read_p/read_ew,read_mv,read_mq,read_ml,read_ckm,read_as
     &     ,read_scheme
      common/alpha_eq1/alphas_eq1,alphaew_eq1
!      common/inh_as/inherit_as,inh_as
      common/ret_poles/ret_poles
