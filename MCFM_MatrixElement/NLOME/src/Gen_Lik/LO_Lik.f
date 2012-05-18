      subroutine LO_LIK(p,mur,muf,ansLO,errLO) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'calc_settings.f' 
      include 'scale.f'
      include 'facscale.f' 
      
      double precision p(n_final,4),mur,muf,ansLO,errLO
      double precision ansit(n_its) 
      integer it,i
      double precision LOME 
      integer int_flag,int_fr 
      integer idum 
      
      common/ranno/idum
      common/int_setup/int_flag,int_fr


!----- set up scale
      scale=mur
      facscale=muf
!------ Vegas parameters 
     
    
      if(int_flag.eq.1) then 
         ansLO=LOME(p) 
         errLO=1d-3*ansLO 
         return 
      endif
     
 111  continue 
      ansLO=0d0 
      errLO=0d0 
    
      if(it.gt.n_itsMX) return 
      do i=1,n_its 
         ansit(i)=0d0 
!         write(6,*) 'Using seed = ',idum 
!===== Calcualte Likelihoods          
         ansit(i)=LOME(p) 
         ansLO=ansLO+ansit(i)
!         write(6,*) 'Found ',i,ansit(i)
      enddo
!====== Calculate average and err  
      ansLO=ansLO/dfloat(n_its) 
      do i=1,n_its 
         errLO=errLO+(ansLO-ansit(i))**2
      enddo
      errLO=dsqrt(errLO/dfloat(n_its))


      return 
      end 
