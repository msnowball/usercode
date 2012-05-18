
      logical function includedipole(nd,p)
      implicit none 
!----- basic version of includdipole from MCFM 
      include 'constants.f' 
      double precision p(mxpart,4),my_pt
      integer j,k,nu,nd
      logical cuts 
      common/cuts/cuts
      
!----- initiailze 
      includedipole=.true. 
      

!------- jet clustering (to be written) 

!------ cuts 
      if(cuts) then 
         call docuts(nd,p,includedipole) 
        
         if(includedipole.eqv..false.) return
      endif 

      return 
      end 
