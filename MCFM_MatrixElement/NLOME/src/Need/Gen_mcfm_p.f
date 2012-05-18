

!----- combines a p_born with p1 p2 and pr the output from the forward branching routine to create 
!----- and MCFM phase space point 

      subroutine Gen_mcfm_p(p1,p2,p3,p) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      double precision pborn(mxpart,4),p(mxpart,4),p1(4),p2(4),p3(4) 
      integer j,nu,ico
      common/pborn/pborn 
 

      do j=1,mxpart 
         do nu=1,4 
            p(j,nu)=0d0 
         enddo
      enddo


      
!---- put in new p1 and p2 
      do nu=1,4
         p(1,nu)=p1(nu) 
         p(2,nu)=p2(nu) 
      enddo
      
!----- copy over p_b 
      do j=3,n_particles
!----- find position of last particle
         do nu=1,4  
            p(j,nu)=pborn(j,nu) 
         enddo
      enddo
     
!---- copy over p_r into final slot 
      do nu=1,4 
         p(n_particles+1,nu)=p3(nu)
      enddo
      
      return 
      end 
      
      
      
