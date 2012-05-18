
      subroutine my_phase6(p) 
      implicit none 
!      include 'constants.f' 
      double precision p(40,4) 
      integer nu,i 
  


!p1  0.0000000000000000E+00  0.0000000000000000E+00 -0.9890272191481677E+02 -0.9890272191481677E+02
! p2  0.0000000000000000E+00  0.0000000000000000E+00  0.1532608498813800E+04 -0.1532608498813800E+04
! p3  0.9453605882857718E+02 -0.7772037264976507E+01 -0.9140569541400263E+03  0.9189654979344058E+03
! p4  0.9024592929740407E+01 -0.3909068861825349E+02 -0.5649453471553933E+03  0.5663680521416483E+03
! p5  0.1045620054917498E+02 -0.2129756387200015E+00 -0.6584333127741990E+01  0.1235843563256839E+02
! p6 -0.1140168523074926E+03  0.4707570152195000E+02  0.5188085752417854E+02  0.1338192350199941E+03

      do i=1,40 
         do nu=1,4
            p(i,nu)=0d0 
         enddo
      enddo

!--- A 6-body phase space point for testing interface 
!--- MCFM notation 
      p(1,1)=0d0
      p(1,2)=0d0 
      p(1,3)=-0.9890272191481677d2 
      p(1,4)=p(1,3)
      
      p(2,1)=0d0
      p(2,2)=0d0 
      p(2,3)= 0.1532608498813800d4
      p(2,4)=-p(2,3)
      
      p(3,1)= 0.9453605882857718d2 
      p(3,2)=  -0.7772037264976507d1 
      p(3,3)=  -0.9140569541400263d3  
      p(3,4)=  0.9189654979344058d3

      p(4,1)=0.9024592929740407d1 
      p(4,2)=-0.3909068861825349d2 
      p(4,3)= -0.5649453471553933d3 
      p(4,4)= 0.5663680521416483d3
      
      p(5,1)=   0.1045620054917498d2 
      p(5,2)=     -0.2129756387200015d0 
      p(5,3)=     -0.6584333127741990d1 
      p(5,4)=         0.1235843563256839d2

      p(6,1)=      -0.1140168523074926d3 
      p(6,2)= 0.4707570152195000d2 
      p(6,3)= 0.5188085752417854d2
      p(6,4)=0.1338192350199941d3
     
      

      !Check 
 !     do nu=1,4 
 !        write(6,*) p(1,nu)+p(2,nu)+p(3,nu)+p(4,nu)+p(5,nu) 
 !     enddo
      ! 
 !     do nu=1,5 
 !        write(6,*) p(nu,1)**2+p(nu,2)**2+p(nu,3)**2-p(nu,4)**2
 !     enddo

      return 
      end subroutine 

      subroutine my_phase6_BLHA(pin,pout)
      implicit none  
      double precision pin(40,4),pout(6,4) 
      integer i,nu
 ! BLHA notation 
 ! MCFM x,y,z,E 
 ! BLH  E,x,y,z
      
      do i=1,6
         pout(i,1)=pin(i,4) 
         pout(i,2)=pin(i,1) 
         pout(i,3)=pin(i,2) 
         pout(i,4)=pin(i,3)
!         pout(i,5)=pin(i,4)**2-pin(i,3)**2-pin(i,2)**2-pin(i,1)**2
      enddo
      
      return 
      end subroutine
      
      
!----- CMS PHASE SPACE POINT (NOTE MCFM notation from email, code expects in the form of BLHA (E,x,y,z) 

      subroutine CMS_Phase(p_cms,st) 
      implicit none 
      double precision p_cms(4,4),swap
      integer st,nu,j 
      logical convert_BLHA 

      convert_BLHA=.false. 

      if(st.eq.1) then 
!---- CMS p_cmsoint 1 (MCFM notation) 
         p_cms(1,1)=35.60d0    
         p_cms(1,2)=42.47d0
         p_cms(1,3)=-55.25d0
         p_cms(1,4)=78.25d0
         p_cms(2,1)=-34.16d0
         p_cms(2,2)= 39.43d0
         p_cms(2,3)=   8.69d0 
         p_cms(2,4)=  52.89d0
         p_cms(3,1)=   -46.17d0
         p_cms(3,2)= -60.85d0
         p_cms(3,3)=  -28.86d0
         p_cms(3,4)=   81.66d0
         p_cms(4,1)=     35.17d0
         p_cms(4,2)=  -14.1d0 
         p_cms(4,3)=  -44.67d0
         p_cms(4,4)=  58.58d0
      elseif(st.eq.2) then 
!----- CMS p_cmsoint 2 (MCFM notation) 
         p_cms(1,1)= -8.141d0    
         p_cms(1,2)=57.81d0  
         p_cms(1,3)=-28.57d0   
         p_cms(1,4)= 65.00d0
         p_cms(2,1)=    0.226d0 
         p_cms(2,2)=  -34.57d0 
         p_cms(2,3)=  -8.60d0  
         p_cms(2,4)=  35.62d0
         p_cms(3,1)=     4.063d0  
         p_cms(3,2)= -3.929d0 
         p_cms(3,3)=  -14.41d0  
         p_cms(3,4)=  15.48d0
         p_cms(4,1)=  -4.43d0    
         p_cms(4,2)=-3.171d0   
         p_cms(4,3)= 9.84d0  
         p_cms(4,4)=  11.25d0
      else
         write(6,*) 'UNKNOWN CMS TEST P_CMSOINT' 
         stop
      endif

      if(convert_BLHA) then 
! MCFM x,y,z,E 
! BLH  E,x,y,z
         do j=1,4  
            swap=p_cms(j,1) 
            p_cms(j,1)=p_cms(j,4) 
            p_cms(j,4)=p_cms(j,3) 
            p_cms(j,3)=p_cms(j,2) 
            p_cms(j,2)=swap 
         enddo
         
      endif
           
      return 
      end subroutine 
