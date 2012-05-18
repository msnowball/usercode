
      double precision function my_simp(p,a,b,n_in,flag) 
!-----simpsons rule, a b are integration limits divide into n pieces 
!-----integrate function f (subroutine in our case)
      implicit none
      include 'constants.f' 
      integer n,flag,i,n_in
      double precision h,fa,fb
      double precision z,p(mxpart,4)
      double precision int_dips
      double precision a,b
     
      my_simp=0d0 
     
     
      n=2*n_in
      h=(b-a)/dfloat(n)
      my_simp=int_dips(p,a,flag)
     

      do i=1,n,2 
         z=a+h*dfloat(i) 
         my_simp=my_simp+4d0*int_dips(p,z,flag) 
      enddo
      
      do i=2,n-1,2
         z=a+h*dfloat(i) 
         my_simp=my_simp+2d0*int_dips(p,z,flag) 
      enddo
      
      my_simp=my_simp+int_dips(p,b,flag) 

      my_simp=h*my_simp/3d0 
      
      return 
      end 


         
         

