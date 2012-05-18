
      double precision function my_simp_born(a,b,n_in,Qsq,log_mode) 
!-----simpsons rule, a b are integration limits divide into n pieces 
!-----integrate function f (subroutine in our case)
      implicit none
      include 'constants.f' 
      integer n,flag,i,n_in
      double precision h,fa,fb
      double precision z
      double precision int_dips
      double precision a,b,Qsq
      double precision my_LO 
      logical log_mode 
      double precision xjac 

     
      
      my_simp_born=0d0 
     
     
      n=2*n_in
      h=(b-a)/dfloat(n)
   
      if(log_mode) then 
         my_simp_born=my_LO(Qsq,dexp(a))*dexp(a)
      else
         my_simp_born=my_LO(Qsq,a)
      endif
      
      

     

      do i=1,n,2 
         z=a+h*dfloat(i) 
         if(log_mode) then             
            my_simp_born=my_simp_born+4d0*my_LO(Qsq,dexp(z))*dexp(z) 
         else
            my_simp_born=my_simp_born+4d0*my_LO(Qsq,z)
      endif
      enddo
      
      do i=2,n-1,2
         z=a+h*dfloat(i)
         if(log_mode) then 
            my_simp_born=my_simp_born+2d0*my_LO(Qsq,dexp(z))*dexp(z) 
         else
            my_simp_born=my_simp_born+2d0*my_LO(Qsq,z)
         endif
         
       
      enddo
      
      if(log_mode) then 
           my_simp_born=my_simp_born+my_LO(Qsq,dexp(b))*dexp(b) 
         else
            my_simp_born=my_simp_born+my_LO(Qsq,b)
         endif
         
         
         
         my_simp_born=h*my_simp_born/3d0 
      
         return 
      end 
      
      
      double precision function my_LO(Qsq,z) 
      implicit none 
      include 'constants.f' 
      include 'facscale.f' 
      double precision lo_me(-nf:nf,-nf:nf)
      common/lo_me/lo_me 
      integer ih1,ih2,j,k,i
      double precision xx(2),Qsq,Q_f(4),sqrts,xzz,z
      double precision fx1(-nf:nf),fx2(-nf:nf),W,flux
      common/energy/sqrts 
      common/density/ih1,ih2

      my_LO=0d0 
      xx(1)=z 
      xx(2)=Qsq/(xx(1)*sqrts**2) 


      if(xx(1).gt.0.99d0) return 
      if(xx(2).gt.0.99d0) return 

     
      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)

     
      do i=-nf,nf 
         fx1(i)=0d0 
         fx2(i)=0d0 
      enddo

    
      
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
       
       xzz=0d0 
       
!------------------ ZZ loop over initial state flavors
       do j=-nf,nf
          do k=-nf,nf
             xzz=fx1(j)*fx2(k)*(lo_me(j,k))*flux             
            
             my_LO=my_LO+xzz
          enddo
       enddo
      
     
      return 
      end 

         

