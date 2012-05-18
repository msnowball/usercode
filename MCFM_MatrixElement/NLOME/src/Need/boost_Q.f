

      subroutine boost_Q(p_in,p_out) 
      implicit none
      include 'constants.f' 
      include 'Gen_Q.f' 
      include 'interface_settings.f' 
      double precision p_in(n_final,4),p_out(n_final,4) 
      integer nu,i 
      double precision vec1(4),vec2(4),pT(4)
      double precision beta_z,gam_z
      double precision beta_st(4),root_bz
      common/beta_st/beta_st
      
!      write(6,*) n_final
!      write(6,*) 'in boost ',p_in
!------ calculate the missing pT vector we would like to boost in, 
      do nu=1,4
         pT(nu)=0d0
    
         do i=1,n_final
            pT(nu)=pT(nu)-p_in(i,nu) 
         enddo 
      enddo

!---- booost x 
      do i=1,n_final 
         do nu=1,4 
            vec1(nu)=p_in(i,nu) 
            vec2(nu)=0d0
         enddo      
         call my_boost(1,pT,vec1,vec2)
         do nu=1,4 
            p_out(i,nu)=vec2(nu)
         enddo
      enddo
      
!--- recast PT under x boost 
      do nu=1,4 
         vec1(nu)=pT(nu) 
         vec2(nu)=0d0 
      enddo
      call my_boost(1,pT,vec1,vec2) 
      do nu=1,4 
         pT(nu)=vec2(nu) 
      enddo

    
!---- boost y 
      do i=1,n_final
         do nu=1,4 
            vec1(nu)=p_out(i,nu) 
            vec2(nu)=0d0
         enddo      
         call my_boost(2,pT,vec1,vec2)
         do nu=1,4 
            p_out(i,nu)=vec2(nu)
         enddo
      enddo

!----- Recast PT under y boost 
      do nu=1,4 
         vec1(nu)=pT(nu) 
         vec2(nu)=0d0 
      enddo
      call my_boost(2,pT,vec1,vec2) 
      do nu=1,4 
         pT(nu)=vec2(nu) 
      enddo
      
!----- Final vector is given by sqrt(Q**2)/2 (0,0,x1-x2,x1+x2) 
!----- Currently have it in the form  (0,0,Qz,EQ), Need final boost to 
!----- Setup Q given x1-x2
      
!      write(6,*) Q_in 
!      write(6,*) pT

!      root_bz=dsqrt(pT(3)**2*pT(4)**2-(Q_in(3)**2
!     &     +pT(4)**2)*(pT(3)**2-Q_in(3)**2))
!      beta_z=(pT(3)*pT(4)+root_bz)/((Q_in(3)**2+pT(4)**2))
!      write(6,*) beta_z 
!      beta_z=(pT(3)*pT(4)-root_bz)/((Q_in(3)**2+pT(4)**2))
!      write(6,*) beta_z 

     
!      gam_z=one/dsqrt(one-beta_z**2) 
      
    
!      do i=1,4    
!         vec1(3)=p_out(i,3) 
!         vec1(4)=p_out(i,4) 
!         p_out(i,3)=gam_z*(vec1(3)-beta_z*vec1(4))
!         p_out(i,4)=gam_z*(vec1(4)-beta_z*vec1(3))
!      enddo

!      do i=1,4 
!         do nu=1,4
!            vec1(nu)=p_out(i,nu) 
!            vec2(nu)=0d0 
!         enddo 
!         call my_boost(3,Q_in,vec1,vec2) 
!         do nu=1,4 
!            p_out(i,nu)=vec2(nu) 
!         enddo
!      enddo
       
!--- writeout check
!      do i=1,4
!         write(6,*) 'p_in (',i,') = ',p_in(i,1),p_in(i,2)
!     &        ,p_in(i,3),p_in(i,4)
!         write(6,*) 'p**2 = ',p_in(i,4)**2-(p_in(i,1)**2
!     &        +p_in(i,2)**2+p_in(i,3)**2)
!      enddo

!      write(6,*) '*****************' 
!      do i=1,4
!         write(6,*) 'p_out (',i,') = ',p_out(i,1),p_out(i,2)
!     &        ,p_out(i,3),p_out(i,4) 
!       write(6,*) 'p**2 = ',p_out(i,4)**2-(p_out(i,1)**2
!     &        +p_out(i,2)**2+p_out(i,3)**2)
!      enddo

!      write(6,*) p_out(1,1)+p_out(2,1)+p_out(3,1)+p_out(4,1)
!      write(6,*) p_out(1,2)+p_out(2,2)+p_out(3,2)+p_out(4,2)

      return 
      end




      subroutine my_boost(i,p1,p2,p3) 
      implicit none 
      include 'store_boost.f' 
      !--- boost p2 by p1 in direction and store in p3 
      double precision p1(4),p2(4),p3(4) 
      integer i,j,nu
      double precision gamma,beta 
      double precision beta_st(4) 
      common/beta_st/beta_st

      if((i.gt.4).or.(i.lt.1)) then 
         write(6,*) 'Invalid boost param i=',i 
         stop 
      endif

      

      beta=p1(i)/p1(4) 
      beta_st(i)=beta
      gamma=1d0/dsqrt(1d0-beta**2) 
      gamma_st(i)=gamma 

      p3(4)=gamma*(p2(4)-beta*p2(i))

      do nu=1,3
         if(nu.eq.i) then
            p3(nu)=gamma*(-beta*p2(4)+p2(nu))
         else
            p3(nu)=p2(nu) 
            
         endif
      enddo

      return 
      end subroutine 

      subroutine my_invboost(pin,pout) 
      implicit none 
      include 'constants.f' 
!----- subroutine takes an mcfm input phase space pin, and returns the lab frame pout 
      double precision pin(mxpart,4),pout(mxpart,4)  
      integer nu,i,j,k 
      double precision Inv_Lam(4,4),Inv_LamX(4,4),Inv_LamY(4,4)
      double precision vec1(4),vec2(4),pT(4),pref_x,pref_y
      double precision beta_st(4),gamma(2)
      common/beta_st/beta_st
      
!----- beta_st contains the beta_x and beta_y parameters, calcualtes in the form  X'  = Lam(beta_y)Lam(beta_x)*X      
!----- new matrix is     Inv_Lam(beta_x)*Inv_Lam(beta_y) X' = X 
      do i=1,2
         gamma(i)=one/dsqrt(one-beta_st(i)**2)
      enddo
      
!----- Calculate Inv_Lam(beta_y) 

      do i=1,4 
         do j=1,4 
            Inv_LamX(i,j)=0d0 
            Inv_LamY(i,j)=0d0 
            Inv_Lam(i,j)=0d0 
         enddo
      enddo

!--------- Diagaonals unchanged 
      Inv_LamX(1,1)=gamma(1) 
      Inv_LamY(1,1)=one 
      do i=2,3
         Inv_LamX(i,i)=one 
      enddo
      Inv_LamY(2,2)=gamma(2) 
      Inv_LamY(3,3)=one 
      Inv_LamX(4,4)=gamma(1)
      Inv_LamY(4,4)=gamma(2) 
      
!------------Off diagonals (minus sign)  
      i=1
      Inv_LamX(i,4)=beta_st(i)*gamma(i) 
      Inv_LamX(4,i)=beta_st(i)*gamma(i)
      i=2
      Inv_LamY(i,4)=beta_st(i)*gamma(i) 
      Inv_LamY(4,i)=beta_st(i)*gamma(i)
      

!------ prefactors 
      pref_x=gamma(1)**2*(one-beta_st(1)**2)
      pref_y=gamma(2)**2*(one-beta_st(2)**2)
      pref_x=one/pref_x
      pref_y=one/pref_y 

!----- Total Matrix 
      do i=1,4 
         do j=1,4
            Inv_Lam(i,j)=0d0 
            Inv_LamX(i,j)=pref_x*Inv_LamX(i,j) 
            Inv_LamY(i,j)=pref_y*Inv_LamY(i,j) 
           
            do k=1,4   
               Inv_Lam(i,j)=Inv_Lam(i,j)+Inv_LamX(i,k)*Inv_LamY(k,j) 
            enddo
         enddo
      enddo

!------ Make new phase space. 
    
      do j=1,mxpart
         do nu=1,4 
            pout(j,nu)=0d0 
            do k=1,4
               pout(j,nu)=pout(j,nu)+Inv_Lam(nu,k)*pin(j,k)
            enddo
         enddo
      enddo

      !--- writeout check
      do i=1,7
         write(6,*) 'pin (',i,') = ',pin(i,1),pin(i,2)
     &        ,pin(i,3),pin(i,4)
         write(6,*) 'p**2 = ',pin(i,4)**2-(pin(i,1)**2
     &        +pin(i,2)**2+pin(i,3)**2)
      enddo
 
      write(6,*) '**********************************' 
      write(6,*) 

      do i=1,7
         write(6,*) 'pout (',i,') = ',pout(i,1),pout(i,2)
     &        ,pout(i,3),pout(i,4) 
         write(6,*) 'p**2 = ',pout(i,4)**2-(pout(i,1)**2
     &        +pout(i,2)**2+pout(i,3)**2)
      enddo

      write(6,*) pout(1,1)+pout(2,1)+pout(3,1)+pout(4,1)
     &    +pout(5,1)+pout(6,1)+pout(7,1)
      write(6,*) pout(1,2)+pout(2,2)+pout(3,2)+pout(4,2)
     &     +pout(5,2)+pout(6,2)+pout(7,2)

      
      return 
      end subroutine 
      
