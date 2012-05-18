
!----- routine for calcualting the phase space weight associated 
!----- with n massless final state partons, 
!----- if var = true then the final states are expressed in terms of 
!----- p_T,eta and phi 
!----- else  px,py and pz are used
      double precision function phase_weight(n,switch) 
      implicit none 
      include 'constants.f' 
      double precision pborn(mxpart,4)
      double precision sqrts
      double precision E_i(3:mxpart),Cal_E
      integer n,j 
      double precision E_sum
      logical switch 
      common/pborn/pborn 
      common/energy/sqrts

      do j=3,mxpart
         E_i=0d0 
      enddo

      phase_weight=1d0 
   

!----- Modified to correct for x dependence of Energy 
!----- Which comes from boost 
      do j=3,n+1 
         E_i(j)=Cal_E(j,n) 
      enddo

      if(switch) then 
         phase_weight=0d0 
!----- to be written 
         stop 
      else
         do j=3,n+1
            phase_weight=phase_weight/(2d0*(twopi)**3*E_i(j))
         enddo
      endif
    
      phase_weight=phase_weight*(four*pi)
      
      return 
      end 

      double precision function Cal_E(j,n) 
      implicit none 
      include 'constants.f' 
      integer n,j
      double precision pborn(mxpart,4),x1,x2
      common/x1x2/x1,x2 
      double precision xt1,xt2,p1x,p2x,p1y,p2y,E1,E2 
      double precision sqrts 
      common/energy/sqrts
      double precision root1,root2
      common/pborn/pborn
      double precision d_root

      root1=0d0 
      root2=0d0 
      xt1=sqrts/2d0*x1 
      xt2=sqrts/2d0*x2 

      p1x=pborn(1,1)
      p2x=pborn(2,1) 
      
      p1y=pborn(1,2)
      p2y=pborn(2,2) 
     
      d_root=Sqrt(xt1*(xt1 - xt2)**2*xt2*(-p1x**2 - p1y**2 + xt1*xt2))
      
      if((j.eq.3).and.(n.eq.2)) then
   
        root1=(xt1**2*xt2 + xt1*xt2**2 - d_root)/(2d0*xt1*xt2)
        root2=(xt1**2*xt2 + xt1*xt2**2 + d_root)/(2d0*xt1*xt2)
        
       
        
!-----ensure E > 0 
         Cal_E=root1
         
 
         write(6,*) 'root1 ',root1
         write(6,*) 'root2 ',root2
         write(6,*) ' x 1 = ',x1
         write(6,*) ' x 2 = ',x2
         write(6,*) 'pborn(3,4) = ',pborn(3,4) 
      else
         stop
      endif

      end

      double precision function pswt_jac(n) 
      implicit none 
      include 'constants.f' 
      double precision pborn(mxpart,4) 
      integer n
      double precision sqrts 
      common/pborn/pborn 
      common/energy/sqrts

      pswt_jac=one

      if(n.eq.2) then 
         pswt_jac=pborn(3,4)*pi/2d0 
      endif

      return 
      end 
