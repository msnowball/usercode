
      subroutine determine_xminmax(xmin,xmax,xl,p,*) 
      implicit none
      include 'constants.f' 
      include 'user_cuts.f' 
      include 'interface_settings.f' 
      include 'process.f' 
     
      double precision xmin,xmax,xl 
      double precision xlow(n_final),xup(n_final)
      double precision p(mxpart,4),eta_mx(n_final) 
      double precision taj,tjb,temp
      double precision ftild(n_final)
      double precision pa(4),pb(4),dot
      common/papb/pa,pb 
      character*2 plabel(mxpart)
      common/plabel/plabel
      integer i,j,ij
!===== xl = Q**2/sqrts**2
!===== x2 = xl/x1
!======== SETUP eta_mx array 
      
      xmin=xl
      xmax=0.99999999d0 
      temp=0d0 

     
  
      do i=3,n_final+2 
         if((plabel(i).eq.'ea').or.(plabel(i).eq.'el').or.
     &      (plabel(i).eq.'ma').or.(plabel(i).eq.'ml')) then 
            eta_mx(i-2)=eta_lept 
         elseif(plabel(i).eq.'ga') then 
            eta_mx(i-2)=eta_phot 
         else
            eta_mx(i-2)=1000d0 
         endif
      enddo

    
      
!======= For each final state particle calcualte constant factor and resutling f
      
      do j=1,n_final 
         ij=j+2
         taj=0d0 
         tjb=0d0 
         ftild(j)=0d0 
         do i=1,3 
            taj=taj-(p(1,i)+p(ij,i))**2
            tjb=tjb-(p(2,i)+p(ij,i))**2 
         enddo
         taj=taj+(p(1,4)+p(ij,4))**2
         tjb=tjb+(p(2,4)+p(ij,4))**2
         
        
      
         ftild(j)=-0.5d0*dlog(taj/tjb)-0.5d0*dlog(xl)        
   
!========= Calculate  xlow and xup 

      
         xlow(j)=dexp(-eta_mx(j)-ftild(j)) 
         xup(j) =dexp(+eta_mx(j)-ftild(j))
            

       

         if(xlow(j).gt.xmin) xmin=xlow(j) 
         if(xup(j).lt.xmax)  xmax = xup(j) 


      enddo
!      write(6,*) xmin,xmax
!      pause
      
      if(xmin.gt.xmax) return 1 
!      write(6,*) xmin,xmax
!      pause
 
      return 
      end 
         
         
      subroutine determine_xnd(nd,xmin,xmax,xl,p,*) 
      implicit none
      include 'constants.f' 
      include 'user_cuts.f' 
      include 'interface_settings.f' 
      include 'process.f' 
      include 'xstore.f'
      double precision xmin,xmax,xl 
      double precision xlow(n_final),xup(n_final)
      double precision p(mxpart,4),eta_mx(n_final) 
      double precision taj,tjb,temp
      double precision ftild(n_final)
      integer nd
     
      character*2 plabel(mxpart)
      common/plabel/plabel
      integer i,j,ij
!===== xl = Q**2/sqrts**2
!===== x2 = xl/x1
!======== SETUP eta_mx array 
      
      xmin=xl
      xmax=0.99999999d0 
      temp=0d0 

     
  
      do i=3,n_final+2 
         if((plabel(i).eq.'ea').or.(plabel(i).eq.'el').or.
     &      (plabel(i).eq.'ma').or.(plabel(i).eq.'ml')) then 
            eta_mx(i-2)=eta_lept
         elseif(plabel(i).eq.'ga') then 
            eta_mx(i-2)=eta_phot+2.5d0
         else
            eta_mx(i-2)=1000d0 
         endif
      enddo
      
!======= For each final state particle calcualte constant factor and resutling f
      
      do j=1,n_final 
         ij=j+2
         taj=0d0 
         tjb=0d0 
         ftild(j)=0d0 
         do i=1,3 
            taj=taj-(p(1,i)+p(ij,i))**2
            tjb=tjb-(p(2,i)+p(ij,i))**2 
         enddo
         taj=taj+(p(1,4)+p(ij,4))**2
         tjb=tjb+(p(2,4)+p(ij,4))**2
       

         ftild(j)=-0.5d0*dlog(taj/tjb)-0.5d0*dlog(xl)
     &        -0.5d0*dlog(xstore(nd,2)/xstore(nd,1))        
   
!========= Calculate  xlow and xup 

      
         xlow(j)=dexp(-eta_mx(j)-ftild(j)) 
         xup(j) =dexp(+eta_mx(j)-ftild(j))
            
         if(xlow(j).gt.xmin) xmin=xlow(j) 
         if(xup(j).lt.xmax)  xmax = xup(j) 

      enddo
      
      if(xmin.gt.xmax) return 1 
      
 
      return 
      end 
         
         
