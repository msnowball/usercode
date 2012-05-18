!----- implementation of the midpoint algorithim for integration 

      subroutine midpnt(func,a,b,s,n) 
      implicit none 
      integer n 
      double precision a,b,s,func 
      external func 
      integer it,j 
      double precision ddel,del,sum,tnm,x 
      
      if(n.eq.1) then 
         s=(b-a)*func(0.5d0*(a+b))
      else
         it=3**(n-2) 
         tnm=it 
         del=(b-a)/(3d0*tnm) 
         ddel=del+del
         x=a+0.5*del 
         sum=0d0 
         do j=1,it 
            sum=sum+func(x) 
            x=x+ddel
            sum=sum+func(x) 
            x=x+del 
         enddo
         s=(s+(b-a)*sum/tnm)/3d0 
      endif
      return 
      end
      
