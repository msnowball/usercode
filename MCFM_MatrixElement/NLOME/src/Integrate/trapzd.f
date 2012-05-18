


      subroutine trapzd(func,a,b,s,n) 
      implicit none 
      integer n 
      double precision a,b,s,func 
      external func 
      integer it,j 
      double precision del,sum,tnm,x 
      
      if(n.eq.1) then 
         s=0.5d0*(b-a)*(func(a)+func(b))
      else
         it=2**(n-2) 
         tnm=it 
         del=(b-a)/tnm 
         x=a+0.5d0*del 
         sum=0d0 
         do j=1,it
            sum=sum+func(x) 
            x=x+del 
         enddo
         s=0.5d0*(s+(b-a)*sum/tnm) 
      endif
      return 
      end
