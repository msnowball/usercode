      subroutine dromo(func,a,b,ss,choose)
      implicit none 
      double precision a,b,func,ss,eps 
      external func,chooser 
      integer jmax,jmaxp,k,km
      parameter(eps=1d-3,JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1)
      integer j 
      double precision dss,h(jmaxp),s(jmaxp) 

      h(1)=1d0 
      do j=1,JMAX
         call choose(func,a,b,s(j),j) 
         if( j.ge.K) then 
            call dpolint(h(j-KM),s(j-KM),K,0d0,ss,dss) 
            if(dabs(dss).le.eps*dabs(ss)) return 
         endif
         s(j+1)=s(j) 
         h(j+1)=h(j)/9d0 
      enddo
      write(6,*) 'Too many steps in dromo'
      stop
      end 
 3    
      
