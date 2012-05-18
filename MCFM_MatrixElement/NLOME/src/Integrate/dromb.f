
      subroutine dromb(func,a,b,ss) 
      implicit none 
      double precision func,a,b,ss,eps 
      external func 
      integer jmax,jmaxp,k,km
      parameter(eps=1e-5,jmax=40,jmaxp=jmax+1,K=5,KM=K-1) 
      integer j 
      double precision dss,h(jmaxp),s(jmaxp) 
      h(1)=1
      do j=1,jmax
         call trapzd(func,a,b,s(j),j) 
         if(j.ge.K) then 
            call polint(h(j-KM),s(j-KM),K,0d0,ss,dss) 
            if(dabs(dss).le.eps*dabs(ss)) return 
         endif
         s(j+1)=s(j) 
         h(j+1)=0.25d0*h(j) 
      enddo
      write(6,*) 'too many steps in dromb' 
      stop 
      end 
