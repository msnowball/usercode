      subroutine DPOLINT(XA,YA,N,X,Y,DY) 
      implicit none
      integer nmax
      parameter(nmax=10)
      integer N,M
      double precision XA(N),YA(N),C(nmax),D(nmax) 
      double precision X,Y,DY,W,HO,HP,den
      integer ns,i
      double precision dif,difT
      

      ns=1
      dif=abs(X-XA(1))
      
      do 11 i=1,N
         difT = abs(X-XA(i))
         if (difT .lt. dif) then 
            ns=i
            dif=difT
            endif
            c(i)=YA(i) 
            d(i)=YA(i) 
 11      continue
         Y=YA(ns) 
         ns=ns-1
         do 13 m=1,N-1
            do 12 I=1,N-M
               HO=XA(i)-X
               HP=XA(i+m)-X
               W=C(i+1)-d(i)
               den=HO-HP
               if (abs(den) .lt. 1d-10) then
                  write(6,*) ' stop in polint ',den
                  stop
               endif
               den=W/den
               d(i)=HP*den
               c(i)=HO*den
 12            continue
               if (2*ns .lt. n-m) then 
                  DY=C(ns+1) 
               else
                  DY=D(ns) 
                  ns=ns-1
               endif
               Y=Y+DY
 13            continue
               return 
               end
