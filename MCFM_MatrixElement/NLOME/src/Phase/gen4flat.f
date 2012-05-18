      subroutine gen4flat(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common_cr.f' 
      include 'limits.f' 
      include 'process.f' 
      double precision p(n_final,4),r(mxdim_cr),k(4),p12(4),p34(4) 
      double precision sqrts,wt,smax,pswt,Qsq,mmax,m12,m34,s,c,sp,cp,phi 
      integer i
      common/energy/sqrts
      common/Qsq/Qsq

      wt=1d0
      mmax=110d0
      smax=sqrts/2d0
      call genP4(k,smax,pswt,r(1),r(2),r(3))
      wt=wt*pswt
      do i=1,3
         p12(i)=k(i)
         p34(i)=-k(i)
      enddo
      m12=r(4)*mmax
      wt=wt*mmax
      wt=wt*m12**2/pi
      p12(4)=sqrt(m12**2+p12(1)**2+p12(2)**2+p12(3)**2)
      m34=r(5)*mmax
      wt=wt*mmax
      wt=wt*m34**2/pi
      p34(4)=sqrt(m34**2+p34(1)**2+p34(2)**2+p34(3)**2)
!--- Generate p1 and p2 in restframe of p12
      c=2d0*r(6)-1d0
      s=sqrt(abs(1d0-c**2))
      wt=wt*2d0
      phi=twopi*r(7)
      sp=sin(phi)
      cp=cos(phi)
      wt=wt*twopi
      k(4)=0.5*m12
      k(1)=k(4)*s*sp
      k(2)=k(4)*s*cp
      k(3)=k(4)*c
      call booster(0,p12,k,k)
      do i=1,4
         p(1,i)=k(i)
         p(2,i)=p12(i)-k(i)
      enddo
!--- Generate p3 and p4 in restframe of p34
      c=2d0*r(8)-1d0
      s=sqrt(abs(1d0-c**2))
      wt=wt*2d0
      phi=twopi*r(9)
      sp=sin(phi)
      cp=cos(phi)
      wt=wt*twopi
      k(4)=0.5*m34
      k(1)=k(4)*s*sp
      k(2)=k(4)*s*cp
      k(3)=k(4)*c
      call booster(0,p34,k,k)
      do i=1,4
         p(3,i)=k(i)
         p(4,i)=p34(i)-k(i)
      enddo

!      write(*,*) p(1,1),p(1,2),p(1,3),p(1,4)
!      write(*,*) p(2,1),p(2,2),p(2,3),p(2,4)
!      write(*,*) p(3,1),p(3,2),p(3,3),p(3,4)
!      write(*,*) p(4,1),p(4,2),p(4,3),p(4,4)
!      write(*,*) sqrt((p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
!     .     -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2),", ",
!     .     sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
!     .     -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
!      write(*,*) '--------------------------------------------'

      Qsq=(p(1,4)+p(2,4)+p(3,4)+p(4,4))**2
      do i=1,3
         Qsq=Qsq-(p(1,i)+p(2,i)+p(3,i)+p(4,i))**2
      enddo

      return
      end subroutine
