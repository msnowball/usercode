      subroutine gen3flat(p,r,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'vegas_common_cr.f' 
      include 'limits.f' 
      include 'process.f' 
      double precision p(n_final,4),r(mxdim_cr),k(4),p12(4)
      double precision sqrts,wt,pswt,smax,Qsq,mmax,m12,s,c,sp,cp,phi
      double precision mmin
      integer i
      common/energy/sqrts
      common/Qsq/Qsq

      mmin=80d0
      mmax=100d0
      wt=1d0
!--- Generate photon and p12=p1+p2 
      smax=sqrts
      call genP4(k,smax,pswt,r(1),r(2),r(3))
      wt=wt*pswt
      do i=1,4
         p(3,i)=k(i)
         p12(i)=-k(i)
      enddo
      m12=mmin+r(4)*(mmax-mmin)
      wt=wt*(mmax-mmin)
      wt=wt*m12**2/pi
      p12(4)=sqrt(m12**2+p12(1)**2+p12(2)**2+p12(3)**2)
!--- Generate p1 and p2 in restframe of p12
      c=2d0*r(5)-1d0
      s=sqrt(abs(1d0-c**2))
      wt=wt*2d0
      phi=twopi*r(6)
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
      Qsq=(p(1,4)+p(2,4)+p(3,4))**2
      do i=1,3
         Qsq=Qsq-(p(1,i)+p(2,i)+p(3,i))**2
      enddo
!      write(*,*) p(1,1),p(1,2),p(1,3),p(1,4)
!      write(*,*) p(2,1),p(2,2),p(2,3),p(2,4)
!      write(*,*) p(3,1),p(3,2),p(3,3),p(3,4)
!      write(*,*) sqrt((p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
!     .     -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2)
!      write(*,*) '--------------------------------------------'

      return
      end subroutine
