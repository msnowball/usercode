

      subroutine genW(p,r,wt,*) 
      implicit none 
!------- integrate the following 
!   1/(2pi)^2 d p_T = 1/(2pi)^2 \int_{ET MIN}^{sqrts} ET dET dphi
      include 'user_cuts.f' 
      include 'constants.f' 
      include 'vegas_common_cr.f' 
      include 'interface_settings.f'

      double precision p(n_final,4),Et_gen,r(mxdim_cr)
      double precision wt0,wt
      parameter(wt0=one/(8d0*pi**2))
      double precision sqrts 
      common/energy/sqrts
      double precision et_miss_max,etr,lnter
      double precision sphi,cphi,phi
      integer nproc,nu 
      common/nproc/nproc
      double precision tsoft 
      parameter(tsoft=1d-1)

      wt=wt0 
      
!====== Generate ET_miss logaritimihcally 
      et_miss_min=0d0
      if(et_miss_min.lt.tsoft) et_miss_min=tsoft
      et_miss_max=sqrts
      etr=et_miss_min/et_miss_max 
      lnter=dlog(etr)
      et_gen=et_miss_max**(r(1))*et_miss_min**(one-r(1))
      wt=wt*et_gen**2*lnter 
 
      phi=twopi*r(2) 
      sphi=dsin(phi) 
      cphi=dcos(phi)
      wt=wt*twopi 

!------ check nproc to see whchi is neutrino 
!1  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))' 'N'
!6  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))' 'N'
!      
      if(nproc.eq.1) then 
         p(1,1)=et_gen*cphi 
         p(1,2)=et_gen*sphi 
         p(1,3)=0d0 
         p(1,4)=0d0 
         do nu=1,4 
            p(2,nu)=-p(1,nu) 
         enddo

      elseif(nproc.eq.6) then 
         p(2,1)=et_gen*cphi 
         p(2,2)=et_gen*sphi 
         p(2,3)=0d0 
         p(2,4)=0d0 
         do nu=1,4 
            p(1,nu)=-p(2,nu) 
         enddo
      else
         write(6,*) 'Not a W production process ',nproc
         return 1 
      endif

      return  
      end 

         
         
