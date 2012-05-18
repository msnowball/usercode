************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 1999        
!===== MODIFIED FOR NLOME 
!===== Need to replace initial-inital with something intial-final 
!====== for correct momentum sturcture                                 *
*     Calculates the nj-jet subtraction term corresponding to dipole   *
*     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
*     Automatically chooses dipole kind                                *
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter parton                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************

      subroutine dips_fb_if(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     . subr_born,subr_corr)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      double precision x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),vtilde
      integer nd,ip,jp,kp,nu,j,k,ipt,i
      integer ii,ij,ik
      logical includedipole
      logical incldip(0:maxd)
      common/incldip/incldip
      external subr_born,subr_corr
            
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      call zeromsq(msq,msqv)
      incldip(nd)=.true.
      
      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

    
    
***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************
      if ((ip .le. 2) .and. (kp .le. 2)) then
         
         x=(sij+sik+sjk)/(sij+sik)
         omx=one-x
         
         z=sij/(sij+sik) 
         omz=one-z
C---npart is the number of particles in the final state
C---transform the momenta so that only the first npart+1 are filled
         call transform_fb_if(p,ptrans,x,ip,jp,kp)
         call storeptilde(nd,ptrans)
       
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
C--   if not return
c        if (incldip(nd) .eqv. .false.) return

c--- if using a dynamic scale, set that scale with dipole kinematics	
!	if (dynamicscale) then
!	  call scaleset(initscale,initfacscale,ptrans)
!	  dipscale(nd)=facscale
!	endif
	
c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq)
C---Modification so that only close to singular subtracted
C---Do not set incldip because initial-final can fail 
C---but final initial needs still to be tested

!        if ((case .eq. 'H_1jet') .and. (ip.eq.1) .and. (jp.eq.6)) then
c--- do nothing
!        else
!          if (u .gt. aif) return
!        endif
        
!        do nu=1,4
!           vec(nu)=(p(jp,nu)/u-p(kp,nu)/omu)/dsqrt(sjk)
!        enddo
        
!        call subr_corr(ptrans,vec,ip,msqv) 

        sub(qq)=-gsq/x/sij*(two/(omx-z)-one-x)
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
  

        return 

!===== NEEDS TO BE RE-IMPLEMENTED .... 
      endif
        
      return
      end
      
