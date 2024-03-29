!---------- Generic Frixione Routine, will take eps and delta_0 from input.DAT 
!--------   n will be hard-coded but can be changed below. 
!--------- C. Williams July 2011 


!----- p -momentum array 
!----- passed, should be obvious!
!----- j - photon identification in p i.e. p(j,i) = photon(i) 
!----- isub wether we are working with a dipole or not 
      subroutine frix(p,passed,j,isub) 
      implicit none 
      include 'constants.f' 
      include 'leptcuts.f' 
      include 'frag.f' 
      include 'npart.f'
      double precision p(mxpart,4),my_R,pref,ret_ET
      integer j,isub
      logical passed,first,is_hadronic,in_cone_n
      integer i,n_pow
      

      data first /.true. /

!---- use this to change n in Frixione isolation 
      n_pow = 1
      passed=.true. 


      if(first) then 
         first=.false. 
       write(6,*)'************** Frixione Isolation    ***************'
       write(6,*)'*                                                  *'
      write(6,99)'*  eps_phot = ',epsilon_h,' delta_0 = ',cone_ang,  '*'
       write(6,97)'*  n = ',n_pow,'                                   *'     
       write(6,*)'****************************************************'
      endif

 99   format (1x,a14,f5.3,a12,f5.3,a16)
 97   format (1x,a7,i1,a44)

!----- Cycle over final state particles, if it is hadronic and inside 
!----- photon cone then check its energy passes frixione requirement 
!----- else fail and return 
     
      
      
      pref=ret_ET(p,j)*epsilon_h/((1d0-dcos(cone_ang))**n_pow)
      do i=3,3+npart-isub
         if(is_hadronic(i).and.(my_R(i,j,p).lt.cone_ang)) then 
            passed=in_cone_n(my_R(i,j,p),ret_ET(p,i),pref,n_pow)
            if(passed.eqv..false.) return            
         endif
      enddo

      return 
      end

            

      
      double precision function ret_ET(p,j) 
      implicit none 
      include 'constants.f' 
      double precision p(mxpart,4)
      integer j 
      double precision ptsq,my_pt
      
!      ptsq=p(j,1)**2+p(j,2)**2
       ret_ET=my_pt(j,p) 
!      ret_ET=p(j,4)*dsqrt(ptsq)/(dsqrt(ptsq+p(j,3)**2))
      
      return 
      end
      
      logical function in_cone_n(Rij,Ejet,pref,n)
      implicit none
      double precision Rij,Ejet,pref
      integer n

      if ( Ejet .lt. pref*(1d0-dcos(Rij))**n) then 
         in_cone_n = .true. 
      else
         in_cone_n = .false. 
      endif

      return 
      end





      

      
      
      
      

      
