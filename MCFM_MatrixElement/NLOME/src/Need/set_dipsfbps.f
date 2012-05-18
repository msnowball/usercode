!===== routine for setting up FBPS style dipoles 
!===== currently implementd for the branching of an initial state quark => gluon 
!===== dipole will be stored in D_sub
      subroutine set_dipsfbps(x,ib,pa,pb,pr,Dsub) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f'
      include 'ptilde.f' 
      include 'qcdcouple.f' 
      include 'qqgg.f'
      
      double precision pa(4),pb(4),pr(4),pborn(mxpart,4) 
      common/pborn/pborn
      double precision Dsub(4),dot_4v,sij
      logical fbps_incldip 
      common/fbps_incldip/fbps_incldip
      integer ib 
      integer j,k 
      logical cuts 
      common/cuts/cuts
      logical passed
      double precision x,z,omz,omx

      fbps_incldip=.true. 

!==== if this point fails the cuts then no dipoles will be needed (real will be finite) 
      if(cuts) call docuts(0,pborn,passed) 
      if(passed.eqv..false.) then 
         fbps_incldip=.false.       
         Dsub=0d0 
         return 
      endif
      
!=====- Now ib will tell us who branched and who is the spectator, determine zi and zj 
      
      if(ib.eq.1) then 
!==== b is the brancher
         z=-dot_4v(pa,pr)/(-dot_4v(pb,pa)+dot_4v(pa,pr))
         sij=2d0*dot_4v(pb,pr)
      elseif(ib.eq.2) then 
         z=-dot_4v(pb,pr)/(dot_4v(pb,pr)-dot_4v(pa,pb))
         sij=2d0*dot_4v(pa,pr)
      endif
      

      omx=one-x
      omz=one-z 

      call storeptilde(1,pborn) 
      
 !     write(6,*) 'in dips',sij 
!      Dsub(qq)=-gsq/x/sij*(two/omx-one-x)

      write(6,*) z,x
      Dsub(qq)=-gsq/sij*(two/(one-z*omx)-one-z)
!      write(6,*) 'dusb', Dsub(qq)
!      Dsub(qq)=-gsq/x/sij*(two/(omz+omx)-one-z)
      write(6,*) 'dusb', Dsub(qq)
      Dsub(gq)=-gsq/sij*(one+z)
      return 
 
      end 

      double precision function dot_4v(pa,pb) 
      implicit none 
      double precision pa(4),pb(4) 

      dot_4v=0d0 
      dot_4v=pa(4)*pb(4)-pa(3)*pb(3)-pa(2)*pb(2)-pa(1)*pb(1)
      return 
      end 
