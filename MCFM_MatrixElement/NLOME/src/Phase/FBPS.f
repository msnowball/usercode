      function fbps(E1,E2,r,pa,pb,pr,pswgt) 
!--- The forward branching phase space generator: 
!--- Input : the born configuration           k_a+k_b --> k_1+...+k_n
!--- Output: the bremsstrahlung configuration p_a+p_b --> k_1+...+k_n+p_r 
      implicit none
!      include 'mxdim.f'
      logical fbps
      integer i,j,k,iab,sample,tag
      double precision E1,E2,Ecm,Scm,PI,tmin,Eb2,s,maxrb,fsoft,t0
      double precision shab,Ea,Eb,phi,sab,tar,trb,pswgt,z,c,beta,CM
      double precision xa,xb,pt_fb,eta_fb,TWOPI,TWOPI3
      double precision r(2),pa(4),pb(4),pr(4),ran2
      parameter(sample=1,t0=1.0d-1,fsoft=0.01d0)
      parameter(PI=3.141592653589793d0)
      parameter(TWOPI=2.0*PI)
      parameter(TWOPI3=TWOPI*TWOPI*TWOPI)
      common/energy/CM
      common/tag/tag
      common/branched_p/pt_fb,eta_fb,xa,xb
      common/iab/iab

!--- initialize variables
      if (ran2().lt.0.5d0) then
         iab=1
         Ea=E1
         Eb=E2
         tag=1
      else
         iab=2
         Ea=E2
         Eb=E1
         tag=2
      endif
      Ecm=CM/2d0
      Scm=4d0*Ecm*Ecm
      shab=4d0*Ea*Eb
      Eb2=4d0*Eb*Eb
      fbps=.true.
      pswgt=1.0d0
!--- generate integration variables
!---
!--- sampling over 1/s_{ab} which sits in the measure
!---
      sab=Scm**r(1)*shab**(1d0-r(1))
      pswgt=pswgt*log(Scm/shab)
!---      
      tmin=sab-shab
      tmin=min(tmin,shab*(Ecm-Ea)/Ea)
      if (sample.eq.1) then
         if (ran2().lt.fsoft) then
            tmin=t0
            trb=-r(2)*tmin
            tmin=tmin/fsoft
         else
            trb=-tmin**r(2)*t0**(1d0-r(2))
            tmin=dlog(tmin/t0)*dabs(trb)
            tmin=tmin/(1d0-fsoft)
         endif
      elseif (sample.eq.0) then
         trb=-r(2)*tmin
      endif
      pswgt=pswgt*tmin
      tmin=TWOPI
      phi=twopi*ran2() 
!----- original line 
!      phi=tmin*r(3)
      pswgt=tmin*pswgt
      tar=shab-sab-trb
      pswgt=pswgt*shab/4.0/sab/TWOPI3
!--- do no include flux-factor in the phase space weight
!     pswgt=pswgt/2d0/sab
!-- generate p_a and p_b
      beta=-trb/shab
      z=(Eb2*sab+tar*trb)/Eb2/(shab-trb)
      c=(Eb2*sab-tar*trb)/(Eb2*sab+tar*trb)
      if ((c.gt.1.0).or.(c.lt.-1.0)) then 
         fbps=.false.    
         return 
      endif
      s=sqrt(1d0-c*c)
      c=-c
      if (fbps) then
         if (iab.eq.1) then
            pa(4)=-(1d0+beta)*Ea
            pa(1)=0d0
            pa(2)=0d0
            pa(3)=pa(4)
            pb(4)=-z*Eb
            pb(1)=pb(4)*s*sin(phi)
            pb(2)=pb(4)*s*cos(phi)
            pb(3)=pb(4)*c
            pr(3)=Eb-Ea
            xa=-pa(4)/Ecm
            xb=sab/Scm/xa
         else
            pa(4)=-z*Eb
            pa(1)=pa(4)*s*sin(phi)
            pa(2)=pa(4)*s*cos(phi)
            pa(3)=-pa(4)*c
            pb(4)=-(1d0+beta)*Ea
            pb(1)=0d0
            pb(2)=0d0
            pb(3)=-pb(4)
            pr(3)=Ea-Eb
            xb=-pb(4)/Ecm
            xa=sab/Scm/xb
         endif

         pr(4)=-Ea-Eb-pa(4)-pb(4)
         pr(1)=-pa(1)-pb(1)
         pr(2)=-pa(2)-pb(2)
         pr(3)=pr(3)-pa(3)-pb(3)
         if(tar*trb/sab.lt.0d0) then 
            fbps=.false. 
            return 
         endif
         pt_fb=sqrt(tar*trb/sab)
         eta_fb=-0.5*log(xb*tar/xa/trb)
      endif
!      if ((xa.gt.1d0).or.(xb.gt.1d0)) fbps=.false.
    
      return
      end

