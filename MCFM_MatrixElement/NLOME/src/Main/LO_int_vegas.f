

      double precision function LO_int_vegas(vector,wgt) 
!----- routine for returning the LO pieces, generating different x1 and x2 for each event 
      implicit none 
      include 'constants.f' 
      include 'Gen_Q.f' 
      include 'n_born.f' 
      include 'vegas_common.f'
      include 'npart.f' 
      include 'interface_settings.f' 
      include 'facscale.f' 
      double precision vector(mxdim),wgt 
      integer j,nu,ih1,ih2,k
      double precision pborn(mxpart,4)
      double precision xx(2),Qsq,Q_f(4),sqrts
      double precision fx1(-nf:nf),fx2(-nf:nf),W,flux
      double precision lo_me(-nf:nf,-nf:nf),xzz
      double precision jac
      double precision x1,x2 
      double precision fpswt,phase_weight
      double precision pa(4),pb(4)
      logical flat 
      common/flat/flat
      common/papb/pa,pb
      common/x1x2/x1,x2
      logical cuts,passed_cuts
      common/cuts/cuts
      common/lo_me/lo_me 
      common/pborn/pborn 
      common/energy/sqrts
      common/density/ih1,ih2
      logical first,logxsam
      double precision xmin,lnxmin,xjac,xmax
      common/xminmax/xmin,xmax
      data first/.true./
      save first
      integer n_p,n_f
      common/npnf/n_p,n_f
      common/logxsam/logxsam
      logical passed_rap
      LO_int_vegas=0d0 
    
     
    

!----- Then calculate Q (final state partons 
      do nu=1,4
         Q_f(nu)=0d0 
         pa(nu)=pborn(1,nu) 
         pb(nu)=pborn(2,nu)
         do j=3,n_particles            
            Q_f(nu)=Q_f(nu)-pborn(j,nu) 
         enddo 
      enddo
      Qsq=Q_f(4)**2-Q_f(3)**2-Q_f(2)**2-Q_f(1)**2

!---- Now use that (x1x2*sqrts**2)=Q**2 to set x2 given that x1 = r1
      if(first) then          
         if(xmin/xmax.lt.1d-2) logxsam=.true.
         write(6,*) '***********************************'
         write(6,*) '* logxsam  = ',logxsam 
         write(6,*) '***********************************'
         first=.false. 
      endif

      if(logxsam.eqv..false.) then          
!---- Generate xx(1) linearly 
         xx(1)=xmax*(one-vector(1))+vector(1)*xmin 
         xjac=dabs(xmin-xmax)
         xx(2)=Qsq/(xx(1)*sqrts**2) 
         x1=xx(1)
         x2=xx(2)
         if(xx(2).gt.1d0) goto 999          
      else  
!---- Generate xx(1) log     
         lnxmin=dlog(xmin/xmax)
         xx(1)=xmax*dexp(lnxmin*(one-vector(1)))
         xjac=-lnxmin*xx(1)
         xx(2)=Qsq/(xx(1)*sqrts**2) 
         x1=xx(1)
         x2=xx(2)
         
         if(xx(2).gt.0.99d0) goto 999 

      endif

!----- Calculate PDFS
      W=sqrts**2 
      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)*xjac
      fpswt=1d0
      if(flat) fpswt=phase_weight(n_final,.false.) 
      flux=flux*fpswt

    
      if(passed_rap(0,pborn).eqv..false.) goto 999 

      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
             
      xzz=0d0 
      
      

!------------------ loop over initial state flavors
       do j=-nf,nf
          do k=-nf,nf
             xzz=fx1(j)*fx2(k)*(lo_me(j,k))*flux                        
             LO_int_vegas=LO_int_vegas+xzz
          enddo
       enddo
       
      
 
      
  

 999   continue
       return 
       end 
      
      
      
