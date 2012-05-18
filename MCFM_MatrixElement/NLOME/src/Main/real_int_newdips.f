
       double precision function real_int_newdips(vector,wgt)
!----- routine for returning the real pieces 
   
      implicit none 
      include 'constants.f' 
      include 'facscale.f'
      include 'interface_settings.f' 
      include 'scale.f' 
      include 'nlooprun.f' 
      include 'qcdcouple.f'
      include 'ptilde.f'
      include 'vegas_common.f'
      include 'user_cuts.f' 
      include 'qqgg.f'
      double precision vector(mxdim),wgt
      double precision sab_born
      
      
     
      double precision p(mxpart,4),ptemp(mxpart,4),pborn(mxpart,4)
      double precision pa(4),pb(4),pr(4),r(3),wgps
      double precision msq_r(-nf:nf,-nf:nf)
      double precision msq_c(-nf:nf,-nf:nf)
      double precision xmsq_zz(0:maxd),xmsq_tzz   
      double precision xmsq_czz
      double precision pat(4),pbt(4)
      double precision sal(0:maxd),sbl(0:maxd)
      double precision sal2(0:maxd),sbl2(0:maxd)
      integer j,k,ih1,ih2,nd,i,flag,ndmax_2
      logical FBPS,inc_real,passed_rap
      logical first 
      data first /.true. /
      save first 
      logical cuts,passed_cuts,includedipole
      common/cuts/cuts
      double precision sqrts,W,flux,my_pt
      double precision sab,s(mxpart,mxpart)
      double precision xx(2),fx1(-nf:nf),fx2(-nf:nf) 
      double precision xxt(2),fx1t(-nf:nf),fx2t(-nf:nf),fluxt
      double precision xa,xb,pt_fb,eta_fb,q(mxpart,4)
      double precision test_p(n_final,4),rte(n_final,4)
      double precision po(mxpart,4),my_etarap,my_dot
      double precision lnxmin,xjac,x1,x2,saj,sbj
      double precision xmin(0:maxd),xmax(0:maxd),xmin_st,xmax_st
      logical incldip(0:maxd)
      common/incldip/incldip
      common/density/ih1,ih2
      common/energy/sqrts
      common/branched_p/pt_fb,eta_fb,xa,xb
      common/x1x2/x1,x2
      common/papb/pa,pb
      common/pborn/pborn
      integer n_fail
      common/n_fail/n_fail 
      double precision Dsub(4)
      common/Dsub/Dsub 
      logical fbps_incldip 
      common/fbps_incldip/fbps_incldip
      double precision dot 

!      do i=1,4 
!         pa(i)=pborn(1,i)
!         pb(i)=pborn(2,i)
 !     enddo
 !     call writeout(pborn) 
 !     write(6,*) my_pt(3,pborn),my_pt(4,pborn) 
   

!----------PHASE SPACE -------------------------------------

!-------BORN PHASE SPACE ------------------------------------
      
      do i=-nf,nf 
         fx1(i)=0d0 
         fx2(i)=0d0 
         
         do j=-nf,nf 
            msq_r(i,j)=0d0 
         
               msq_c(i,j)=0d0 
        
         enddo
      enddo
      xmin_st=1d0 
      xmax_st=0d0 
      do i=1,4 
         pa(i)=0d0 
         pb(i)=0d0 
         pr(i)=0d0 
      enddo
      sab=0d0 
!---------- GENERATE REAL PHASE SPACE ---------------------------
      
      wgps=0d0 
      real_int_newdips=0d0 
      sab_born=0d0
      do i=1,3 
         sab_born=sab_born-(pborn(1,i)+pborn(2,i))**2
      enddo
      sab_born=sab_born+(pborn(1,4)+pborn(2,4))**2
      inc_real=FBPS(-pborn(1,4),-pborn(2,4),vector(1:2),pa,pb,pr,wgps)
      if(inc_real.eqv..false.) goto 999
      wgps=wgps*xa*xb*sqrts**2/(sab_born)

      call Gen_mcfm_p(pa,pb,pr,p) 
      call dotem(n_final+3,p,s)
      call smalls(s,n_final+1,*999)

      incldip(0)=.true. 
      if(cuts) call docuts(0,p,incldip(0))
!----- Calculate RealME given this phase space point 
      call RealME_newdips(p,msq_r,msq_c,incldip(0)) 
      ndmax=1 
      incldip(1)=fbps_incldip
     
!----- Generate PDFS 
      sab=0d0 

      do i=1,3 
         sab=sab-(pa(i)+pb(i))**2
      enddo
      sab=sab+(pa(4)+pb(4))**2 

!----- Determine xmin and xmax
      do nd=0,ndmax
          if(nd.eq.0) then 
            do j=1,mxpart 
                do k=1,4 
                   q(j,k)=p(j,k) 
                enddo
             enddo             
          else
             do j=1,mxpart 
                do k=1,4 
                   q(j,k)=ptilde(nd,j,k) 
                enddo
             enddo
          endif
          do k=1,4 
             pa(k)=q(1,k) 
             pb(k)=q(2,k) 
          enddo
          if(incldip(nd)) then 
           Call determine_xminmax(xmin(nd),xmax(nd),sab/sqrts**2,q,*998)
             if(xmin(nd).lt.xmin_st) xmin_st=xmin(nd)
             if(xmax(nd).gt.xmax_st) xmax_st=xmax(nd)
             goto 112
 998         continue
             incldip(nd)=.false.
 112         continue 
          endif
       enddo
     
!----- if all failed then exit 
       do nd=0,ndmax 
          if(incldip(nd)) goto 113 
       enddo
       goto 999 

 113   continue 
     
!      if(xmin_st/xmax_st.gt.1d-2) then 
         
!         xx(1)=xmin_st*(one-vector(3))+xmax_st*vector(3) 
!         xjac=xmax_st-xmin_st
!         xx(2)=sab/(xx(1)*sqrts**2) 
!         if(xx(2).gt.0.99d0) goto 999 
         
         
!      else

!---- Generate xx(1) log
      xmin_st=sab/sqrts**2
!      xmin_st=xmin(0)
      lnxmin=dlog(xmin_st)
      xx(1)=dexp(lnxmin*(one-vector(3)))
      xjac=-lnxmin*xx(1)
      xx(2)=sab/(xx(1)*sqrts**2) 
      if(xx(1).gt.0.99d0) goto 999
      if(xx(2).gt.0.99d0) goto 999 
      xxt(1)=xx(1) 
      xxt(2)=sab_born/(xxt(1)*sqrts**2) 

!      endif

      W=sqrts**2
      
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W**2)*wgps*xjac
      flux=flux/xx(1)      
      x1=xx(1)
      x2=xx(2)

      

!----- jet veto / require jet 
      if(jet_veto) then 
         if(first) then 
            first=.false.
            write(6,*) '*****************************************' 
            write(6,121) '* Vetoing jets with pt > ',jet_pt_max,' *'
            write(6,*) '*****************************************'
         endif
        if(pt_fb.gt.jet_pt_max) incldip(0)=.false.       
      endif
      
     
     
  
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
      
     
      xmsq_zz(0)=0d0 
!------Combine events and counter events with pdfs 
       do j=-nf,nf
          do k=-nf,nf
             xmsq_tzz=fx1(j)*fx2(k)*msq_r(j,k)
             xmsq_zz(0)=xmsq_zz(0)+xmsq_tzz
          enddo
       enddo

    
       do nd=1,ndmax 
          xmsq_zz(nd)=0d0 
          do j=-nf,nf
             do k=-nf,nf
!                if((j.ne.0).and.(k.eq.0)) then 
!!!                   xmsq_czz=fx1(j)*fx2(k)*(-2d0*tr*Dsub(gq)*msq_c(-j,j))
!!                elseif((j.eq.0).and.(k.ne.0)) then                   !
!                   xmsq_czz=fx1(j)*fx2(k)*(-2d0*tr*Dsub(gq)*msq_c(-k,k))
!                xmsq_zz(nd)=xmsq_zz(nd)+xmsq_czz
             enddo
          enddo
       enddo
       

       call writeout(p) 
       call writeout(pborn) 
       write(6,*) dot(p,1,5),dot(p,2,5) 
       write(6,*) msq_r(0,1),msq_r(-1,0) 
       write(6,*)  2d0*cf*Dsub(qq)*msq_c(-1,1) 
       write(6,*)   +4d0*tr*Dsub(gq)*msq_c(-1,1)
    

       if((dabs(dot(p,1,5)).lt.5d-4).or.(-dot(p,2,5).lt.5d-4)) pause

       real_int_newdips=0d0
     
     
       do nd=0,ndmax
          if(nd.eq.0) then 
             do j=1,mxpart 
                do k=1,4 
                   q(j,k)=p(j,k) 
                enddo
             enddo             
          else
             do j=1,mxpart 
                do k=1,4 
                   q(j,k)=ptilde(nd,j,k) 
                enddo
             enddo
          endif
          do k=1,4 
             pa(k)=pat(k) 
             pb(k)=pbt(k) 
          enddo
!--------
          
        
!          if(incldip(nd)) incldip(nd)=passed_rap(q) 
!          do k=1,4
!             pa(k)=q(1,k) 
!             pb(k)=q(2,k) 
!          enddo

!        call writeout(q)
!          write(6,*) nd,my_pt(3,q),my_pt(4,q) 
!---------if the event passes the cuts we include it
!          write(6,*) nd,includedipole(nd,q),incldip(nd)
!          if((nd.gt.0).and.incldip(nd).and.includedipole(nd,q)) pause
          if(incldip(nd)) then      
        
             real_int_newdips=real_int_newdips+xmsq_zz(nd)*flux

          else 
             goto 997
          endif
          
 997      continue
       enddo
       
       return 
   
 999   continue 
!----- nb n_fail will be zeroed outside of this routine, hence 
!----- no initializtion here. 
       n_fail=n_fail+1 
       return 
  
       
 121  format(1x,a30,f8.2,a4) 
   
       end 


      logical function passed_rap(p) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'user_cuts.f' 
      double precision p(mxpart,4),my_etarap 
      character*2 plabel
      common/plabel/plabel(mxpart)
      integer i 
      double precision eta_mx(3:n_final+2) 

      passed_rap=.true. 

      do i=3,n_final+2 
         if((plabel(i).eq.'ea').or.(plabel(i).eq.'el').or.
     &        (plabel(i).eq.'ma').or.(plabel(i).eq.'ml')) then 
            eta_mx(i)=eta_lept 
          
         elseif(plabel(i).eq.'ga') then 
            eta_mx(i)=eta_phot 
         else
            eta_mx(i)=1000d0 
         endif
      enddo

     

      do i =3,n_final+2
         if(dabs(my_etarap(i,p)).gt.eta_mx(i)) then 
           
            passed_rap=.false. 
            return 
         endif
      enddo

      return 
      end 
        
      
