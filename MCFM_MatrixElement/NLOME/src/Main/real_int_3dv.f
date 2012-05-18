
       double precision function real_int_3dv(vector,wgt)
!----- routine for returning the real pieces 
   !----- 3dim vegas integration 
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
      double precision vector(mxdim),wgt
      double precision sab_born
      
      
     
      double precision p(mxpart,4),ptemp(mxpart,4),pborn(mxpart,4)
      double precision pa(4),pb(4),pr(4),r(3),wgps
      double precision msq_r(-nf:nf,-nf:nf)
      double precision msq_c(maxd,-nf:nf,-nf:nf)
      double precision xmsq_zz(0:maxd),xmsq_tzz   
      double precision xmsq_czz
      double precision sal(0:maxd),sbl(0:maxd)
      double precision sal2(0:maxd),sbl2(0:maxd)
      integer j,k,ih1,ih2,nd,i,flag,ndmax_2
      logical FBPS,inc_real
      logical first 
      data first /.true. /
      save first 
      logical cuts,passed_cuts,includedipole
      common/cuts/cuts
      double precision sqrts,W,flux
      double precision sab,s(mxpart,mxpart)
      double precision xx(2),fx1(-nf:nf),fx2(-nf:nf) 
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
      
    
   
    
      

!----------PHASE SPACE -------------------------------------

!-------BORN PHASE SPACE ------------------------------------
      
      do i=-nf,nf 
         fx1(i)=0d0 
         fx2(i)=0d0 
         
         do j=-nf,nf 
            msq_r(i,j)=0d0 
            do nd=1,ndmax 
               msq_c(nd,i,j)=0d0 
            enddo
         enddo
      enddo
      
      do nd=1,maxd 
         xmin(nd)=1d0 
         xmax(nd)=0d0 
      enddo
      xmin_st=1d0
      xmax_st=0d0 
      
   
!---------- GENERATE REAL PHASE SPACE ---------------------------
      
      wgps=0d0 
      real_int_3dv=0d0 
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
  
!----- Generate PDFS 
      sab=0d0 
      do i=1,3 
         sab=sab-(pa(i)+pb(i))**2
      enddo
      sab=sab+(pa(4)+pb(4))**2 
      if(cuts) call docuts(0,p,incldip(0))
!-----jet veto / require jet 
      if(jet_veto) then 
         if(first) then 
            first=.false.
            write(6,*) '*****************************************' 
            write(6,121) '* Vetoing jets with pt > ',jet_pt_max,' *'
            write(6,*) '*****************************************'
         endif
         if(pt_fb.gt.jet_pt_max) incldip(0)=.false.       
      endif
      

     

!-----Calculate RealME given this phase space point 
      call RealME(p,msq_r,msq_c,incldip(0)) 

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
     
          call determine_xminmax(xmin(nd),xmax(nd),sab/sqrts**2,q,*998)
          if(xmin(nd).lt.xmin_st) xmin_st=xmin(nd)
          if(xmax(nd).gt.xmax_st) xmax_st=xmax(nd)
          goto 112
 998      continue
          incldip(nd)=.false.
 112      continue 
       enddo
      
  
      if(xmin_st/xmax_st.gt.1d-2) then 

!---- Generate xx(1) linearly 
      xx(1)=xmax_st*(one-vector(3))+vector(3)*xmin_st 
      xjac=dabs(xmin_st-xmax_st)     
      xx(2)=sab/(xx(1)*sqrts**2) 
      if(xx(2).gt.0.99d0) goto 999 
      x1=xx(1)
      x2=xx(2)
     
      else

!---- Generate xx(1) log     
      lnxmin=dlog(xmin_st/xmax_st)
      xx(1)=xmax_st*dexp(lnxmin*(one-vector(3)))
      xjac=-lnxmin*xx(1)*xmax_st
      xx(2)=sab/(xx(1)*sqrts**2) 
      x1=xx(1)
      x2=xx(2)
     
      if(xx(2).gt.0.99d0) goto 999 

      endif
    
!      write(6,*) xmin(0:2) 
!      write(6,*) xmax(0:2)
!      write(6,*) xmin_st,xmax_st 
!      write(6,*) x1,x2
!      pause


       W=sqrts**2
       
       flux=fbGeV2/(2d0*xx(1)*xx(2)*W**2)*wgps*xjac
       flux=flux/xx(1)      


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
                xmsq_czz=fx1(j)*fx2(k)*(-msq_c(nd,j,k))
                xmsq_zz(nd)=xmsq_zz(nd)+xmsq_czz
             enddo
          enddo
       enddo
       



       real_int_3dv=0d0
     
     
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

       

!          call writeout(q) 
!          write(6,*) nd 
!          pause 
          if((x1.gt.xmax(nd)).or.(x1.lt.xmin(nd))) then 
             incldip(nd)=.false. 
          endif

          if(incldip(nd).and.includedipole(nd,q)) then      
!---include 
             real_int_3dv=real_int_3dv+xmsq_zz(nd)*flux

          else 
             goto 997
          endif
          
 997      continue
       enddo
       


 999   continue 
    
       return 
       
 121  format(1x,a30,f8.2,a4) 
   
       end 
