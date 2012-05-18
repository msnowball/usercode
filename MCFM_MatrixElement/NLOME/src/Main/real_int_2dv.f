

      double precision function real_int_2dv(vector,wgt) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f' 
      include 'mxdim.f' 
      include 'ptilde.f' 
      include 'user_cuts.f'
      double precision vector(mxdim),wgt 
      double precision xmin(0:maxd),xmax(0:maxd) 
      double precision p(mxpart,4),pborn(mxpart,4),q(mxpart,4) 
      double precision pa(4),pb(4),pr(4),sab,sab_born 
      double precision wgps,sqrts
      double precision pt_fb,eta_fb,xa,xb
      logical inc_real,FBPS,cuts
      common/cuts/cuts
      integer i,j,k,nd
      double precision s(mxpart,mxpart)
      double precision cont_nd,my_simp_real
      logical incldip(0:maxd),includedipole
      double precision msq_r(-nf:nf,-nf:nf),msq_c(1:maxd,-nf:nf,-nf:nf)
      double precision msq_t(-nf:nf,-nf:nf)
      double precision real_int_rb 
      external real_int_rb
      logical first 
      data first /.true./
      save first
      common/incldip/incldip
      common/pborn/pborn
      common/msq_t/msq_t
      common/sab/sab
      common/papb/pa,pb
      common/energy/sqrts
      common/branched_p/pt_fb,eta_fb,xa,xb 
  
      real_int_2dv=0d0 
   

      do j=-nf,nf
         do k=-nf,nf 
            msq_r(j,k)=0d0 
            do i=1,ndmax 
               msq_c(i,j,k)=0d0 
            enddo
         enddo
      enddo

      do i=1,4 
         pa(i)=0d0 
         pb(i)=0d0 
         pr(i)=0d0 
      enddo
      inc_real=.true.

!---------- GENERATE REAL PHASE SPACE ---------------------------
      
      wgps=0d0 
      sab_born=0d0
      do i=1,3 
         sab_born=sab_born-(pborn(1,i)+pborn(2,i))**2
      enddo
      sab_born=sab_born+(pborn(1,4)+pborn(2,4))**2
     
      inc_real=FBPS(-pborn(1,4),-pborn(2,4),
     &     vector(1:2),pa,pb,pr,wgps)
      if(inc_real.eqv..false.) goto 999 
      wgps=wgps*xa*xb*sqrts**2/(sab_born)
      incldip(0)=.true. 
      call Gen_mcfm_p(pa,pb,pr,p) 

      if(cuts) call docuts(0,p,incldip(0))
  
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
    
      if(pt_fb.lt.jet_pt_max) goto 999 
        
!-----Calculate RealME given this phase space point 
      call RealME(p,msq_r,msq_c,incldip(0),.true.) 



!----- Generate PDFS 
      sab=0d0 
      do i=1,3 
         sab=sab-(pa(i)+pb(i))**2
      enddo
      sab=sab+(pa(4)+pb(4))**2 
      
!----- Loop over dipoles 
     
      do i=0,ndmax 
         xmin(i)=0d0 
         xmax(i)=0d0
         cont_nd=0d0 
         do j=1,mxpart 
            do k=1,4 
               q(j,k)=0d0 
            enddo
         enddo
         if(i.eq.0) then 
            do j=1,mxpart 
               do k=1,4 
                  q(j,k)=p(j,k) 
               enddo
             enddo             
          else
             do j=1,mxpart 
                do k=1,4 
                   q(j,k)=ptilde(i,j,k) 
                enddo
             enddo
        endif
               
     


!-------- check that this piece passes cuts 
          if(incldip(i).and.includedipole(i,q).and.(i.eq.0)) then         
             
             call determine_xminmax(xmin(i),xmax(i),sab/sqrts**2
     &            ,q,*997)
          
             do j=-nf,nf 
                do k=-nf,nf 
                   msq_t(j,k)=0d0 
                   if(i.eq.0) then 
                      msq_t(j,k)=msq_r(j,k)                      
                   else
                      msq_t(j,k)=msq_c(i,j,k)
                   endif
                enddo
             enddo
             cont_nd=0d0 
             call dromb(real_int_rb,xmin(i),xmax(i),cont_nd)
             if(i.ne.0) cont_nd=-cont_nd            
             real_int_2dv=real_int_2dv+cont_nd*wgps 
          endif
 997      continue 
       enddo
       
   
    

      return 
 999  continue
 121  format(1x,a30,f8.2,a4)  
      return 
      end 

        double precision function real_int_rb(x)
      implicit none 
      include 'constants.f' 
      include 'facscale.f'
      double precision x,msq(-nf:nf,-nf:nf),fx1(-nf:nf),fx2(-nf:nf) 
      double precision Qsq,sqrts,xjac,xmsq 
      logical log_mode_x
      common/log_mode_x/log_mode_x
      common/sab/Qsq
      common/energy/sqrts 
      common/msq_t/msq
      integer i,j,k,ih1,ih2
      common/density/ih1,ih2
      double precision xx(2),flux,t


     
      real_int_rb=0d0 

      
      xx(1)=x
      
     
      do i=-nf,nf 
         fx1(i)=0d0 
         fx2(i)=0d0 
      enddo 

   
      xx(2)=Qsq/(xx(1)*sqrts**2) 
      if(xx(2).gt.0.999d0) return 
      if(xx(1).gt.0.999d0) return 
      
    
      flux=fbGeV2/(2d0*Qsq*xx(1)*sqrts**2)
      
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
      
      xmsq=0d0 

      do j=-nf,nf
         do k=-nf,nf
            xmsq=fx1(j)*fx2(k)*msq(j,k)
            real_int_rb=real_int_rb+xmsq*flux
         enddo
      enddo
      
      

      return 
      end
