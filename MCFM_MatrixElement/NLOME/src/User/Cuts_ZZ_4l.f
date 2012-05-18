
!------- This subroutine implements cuts on ZZ-> 4l 

!------- pt_mu_h > 18 GeV , pt_mu_s > 8 GeV  80 < m_mumu 100, |eta_mu| < 2.4 

      subroutine cuts_ZZ_4l(p,passed) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f'
      include 'user_cuts.f' 
      double precision p(mxpart,4)
      logical passed 
      double precision Qsq
      integer n_s1,n_s2,n_h1,n_h2,j
      double precision Qsq_34,Qsq_56
      double precision my_pt,my_etarap,my_dot
      double precision pt_h1,pt_h2,pt_s1,pt_s2
      double precision mll_min,mll_max,eta_max 
      
      logical first 
      data first/.true./
      save first
      
     

!----- For Z->4l n_final = 4 
!----- In general Drell-Yan cuts will be on 3,4 

      passed=.true. 

      
!------ Cut parameters 
!------pt_h1 = hardest lepton 
!----- pt_h2 = second hardest lepton 
!----- pt_s2 = second softest lepton 
!----- pt_s1 = softest lepton 

      pt_h1 = lept_pt_1 
      pt_h2 =  lept_pt_2
      pt_s2 =  lept_pt_3
      pt_s1 =  lept_pt_4
   
      eta_max = eta_lept

      

!----- Writeout 

      if(first) then 
         first=.false. 
        write(6,99) '*****************  ZZ => 4l Cuts *****************' 
        write(6,59) '*          pt_mu(hard) > ', pt_h1, '             *'
        write(6,59) '*       pt_mu(2nd hard) > ',pt_h2,'              *'
        write(6,59) '*          pt_mu(soft) > ', pt_s1, '             *'
        write(6,59) '*       pt_mu(2nd soft) > ',pt_s2,'              *'
!        write(6,79) '*       ',mll_min,' < m_{mu,mu} < ',mll_max,'    *' 
!        write(6,59) '*           |eta_mu | < ',eta_max, '             *'
        write(6,99) '**************************************************'
      endif

 99   format(1x,a52) 
 89   format(1x,a20,f5.1,a20,f5.1,a6)
 79   format(1x,a10,f5.1,a20,f5.1,a10) 
 59   format(1x,a30,f5.1,a14)


      n_s1=0 
      n_s2=0 
      n_h1=0
      n_h2=0

!------ count numbers greater than each cut 
      do j=3,6 
         if(my_pt(j,p).gt.pt_s1) then 
            n_s1=n_s1+1 
         endif
         if(my_pt(j,p).gt.pt_s2) then 
            n_s2=n_s2+1 
         endif
         if(my_pt(j,p).gt.pt_h2) then 
            n_h2=n_h2+1 
         endif
         if(my_pt(j,p).gt.pt_h1) then 
            n_h1=n_h1+1 
         endif
      enddo
         
!---- Cuts 
      if(n_s1.lt.4) passed = .false. 
      if(n_s2.lt.3) passed = .false. 
      if(n_h2.lt.2) passed = .false. 
      if(n_h1.lt.1) passed = .false. 
  
      return 
!------------ CUTS BELOW ARE NO LONGER DONE HERE
      if(passed.eqv..false.) return 
!------ eta cuts 
      if((dabs(my_etarap(3,p)).gt.eta_max)
     &     .or.(dabs(my_etarap(4,p)).gt.eta_max)) passed =.false. 
      
      if((dabs(my_etarap(5,p)).gt.eta_max)
     &     .or.(dabs(my_etarap(6,p)).gt.eta_max)) passed =.false. 
   

      return 
      !------------ Cuts -----------------------------------------------

      Qsq_34=(p(3,4)+p(4,4))**2-(p(3,3)+p(4,3))**2
     &-(p(3,2)+p(4,2))**2-(p(3,1)+p(4,1))**2
     
      Qsq_56=(p(5,4)+p(6,4))**2-(p(5,3)+p(6,3))**2
     &-(p(5,2)+p(6,2))**2-(p(5,1)+p(6,1))**2


!------ invariant mass (inhenrit Qsq from phase space) 
      if((dsqrt(Qsq_34).lt.mll_min).or.
     &     (dsqrt(Qsq_34).gt.mll_max)) passed =.false. 

      if(passed.eqv..false.) return 

      if((dsqrt(Qsq_56).lt.mll_min).or.
     &     (dsqrt(Qsq_56).gt.mll_max)) passed =.false. 

      if(passed.eqv..false.) return 

 
      return 
      end subroutine 



      

