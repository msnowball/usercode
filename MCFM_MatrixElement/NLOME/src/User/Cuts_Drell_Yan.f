
!------- This subroutine implements the Drell-Yan cuts used by CMS in arXiv:1110.2682.
!------- Namely 

!------- pt_mu_h > 18 GeV , pt_mu_s > 8 GeV  80 < m_mumu 100, |eta_mu| < 2.4 

      subroutine cuts_drell_yan(p,passed) 
      implicit none 
      include 'constants.f' 
      include 'interface_settings.f'
      include 'user_cuts.f' 
      double precision p(mxpart,4)
      logical passed 
      double precision Qsq
      double precision my_pt,my_etarap,my_dot
      double precision pt_h,pt_s,m_min,m_max
      logical first 
      data first/.true./
      save first
      logical calc_real
      common/calc_real/calc_real
!==== if both leptons are too soft then no FBPS will save them 
!==== (since one will get harder and the other will get softer) 
!==== therefore we need to store the outcome of the soft test in 
!==== a logical variable calc_real
     

!----- For Z->ll n_final = 2 
!----- In general Drell-Yan cuts will be on 3,4 
      calc_real=.true. 
      passed=.true. 
   
!------ Cut parameters 
      pt_h = lept_pt_1 
      pt_s = lept_pt_2 
      m_min = mll_s 
      m_max = mll_h 

      

      
!----- Writeout 

      if(first) then 
         first=.false. 
        write(6,99) '***************** Drell Yan Cuts *****************' 
        write(6,89) '* pt_mu(hard) > ', pt_h,'pt_mu(soft) > ',pt_s,'  *'
        write(6,79) '*       ',m_min,' < m_{mu,mu} < ',m_max,'        *' 
        write(6,59) '*           |eta_mu | < ',eta_lept, '            *'
        write(6,99) '**************************************************'
      endif

 99   format(1x,a52) 
 89   format(1x,a20,f5.1,a20,f5.1,a6)
 79   format(1x,a10,f5.1,a20,f5.1,a10) 
 59   format(1x,a30,f5.1,a14)
!------------ Cuts -----------------------------------------------

      Qsq=(p(3,4)+p(4,4))**2-(p(3,3)+p(4,3))**2
     &-(p(3,2)+p(4,2))**2-(p(3,1)+p(4,1))**2
     
      if((my_pt(3,p).lt.pt_s).and.(my_pt(4,p).lt.pt_s))calc_real=.false. 
         
     
         

!------ invariant mass (inhenrit Qsq from phase space) 
      if((dsqrt(Qsq).lt.m_min).or.(dsqrt(Qsq).gt.m_max)) passed =.false. 

      if(passed.eqv..false.) return 
     
      
     
!--   Make sure both are hard enough to pass soft 
      if((my_pt(3,p).lt.pt_s).or.(my_pt(4,p).lt.pt_s)) then 
         passed = .false. 
!            write(6,*) 'soft',my_pt(3,p),my_pt(4,p) 
         return 
      endif

!--   Make sure at least one is hard enough to pass hard
      if((my_pt(3,p).lt.pt_h).and.(my_pt(4,p).lt.pt_h)) then 
         passed = .false. 
!         write(6,*) 'hard',my_pt(3,p),my_pt(4,p) 
         return 
      endif
      
         
      return 
      end subroutine 



      

