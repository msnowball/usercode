
!----- pt_lx = pt minmum of x hardest lepton (up to four) 
      double precision pt_l1,pt_l2,pt_l3,pt_l4
      double precision mll_1,mll_2 
      common/leptcuts/pt_l1,pt_l2,pt_l3,pt_l4,mll_1,mll_2
!------ gampt_x is minmum pt of x hardest photon (up to 2) 
      double precision gampt_1,gampt_2,mgamgam
      common/gamcuts/gampt_1,gampt_2,mgamgam 
!------ pt_jx = pt minimum of x hardest jet
      double precision pt_j1,pt_j2,pt_veto 
      common/jetcuts/pt_j1,pt_j2,pt_veto
