#include <iostream>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"



extern "C"{
  
  void i2mcfm_fill_alphas_(double *my_as); 
  void i2mcfm_fill_ew_(int *tmp1,double *my_gf,double *my_ew,double *my_xw,int *tmp2);
  void i2mcfm_fill_mv_(double *my_mw,double *my_ww,double *my_mz,double *my_wz);
  void i2mcfm_fill_mq_(double mq[]);
  void i2mcfm_fill_ml_(double ml[],int *tmp3);
  void i2mcfm_fill_ckm_(double my_V[]);
  void i2mcfm_fill_higgs_(double *mH,bool *bltmp1);
  void i2mcfm_fill_sqrts_(double *sqrts);
  void i2mcfm_pdfset_(char pdfLabel[]); 
  void init_nlome_(double *mu_r,double *mu_f);
  void nlo_lik_(double p[4][4],double *mu_r,double *mu_f,double *nlowt,double *nloer);
  void swap_process_(int *proc1,int *proc2); 
  void lo_lik_(double p[4][4],double *mu_r,double *mu_f,double *nlowt,double *nloer);
  void chooser_(void);

}


int main(int argc,char *argv[])
{
  using namespace std;
  
  if(argc != 5)
    {
      cout << "Usage: ./nlome_root.exe <mH> <+/- window> <inputFile> <tree number>" << endl;
      return -1;
    }
  
  char treeName[100];
  int treeNum = atoi(argv[4]);
  sprintf(treeName,"passedEvents_%i",treeNum);
  string inputFile = argv[3];
  double window = atof(argv[2]);
  double mH = atof(argv[1]);
  cout << "Higgs mass set to " << mH << " GeV" << endl;
  cout << "Input file is " << inputFile << endl;
  cout << "Tree is " << treeName << endl;

  int proc1 = 81, proc2 = 114;
  double mu_r,mu_f;
  double pte[4][4];
  double nlowt, nloer;
  double mHLow = mH - window;
  double mHHigh = mH + window;

  TFile* sigFile = new TFile(inputFile.c_str(),"UPDATE");
  TTree* sigTree;
  if(sigFile) sigTree = (TTree*) sigFile->Get(treeName);
  if(!sigTree)
    {
      std::cout << "Cannot find tree!" << std::endl;
      return -1;
    }

  double m4l,e1,e2,e3,e4,px1,px2,px3,px4,D;
  double py1,py2,py3,py4,pz1,pz2,pz3,pz4;
  double massZ1; 
  int idL1, idL2, idL3, idL4;

  sigTree->SetBranchAddress("mass4l",&m4l);
  sigTree->SetBranchAddress("EL1",&e1);
  sigTree->SetBranchAddress("EL2",&e2);
  sigTree->SetBranchAddress("EL3",&e3);
  sigTree->SetBranchAddress("EL4",&e4);
  sigTree->SetBranchAddress("pXL1",&px1);
  sigTree->SetBranchAddress("pXL2",&px2);
  sigTree->SetBranchAddress("pXL3",&px3);
  sigTree->SetBranchAddress("pXL4",&px4);
  sigTree->SetBranchAddress("pYL1",&py1);
  sigTree->SetBranchAddress("pYL2",&py2);
  sigTree->SetBranchAddress("pYL3",&py3);
  sigTree->SetBranchAddress("pYL4",&py4);
  sigTree->SetBranchAddress("pZL1",&pz1);
  sigTree->SetBranchAddress("pZL2",&pz2);
  sigTree->SetBranchAddress("pZL3",&pz3);
  sigTree->SetBranchAddress("pZL4",&pz4);
  sigTree->SetBranchAddress("massZ1",&massZ1);
  sigTree->SetBranchAddress("idL1",&idL1);
  sigTree->SetBranchAddress("idL2",&idL2);
  sigTree->SetBranchAddress("idL3",&idL3);
  sigTree->SetBranchAddress("idL4",&idL4);

  double HZZ, HZZerr, ZZ, ZZerr;
  TBranch *newBranchHZZ = sigTree->Branch("mcfm_HZZ",&HZZ,"mcfm_HZZ/D");
  TBranch *newBranchHZZerr = sigTree->Branch("mcfm_HZZerr",&HZZerr,"mcfm_HZZerr/D");
  TBranch *newBranchZZ = sigTree->Branch("mcfm_ZZ",&ZZ,"mcfm_ZZ/D");
  TBranch *newBranchZZerr = sigTree->Branch("mcfm_ZZerr",&ZZerr,"mcfm_ZZerr/D");
  TBranch *newBranchD = sigTree->Branch("D",&D,"D/D"); 
 
  Long64_t nentries = sigTree->GetEntries();

  for(int iEvt=0; iEvt < nentries; iEvt++){

    sigTree->GetEntry(iEvt);

    nlowt = -1; nloer = -1;
    HZZ = -1; HZZerr = -1;
    ZZ = -1;  ZZerr = -1;
        
    if( m4l >= mHLow && m4l <= mHHigh  && massZ1 >50 )
      {
	
	pte[3][0] = e1;
	pte[2][0] = pz1;
	pte[1][0] = py1;
	pte[0][0] = px1;
	
	pte[3][1] = e2;
	pte[2][1] = pz2;
	pte[1][1] = py2;
	pte[0][1] = px2;
	

	pte[3][2] = e3;
	pte[2][2] = pz3;
	pte[1][2] = py3;
	pte[0][2] = px3;
        	
	pte[3][3] = e4;
	pte[2][3] = pz4;
	pte[1][3] = py4;
	pte[0][3] = px4;

        if(idL1<0 && idL2 > 0 ) 
	  {
	    pte[3][0] = e2;
	    pte[2][0] = pz2;
	    pte[1][0] = py2;
	    pte[0][0] = px2;
	    
	    pte[3][1] = e1;
	    pte[2][1] = pz1;
	    pte[1][1] = py1;
	    pte[0][1] = px1;
	    
	  }  //keep first pair order -+
  	
        if(idL3<0 && idL4 > 0 ) 
	  { 
	    pte[3][2] = e4;
	    pte[2][2] = pz4;
	    pte[1][2] = py4;
	    pte[0][2] = px4;
	    
	    pte[3][3] = e3;
	    pte[2][3] = pz3;
	    pte[1][3] = py3;
	    pte[0][3] = px3;
	    
	  } //keep second pair order -+
	
	mu_r = m4l;
	mu_f = m4l;
	
        bool bltmp1 = true; 

        init_nlome_(&mu_r,&mu_f); // <--- reads in all the .DAT files in NLOME, sets proc = Higgs

        i2mcfm_fill_higgs_(&m4l, &bltmp1); // <--- sets mH = m4l
 
	chooser_(); // <--- makes sure the settings are correct

        lo_lik_(pte,&mu_r,&mu_f,&nlowt,&nloer); // <--- gets event weight  

	HZZ = nlowt;
	HZZerr = nloer;
	nlowt = 0; nloer = 0;
	
	mu_r = 2*91.1876;
	mu_f = mu_r;
	
        i2mcfm_fill_higgs_(&m4l, &bltmp1);

	init_nlome_(&mu_r,&mu_f);

        i2mcfm_fill_higgs_(&m4l, &bltmp1); 

	swap_process_(&proc1,&proc2); // <--- change to ZZ since init_nlome sets proc to Higgs
        
	chooser_();
        
	lo_lik_(pte,&mu_r,&mu_f,&nlowt,&nloer);
        
	ZZ = nlowt;
	ZZerr = nloer;
	
        D = TMath::Log(HZZ/ZZ);   

	cout << "mH is " << m4l << "  " << HZZ << "  " << ZZ << endl;

      }

    newBranchHZZ->Fill();
    newBranchHZZerr->Fill();
    newBranchZZ->Fill();
    newBranchZZerr->Fill();
    newBranchD->Fill();
  }

   sigFile->cd();

   sigTree->Write("passedEvents");

   sigFile->Write();

   sigFile->Close(); 

  return 0;

} 
             
