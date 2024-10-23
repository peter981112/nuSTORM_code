#include "stdlib.h"

#include "TMatrixD.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>
#include <string>

#include "/home/wonjong/GiBUU_2024/NuGenTKI/include/HistIO.h"
#include "/home/wonjong/GiBUU_2024/NuGenTKI/include/GeneratorUtils.h"
#include "/home/wonjong/GiBUU_2024/NuGenTKI/style/style.h"

using namespace HistIO;
using namespace std;


void test_tree_rewrite()
{
   /*
  TString fn_lst = outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"rewighting outana9 - case with no pi 0"<<endl;

*/
  
  TString fn = "outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{2};
  int noff{1};
  auto anatag{"outana7"};
  cout<<"rewighting outana7 - case with pi 0"<<endl;

  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot_outana7.root","recreate");

  TTree * theader = (TTree*) f->Get("header");
  TTree *t1 = (TTree*) f->Get("tree");
  TTree *new_t1, *new_theader;
  new_t1 = t1->CloneTree(0); 
  new_theader = theader->CloneTree(0);
  
  Double_t perweight;
  Int_t nrun;
  //Int_t prod,evtMode;

  t1->SetBranchAddress("perweight",&perweight);
  theader->SetBranchAddress("nrun",&nrun);

  Int_t nentries = (Int_t)t1->GetEntries();
  cout<<"number of events : "<<nentries<<endl;
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i); 
    if(i%200000==0){cout<<i<<"th weight : "<<perweight<<endl;}
    perweight = 1.0; // Make changes here 
    new_t1->Fill();
  }
  theader->GetEntry(0);
  new_theader->Fill();
  if (t1->GetEntries() == new_t1->GetEntries()) {
      t1->Delete();
      theader->Delete();
  };

  new_t1->Write();
  new_theader->Write();
  fout->Write();

  return;

}