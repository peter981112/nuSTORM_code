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

#include "/home/wjchang/Desktop/GiBUU/NuGenTKI/include/HistIO.h"
#include "/home/wjchang/Desktop/GiBUU/NuGenTKI/include/GeneratorUtils.h"
#include "/home/wjchang/Desktop/GiBUU/NuGenTKI/style/style.h"

void flux_plot()
{
  TString fn3{"E3SpectraMuSig560.root"};
  TString fn4{"E4SpectraMuSig566.root"}; 
  TString fn5{"E5SpectraMuSig558.root"};
  TString fn6{"E6SpectraMuSig557.root"}; 
  TString fn7{"E7SpectraMuBack561.root"};

  TFile *f3 = new TFile(fn3);
  TFile *f4 = new TFile(fn4);
  TFile *f5 = new TFile(fn5);
  TFile *f6 = new TFile(fn6);
  TFile *f7 = new TFile(fn7);

  TFile *fout = new TFile("tmpplot.root","recreate");

//const TString varname[]={"muonmomentum","muontheta","protonmomentum","protontheta","dalphat","dphit","dpt", "IApN", "dpTx", "dpTy"};
  TString varname = "Enumu Energy";
  TH1D * hh3  = (TH1D*) f3->Get(varname);
  TH1D * hh4  = (TH1D*) f4->Get(varname);   
  TH1D * hh5  = (TH1D*) f5->Get(varname);
  TH1D * hh6  = (TH1D*) f6->Get(varname);  
  TH1D * hh7  = (TH1D*) f7->Get(varname);
 
  TH1D *hh_flux[]={hh3,hh4,hh5,hh6,hh7};

  const int nh_flux = sizeof(hh_flux)/sizeof(TH1D*);

  TCanvas *c1 = new TCanvas("c","",900,700);
  c1->Print("flux_plot.pdf[");

  THStack *stk_flux = new THStack;
  TLegend *leg_flux = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_flux->SetFillStyle(0);

  leg_flux->SetTextFont(42);
  leg_flux->SetTextSize(0.023);
  gStyle-> SetOptStat(0); // turn off stat box
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.1);
  gStyle-> SetPalette(kRainbow);

  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange, kAzure+1};

  int aa{0};

  for(int ii=0; ii<nh_flux; ii++){
    
    //leg->SetFillStyle(0);
    hh_flux[ii]->Scale(1/hh_flux[ii]->Integral());
    hh_flux[ii]->SetLineStyle(kSolid);
    hh_flux[ii]->SetLineColor(cols[ii]);
    hh_flux[ii]->SetFillStyle(4050);
    hh_flux[ii]->SetFillColorAlpha(cols[ii],0.35);

    //TString temp_tit="";
    //temp_tit = tit[ii*4];
    stk_flux->Add(hh_flux[ii]);
    TString title = "#nu_{#mu} from ";
    aa = ii+3;
    title = title + aa + " GeV #pi";
    leg_flux->AddEntry(hh_flux[ii],title,"fl");
    
  }

  const TString opt_nostack="hist nostack";

  stk_flux->Draw(opt_nostack);
  stk_flux->GetXaxis()->SetTitle("E of #nu_{#mu}(GeV)");
  stk_flux->GetYaxis()->SetTitle("flux");
  stk_flux->SetTitle("#nu_{#mu} flux");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_flux->Draw();
  c1->Print("flux_plot.pdf","normalised #nu_{#mu} flux");
  //c1->Print(Form("png/%s_%s.png", anatag, "#theta of #pi of 4 event modes"));

  c1->Print("flux_plot.pdf]","nuSTORM neutrino flux for muSig");
  fout->Close();
}