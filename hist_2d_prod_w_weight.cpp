#include "stdlib.h"

#include "TMatrixD.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH3D.h"
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

using namespace HistIO;
using namespace std;

void hist_2d_prod_w_weight()
{
  TString fn{"outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  int anaid{1};
  cout<<"plotting outana9 - case with no pi 0"<<endl;
  //TString fn{"outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  //int anaid{2};
  //cout<<"plotting outana7 - case with pi 0"<<endl;
  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot.root","recreate");
  TTree *t1 = (TTree*) f->Get("tree");
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("plot.pdf[");

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*1e3*13;//Scaling to per nucleon in CH 13=1+12
  
  TList *lout=new TList;
  TLegend *leg = new TLegend(0.7, 0.5, 0.9,0.9);
  leg->SetFillStyle(0);
  //leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.026);
  gStyle-> SetOptStat(0); // turn off stat box
  gStyle->SetPadRightMargin(0.2);
  gStyle-> SetPalette(kRainbow);
  Double_t xBj, Q2, Wrest;
  Int_t prod;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("Wrest",&Wrest);
  t1->SetBranchAddress("prod",&prod);
  TH2D *hxQ2_QE = new TH2D("hxQ2_QE","distribution of x_{Bj} and Q^{2} mode : QE",80,0.1,2,80,0.4,10);
  TH1D *hxQ2_QE_x = new TH1D("hxQ2_QE_x","x_{Bj} mode : QE",80,0,2);
  TH1D *hxQ2_QE_Q2 = new TH1D("hxQ2_QE_Q2","Q^{2} mode : QE",80,0.1,5);
  TH1D *hxQ2_QE_W = new TH1D("hxQ2_QE_W","W mode : QE",80,0,4);
  TH2D *hxQ2_RES = new TH2D("hxQ2_RES","distribution of x_{Bj} and Q^{2} mode : RES",80,0.1,1.2,80,0.1,10);
  TH1D *hxQ2_RES_x = new TH1D("hxQ2_RES_x","x_{Bj} mode : RES",80,0,2);
  TH1D *hxQ2_RES_Q2 = new TH1D("hxQ2_RES_Q2","Q^{2} mode : RES",80,0.1,5);
  TH1D *hxQ2_RES_W = new TH1D("hxQ2_RES_W","W mode : RES",80,0,4);
  TH2D *hxQ2_DIS = new TH2D("hxQ2_DIS","distribution of x_{Bj} and Q^{2} mode : DIS",80,0.1,1.2,80,0.4,10);
  TH1D *hxQ2_DIS_x = new TH1D("hxQ2_DIS_x","x_{Bj} mode : DIS",80,0,2);
  TH1D *hxQ2_DIS_Q2 = new TH1D("hxQ2_DIS_Q2","^{2} mode : DIS",80,0.1,5);
  TH1D *hxQ2_DIS_W = new TH1D("hxQ2_DIS_W","W mode : DIS",80,0,4);
  TH2D *hxQ2_2p2h = new TH2D("hxQ2_2p2h","distribution of x_{Bj} and Q^{2} mode : 2p2h",80,0.1,2,80,0.5,10);
  TH1D *hxQ2_2p2h_x = new TH1D("hxQ2_2p2h_x","x_{Bj} mode : 2p2h",80,0,2);
  TH1D *hxQ2_2p2h_Q2 = new TH1D("hxQ2_2p2h_Q2","Q^{2} mode : 2p2h",80,0.1,5);
  TH1D *hxQ2_2p2h_W = new TH1D("hxQ2_2p2h_W","W mode : 2p2h",80,0,4);
  
  TString tit[]={"evtMode QE x_{Bj}-Q^{2}", "QE x_{Bj}", "QE Q^{2}(GeV^{2})", "QE W(GeV)", "evtMode RES x_{Bj}-Q^{2}", "RES x_{Bj}", "RES Q^{2}(GeV^{2})", "RES W(GeV)", "evtMode DIS x_{Bj}-Q^{2}", "DIS x_{Bj}", "DIS Q^{2}(GeV^{2})", "DIS W(GeV)", "evtMode 2p2h x_{Bj}-Q^{2}","2p2h x_{Bj}", "2p2h Q^{2}(GeV^{2})", "2p2h W(GeV)"};
  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange};
  TH2D *hh2d[]={hxQ2_QE, hxQ2_RES, hxQ2_DIS, hxQ2_2p2h}; //list of 2d hist
  TH1D *hh1d[]={hxQ2_QE_x, hxQ2_QE_Q2, hxQ2_QE_W, hxQ2_RES_x, hxQ2_RES_Q2, hxQ2_RES_W, hxQ2_DIS_x, hxQ2_DIS_Q2, hxQ2_DIS_W, 
  	hxQ2_2p2h_x, hxQ2_2p2h_Q2, hxQ2_2p2h_W}; //list of 1d hist
  const Int_t aln_by_var[]={0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15};//list for index of hh1d in order of QE x,RES x...Q2,Q2,..,W,W
  
  THStack *stk_x = new THStack;
  TLegend *leg_x = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_x->SetFillStyle(0);
  
  THStack *stk_Q2 = new THStack;
  TLegend *leg_Q2 = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_Q2->SetFillStyle(0);
  
  THStack *stk_W = new THStack;
  TLegend *leg_W = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_W->SetFillStyle(0);
  
  Int_t nentries = (Int_t)t1->GetEntries();
  const int nh1d = sizeof(hh1d)/sizeof(TH1D*);
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    if(prod == 1) //fil histogram bins for quasi elastic case
    {
      hxQ2_QE->Fill(xBj,Q2);
      hxQ2_QE_x->Fill(xBj);
      hxQ2_QE_Q2->Fill(Q2);
      hxQ2_QE_W->Fill(Wrest);
    }
    else if(prod>1 && prod<34) //fil histogram bins for resonance case
    {
      hxQ2_RES->Fill(xBj,Q2);
      hxQ2_RES_x->Fill(xBj);
      hxQ2_RES_Q2->Fill(Q2);
      hxQ2_RES_W->Fill(Wrest);
    }
    else if(prod == 34) //fil histogram bins for deep inelastic case
    {
      hxQ2_DIS->Fill(xBj,Q2);
      hxQ2_DIS_x->Fill(xBj);
      hxQ2_DIS_Q2->Fill(Q2);
      hxQ2_DIS_W->Fill(Wrest);
    }
    else if(prod == 35 || prod == 36) //fil histogram bins for 2 part 2 hole case
    {
      hxQ2_2p2h->Fill(xBj,Q2);
      hxQ2_2p2h_x->Fill(xBj);
      hxQ2_2p2h_Q2->Fill(Q2);
      hxQ2_2p2h_W->Fill(Wrest);
    }
   
  }

  hxQ2_QE->Scale(1/hxQ2_QE->GetMaximum());
  hxQ2_RES->Scale(1/hxQ2_RES->GetMaximum());
  hxQ2_DIS->Scale(1/hxQ2_DIS->GetMaximum());
  hxQ2_2p2h->Scale(1/hxQ2_2p2h->GetMaximum());
  
  hxQ2_QE_x->Scale(1/histnormFactor);
  hxQ2_QE_Q2->Scale(1/histnormFactor);
  hxQ2_QE_W->Scale(1/histnormFactor);

  hxQ2_RES_x->Scale(1/histnormFactor);
  hxQ2_RES_Q2->Scale(1/histnormFactor);
  hxQ2_RES_W->Scale(1/histnormFactor);

  hxQ2_DIS_x->Scale(1/histnormFactor);
  hxQ2_DIS_Q2->Scale(1/histnormFactor);
  hxQ2_DIS_W->Scale(1/histnormFactor);

  hxQ2_2p2h_x->Scale(1/histnormFactor);
  hxQ2_2p2h_Q2->Scale(1/histnormFactor);
  hxQ2_2p2h_W->Scale(1/histnormFactor);
  

  for(int ii=0; ii<nh2d; ii++){
    hh2d[ii]->GetXaxis()->SetTitle("x_{Bj}");
    hh2d[ii]->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

    hh2d[ii]->SetContour(1000);
    hh2d[ii]->Draw("colz");

    gStyle->SetPadRightMargin(0.2);
    //redraw legend and update bottom margin dependent on number of entries

    gPad->SetLogy(1); //set y axis to log scale
    gPad->Update();

    // add this plot to plotbook
    TString temp_tit="";
    temp_tit = "ratio of cross section for "+tit[ii*4];
    c1->Print("plot.pdf",temp_tit);
    c1->Print(Form("outana9_%s.png", temp_tit.Data()));

    leg->Clear();
  }
  
  int aa{0};
  for(int ii=0; ii<nh1d; ii++){
    if(aa%4==0){++aa;};
    hh1d[ii]->GetXaxis()->SetTitle(tit[aa]);
    hh1d[ii]->GetYaxis()->SetTitle("cross section(cm^{2})");
    
    //leg->SetFillStyle(0);
    hh1d[ii]->SetLineStyle(kSolid);
    hh1d[ii]->SetLineColor(cols[ii/3]);
    hh1d[ii]->SetFillStyle(4050);
    hh1d[ii]->SetFillColorAlpha(cols[ii/3],0.35);

    //TString temp_tit="";
    //temp_tit = tit[ii*4];
    if(ii%3==0){
      stk_x->Add(hh1d[ii]);
      leg_x->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }
    
    else if(ii%3==1){
      stk_Q2->Add(hh1d[ii]);
      leg_Q2->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }

    else{
      stk_W->Add(hh1d[ii]);
      leg_W->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }
    ++aa;
  }
  
  gPad->SetLogy(1);//set y axis to linear scale
  gPad->Update();
  
  stk_x->Draw("hist, nostack");
  stk_x->GetXaxis()->SetTitle("x_{Bj}");
  stk_x->GetYaxis()->SetTitle("#sigma per nucleon(cm^{2})");
  stk_x->SetTitle("x_{Bj} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  leg_x->Draw();
  c1->Print("plot.pdf","x_{Bj} of 4 event modes");
  c1->Print(Form("outana9_%s.png", "x_{Bj} of 4 event modes"));

  stk_Q2->Draw("hist, nostack");
  stk_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  stk_Q2->GetYaxis()->SetTitle("#sigma per nucleon(cm^{2})");
  stk_Q2->SetTitle("Q^{2} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  leg_Q2->Draw();
  c1->Print("plot.pdf","Q^{2} of 4 event modes");
  c1->Print(Form("outana9_%s.png", "Q^{2} of 4 event modes"));

  stk_W->Draw("hist, nostack");
  stk_W->GetXaxis()->SetTitle("W(GeV)");
  stk_W->GetYaxis()->SetTitle("#sigma per nucleon(cm^{2})");
  stk_W->SetTitle("W of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  leg_W->Draw();
  c1->Print("plot.pdf","W of 4 event modes");
  c1->Print(Form("outana9_%s.png", "W of 4 event modes"));


  c1->Print("plot.pdf]","W of 4 event modes");
  fout->Close();
}
