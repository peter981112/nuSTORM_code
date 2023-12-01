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

#include "/home/wjchang/Desktop/GiBUU/NuGenTKI/include/HistIO.h"
#include "/home/wjchang/Desktop/GiBUU/NuGenTKI/include/GeneratorUtils.h"

using namespace HistIO;
using namespace std;

void hist_2d_prod_w()
{
  //TString fn{"outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  TString fn{"outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot.root","recreate");
  TTree *t1 = (TTree*) f->Get("tree");
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("plot.pdf[");

  const int anaid{1};
  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid);
  
  TLegend *leg = new TLegend(0.4,0.5,1.-gPad->GetRightMargin()-0.04,1.-gPad->GetTopMargin()-0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.026);
  gStyle -> SetOptStat(0); // turn off stat box
  gStyle -> SetPalette(kRainbow);
  Double_t xBj, Q2;
  Int_t prod;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("prod",&prod);
  TH2F *hxQ2_QE = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : QE",80,0.1,2,80,0.4,10);
  TH1F *hxQ2_QE_x = new TH1F("hxQ2","distribution of x_{Bj} mode : QE",80,0.1,2);
  TH1F *hxQ2_QE_Q2 = new TH1F("hxQ2","distribution of Q^{2} mode : QE",80,0.1,5);
  TH2F *hxQ2_RES = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : RES",80,0.1,1.2,80,0.1,10);
  TH1F *hxQ2_RES_x = new TH1F("hxQ2","distribution of x_{Bj} mode : RES",80,0.1,1.2);
  TH1F *hxQ2_RES_Q2 = new TH1F("hxQ2","distribution of Q^{2} mode : RES",80,0.1,5);
  TH2F *hxQ2_DIS = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : DIS",80,0.1,1.2,80,0.4,10);
  TH1F *hxQ2_DIS_x = new TH1F("hxQ2","distribution of x_{Bj} mode : DIS",80,0.1,1.2);
  TH1F *hxQ2_DIS_Q2 = new TH1F("hxQ2","distribution of Q^{2} mode : DIS",80,0.1,5);
  TH2F *hxQ2_2p2h = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : 2p2h",80,0.1,2,80,0.5,10);
  TH1F *hxQ2_2p2h_x = new TH1F("hxQ2","distribution of x_{Bj} mode : 2p2h",80,0.1,2);
  TH1F *hxQ2_2p2h_Q2 = new TH1F("hxQ2","distribution of Q^{2} mode : 2p2h",80,0.1,5);
  Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    if(prod == 1)
    {
      hxQ2_QE->Fill(xBj,Q2);
      hxQ2_QE_x->Fill(xBj);
      hxQ2_QE_Q2->Fill(Q2);
    }
    else if(prod>1 && prod<34)
    {
      hxQ2_RES->Fill(xBj,Q2);
      hxQ2_RES_x->Fill(xBj);
      hxQ2_RES_Q2->Fill(Q2);
    }
    else if(prod == 34)
    {
      hxQ2_DIS->Fill(xBj,Q2);
      hxQ2_DIS_x->Fill(xBj);
      hxQ2_DIS_Q2->Fill(Q2);
    }
    else if(prod == 35 || prod == 36)
    {
      hxQ2_2p2h->Fill(xBj,Q2);
      hxQ2_2p2h_x->Fill(xBj);
      hxQ2_2p2h_Q2->Fill(Q2);
    }
   
  }
  
  hxQ2_QE->Scale(1/histnormFactor);
  hxQ2_QE_x->Scale(1/histnormFactor);
  hxQ2_QE_Q2->Scale(1/histnormFactor);
  hxQ2_RES->Scale(1/histnormFactor);
  hxQ2_RES_x->Scale(1/histnormFactor);
  hxQ2_RES_Q2->Scale(1/histnormFactor);
  hxQ2_DIS->Scale(1/histnormFactor);
  hxQ2_DIS_x->Scale(1/histnormFactor);
  hxQ2_DIS_Q2->Scale(1/histnormFactor);
  hxQ2_2p2h->Scale(1/histnormFactor);
  hxQ2_2p2h_x->Scale(1/histnormFactor);
  hxQ2_2p2h_Q2->Scale(1/histnormFactor);
  

  hxQ2_QE->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_QE->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_QE->SetContour(1000);
  hxQ2_QE->Draw("colz");

  //leg->AddEntry(hxQ2_QE,"evtMode 1 x_{Bj}-Q^{2}(GeV^{2})","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->SetLogy(1);
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 1 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_QE_x->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_QE_x->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_QE_x->SetLineStyle(kSolid);
  hxQ2_QE_x->SetLineColor(kRed-3);
  hxQ2_QE_x->SetFillColor(kRed-3);
  hxQ2_QE_x->SetFillStyle(1001);
  hxQ2_QE_x->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 1 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_QE_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  hxQ2_QE_Q2->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_QE_Q2->SetLineStyle(kSolid);
  hxQ2_QE_Q2->SetLineColor(kBlue);
  hxQ2_QE_Q2->SetFillColor(kBlue);
  hxQ2_QE_Q2->SetFillStyle(1001);
  gPad->SetLogy(0);
  gPad->Update();
  hxQ2_QE_Q2->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();

  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 1 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  gPad->SetLogy(1);
  gPad->Update();
  hxQ2_RES->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_RES->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_RES->SetContour(1000);
  hxQ2_RES->Draw("colz");

  //leg->AddEntry(hxQ2_RES,"evtMode 2 x_{Bj}-Q^{2}(GeV^{2})","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 2 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
    
  hxQ2_RES_x->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_RES_x->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_RES_x->SetLineStyle(kSolid);
  hxQ2_RES_x->SetLineColor(kRed-3);
  hxQ2_RES_x->SetFillColor(kRed-3);
  hxQ2_RES_x->SetFillStyle(1001);
  hxQ2_RES_x->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 2 x_{Bj}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_RES_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  hxQ2_RES_Q2->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_RES_Q2->SetLineStyle(kSolid);
  hxQ2_RES_Q2->SetLineColor(kBlue);
  hxQ2_RES_Q2->SetFillColor(kBlue);
  hxQ2_RES_Q2->SetFillStyle(1001);
  gPad->SetLogy(0);
  gPad->Update();
  hxQ2_RES_Q2->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 2 Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  gPad->SetLogy(1);
  gPad->Update();
  
  hxQ2_DIS->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_DIS->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_DIS->SetContour(1000);
  hxQ2_DIS->Draw("colz");

  //leg->AddEntry(hxQ2_DIS,"evtMode 3 x_{Bj}-Q^{2}(GeV^{2})","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 3 x_{Bj}-Q^{2}(GeV^{2})");
  //gPad->SetLogy(1);

  leg->Clear();
  
    
  hxQ2_DIS_x->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_DIS_x->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_DIS_x->SetLineStyle(kSolid);
  hxQ2_DIS_x->SetLineColor(kRed-3);
  hxQ2_DIS_x->SetFillColor(kRed-3);
  hxQ2_DIS_x->SetFillStyle(1001);
  hxQ2_DIS_x->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 3 x_{Bj}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_DIS_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  hxQ2_DIS_Q2->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_DIS_Q2->SetLineStyle(kSolid);
  hxQ2_DIS_Q2->SetLineColor(kBlue);
  hxQ2_DIS_Q2->SetFillColor(kBlue);
  hxQ2_DIS_Q2->SetFillStyle(1001);
  gPad->SetLogy(0);
  gPad->Update();
  hxQ2_DIS_Q2->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 3 Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  gPad->SetLogy(1);
  gPad->Update();

  hxQ2_2p2h->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_2p2h->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_2p2h->SetContour(1000);
  hxQ2_2p2h->Draw("colz");

  //leg->AddEntry(hxQ2_2p2h,"evtMode 4 x_{Bj}-Q^{2}(GeV^{2})","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 4 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);
  leg->Clear();
   
  hxQ2_2p2h_x->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_2p2h_x->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_2p2h_x->SetLineStyle(kSolid);
  hxQ2_2p2h_x->SetLineColor(kRed-3);
  hxQ2_2p2h_x->SetFillColor(kRed-3);
  hxQ2_2p2h_x->SetFillStyle(1001);
  hxQ2_2p2h_x->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 4 x_{Bj}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_2p2h_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  hxQ2_2p2h_Q2->GetYaxis()->SetTitle("cross section(cm^{2})");

  //leg->SetFillStyle(0);
  hxQ2_2p2h_Q2->SetLineStyle(kSolid);
  hxQ2_2p2h_Q2->SetLineColor(kBlue);
  hxQ2_2p2h_Q2->SetFillColor(kBlue);
  hxQ2_2p2h_Q2->SetFillStyle(1001);
  
  gPad->SetLogy(0);
  gPad->Update();
  hxQ2_2p2h_Q2->Draw();
    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 4 Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  c1->Print("plot.pdf]","Title: evtMode 4 Q^{2}(GeV^{2})");
  fout->Close();
}
