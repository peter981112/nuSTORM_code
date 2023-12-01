#include "iostream"
#include "TCanvas.h"
#include "TClass.h"
#include "TGraph.h"

void hist_2d_ch()
{
  //TFile *f = new TFile("outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root");
  TFile *f = new TFile("outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root");
  TFile *fout = new TFile("tmpplot.root","recreate");
  TTree *t1 = (TTree*)f->Get("tree");
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("plot.pdf[");

  TLegend *leg = new TLegend(0.4,0.5,1.-gPad->GetRightMargin()-0.04,1.-gPad->GetTopMargin()-0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.026);
  gStyle -> SetPalette(kRainbow);
  Double_t xBj, Q2;
  Int_t evtMode;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("evtMode",&evtMode);
  TH2F *hxQ2_1 = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : 1",80,0.1,2,80,0,3.5);
  TH2F *hxQ2_2 = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : 2",80,0.1,2,80,0,3.5);
  TH2F *hxQ2_3 = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : 3",80,0.1,2,80,0,3.5);
  TH2F *hxQ2_4 = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2} mode : 4",80,0.1,2,80,0,3.5);
  Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    if(evtMode == 1)
    {
      hxQ2_1->Fill(xBj,Q2);
    }
    else if(evtMode == 2)
    {
      hxQ2_2->Fill(xBj,Q2);
    }
    else if(evtMode == 3)
    {
      hxQ2_3->Fill(xBj,Q2);
    }
    else if(evtMode == 4)
    {
      hxQ2_4->Fill(xBj,Q2);
    }
   
  }

  hxQ2_1->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_1->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_1->SetContour(1000);
  hxQ2_1->Draw("colz");

  //leg->AddEntry(hxQ2_1,"evtMode 1 x_{Bj}-Q^{2}","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 1 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();

  hxQ2_2->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_2->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_2->SetContour(1000);
  hxQ2_2->Draw("colz");

  //leg->AddEntry(hxQ2_2,"evtMode 2 x_{Bj}-Q^{2}","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 2 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();
  
  hxQ2_3->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_3->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_3->SetContour(1000);
  hxQ2_3->Draw("colz");

  //leg->AddEntry(hxQ2_3,"evtMode 3 x_{Bj}-Q^{2}","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 3 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);

  leg->Clear();

  hxQ2_4->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2_4->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2_4->SetContour(1000);
  hxQ2_4->Draw("colz");

  //leg->AddEntry(hxQ2_4,"evtMode 4 x_{Bj}-Q^{2}","L");

    // redraw legend and update bottom margin dependent on number of entries
  leg->Draw();
  gPad->Update();
  //leg->SetY1NDC(leg.GetY2NDC()-leg.GetNRows()*leg.GetTextSize()*1.2);
  // add this plot to plotbook
  c1->Print("plot.pdf","Title: evtMode 4 x_{Bj}-Q^{2}");
  //gPad->SetLogy(1);
  leg->Clear();
  c1->Print("plot.pdf]","Title: evtMode 4 x_{Bj}-Q^{2}");
  fout->Close();
}
