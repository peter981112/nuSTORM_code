#include "iostream"
#include "TCanvas.h"
#include "TClass.h"
#include "TGraph.h"

void BinLogY(TH2F*h)
{

   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
}

void hist_2d()
{
  TFile *f = new TFile("outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root");
  TTree *t1 = (TTree*)f->Get("tree");
  TCanvas *c1 = new TCanvas();
  gStyle -> SetPalette(kRainbow);
  Double_t xBj, Q2;
  Int_t ev;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("ev",&ev);
  TH2F *hxQ2 = new TH2F("hxQ2","distribution of x_{Bj} and Q^{2}",80,0.1,2,80,0,3.5);
  Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    hxQ2->Fill(xBj,Q2);
  }
  //BinLogY(hxQ2);
  //gPad->SetLogy(1);
  hxQ2->GetXaxis()->SetTitle("x_{Bj}");
  hxQ2->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");

  hxQ2->SetContour(1000);
  hxQ2->Draw("colz");
}
