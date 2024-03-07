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

void data_test_plt()
{
  TString fn1{"PhysRevD.101.092001.root"}; //getMINERvA0PI
  TString fn2{"PIZEROTKI_MINERvA.root"};  //getMINERvAPIZERO
  TFile *f1 = new TFile(fn1);
  TFile *f2 = new TFile(fn2);
  TFile *fout = new TFile("tmpplot.root","recreate");

//const TString varname[]={"muonmomentum","muontheta","protonmomentum","protontheta","dalphat","dphit","dpt", "IApN", "dpTx", "dpTy"};
  auto varname = "dpt"; 
  TList * tmplist1 = (TList*)f1->Get(varname);
  auto varname2 = "neutronmomentum"; 
  TList * tmplist2 = (TList*)f2->Get(varname2);

  TH1D * hh1  = (TH1D*) tmplist1->At(0);
  TH1D * hh2  = (TH1D*) tmplist2->At(0);   
  style::ScaleXaxis(hh2, 1E-3);
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("data_test_plot.pdf[");

  Double_t xsec;
  xsec = hh1->Integral(0,100000,"width");
  cout<<"xesc for 0 pi : "<<xsec<<" cm^{2}"<<" from integrating "<<varname<<endl;
  xsec = hh2->Integral(0,100000,"width");
  cout<<"xesc for with pi : "<<xsec<<" cm^{2}"<<" from integrating "<<varname2<<endl;
  TLegend *leg = new TLegend(0.65, 0.45, 0.85, 0.85);
  leg->SetFillStyle(0);
  //leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.023);
  gStyle-> SetOptStat(0); // turn off stat box
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.1);
  gStyle-> SetPalette(kRainbow);

  hh1->GetXaxis()->SetTitle("muonmoentum");
  hh1->GetYaxis()->SetTitle("diff xsec");
  hh1->Draw("hist");
  hh1->SetTitle("xsec for 0 pi");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg->Draw();
  c1->Print("data_test_plot.pdf","xsec for 0 pi");

  hh2->GetXaxis()->SetTitle("neutronmomentum");
  hh2->GetYaxis()->SetTitle("diff xsec");
  hh2->Draw("hist");
  hh2->SetTitle("xsec for pi");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg->Draw();
  c1->Print("data_test_plot.pdf","xsec for 0 pi");

  c1->Print("data_test_plot.pdf]","xsec for 0 pi");
  fout->Close();
}