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

using namespace HistIO;
using namespace std;

void col_normal(TH2D* &hist) //normalize 2 hist by max value for each column
{
  Int_t x_bin_num,y_bin_num;
  x_bin_num = hist->GetNbinsX();
  y_bin_num = hist->GetNbinsY();
  //cout<<"# of x bin : "<<x_bin_num<<" # of y bin : "<<y_bin_num<<endl;
  double col_max{0};
  double bin_content{0};
  for (Int_t xx=1; xx<=x_bin_num; xx++)
  {
    col_max = 0;
    for (Int_t yy=1; yy<=y_bin_num; yy++) //loop to get max bin content within column
    {
      bin_content = hist->GetBinContent(xx,yy); //get bin content for 2d histogram
      if(bin_content<0){hist->SetBinContent(xx,yy,0);} //if bin content is <0 delete it
      if(bin_content>col_max){col_max=bin_content;}
    }
    //if(xx%3==0){cout<<"col max : "<<col_max<<endl;}
    if(col_max==0){col_max=1;} //if there is nothing in the column avoide deviding by 0
    for (Int_t yy=1; yy<=y_bin_num; yy++) //loop to normalize column by the max within col
    {
      bin_content = hist->GetBinContent(xx,yy); //get bin content for 2d histogram
      hist->SetBinContent(xx,yy,bin_content/col_max);
    }
  }
  return;
}

void remove_neg(TH2D* &hist) //remove negative cross section values from 2d hist
{
  Int_t x_bin_num,y_bin_num;
  x_bin_num = hist->GetNbinsX();
  y_bin_num = hist->GetNbinsY();
  for (Int_t xx=1; xx<=x_bin_num; xx++)
  {
    for (Int_t yy=1; yy<=y_bin_num; yy++) //loop to get max bin content within column
    {
      if(hist->GetBinContent(xx,yy)<0){hist->SetBinContent(xx,yy,0);} //if bin content is <0 delete it
    }
  }
  return;
}

void wipe_hist(TH1D* &hist) //makes all bin contents to 0
{
  Int_t x_bin_num,y_bin_num;
  x_bin_num = hist->GetNbinsX();
  for (Int_t xx=1; xx<=x_bin_num; xx++){hist->SetBinContent(xx,0);}
  return;
}

void diff_xsec(TH1D* &hist, const double histnormFactor) //makes all bin contents to 0
{
  Int_t x_bin_num,y_bin_num;
  Double_t binwidth, bincontents;
  x_bin_num = hist->GetNbinsX();
  for (Int_t xx=1; xx<=x_bin_num; xx++)
  {
    binwidth = hist->GetXaxis()->GetBinWidth(xx);
    bincontents = hist->GetBinContent(xx);
    hist->SetBinContent(xx, bincontents/(histnormFactor*binwidth)); //to get diff cross section
  }
  return;
}

double GetChi2_stack(THStack *sim_stack, TH1D * data_hist, 
  const TMatrixD *rawcov, const int noff, const double dummyunit)
{
  double chi2{0};
  cout<<"stack pointer : "<<sim_stack<<endl;
  TList * hist_lst = sim_stack->GetHists();
  TAxis * axis = data_hist->GetXaxis();
  //cout<<"print axis N bin : "<<axis->GetNbins()<<" axis min : "<<
  //axis->GetXmin()<<" axis x max : "<<axis->GetXmax()<<endl;
  TH1D * sum_hist = new TH1D(*(TH1D*) hist_lst->At(0));
  wipe_hist(sum_hist);

  TIter next(hist_lst);
	TObject* object = 0;
	while ((object = next()))
  {
    sum_hist->Add((TH1D*)object,1);
  }

/*
  TCanvas *c1 = new TCanvas("c","",800,600);
  data_hist->Draw("X1 E1");
  data_hist->SetTitle("test sum hist");
  sum_hist->Draw("same, hist");
  sum_hist->GetXaxis()->SetTitle("x");
  sum_hist->GetYaxis()->SetTitle("y-#sigma");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  c1->Print(Form("png/%s_%s.png", data_hist->GetTitle(), "test sum hist"));
*/
  chi2 = style::GetChi2(data_hist, *rawcov, 
      noff, data_hist->GetNbinsX()-1, sum_hist, dummyunit);

  if(fabs(chi2)>500){chi2 = style::GetChi2(data_hist, *rawcov, 
      noff, data_hist->GetNbinsX()-2, sum_hist, dummyunit);}

  else if(fabs(chi2)>500){chi2 = style::GetChi2(data_hist, *rawcov, 
      noff, data_hist->GetNbinsX(), sum_hist, dummyunit);}

  cout<<"chi2 : "<<chi2<<" for "<<data_hist->GetTitle()<<endl;
  return chi2;
}

void hist_2d_prod_w_weight()
{
  TString fn{"outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"plotting outana9 - case with no pi 0"<<endl;
  //TString fn{"outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  //int anaid{2};
  //int noff{1};
  //cout<<"plotting outana7 - case with pi 0"<<endl;
  //auto anatag{"outana7"};

  int nodata = 1; //nodata 1 makes sure no data is polotted on the canvas

  TString fn1{"PhysRevD.101.092001.root"}; //getMINERvA0PI
  TString fn2{"PIZEROTKI_MINERvA.root"};  //getMINERvAPIZERO
  TFile *f1 = new TFile(fn1);
  TFile *f2 = new TFile(fn2);


  const double dummyunit = 1E-42; //data
  TList * tmplist1_scattermomentum = (TList*)f1->Get("muonmomentum");
  TList * tmplist1_scattertheta = (TList*)f1->Get("muontheta");
  TList * tmplist1_recoilmomentum = (TList*)f1->Get("protonmomentum");
  TList * tmplist1_recoiltheta = (TList*)f1->Get("protontheta");
  TList * tmplist1_neutronmomentum = (TList*)f1->Get("neutronmomentum");
  TList * tmplist1_dpt = (TList*)f1->Get("dpt");
  TList * tmplist1_dphit = (TList*)f1->Get("dphit");
  TList * tmplist1_dalphat = (TList*)f1->Get("dalphat");

  TList * tmplist2_neutronmomentum = (TList*)f2->Get("neutronmomentum");
  TList * tmplist2_dalphat = (TList*)f2->Get("dalphat");   
  TList * tmplist2_dpTT = (TList*)f2->Get("dpTT");

  TH1D * hh_0pi_scattermomentum  = (TH1D*) tmplist1_scattermomentum->At(0);
  TMatrixD * hh_0pi_scattermomentum_cov = (TMatrixD*) (tmplist1_scattermomentum->At(2));
  hh_0pi_scattermomentum_cov = (TMatrixD*) hh_0pi_scattermomentum_cov->Clone("cov");
  Int_t a0pi_scattermomentum_binN{hh_0pi_scattermomentum->GetNbinsX()};
  //Double_t a0pi_scattermomentum_binmin{hh_0pi_scattermomentum->GetXaxis()->GetXmin()},
     //a0pi_scattermomentum_binmax{hh_0pi_scattermomentum->GetXaxis()->GetXmax()};

  TH1D * hh_0pi_scattertheta  = (TH1D*) tmplist1_scattertheta->At(0);
  TMatrixD * hh_0pi_scattertheta_cov = (TMatrixD*) (tmplist1_scattermomentum->At(2));
  hh_0pi_scattertheta_cov = (TMatrixD*) hh_0pi_scattertheta_cov->Clone("cov");
  Int_t a0pi_scattertheta_binN{hh_0pi_scattertheta->GetNbinsX()};
  //Double_t  a0pi_scattertheta_binmin{hh_0pi_scattertheta->GetXaxis()->GetXmin()},
     //a0pi_scattertheta_binmax{hh_0pi_scattertheta->GetXaxis()->GetXmax()};

  TH1D * hh_0pi_recoilmomentum = (TH1D*) tmplist1_recoilmomentum->At(0);
  TMatrixD * hh_0pi_recoilmomentum_cov = (TMatrixD*) (tmplist1_recoilmomentum->At(2));
  hh_0pi_recoilmomentum_cov = (TMatrixD*) hh_0pi_recoilmomentum_cov->Clone("cov");
  Int_t a0pi_recoilmomentum_binN{hh_0pi_recoilmomentum->GetNbinsX()};
  //Double_t  a0pi_recoilmomentum_binmin{hh_0pi_recoilmomentum->GetXaxis()->GetXmin()},
     //a0pi_recoilmomentum_binmax{hh_0pi_recoilmomentum->GetXaxis()->GetXmax()};


  TH1D * hh_0pi_recoiltheta  = (TH1D*) tmplist1_recoiltheta->At(0);
  TMatrixD * hh_0pi_recoiltheta_cov = (TMatrixD*) (tmplist1_recoiltheta->At(2));
  hh_0pi_recoiltheta_cov = (TMatrixD*) hh_0pi_recoiltheta_cov->Clone("cov");
  Int_t a0pi_recoiltheta_binN{hh_0pi_recoiltheta->GetNbinsX()};
  //Double_t a0pi_recoiltheta_binmin{hh_0pi_recoiltheta->GetXaxis()->GetXmin()},
     //a0pi_recoiltheta_binmax{hh_0pi_recoiltheta->GetXaxis()->GetXmax()};

  TH1D * hh_0pi_neutronmomentum  = (TH1D*) tmplist1_neutronmomentum->At(0);   
  TMatrixD * hh_0pi_neutronmomentum_cov = (TMatrixD*) (tmplist1_neutronmomentum->At(2));
  hh_0pi_neutronmomentum_cov = (TMatrixD*) hh_0pi_neutronmomentum_cov->Clone("cov");
  //(*hh_0pi_neutronmomentum_cov) *= 1E6;
  //style::ScaleXaxis(hh_0pi_neutronmomentum, 1E-3);

  TH1D * hh_0pi_dpt  = (TH1D*) tmplist1_dpt->At(0);
  TMatrixD * hh_0pi_dpt_cov = (TMatrixD*) (tmplist1_dpt->At(2));
  hh_0pi_dpt_cov = (TMatrixD*) hh_0pi_dpt_cov->Clone("cov");
  Int_t a0pi_dpt_binN{hh_0pi_dpt->GetNbinsX()};
  //Double_t a0pi_dpt_binmin{hh_0pi_dpt->GetXaxis()->GetXmin()},
     //a0pi_dpt_binmax{hh_0pi_dpt->GetXaxis()->GetXmax()};


  TH1D * hh_0pi_dphit = (TH1D*) tmplist1_dphit->At(0);
  TMatrixD * hh_0pi_dphit_cov = (TMatrixD*) (tmplist1_dphit->At(2));
  hh_0pi_dphit_cov = (TMatrixD*) hh_0pi_dphit_cov->Clone("cov");
  Int_t a0pi_dphit_binN{hh_0pi_dphit->GetNbinsX()};
  //Double_t a0pi_dphit_binmin{hh_0pi_dphit->GetXaxis()->GetXmin()},
     //a0pi_dphit_binmax{hh_0pi_dphit->GetXaxis()->GetXmax()};


  TH1D * hh_0pi_dalphat  = (TH1D*) tmplist1_dalphat->At(0);
  TMatrixD * hh_0pi_dalphat_cov = (TMatrixD*) (tmplist1_dalphat->At(2));
  hh_0pi_dalphat_cov = (TMatrixD*) hh_0pi_dalphat_cov->Clone("cov");
  Int_t a0pi_dalphat_binN{hh_0pi_dalphat->GetNbinsX()};
  //Double_t a0pi_dalphat_binmin{hh_0pi_dalphat->GetXaxis()->GetXmin()},
     //a0pi_dalphat_binmax{hh_0pi_dalphat->GetXaxis()->GetXmax()};



  TH1D * hh_pi_neutronmomentum  = (TH1D*) tmplist2_neutronmomentum->At(0);   
  TMatrixD * hh_pi_neutronmomentum_cov = (TMatrixD*) (tmplist2_neutronmomentum->At(2));
  hh_pi_neutronmomentum_cov = (TMatrixD*) hh_pi_neutronmomentum_cov->Clone("cov");
  (*hh_pi_neutronmomentum_cov) *= 1E6;
  style::ScaleXaxis(hh_pi_neutronmomentum, 1E-3);

  TH1D * hh_pi_dalphat  = (TH1D*) tmplist2_dalphat->At(0);   
  TMatrixD * hh_pi_dalphat_cov = (TMatrixD*) (tmplist2_dalphat->At(2));
  hh_pi_dalphat_cov = (TMatrixD*) hh_pi_dalphat_cov->Clone("cov");

  TH1D * hh_pi_dpTT  = (TH1D*) tmplist2_dpTT->At(0);   
  TMatrixD * hh_pi_dpTT_cov = (TMatrixD*) (tmplist2_dpTT->At(2));
  hh_pi_dpTT_cov = (TMatrixD*) hh_pi_dpTT_cov->Clone("cov");

  TH1D * hh_0pi_data[]={hh_0pi_scattermomentum, hh_0pi_scattertheta, 
    hh_0pi_recoilmomentum, hh_0pi_recoiltheta, hh_0pi_neutronmomentum, 
    hh_0pi_dpt, hh_0pi_dphit, hh_0pi_dalphat}; 

  TH1D * hh_pi_data[]={hh_pi_neutronmomentum, hh_pi_dalphat, hh_pi_dpTT}; 

  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot.root","recreate");
  TTree *t1 = (TTree*) f->Get("tree");
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("plot.pdf[");

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*13;//Scaling to per nucleon in CH 13=1+12
  
  const TString modes[]={"_all","_qe","_res","_dis","_2p2h", "_other"};
  const int pdsmin[]={0,     1, 2, 3, 4, 5};
  const int pdsmax[]={10000, 1, 2, 3, 4, 5};
  const Int_t nmode = sizeof(modes)/sizeof(TString);

  ULong64_t minnp, maxnp;
  const bool isTopoTask = GeneratorUtils::SetTopoCut(anaid, minnp, maxnp);

  const TString npiDecomp[]={"","_0pi","_1pi","_Mpi","_other"};
  //only gointo sub-category for TopoTask
  const int nNpiCut = isTopoTask? (sizeof(npiDecomp)/sizeof(TString)) : 1;
  cout<<"n pi cut : "<<nNpiCut<<endl;

  TList *lout=new TList;
  TLegend *leg = new TLegend(0.65, 0.45, 0.85, 0.85);
  leg->SetFillStyle(0);
  //leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.023);
  gStyle-> SetOptStat(0); // turn off stat box
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.1);
  gStyle-> SetPalette(kRainbow);
  Double_t xBj, Q2, Wtrue, scattertheta, scattermomentum, mesontheta, mesonmomentum, 
    recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat, dpTT;
  Int_t prod,evtMode;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("Wtrue",&Wtrue);
  t1->SetBranchAddress("scattertheta",&scattertheta);
  t1->SetBranchAddress("scattermomentum",&scattermomentum);
  t1->SetBranchAddress("mesontheta",&mesontheta);
  t1->SetBranchAddress("mesonmomentum",&mesonmomentum);
  t1->SetBranchAddress("recoiltheta",&recoiltheta);
  t1->SetBranchAddress("recoilmomentum",&recoilmomentum);
  t1->SetBranchAddress("baryontheta",&baryontheta);
  t1->SetBranchAddress("baryonmomentum",&baryonmomentum);
  t1->SetBranchAddress("IApN",&IApN); //p_n
  t1->SetBranchAddress("dpt",&dpt);
  t1->SetBranchAddress("dphit",&dphit);
  t1->SetBranchAddress("dalphat",&dalphat);
  t1->SetBranchAddress("dpTT",&dpTT);
  t1->SetBranchAddress("prod",&prod);
  t1->SetBranchAddress("evtMode",&evtMode);
  t1->SetBranchAddress("perweight",&perweight);

  TH2D *h_weight_QE_scattertheta = new TH2D("h_weight_QE_scattertheta","distribution of weight and scattertheta mode : QE",80,0,26,100,0.0000001,0.001);

  TH2D *hxQ2_QE = new TH2D("hxQ2_QE","distribution of x_{Bj} and Q^{2} mode : QE",60,0,2,80,0.1,10);
  TH2D *hxW_QE = new TH2D("hxW_QE","distribution of x_{Bj} and W mode : QE",60,0,2,80,0.6,3);
  TH1D *hxQ2_QE_x = new TH1D("hxQ2_QE_x","x_{Bj} mode : QE",80,0,2);
  TH1D *hxQ2_QE_Q2 = new TH1D("hxQ2_QE_Q2","Q^{2} mode : QE",80,0,5);
  TH1D *hxQ2_QE_W = new TH1D("hxQ2_QE_W","W mode : QE",80,0,4);
  
    cout<<"init hist in process"<<endl;
    TH1D *hxQ2_QE_scattertheta = new TH1D(*hh_0pi_scattertheta); //muon theta
    hxQ2_QE_scattertheta->SetNameTitle("hxQ2_QE_scattertheta","#mu #theta mode : QE");
    wipe_hist(hxQ2_QE_scattertheta); //make all the bin content to 0 before filling
    TH1D *hxQ2_QE_scattermomentum = new TH1D(*hh_0pi_scattermomentum);
    hxQ2_QE_scattermomentum->SetNameTitle("hxQ2_QE_scattermomentum","#mu p mode : QE");
    wipe_hist(hxQ2_QE_scattermomentum);
    TH1D *hxQ2_QE_mesontheta = new TH1D("hxQ2_QE_mesontheta","#pi #theta mode : QE",
      80,0,180); //muon theta
    TH1D *hxQ2_QE_mesonmomentum = new TH1D("hxQ2_QE_mesonmomentum","#pi p mode : QE",
      80,0,7);
    TH1D *hxQ2_QE_recoiltheta = new TH1D(*hh_0pi_recoiltheta);
    hxQ2_QE_recoiltheta->SetNameTitle("hxQ2_QE_recoiltheta","recoil #theta mode : QE");
    wipe_hist(hxQ2_QE_recoiltheta);
    TH1D *hxQ2_QE_recoilmomentum = new TH1D(*hh_0pi_recoilmomentum);
    hxQ2_QE_recoilmomentum->SetNameTitle("hxQ2_QE_recoilmomentum","recoil p mode : QE");
    wipe_hist(hxQ2_QE_recoilmomentum);
    TH1D *hxQ2_QE_baryontheta = new TH1D("hxQ2_QE_baryontheta","baryon #theta mode : QE",
      80,0,180); //muon theta
    TH1D *hxQ2_QE_baryonmomentum = new TH1D("hxQ2_QE_baryonmomentum","baryon p mode : QE",
      80,0,10);
    TH1D *hxQ2_QE_neutronmomentum = new TH1D(*hh_0pi_neutronmomentum); 
    hxQ2_QE_neutronmomentum ->SetNameTitle("hxQ2_QE_neutronmomentum","neutron p mode : QE");
    wipe_hist(hxQ2_QE_neutronmomentum);
    TH1D *hxQ2_QE_dpt = new TH1D(*hh_0pi_dpt);
    hxQ2_QE_dpt->SetNameTitle("hxQ2_QE_dpt","dpt mode : QE");
    wipe_hist(hxQ2_QE_dpt);
    TH1D *hxQ2_QE_dphit = new TH1D(*hh_0pi_dphit);
    hxQ2_QE_dphit->SetNameTitle("hxQ2_QE_dphit","#delta#phi_{t} mode : QE");
    wipe_hist(hxQ2_QE_dphit);
    TH1D *hxQ2_QE_dalphat = new TH1D(*hh_0pi_dalphat);
    hxQ2_QE_dalphat->SetNameTitle("hxQ2_QE_dalphat","#delta#alpha_{t} mode : QE");
    wipe_hist(hxQ2_QE_dalphat);
  
  if(anaid == 2)
  {
    TH1D *hxQ2_QE_scattertheta = new TH1D("hxQ2_QE_scattertheta","#mu #theta mode : QE",80,0,26); //muon theta
    TH1D *hxQ2_QE_scattermomentum = new TH1D("hxQ2_QE_scattermomentum","#mu p mode : QE",80,1,15);
    TH1D *hxQ2_QE_mesontheta = new TH1D("hxQ2_QE_mesontheta","#pi #theta mode : QE",80,0,180); //muon theta
    TH1D *hxQ2_QE_mesonmomentum = new TH1D("hxQ2_QE_mesonmomentum","#pi p mode : QE",80,0,7);
    TH1D *hxQ2_QE_recoiltheta = new TH1D("hxQ2_QE_recoiltheta","recoil #theta mode : QE",80,0,180); //muon theta
    TH1D *hxQ2_QE_recoilmomentum = new TH1D("hxQ2_QE_recoilmomentum","recoil p mode : QE",80,0,6);
    TH1D *hxQ2_QE_baryontheta = new TH1D("hxQ2_QE_baryontheta","baryon #theta mode : QE",80,0,180); //muon theta
    TH1D *hxQ2_QE_baryonmomentum = new TH1D("hxQ2_QE_baryonmomentum","baryon p mode : QE",80,0,10);
    TH1D *hxQ2_QE_neutronmomentum = new TH1D("hxQ2_QE_neutronmomentum","neutron p mode : QE",80,0,2);
    TH1D *hxQ2_QE_dpt = new TH1D("hxQ2_QE_dpt","dpt mode : QE",80,0,2.5);
    TH1D *hxQ2_QE_dphit = new TH1D("hxQ2_QE_dphit","#delta#phi_{t} mode : QE",80,0,180);
    TH1D *hxQ2_QE_dalphat = new TH1D("hxQ2_QE_dalphat","#delta#alpha_{t} mode : QE",80,0,180);
  }

  TH2D *hxQ2_RES = new TH2D("hxQ2_RES","distribution of x_{Bj} and Q^{2} mode : RES",60,0,1.1,80,0.1,10);
  TH2D *hxW_RES = new TH2D("hxW_RES","distribution of x_{Bj} and W mode : RES",60,0,1.1,80,0.9,3);
  TH1D *hxQ2_RES_x = new TH1D("hxQ2_RES_x","x_{Bj} mode : RES",80,0,2);
  TH1D *hxQ2_RES_Q2 = new TH1D("hxQ2_RES_Q2","Q^{2} mode : RES",80,0,5);
  TH1D *hxQ2_RES_W = new TH1D("hxQ2_RES_W","W mode : RES",80,0,4);
  
    TH1D *hxQ2_RES_scattertheta = new TH1D(*hh_0pi_scattertheta); //muon theta
    hxQ2_RES_scattertheta->SetNameTitle("hxQ2_RES_scattertheta","#mu #theta mode : RES");
    wipe_hist(hxQ2_RES_scattertheta); //make all the bin content to 0 before filling
    TH1D *hxQ2_RES_scattermomentum = new TH1D(*hh_0pi_scattermomentum);
    hxQ2_RES_scattermomentum->SetNameTitle("hxQ2_RES_scattermomentum","#mu p mode : RES");
    wipe_hist(hxQ2_RES_scattermomentum);
    TH1D *hxQ2_RES_mesontheta = new TH1D("hxQ2_RES_mesontheta","#pi #theta mode : RES",
      80,0,180); //muon theta
    TH1D *hxQ2_RES_mesonmomentum = new TH1D("hxQ2_RES_mesonmomentum","#pi p mode : RES",
      80,0,7);
    TH1D *hxQ2_RES_recoiltheta = new TH1D(*hh_0pi_recoiltheta);
    hxQ2_RES_recoiltheta->SetNameTitle("hxQ2_RES_recoiltheta","recoil #theta mode : RES");
    wipe_hist(hxQ2_RES_recoiltheta);
    TH1D *hxQ2_RES_recoilmomentum = new TH1D(*hh_0pi_recoilmomentum);
    hxQ2_RES_recoilmomentum->SetNameTitle("hxQ2_RES_recoilmomentum","recoil p mode : RES");
    wipe_hist(hxQ2_RES_recoilmomentum);
    TH1D *hxQ2_RES_baryontheta = new TH1D("hxQ2_RES_baryontheta","baryon #theta mode : RES",
      80,0,180); //muon theta
    TH1D *hxQ2_RES_baryonmomentum = new TH1D("hxQ2_RES_baryonmomentum","baryon p mode : RES",
      80,0,10);
    TH1D *hxQ2_RES_neutronmomentum = new TH1D(*hh_0pi_neutronmomentum); 
    hxQ2_RES_neutronmomentum ->SetNameTitle("hxQ2_RES_neutronmomentum","neutron p mode : RES");
    wipe_hist(hxQ2_RES_neutronmomentum);
    TH1D *hxQ2_RES_dpt = new TH1D(*hh_0pi_dpt);
    hxQ2_RES_dpt->SetNameTitle("hxQ2_RES_dpt","dpt mode : RES");
    wipe_hist(hxQ2_RES_dpt);
    TH1D *hxQ2_RES_dphit = new TH1D(*hh_0pi_dphit);
    hxQ2_RES_dphit->SetNameTitle("hxQ2_RES_dphit","#delta#phi_{t} mode : RES");
    wipe_hist(hxQ2_RES_dphit);
    TH1D *hxQ2_RES_dalphat = new TH1D(*hh_0pi_dalphat);
    hxQ2_RES_dalphat->SetNameTitle("hxQ2_RES_dalphat","#delta#alpha_{t} mode : RES");
    wipe_hist(hxQ2_RES_dalphat);
  
  if(anaid == 2)
  {
    TH1D *hxQ2_RES_scattertheta = new TH1D("hxQ2_RES_scattertheta","#mu #theta mode : RES",80,0,26); //muon theta
    TH1D *hxQ2_RES_scattermomentum = new TH1D("hxQ2_RES_scattermomentum","#mu p mode : RES",80,1,15);
    TH1D *hxQ2_RES_mesontheta = new TH1D("hxQ2_RES_mesontheta","#pi #theta mode : RES",80,0,180); //muon theta
    TH1D *hxQ2_RES_mesonmomentum = new TH1D("hxQ2_RES_mesonmomentum","#pi p mode : RES",80,0,7);
    TH1D *hxQ2_RES_recoiltheta = new TH1D("hxQ2_RES_recoiltheta","recoil #theta mode : RES",80,0,180); //muon theta
    TH1D *hxQ2_RES_recoilmomentum = new TH1D("hxQ2_RES_recoilmomentum","recoil p mode : RES",80,0,6);
    TH1D *hxQ2_RES_baryontheta = new TH1D("hxQ2_RES_baryontheta","baryon #theta mode : RES",80,0,180); //muon theta
    TH1D *hxQ2_RES_baryonmomentum = new TH1D("hxQ2_RES_baryonmomentum","baryon p mode : RES",80,0,10);
    TH1D *hxQ2_RES_neutronmomentum = new TH1D("hxQ2_RES_neutronmomentum","neutron p mode : RES",80,0,2);
    TH1D *hxQ2_RES_dpt = new TH1D("hxQ2_RES_dpt","dpt mode : RES",80,0,2.5);
    TH1D *hxQ2_RES_dphit = new TH1D("hxQ2_RES_dphit","#delta#phi_{t} mode : RES",80,0,180);
    TH1D *hxQ2_RES_dalphat = new TH1D("hxQ2_RES_dalphat","#delta#alpha_{t} mode : RES",80,0,180);
  }

  TH2D *hxQ2_DIS = new TH2D("hxQ2_DIS","distribution of x_{Bj} and Q^{2} mode : DIS",60,0,1.2,80,0.1,10);
  TH2D *hxW_DIS = new TH2D("hxW_DIS","distribution of x_{Bj} and W mode : DIS",60,0,2,80,0.8,3);
  TH1D *hxQ2_DIS_x = new TH1D("hxQ2_DIS_x","x_{Bj} mode : DIS",80,0,2);
  TH1D *hxQ2_DIS_Q2 = new TH1D("hxQ2recoiltheta_DIS_Q2","Q^{2} mode : DIS",80,0,5);
  TH1D *hxQ2_DIS_W = new TH1D("hxQ2_DIS_W","W mode : DIS",80,0,4);
 
    TH1D *hxQ2_DIS_scattertheta = new TH1D(*hh_0pi_scattertheta); //muon theta
    hxQ2_DIS_scattertheta->SetNameTitle("hxQ2_DIS_scattertheta","#mu #theta mode : DIS");
    wipe_hist(hxQ2_DIS_scattertheta); //make all the bin content to 0 before filling
    TH1D *hxQ2_DIS_scattermomentum = new TH1D(*hh_0pi_scattermomentum);
    hxQ2_DIS_scattermomentum->SetNameTitle("hxQ2_DIS_scattermomentum","#mu p mode : DIS");
    wipe_hist(hxQ2_DIS_scattermomentum);
    TH1D *hxQ2_DIS_mesontheta = new TH1D("hxQ2_DIS_mesontheta","#pi #theta mode : DIS",
      80,0,180); //muon theta
    TH1D *hxQ2_DIS_mesonmomentum = new TH1D("hxQ2_DIS_mesonmomentum","#pi p mode : DIS",
      80,0,7);
    TH1D *hxQ2_DIS_recoiltheta = new TH1D(*hh_0pi_recoiltheta);
    hxQ2_DIS_recoiltheta->SetNameTitle("hxQ2_DIS_recoiltheta","recoil #theta mode : DIS");
    wipe_hist(hxQ2_DIS_recoiltheta);
    TH1D *hxQ2_DIS_recoilmomentum = new TH1D(*hh_0pi_recoilmomentum);
    hxQ2_DIS_recoilmomentum->SetNameTitle("hxQ2_DIS_recoilmomentum","recoil p mode : DIS");
    wipe_hist(hxQ2_DIS_recoilmomentum);
    TH1D *hxQ2_DIS_baryontheta = new TH1D("hxQ2_DIS_baryontheta","baryon #theta mode : DIS",
      80,0,180); //muon theta
    TH1D *hxQ2_DIS_baryonmomentum = new TH1D("hxQ2_DIS_baryonmomentum","baryon p mode : DIS",
      80,0,10);
    TH1D *hxQ2_DIS_neutronmomentum = new TH1D(*hh_0pi_neutronmomentum); 
    hxQ2_DIS_neutronmomentum ->SetNameTitle("hxQ2_DIS_neutronmomentum","neutron p mode : DIS");
    wipe_hist(hxQ2_DIS_neutronmomentum);
    TH1D *hxQ2_DIS_dpt = new TH1D(*hh_0pi_dpt);
    hxQ2_DIS_dpt->SetNameTitle("hxQ2_DIS_dpt","dpt mode : DIS");
    wipe_hist(hxQ2_DIS_dpt);
    TH1D *hxQ2_DIS_dphit = new TH1D(*hh_0pi_dphit);
    hxQ2_DIS_dphit->SetNameTitle("hxQ2_DIS_dphit","#delta#phi_{t} mode : DIS");
    wipe_hist(hxQ2_DIS_dphit);
    TH1D *hxQ2_DIS_dalphat = new TH1D(*hh_0pi_dalphat);
    hxQ2_DIS_dalphat->SetNameTitle("hxQ2_DIS_dalphat","#delta#alpha_{t} mode : DIS");
    wipe_hist(hxQ2_DIS_dalphat);
  
 
  if(anaid == 2)
  {
    TH1D *hxQ2_DIS_scattertheta = new TH1D("hxQ2_DIS_scattertheta","#mu #theta mode : DIS",80,0,26); //muon theta
    TH1D *hxQ2_DIS_scattermomentum = new TH1D("hxQ2_DIS_scattermomentum","#mu p mode : DIS",80,1,15);
    TH1D *hxQ2_DIS_mesontheta = new TH1D("hxQ2_DIS_mesontheta","#pi #theta mode : DIS",80,0,180); //muon theta
    TH1D *hxQ2_DIS_mesonmomentum = new TH1D("hxQ2_DIS_mesonmomentum","#pi p mode : DIS",80,0,7);
    TH1D *hxQ2_DIS_recoiltheta = new TH1D("hxQ2_DIS_recoiltheta","recoil #theta mode : DIS",80,0,180); //muon theta
    TH1D *hxQ2_DIS_recoilmomentum = new TH1D("hxQ2_DIS_recoilmomentum","recoil p mode : DIS",80,0,6);
    TH1D *hxQ2_DIS_baryontheta = new TH1D("hxQ2_DIS_baryontheta","baryon #theta mode : DIS",80,0,180); //muon theta
    TH1D *hxQ2_DIS_baryonmomentum = new TH1D("hxQ2_DIS_baryonmomentum","baryon p mode : DIS",80,0,10);
    TH1D *hxQ2_DIS_neutronmomentum = new TH1D("hxQ2_DIS_neutronmomentum","neutron p mode : DIS",80,0,2);
    TH1D *hxQ2_DIS_dpt = new TH1D("hxQ2_DIS_dpt","dpt mode : DIS",80,0,2.5);
    TH1D *hxQ2_DIS_dphit = new TH1D("hxQ2_DIS_dphit","#delta#phi_{t} mode : DIS",80,0,180);
    TH1D *hxQ2_DIS_dalphat = new TH1D("hxQ2_DIS_dalphat","#delta#alpha_{t} mode : DIS",80,0,180);
  }
  TH2D *hxQ2_2p2h = new TH2D("hxQ2_2p2h","distribution of x_{Bj} and Q^{2} mode : 2p2h",60,0,2,80,0.2,10);
  TH2D *hxW_2p2h = new TH2D("hxW_2p2h","distribution of x_{Bj} and W mode : 2p2h",60,0,2,80,0.1,1.5);
  TH1D *hxQ2_2p2h_x = new TH1D("hxQ2_2p2h_x","x_{Bj} mode : 2p2h",80,0,2);
  TH1D *hxQ2_2p2h_Q2 = new TH1D("hxQ2_2p2h_Q2","Q^{2} mode : 2p2h",80,0,5);
  TH1D *hxQ2_2p2h_W = new TH1D("hxQ2_2p2h_W","W mode : 2p2h",80,0,4);
  
    TH1D *hxQ2_2p2h_scattertheta = new TH1D(*hh_0pi_scattertheta); //muon theta
    hxQ2_2p2h_scattertheta->SetNameTitle("hxQ2_2p2h_scattertheta","#mu #theta mode : 2p2h");
    wipe_hist(hxQ2_2p2h_scattertheta); //make all the bin content to 0 before filling
    TH1D *hxQ2_2p2h_scattermomentum = new TH1D(*hh_0pi_scattermomentum);
    hxQ2_2p2h_scattermomentum->SetNameTitle("hxQ2_2p2h_scattermomentum","#mu p mode : 2p2h");
    wipe_hist(hxQ2_2p2h_scattermomentum);
    TH1D *hxQ2_2p2h_mesontheta = new TH1D("hxQ2_2p2h_mesontheta","#pi #theta mode : 2p2h",
      80,0,180); //muon theta
    TH1D *hxQ2_2p2h_mesonmomentum = new TH1D("hxQ2_2p2h_mesonmomentum","#pi p mode : 2p2h",
      80,0,7);
    TH1D *hxQ2_2p2h_recoiltheta = new TH1D(*hh_0pi_recoiltheta);
    hxQ2_2p2h_recoiltheta->SetNameTitle("hxQ2_2p2h_recoiltheta","recoil #theta mode : 2p2h");
    wipe_hist(hxQ2_2p2h_recoiltheta);
    TH1D *hxQ2_2p2h_recoilmomentum = new TH1D(*hh_0pi_recoilmomentum);
    hxQ2_2p2h_recoilmomentum->SetNameTitle("hxQ2_2p2h_recoilmomentum","recoil p mode : 2p2h");
    wipe_hist(hxQ2_2p2h_recoilmomentum);
    TH1D *hxQ2_2p2h_baryontheta = new TH1D("hxQ2_2p2h_baryontheta","baryon #theta mode : 2p2h",
      80,0,180); //muon theta
    TH1D *hxQ2_2p2h_baryonmomentum = new TH1D("hxQ2_2p2h_baryonmomentum","baryon p mode : 2p2h",
      80,0,10);
    TH1D *hxQ2_2p2h_neutronmomentum = new TH1D(*hh_0pi_neutronmomentum); 
    hxQ2_2p2h_neutronmomentum ->SetNameTitle("hxQ2_2p2h_neutronmomentum","neutron p mode : 2p2h");
    wipe_hist(hxQ2_2p2h_neutronmomentum);
    TH1D *hxQ2_2p2h_dpt = new TH1D(*hh_0pi_dpt);
    hxQ2_2p2h_dpt->SetNameTitle("hxQ2_2p2h_dpt","dpt mode : 2p2h");
    wipe_hist(hxQ2_2p2h_dpt);
    TH1D *hxQ2_2p2h_dphit = new TH1D(*hh_0pi_dphit);
    hxQ2_2p2h_dphit->SetNameTitle("hxQ2_2p2h_dphit","#delta#phi_{t} mode : 2p2h");
    wipe_hist(hxQ2_2p2h_dphit);
    TH1D *hxQ2_2p2h_dalphat = new TH1D(*hh_0pi_dalphat);
    hxQ2_2p2h_dalphat->SetNameTitle("hxQ2_2p2h_dalphat","#delta#alpha_{t} mode : 2p2h");
    wipe_hist(hxQ2_2p2h_dalphat);
  if(anaid == 2)
  {
    TH1D *hxQ2_2p2h_scattertheta = new TH1D("hxQ2_2p2h_scattertheta","#mu #theta mode : 2p2h",80,0,26); //muon theta
    TH1D *hxQ2_2p2h_scattermomentum = new TH1D("hxQ2_2p2h_scattermomentum","#mu p mode : 2p2h",80,1,15);
    TH1D *hxQ2_2p2h_mesontheta = new TH1D("hxQ2_2p2h_mesontheta","#pi #theta mode : 2p2h",80,0,180); //muon theta
    TH1D *hxQ2_2p2h_mesonmomentum = new TH1D("hxQ2_2p2h_mesonmomentum","#pi p mode : 2p2h",80,0,7);
    TH1D *hxQ2_2p2h_recoiltheta = new TH1D("hxQ2_2p2h_recoiltheta","recoil #theta mode : 2p2h",80,0,180); //muon theta
    TH1D *hxQ2_2p2h_recoilmomentum = new TH1D("hxQ2_2p2h_recoilmomentum","recoil p mode : 2p2h",80,0,6);
    TH1D *hxQ2_2p2h_baryontheta = new TH1D("hxQ2_2p2h_baryontheta","baryon #theta mode : 2p2h",80,0,180); //muon theta
    TH1D *hxQ2_2p2h_baryonmomentum = new TH1D("hxQ2_2p2h_baryonmomentum","baryon p mode : 2p2h",80,0,10);
    TH1D *hxQ2_2p2h_neutronmomentum = new TH1D("hxQ2_2p2h_neutronmomentum","neutron p mode : 2p2h",80,0,2);
    TH1D *hxQ2_2p2h_dpt = new TH1D("hxQ2_2p2h_dpt","dpt mode : 2p2h",80,0,2.5);
    TH1D *hxQ2_2p2h_dphit = new TH1D("hxQ2_2p2h_dphit","#delta#phi_{t} mode : 2p2h",80,0,180);
    TH1D *hxQ2_2p2h_dalphat = new TH1D("hxQ2_2p2h_dalphat","#delta#alpha_{t} mode : 2p2h",80,0,180);
  }

  TString tit[]={"evtMode QE x_{Bj}-Q^{2}", "evtMode QE x_{Bj}-W","QE x_{Bj}", "QE Q^{2}(GeV^{2})", "QE W(GeV)"
  , "evtMode RES x_{Bj}-Q^{2}", "evtMode RES x_{Bj}-W","RES x_{Bj}", "RES Q^{2}(GeV^{2})", "RES W(GeV)"
  , "evtMode DIS x_{Bj}-Q^{2}", "evtMode DIS x_{Bj}-W","DIS x_{Bj}", "DIS Q^{2}(GeV^{2})", "DIS W(GeV)"
  , "evtMode 2p2h x_{Bj}-Q^{2}","evtMode 2p2h x_{Bj}-W","2p2h x_{Bj}", "2p2h Q^{2}(GeV^{2})", "2p2h W(GeV)"};
  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange};
  TH2D *hh2d[]={hxQ2_QE, hxW_QE, hxQ2_RES, hxW_RES, hxQ2_DIS, hxW_DIS, hxQ2_2p2h, hxW_2p2h}; //list of 2d hist
  TH1D *hh1d[]={hxQ2_QE_x, hxQ2_QE_Q2, hxQ2_QE_W, hxQ2_RES_x, hxQ2_RES_Q2, hxQ2_RES_W, hxQ2_DIS_x, hxQ2_DIS_Q2, hxQ2_DIS_W, 
  	hxQ2_2p2h_x, hxQ2_2p2h_Q2, hxQ2_2p2h_W}; //list of 1d hist
    
  TH1D *hh1d_mom[]={hxQ2_QE_scattertheta, hxQ2_QE_scattermomentum,  hxQ2_QE_mesontheta, hxQ2_QE_mesonmomentum, 
    hxQ2_QE_recoiltheta, hxQ2_QE_recoilmomentum,  hxQ2_QE_baryontheta, hxQ2_QE_baryonmomentum, 
    hxQ2_QE_neutronmomentum, hxQ2_QE_dpt, hxQ2_QE_dphit, hxQ2_QE_dalphat, 
    hxQ2_RES_scattertheta, hxQ2_RES_scattermomentum,  hxQ2_RES_mesontheta, hxQ2_RES_mesonmomentum, 
    hxQ2_RES_recoiltheta, hxQ2_RES_recoilmomentum,  hxQ2_RES_baryontheta, hxQ2_RES_baryonmomentum, 
    hxQ2_RES_neutronmomentum, hxQ2_RES_dpt, hxQ2_RES_dphit, hxQ2_RES_dalphat, 
    hxQ2_DIS_scattertheta, hxQ2_DIS_scattermomentum, hxQ2_DIS_mesontheta, hxQ2_DIS_mesonmomentum, 
    hxQ2_DIS_recoiltheta, hxQ2_DIS_recoilmomentum,  hxQ2_DIS_baryontheta, hxQ2_DIS_baryonmomentum, 
    hxQ2_DIS_neutronmomentum, hxQ2_DIS_dpt, hxQ2_DIS_dphit, hxQ2_DIS_dalphat, 
    hxQ2_2p2h_scattertheta, hxQ2_2p2h_scattermomentum, hxQ2_2p2h_mesontheta, hxQ2_2p2h_mesonmomentum,
    hxQ2_2p2h_recoiltheta, hxQ2_2p2h_recoilmomentum,  hxQ2_2p2h_baryontheta, hxQ2_2p2h_baryonmomentum,
    hxQ2_2p2h_neutronmomentum, hxQ2_2p2h_dpt, hxQ2_2p2h_dphit, hxQ2_2p2h_dalphat};

  TH1D *hh1d_TKI[]={hxQ2_QE_dpt, hxQ2_QE_dphit, hxQ2_RES_dalphat, 
    hxQ2_RES_dpt, hxQ2_RES_dphit, hxQ2_RES_dalphat, hxQ2_DIS_dpt, hxQ2_DIS_dphit, hxQ2_DIS_dalphat, 
    hxQ2_DIS_dalphat, hxQ2_2p2h_dphit, hxQ2_2p2h_dalphat};
  //const Int_t aln_by_var[]={0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15};//list for index of hh1d in order of QE x,RES x...Q2,Q2,..,W,W
  string mode[] = {"QE", "RES", "DIS", "2p2h"};
  string var[] = {"x", "Q2", "W", "scattertheta", "scattermomentum", 
    "mesontheta", "mesonmomentum", "recoiltheta", "recoilmomentum", "baryontheta", "baryonmomentum",
    "neutronmomentum", "dpt", "dphit", "dalphat"};
  
  THStack *stk_x = new THStack;
  TLegend *leg_x = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_x->SetFillStyle(0);
  
  THStack *stk_Q2 = new THStack;
  TLegend *leg_Q2 = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_Q2->SetFillStyle(0);
  
  THStack *stk_W = new THStack;
  TLegend *leg_W = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_W->SetFillStyle(0);
  
  THStack *stk_mu_theta = new THStack;
  TLegend *leg_mu_theta = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_mu_theta->SetFillStyle(0);
  
  THStack *stk_mu_mom = new THStack;
  TLegend *leg_mu_mom = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_mu_mom->SetFillStyle(0);

  THStack *stk_pi_theta = new THStack;
  TLegend *leg_pi_theta = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_pi_theta->SetFillStyle(0);
  
  THStack *stk_pi_mom = new THStack;
  TLegend *leg_pi_mom = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_pi_mom->SetFillStyle(0);

  THStack *stk_recoil_theta = new THStack;
  TLegend *leg_recoil_theta = new TLegend(0.13, 0.5, 0.33,0.9);
  leg_recoil_theta->SetFillStyle(0);
  
  THStack *stk_recoil_mom = new THStack;
  TLegend *leg_recoil_mom = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_recoil_mom->SetFillStyle(0);

  THStack *stk_baryon_theta = new THStack;
  TLegend *leg_baryon_theta = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_baryon_theta->SetFillStyle(0);
  
  THStack *stk_baryon_mom = new THStack;
  TLegend *leg_baryon_mom = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_baryon_mom->SetFillStyle(0);

  THStack *stk_neutron_mom = new THStack;
  TLegend *leg_neutron_mom = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_neutron_mom->SetFillStyle(0);
 
  THStack *stk_dpt = new THStack;
  TLegend *leg_dpt = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_dpt->SetFillStyle(0);

  THStack *stk_dphit = new THStack;
  TLegend *leg_dphit = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_dphit->SetFillStyle(0);
 
  THStack *stk_dalphat = new THStack;
  TLegend *leg_dalphat = new TLegend(0.13, 0.5, 0.33,0.9);
  leg_dalphat->SetFillStyle(0);
 
 
  if(anaid == 1)
  {
    
    for(Int_t i=0; i<7; i++)
    {
      hh_0pi_data[i]->SetMarkerSize(1);
      hh_0pi_data[i]->SetMarkerStyle(20);
      hh_0pi_data[i]->SetLineStyle(kSolid);
      hh_0pi_data[i]-> SetLineWidth(2); 
    }

    leg_mu_mom->AddEntry(hh_0pi_scattermomentum,"Data total","fl");
    leg_mu_theta->AddEntry(hh_0pi_scattertheta,"Data total","fl");
    leg_recoil_mom->AddEntry(hh_0pi_recoilmomentum,"Data total","fl");
    leg_recoil_theta->AddEntry(hh_0pi_recoiltheta,"Data total","fl");
    leg_neutron_mom->AddEntry(hh_0pi_neutronmomentum,"Data total","fl");
    leg_dpt->AddEntry(hh_0pi_dpt,"Data total","fl");
    leg_dphit->AddEntry(hh_0pi_dphit,"Data total","fl");
    leg_dalphat->AddEntry(hh_0pi_dalphat,"Data total","fl");
  }
  else
  {
    hh_pi_dalphat->SetMarkerSize(1);
    hh_pi_dalphat->SetMarkerStyle(20);
    hh_pi_dalphat->SetLineStyle(kSolid);
    hh_pi_dalphat-> SetLineWidth(2); 

    leg_dalphat->AddEntry(hh_pi_dalphat,"Data total","fl");
  }
 
  Int_t nentries = (Int_t)t1->GetEntries();
  const int nh1d = sizeof(hh1d)/sizeof(TH1D*);
  const int nh1d_mom = sizeof(hh1d_mom)/sizeof(TH1D*);
  cout<<"mom array len : "<<nh1d_mom<<endl;
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);
  int npi_err{0};
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    //if(perweight<0 && i%500==0){cout<<"weight : "<<perweight<<endl;}
    if(i%10000==0){cout<<"meson(pi) p : "<<mesonmomentum<<endl;}
    if(perweight>4000){
      printf("\n\n\nAlert!!  Filling Super weight!! %d %f skiping...\n\n\n", i, perweight);
    }
    if(anaid == CC1piNpID){
      if(targetZ == 1){//no hydrogen for GFS1
        cout<<"no hydrogen for GFS1"<<endl;
      }
    }
    if(anaid == 1){if(GeneratorUtils::GetNpi(npar)==0){npi_err++;}}
    else{if(!(GeneratorUtils::GetNpi(npar)==0)){npi_err++;}}
    if(prod==1 && evtMode==1) //fil histogram bins for quasi elastic case
    {
      hxQ2_QE->Fill(xBj,Q2,perweight);
      hxW_QE->Fill(xBj,Wtrue,perweight);
      hxQ2_QE_x->Fill(xBj,perweight);
      hxQ2_QE_Q2->Fill(Q2,perweight);
      hxQ2_QE_W->Fill(Wtrue,perweight);
      hxQ2_QE_scattertheta->Fill(scattertheta,perweight);
      hxQ2_QE_scattermomentum->Fill(scattermomentum,perweight);
      hxQ2_QE_mesontheta->Fill(mesontheta,perweight);
      hxQ2_QE_mesonmomentum->Fill(mesonmomentum,perweight);
      hxQ2_QE_recoiltheta->Fill(recoiltheta,perweight);
      hxQ2_QE_recoilmomentum->Fill(recoilmomentum,perweight);
      hxQ2_QE_baryontheta->Fill(baryontheta,perweight);
      hxQ2_QE_baryonmomentum->Fill(baryonmomentum,perweight);
      hxQ2_QE_neutronmomentum->Fill(IApN,perweight);
      hxQ2_QE_dpt->Fill(dpt,perweight);
      hxQ2_QE_dphit->Fill(dphit,perweight);
      hxQ2_QE_dalphat->Fill(dalphat,perweight);

      h_weight_QE_scattertheta->Fill(scattertheta,perweight);
    }
    else if(prod>=2 && prod<=33 && evtMode==2) //fil histogram bins for resonance case
    {
      hxQ2_RES->Fill(xBj,Q2,perweight);
      hxW_RES->Fill(xBj,Wtrue,perweight);
      hxQ2_RES_x->Fill(xBj,perweight);
      hxQ2_RES_Q2->Fill(Q2,perweight);
      hxQ2_RES_W->Fill(Wtrue,perweight);
      hxQ2_RES_scattertheta->Fill(scattertheta,perweight);
      hxQ2_RES_scattermomentum->Fill(scattermomentum,perweight);
      hxQ2_RES_mesontheta->Fill(mesontheta,perweight);
      hxQ2_RES_mesonmomentum->Fill(mesonmomentum,perweight);
      hxQ2_RES_recoiltheta->Fill(recoiltheta,perweight);
      hxQ2_RES_recoilmomentum->Fill(recoilmomentum,perweight);
      hxQ2_RES_baryontheta->Fill(baryontheta,perweight);
      hxQ2_RES_baryonmomentum->Fill(baryonmomentum,perweight);
      hxQ2_RES_neutronmomentum->Fill(IApN,perweight);
      hxQ2_RES_dpt->Fill(dpt,perweight);
      hxQ2_RES_dphit->Fill(dphit,perweight);
      hxQ2_RES_dalphat->Fill(dalphat,perweight);
    }
    else if(prod==34 && evtMode==3) //fil histogram bins for deep inelastic case
    {
      hxQ2_DIS->Fill(xBj,Q2,perweight);
      hxW_DIS->Fill(xBj,Wtrue,perweight);
      hxQ2_DIS_x->Fill(xBj,perweight);
      hxQ2_DIS_Q2->Fill(Q2,perweight);
      hxQ2_DIS_W->Fill(Wtrue,perweight);
      hxQ2_DIS_scattertheta->Fill(scattertheta,perweight);
      hxQ2_DIS_scattermomentum->Fill(scattermomentum,perweight);
      hxQ2_DIS_mesontheta->Fill(mesontheta,perweight);
      hxQ2_DIS_mesonmomentum->Fill(mesonmomentum,perweight);
      hxQ2_DIS_recoiltheta->Fill(recoiltheta,perweight);
      hxQ2_DIS_recoilmomentum->Fill(recoilmomentum,perweight);
      hxQ2_DIS_baryontheta->Fill(baryontheta,perweight);
      hxQ2_DIS_baryonmomentum->Fill(baryonmomentum,perweight);
      hxQ2_DIS_neutronmomentum->Fill(IApN,perweight);
      hxQ2_DIS_dpt->Fill(dpt,perweight);
      hxQ2_DIS_dphit->Fill(dphit,perweight);
      hxQ2_DIS_dalphat->Fill(dalphat,perweight);
    }
    else if(prod==35 && evtMode==4) //fil histogram bins for 2 part 2 hole case
    {
      hxQ2_2p2h->Fill(xBj,Q2,perweight);
      hxW_2p2h->Fill(xBj,Wtrue,perweight);
      hxQ2_2p2h_x->Fill(xBj,perweight);
      hxQ2_2p2h_Q2->Fill(Q2,perweight);
      hxQ2_2p2h_W->Fill(Wtrue,perweight);
      hxQ2_2p2h_scattertheta->Fill(scattertheta,perweight);
      hxQ2_2p2h_scattermomentum->Fill(scattermomentum,perweight);
      hxQ2_2p2h_mesontheta->Fill(mesontheta,perweight);
      hxQ2_2p2h_mesonmomentum->Fill(mesonmomentum,perweight);
      hxQ2_2p2h_recoiltheta->Fill(recoiltheta,perweight);
      hxQ2_2p2h_recoilmomentum->Fill(recoilmomentum,perweight);
      hxQ2_2p2h_baryontheta->Fill(baryontheta,perweight);
      hxQ2_2p2h_baryonmomentum->Fill(baryonmomentum,perweight);
      hxQ2_2p2h_neutronmomentum->Fill(IApN,perweight);
      hxQ2_2p2h_dpt->Fill(dpt,perweight);
      hxQ2_2p2h_dphit->Fill(dphit,perweight);
      hxQ2_2p2h_dalphat->Fill(dalphat,perweight);
    }
   
  }
  cout<<" n pi err # : "<<npi_err<<endl;

  Double_t temp_xsec{0};

  Double_t QE_x_xsec, QE_Q2_xsec, QE_W_xsec, QE_scattertheta_xsec, QE_scattermomentum_xsec, 
    QE_mesontheta_xsec, QE_mesonmomentum_xsec, QE_recoiltheta_xsec, QE_recoilmomentum_xsec, 
    QE_baryontheta_xsec, QE_baryonmomentum_xsec, QE_neutronmomentum_xsec, QE_dpt_xsec, QE_dphit_xsec, QE_dalphat_xsec, 
    RES_x_xsec, RES_Q2_xsec, RES_W_xsec, RES_scattertheta_xsec, RES_scattermomentum_xsec, 
    RES_mesontheta_xsec, RES_mesonmomentum_xsec, RES_recoiltheta_xsec, RES_recoilmomentum_xsec, 
    RES_baryontheta_xsec, RES_baryonmomentum_xsec, RES_neutronmomentum_xsec, RES_dpt_xsec, RES_dphit_xsec, RES_dalphat_xsec, 
    DIS_x_xsec, DIS_Q2_xsec, DIS_W_xsec, DIS_scattertheta_xsec, DIS_scattermomentum_xsec, 
    DIS_mesontheta_xsec, DIS_mesonmomentum_xsec, DIS_recoiltheta_xsec, DIS_recoilmomentum_xsec, 
    DIS_baryontheta_xsec, DIS_baryonmomentum_xsec, DIS_neutronmomentum_xsec, DIS_dpt_xsec, DIS_dphit_xsec, DIS_dalphat_xsec, 
    a2p2h_x_xsec, a2p2h_Q2_xsec, a2p2h_W_xsec, a2p2h_scattertheta_xsec, a2p2h_scattermomentum_xsec, 
    a2p2h_mesontheta_xsec, a2p2h_mesonmomentum_xsec, a2p2h_recoiltheta_xsec, 
    a2p2h_recoilmomentum_xsec, a2p2h_baryontheta_xsec, a2p2h_baryonmomentum_xsec, a2p2h_neutronmomentum_xsec, 
    a2p2h_dpt_xsec, a2p2h_dphit_xsec, a2p2h_dalphat_xsec;

  Double_t xsec_lst[] = {QE_x_xsec, QE_Q2_xsec, QE_W_xsec, QE_scattertheta_xsec, QE_scattermomentum_xsec, 
    QE_mesontheta_xsec, QE_mesonmomentum_xsec, QE_recoiltheta_xsec, QE_recoilmomentum_xsec, 
    QE_baryontheta_xsec, QE_baryonmomentum_xsec, QE_neutronmomentum_xsec, QE_dpt_xsec, QE_dphit_xsec, QE_dalphat_xsec, 
    RES_x_xsec, RES_Q2_xsec, RES_W_xsec, RES_scattertheta_xsec, RES_scattermomentum_xsec, 
    RES_mesontheta_xsec, RES_mesonmomentum_xsec, RES_recoiltheta_xsec, RES_recoilmomentum_xsec, 
    RES_baryontheta_xsec, RES_baryonmomentum_xsec, RES_neutronmomentum_xsec, RES_dpt_xsec, RES_dphit_xsec, RES_dalphat_xsec, 
    DIS_x_xsec, DIS_Q2_xsec, DIS_W_xsec, DIS_scattertheta_xsec, DIS_scattermomentum_xsec, 
    DIS_mesontheta_xsec, DIS_mesonmomentum_xsec, DIS_recoiltheta_xsec, DIS_recoilmomentum_xsec, 
    DIS_baryontheta_xsec, DIS_baryonmomentum_xsec, DIS_neutronmomentum_xsec, DIS_dpt_xsec, DIS_dphit_xsec, DIS_dalphat_xsec, 
    a2p2h_x_xsec, a2p2h_Q2_xsec, a2p2h_W_xsec, a2p2h_scattertheta_xsec, a2p2h_scattermomentum_xsec, 
    a2p2h_mesontheta_xsec, a2p2h_mesonmomentum_xsec, a2p2h_recoiltheta_xsec, 
    a2p2h_recoilmomentum_xsec, a2p2h_baryontheta_xsec, a2p2h_baryonmomentum_xsec, a2p2h_neutronmomentum_xsec, 
    a2p2h_dpt_xsec, a2p2h_dphit_xsec, a2p2h_dalphat_xsec};
  const int nxsec = sizeof(xsec_lst)/sizeof(Double_t);

  int normal_choice{1}; //if 0 just remove negative cross section if 1 do column normalization

  for(int ii=0; ii<nh2d; ii++){
    if(normal_choice==0)
    {
      hh2d[ii]->Scale(1/histnormFactor);
      remove_neg(hh2d[ii]);
    }
    else if(normal_choice==1)
    {
      col_normal(hh2d[ii]);
    }
    hh2d[ii]->GetXaxis()->SetTitle("x_{Bj}");
    if(ii%2==0)
    {
      hh2d[ii]->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");
      gPad->SetLogy(1); //set y axis to log-1 lin-0 scale
      gPad->Update();
    }
    else if(ii%2==1)
    {
      hh2d[ii]->GetYaxis()->SetTitle("W(GeV)");
      gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
      gPad->Update();
    }
    hh2d[ii]->SetContour(1000);
    //hh2d[ii]->Draw("TEXT");
    hh2d[ii]->Draw("colz");

    gStyle->SetPadRightMargin(0.2);
    //redraw legend and update bottom margin dependent on number of entries

    // add this plot to plotbook
    c1->Print("plot.pdf",hh2d[ii]->GetTitle());
    c1->Print(Form("png/%s_%s.png", anatag, hh2d[ii]->GetTitle()));

    leg->Clear();
  }
  /*
  remove_neg(h_weight_QE_scattertheta);
  h_weight_QE_scattertheta->GetXaxis()->SetTitle("QE #theta_{mu}");
  h_weight_QE_scattertheta->GetYaxis()->SetTitle("weight");
  gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
  gPad->Update();
  h_weight_QE_scattertheta->SetContour(1000);
  h_weight_QE_scattertheta->Draw("colz");
  gStyle->SetPadRightMargin(0.2);
  c1->Print("plot.pdf",h_weight_QE_scattertheta->GetTitle());
  c1->Print(Form("png/%s_%s.png", anatag, h_weight_QE_scattertheta->GetTitle()));
  */
  
  int aa{0};
  Double_t binwidth{0};
  for(int ii=0; ii<nh1d; ii++){
    if(aa%4==0){++aa;};
    binwidth = hh1d[ii]->GetXaxis()->GetBinWidth(1);
    hh1d[ii]->Scale(1/(histnormFactor*binwidth)); //to get diff cross section
    hh1d[ii]->GetXaxis()->SetTitle(tit[aa]);
    hh1d[ii]->GetYaxis()->SetTitle("differential cross section(cm^{2}/GeV/nucleon)");
    
    //leg->SetFillStyle(0);
    hh1d[ii]->SetLineStyle(kSolid);
    hh1d[ii]->SetLineColor(cols[ii/3]);
    hh1d[ii]->SetFillStyle(4050);
    hh1d[ii]->SetFillColorAlpha(cols[ii/3],0.35);

    temp_xsec = hh1d[ii]->Integral(0,10000,"width");
    xsec_lst[(ii/3)*14+ii%3] = temp_xsec;

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
  //binwidth = 0;
  int count_num{12}, stk_num{0};
  for(int ii=0; ii<nh1d_mom; ii++){
    if(aa%4==0){++aa;};
    diff_xsec(hh1d_mom[ii], histnormFactor);
    //hh1d_mom[ii]->GetXaxis()->SetTitle(tit[aa]);
    hh1d_mom[ii]->GetYaxis()->SetTitle("differential cross section(cm^{2}/GeV/nucleon)");
    
    //leg->SetFillStyle(0);
    hh1d_mom[ii]->SetLineStyle(kSolid);
    hh1d_mom[ii]->SetLineColor(cols[ii/count_num]);
    hh1d_mom[ii]->SetFillStyle(4050);
    hh1d_mom[ii]->SetFillColorAlpha(cols[ii/count_num],0.35);
    //cout<<"title : "<<hh1d_mom[ii]->GetTitle()<<endl;

    temp_xsec = hh1d_mom[ii]->Integral(0,10000,"width");
    xsec_lst[(ii/count_num)*15+3+ii%count_num] = temp_xsec;
    if(ii%count_num==0){
      stk_mu_theta->Add(hh1d_mom[ii]);
      leg_mu_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==1){
      stk_mu_mom->Add(hh1d_mom[ii]);
      leg_mu_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==2){
      stk_pi_theta->Add(hh1d_mom[ii]);
      leg_pi_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==3){
      stk_pi_mom->Add(hh1d_mom[ii]);
      leg_pi_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==4){
      stk_recoil_theta->Add(hh1d_mom[ii]);
      leg_recoil_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==5){
      stk_recoil_mom->Add(hh1d_mom[ii]);
      leg_recoil_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==6){
      stk_baryon_theta->Add(hh1d_mom[ii]);
      leg_baryon_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==7){
      stk_baryon_mom->Add(hh1d_mom[ii]);
      leg_baryon_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==8){
      stk_neutron_mom->Add(hh1d_mom[ii]);
      leg_neutron_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==9){
      stk_dpt->Add(hh1d_mom[ii]);
      leg_dpt->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==10){
      stk_dphit->Add(hh1d_mom[ii]);
      leg_dphit->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    else if(ii%count_num==11){
      stk_dalphat->Add(hh1d_mom[ii]);
      leg_dalphat->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
    }
    
    ++aa;
  }
  
  gPad->SetLogy(0);//set y axis to linear scale
  gPad->Update();

  aa=0;
  cout<<"xsec for "<<anatag<<endl;
  for(int ii=0; ii<nxsec; ii++)
  {
    cout<<"xsec for mode : "<<mode[ii/(nxsec/4)]<<" calc with : "<<var[ii%(nxsec/4)]
    <<" is "<<xsec_lst[ii]<<" cm^{2}"<<endl;
  }
  cout<<endl;
  Double_t comb_xsec{0};
  for(int ii=0; ii<(nxsec/4); ii++)
  {
    for(int jj=0; jj<4; jj++)
    {
      comb_xsec = comb_xsec+xsec_lst[jj*(nxsec/4)+ii];
    }
    cout<<"combind total xsec cal with "<<var[ii]<<" : "<<comb_xsec<<" cm^{2}"<<endl;
    comb_xsec = 0;
  }

  double chi2_scattertheta, chi2_scattermomentum, chi2_recoiltheta, chi2_recoilmomentum, 
    chi2_neutronmomentum, chi2_dpt, chi2_dphit, chi2_dalphat;
  
  TString chi2_name_lst[] = {"scattertheta", "scattermomentum", "recoiltheta", "recoilmomentum", 
    "neutronmomentum", "dpt", "dphit", "dalphat"};

  TString leghead = "GIBUU", leghead_ori = "GIBUU";
  if(anaid == 1)
  {
    chi2_scattertheta = GetChi2_stack(stk_mu_theta, hh_0pi_scattertheta, 
      hh_0pi_scattertheta_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_scattertheta, hh_0pi_scattertheta->GetNbinsX()-1);
    leg_mu_theta->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_scattermomentum = GetChi2_stack(stk_mu_mom, hh_0pi_scattermomentum, 
      hh_0pi_scattermomentum_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_scattermomentum, hh_0pi_scattermomentum->GetNbinsX()-1);
    leg_mu_mom->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_recoiltheta = GetChi2_stack(stk_recoil_theta, hh_0pi_recoiltheta, 
      hh_0pi_recoiltheta_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_recoiltheta, hh_0pi_recoiltheta->GetNbinsX()-1);
    leg_recoil_theta->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_recoilmomentum = GetChi2_stack(stk_recoil_mom, hh_0pi_recoilmomentum, 
      hh_0pi_recoilmomentum_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_recoilmomentum, hh_0pi_recoilmomentum->GetNbinsX()-1);
    leg_recoil_mom->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_neutronmomentum = GetChi2_stack(stk_neutron_mom, hh_0pi_neutronmomentum, 
      hh_0pi_neutronmomentum_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_neutronmomentum, hh_0pi_neutronmomentum->GetNbinsX()-1);
    leg_neutron_mom->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_dpt = GetChi2_stack(stk_dpt, hh_0pi_dpt, 
      hh_0pi_dpt_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_dpt, hh_0pi_dpt->GetNbinsX()-1);
    leg_dpt->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_dphit = GetChi2_stack(stk_dphit, hh_0pi_dphit, 
      hh_0pi_dphit_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_dphit, hh_0pi_dphit->GetNbinsX()-1);
    leg_dphit->SetHeader(leghead);
    leghead = leghead_ori;
    chi2_dalphat = GetChi2_stack(stk_dalphat, hh_0pi_dalphat, 
      hh_0pi_dalphat_cov, noff, dummyunit);
    leghead += Form("  #chi^{2}/ndf=%.0f/%i", chi2_dalphat, hh_0pi_dalphat->GetNbinsX()-1);
    leg_dalphat->SetHeader(leghead);
    //cout<<"chi2 - "<<chi2_scattertheta<<" chi2 - "<<chi2_dalphat<<endl;
  }
  else
  {
    chi2_dalphat = GetChi2_stack(stk_dalphat, hh_pi_dalphat, 
      hh_pi_dalphat_cov, noff, dummyunit);
  }

  double chi2_lst[] = {chi2_scattertheta, chi2_scattermomentum, chi2_recoiltheta, chi2_recoilmomentum, 
    chi2_neutronmomentum, chi2_dpt, chi2_dphit, chi2_dalphat};
  for(int ii=0; ii<8; ii++){cout<<"chi2 of "<<chi2_name_lst[ii]<<" : "<<chi2_lst[ii]<<endl;}

  const TString opt_same="same hist";
  const TString opt_err_bar="X1 E0";
  const TString opt_err_bar_same="same X1 E0";
  const TString opt_nostack="hist nostack";

  stk_x->Draw(opt_nostack);
  stk_x->GetXaxis()->SetTitle("x_{Bj}");
  stk_x->GetYaxis()->SetTitle("#frac{d#sigma}{dx_{Bj}} (cm^{2}/x_{Bj}/nucleon)");
  stk_x->SetTitle("x_{Bj} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_x->Draw();
  c1->Print("plot.pdf","x_{Bj} of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "x_{Bj} of 4 event modes"));

  stk_Q2->Draw(opt_nostack);
  stk_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  stk_Q2->GetYaxis()->SetTitle("#frac{d#sigma}{dQ^{2}}(cm^{2}/GeV^{2}/nucleon)");
  stk_Q2->SetTitle("Q^{2} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_Q2->Draw();
  c1->Print("plot.pdf","Q^{2} of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "Q^{2} of 4 event modes"));

  stk_W->Draw(opt_nostack);
  stk_W->GetXaxis()->SetTitle("W(GeV)");
  stk_W->GetYaxis()->SetTitle("#frac{d#sigma}{dW}(cm^{2}/GeV/nucleon)");
  stk_W->SetTitle("W of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_W->Draw();
  c1->Print("plot.pdf","W of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "W of 4 event modes"));

  if(nodata == 1){stk_mu_theta->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_scattertheta->SetAxisRange(0,3.0e-40,"Y");
    hh_0pi_scattertheta->Draw(opt_err_bar);
    hh_0pi_scattertheta->SetTitle("#theta of #mu of 4 event modes");
    stk_mu_theta->Draw(opt_same);
    hh_0pi_scattertheta->Draw(opt_err_bar_same);
  }
  else{stk_mu_theta->Draw(opt_nostack);}
  stk_mu_theta->GetXaxis()->SetTitle("#theta of #mu(degree)");
  stk_mu_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_mu_theta->SetTitle("#theta of #mu of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_mu_theta->Draw();
  c1->Print("plot.pdf","#theta of #mu of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of #mu of 4 event modes"));

  if(nodata == 1){stk_mu_mom->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_scattermomentum->Draw(opt_err_bar);
    hh_0pi_scattermomentum->SetTitle("momentum of #mu of 4 event modes");
    stk_mu_mom->Draw(opt_same);
    hh_0pi_scattermomentum->Draw(opt_err_bar_same);
  }
  else{stk_mu_mom->Draw(opt_nostack);}
  stk_mu_mom->GetXaxis()->SetTitle("momentum of #mu(GeV)");
  stk_mu_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}}(cm^{2}/GeV/nucleon)");
  stk_mu_mom->SetTitle("momentum of #mu of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_mu_mom->Draw();
  c1->Print("plot.pdf","momentum of #mu of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of #mu of 4 event modes"));

  stk_pi_theta->Draw(opt_nostack);
  stk_pi_theta->GetXaxis()->SetTitle("#theta of #pi(degree)");
  stk_pi_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_pi_theta->SetTitle("#theta of #pi of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_pi_theta->Draw();
  c1->Print("plot.pdf","#theta of #pi of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of #pi of 4 event modes"));

  stk_pi_mom->Draw(opt_nostack);
  stk_pi_mom->GetXaxis()->SetTitle("momentum of #pi(GeV)");
  stk_pi_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#pi}}(cm^{2}/GeV/nucleon)");
  stk_pi_mom->SetTitle("momentum of #pi of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_pi_mom->Draw();
  c1->Print("plot.pdf","momentum of #pi of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of #pi of 4 event modes"));

  if(nodata == 1){stk_recoil_theta->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_recoiltheta->Draw(opt_err_bar);
    hh_0pi_recoiltheta->SetTitle("#theta of recoil of 4 event modes");
    stk_recoil_theta->Draw(opt_same);
    hh_0pi_recoiltheta->Draw(opt_err_bar_same);
  }
  else{stk_recoil_theta->Draw(opt_nostack);}
  stk_recoil_theta->GetXaxis()->SetTitle("#theta of recoil(degree)");
  stk_recoil_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_recoil_theta->SetTitle("#theta of recoil of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_recoil_theta->Draw();
  c1->Print("plot.pdf","#theta of recoil of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of recoil of 4 event modes"));

  if(nodata == 1){stk_recoil_mom->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_recoilmomentum->Draw(opt_err_bar);
    hh_0pi_recoilmomentum->SetTitle("momentum of recoil of 4 event modes");
    stk_recoil_mom->Draw(opt_same);
    hh_0pi_recoilmomentum->Draw(opt_err_bar_same);
  }
  else{stk_recoil_mom->Draw(opt_nostack);}
  stk_recoil_mom->GetXaxis()->SetTitle("momentum of recoil(GeV)");
  stk_recoil_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{recoil}}(cm^{2}/GeV/nucleon)");
  stk_recoil_mom->SetTitle("momentum of recoil of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_recoil_mom->Draw();
  c1->Print("plot.pdf","momentum of recoil of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of recoil of 4 event modes"));

  stk_baryon_theta->Draw(opt_nostack);
  stk_baryon_theta->GetXaxis()->SetTitle("#theta of baryon(degree)");
  stk_baryon_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_baryon_theta->SetTitle("#theta of baryon of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_baryon_theta->Draw();
  c1->Print("plot.pdf","#theta of baryon of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of baryon of 4 event modes"));

  stk_baryon_mom->Draw(opt_nostack);
  stk_baryon_mom->GetXaxis()->SetTitle("momentum of baryon(GeV)");
  stk_baryon_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{baryon}}(cm^{2}/GeV/nucleon)");
  stk_baryon_mom->SetTitle("momentum of baryon of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_baryon_mom->Draw();
  c1->Print("plot.pdf","momentum of baryon of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of baryon of 4 event modes"));
  
  if(nodata == 1){stk_neutron_mom->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_neutronmomentum->Draw(opt_err_bar);
    hh_0pi_neutronmomentum->SetTitle("p_{n} of 4 event modes");
    stk_neutron_mom->Draw(opt_same);
    hh_0pi_neutronmomentum->Draw(opt_err_bar_same);
  }
  else{stk_neutron_mom->Draw(opt_nostack);}
  stk_neutron_mom->GetXaxis()->SetTitle("p_{n}(GeV)");
  stk_neutron_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{n}}(cm^{2}/GeV/nucleon)");
  stk_neutron_mom->SetTitle("p_{n} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_neutron_mom->Draw();
  c1->Print("plot.pdf","p_{n} of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "p_{n} of 4 event modes"));
  
  if(nodata == 1){stk_dpt->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_dpt->Draw(opt_err_bar);
    hh_0pi_dpt->SetTitle("dpt of 4 event modes");
    stk_dpt->Draw(opt_same);
    hh_0pi_dpt->Draw(opt_err_bar_same);
  }
  else{  stk_dpt->Draw(opt_nostack);}
  stk_dpt->GetXaxis()->SetTitle("dpt(GeV)");
  stk_dpt->GetYaxis()->SetTitle("#frac{d#sigma}{d dpt}(cm^{2}/GeV/nucleon)");
  stk_dpt->SetTitle("dpt of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dpt->Draw();
  c1->Print("plot.pdf","dpt of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dpt of 4 event modes"));

  if(nodata == 1){stk_dphit->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_dphit->Draw(opt_err_bar);
    hh_0pi_dphit->SetTitle("dphit of 4 event modes");
    stk_dphit->Draw(opt_same);
    hh_0pi_dphit->Draw(opt_err_bar_same);
  }
  else{stk_dphit->Draw(opt_nostack);}
  stk_dphit->GetXaxis()->SetTitle("dphit(degree)");
  stk_dphit->GetYaxis()->SetTitle("#frac{d#sigma}{d dphit}(cm^{2}/GeV/nucleon)");
  stk_dphit->SetTitle("dphit of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dphit->Draw();
  c1->Print("plot.pdf","dphit of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dphit of 4 event modes"));
  
  
  if(nodata == 1){stk_dalphat->Draw(opt_nostack);}
  else if(anaid == 1)
  {
    hh_0pi_dalphat->Draw(opt_err_bar);
    hh_0pi_dalphat->SetTitle("dalphat of 4 event modes");
    stk_dalphat->Draw(opt_same);
    hh_0pi_dalphat->Draw(opt_err_bar_same);
  }
  else
  {
    hh_pi_dalphat->Draw(opt_err_bar);
    hh_pi_dalphat->SetTitle("dalphat of 4 event modes");
    stk_dalphat->Draw(opt_same);
    hh_pi_dalphat->Draw(opt_err_bar_same);
  }
  stk_dalphat->GetXaxis()->SetTitle("dalphat(degree)");
  stk_dalphat->GetYaxis()->SetTitle("#frac{d#sigma}{d dalphat}(cm^{2}/GeV/nucleon)");
  stk_dalphat->SetTitle("dalphat of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dalphat->Draw();
  c1->Print("plot.pdf","dalphat of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dalphat of 4 event modes"));
  
  c1->Print("plot.pdf]","dalphat of 4 event modess");
  fout->Close();
}
