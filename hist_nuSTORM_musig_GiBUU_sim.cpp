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

void hist_nuSTORM_musig_GiBUU_sim()
{
/*
  TString fn{"E7_outAna9_GFS0PIa9nuC_Filelist_GiBUU.root"};
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"plotting outana9 - case with no pi 0"<<endl;
  */
  TString fn{"E7_outAna7_GFSPIZEROa7nuC_Filelist_GiBUU.root"};
  int anaid{2};
  int noff{1};
  cout<<"plotting outana7 - case with pi 0"<<endl;
  auto anatag{"outana7"};
  

  const double dummyunit = 1E-42; //data

  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot.root","recreate");
  TTree *t1 = (TTree*) f->Get("tree");
  TCanvas *c1 = new TCanvas("c","",800,600);
  c1->Print("plot.pdf[");

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*13;//Scaling to per nucleon in CH 13=1+12
  
  const TString modes[]={"_all","_E7","_res","_dis","_2p2h", "_other"};
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

  cout<<"init hist in process"<<endl;

  TH2D *h_weight_E7_scattertheta = new TH2D("h_weight_E7_scattertheta","distribution of weight and scattertheta mode : E7",80,0,26,100,0.0000001,0.001);

  TH2D *hxQ2_E7 = new TH2D("hxQ2_E7","distribution of x_{Bj} and Q^{2} mode : E7",60,0,2,80,0.1,10);
  TH2D *hxW_E7 = new TH2D("hxW_E7","distribution of x_{Bj} and W mode : E7",60,0,2,80,0.6,3);
  TH1D *hxQ2_E7_x = new TH1D("hxQ2_E7_x","x_{Bj} mode : E7",80,0,2);
  TH1D *hxQ2_E7_Q2 = new TH1D("hxQ2_E7_Q2","Q^{2} mode : E7",80,0,5);
  TH1D *hxQ2_E7_W = new TH1D("hxQ2_E7_W","W mode : E7",80,0,4);


  TH1D *hxQ2_E7_scattertheta = new TH1D("hxQ2_E7_scattertheta","#mu #theta mode : E7",80,0,26); //muon theta
  TH1D *hxQ2_E7_scattermomentum = new TH1D("hxQ2_E7_scattermomentum","#mu p mode : E7",80,1,15);
  TH1D *hxQ2_E7_mesontheta = new TH1D("hxQ2_E7_mesontheta","#pi #theta mode : E7",80,0,180); //muon theta
  TH1D *hxQ2_E7_mesonmomentum = new TH1D("hxQ2_E7_mesonmomentum","#pi p mode : E7",80,0,7);
  TH1D *hxQ2_E7_recoiltheta = new TH1D("hxQ2_E7_recoiltheta","recoil #theta mode : E7",80,0,180); //muon theta
  TH1D *hxQ2_E7_recoilmomentum = new TH1D("hxQ2_E7_recoilmomentum","recoil p mode : E7",80,0,6);
  TH1D *hxQ2_E7_baryontheta = new TH1D("hxQ2_E7_baryontheta","baryon #theta mode : E7",80,0,180); //muon theta
  TH1D *hxQ2_E7_baryonmomentum = new TH1D("hxQ2_E7_baryonmomentum","baryon p mode : E7",80,0,10);
  TH1D *hxQ2_E7_neutronmomentum = new TH1D("hxQ2_E7_neutronmomentum","neutron p mode : E7",80,0,2);
  TH1D *hxQ2_E7_dpt = new TH1D("hxQ2_E7_dpt","dpt mode : E7",80,0,2.5);
  TH1D *hxQ2_E7_dphit = new TH1D("hxQ2_E7_dphit","#delta#phi_{t} mode : E7",80,0,180);
  TH1D *hxQ2_E7_dalphat = new TH1D("hxQ2_E7_dalphat","#delta#alpha_{t} mode : E7",80,0,180);

   const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange};
  TH2D *hh2d[]={hxQ2_E7, hxW_E7}; //list of 2d hist
  TH1D *hh1d[]={hxQ2_E7_x, hxQ2_E7_Q2, hxQ2_E7_W}; //list of 1d hist
    
  TH1D *hh1d_mom[]={hxQ2_E7_scattertheta, hxQ2_E7_scattermomentum,  hxQ2_E7_mesontheta, hxQ2_E7_mesonmomentum, 
    hxQ2_E7_recoiltheta, hxQ2_E7_recoilmomentum,  hxQ2_E7_baryontheta, hxQ2_E7_baryonmomentum, 
    hxQ2_E7_neutronmomentum, hxQ2_E7_dpt, hxQ2_E7_dphit, hxQ2_E7_dalphat};

  TH1D *hh1d_TKI[]={hxQ2_E7_dpt, hxQ2_E7_dphit, hxQ2_E7_dalphat};
  string mode[] = {"E7", "RES", "DIS", "2p2h"};
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
 
 
  Int_t nentries = (Int_t)t1->GetEntries();
  const int nh1d = sizeof(hh1d)/sizeof(TH1D*);
  const int nh1d_mom = sizeof(hh1d_mom)/sizeof(TH1D*);
  cout<<"mom array len : "<<nh1d_mom<<endl;
  cout<<"events : "<<nentries<<endl;
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);
  int npi_err{0};
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    //if(perweight<0 && i%500==0){cout<<"weight : "<<perweight<<endl;}
    if(i%50000==0){cout<<"meson(pi) p : "<<mesonmomentum<<endl;}
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

    hxQ2_E7->Fill(xBj,Q2,perweight);
    hxW_E7->Fill(xBj,Wtrue,perweight);
    hxQ2_E7_x->Fill(xBj,perweight);
    hxQ2_E7_Q2->Fill(Q2,perweight);
    hxQ2_E7_W->Fill(Wtrue,perweight);
    hxQ2_E7_scattertheta->Fill(scattertheta,perweight);
    hxQ2_E7_scattermomentum->Fill(scattermomentum,perweight);
    hxQ2_E7_mesontheta->Fill(mesontheta,perweight);
    hxQ2_E7_mesonmomentum->Fill(mesonmomentum,perweight);
    hxQ2_E7_recoiltheta->Fill(recoiltheta,perweight);
    hxQ2_E7_recoilmomentum->Fill(recoilmomentum,perweight);
    hxQ2_E7_baryontheta->Fill(baryontheta,perweight);
    hxQ2_E7_baryonmomentum->Fill(baryonmomentum,perweight);
    hxQ2_E7_neutronmomentum->Fill(IApN,perweight);
    hxQ2_E7_dpt->Fill(dpt,perweight);
    hxQ2_E7_dphit->Fill(dphit,perweight);
    hxQ2_E7_dalphat->Fill(dalphat,perweight);

    h_weight_E7_scattertheta->Fill(scattertheta,perweight);
  }
  cout<<" n pi err # : "<<npi_err<<endl;

  Double_t temp_xsec{0};

  Double_t E7_x_xsec, E7_Q2_xsec, E7_W_xsec, E7_scattertheta_xsec, E7_scattermomentum_xsec, 
    E7_mesontheta_xsec, E7_mesonmomentum_xsec, E7_recoiltheta_xsec, E7_recoilmomentum_xsec, 
    E7_baryontheta_xsec, E7_baryonmomentum_xsec, E7_neutronmomentum_xsec, E7_dpt_xsec, E7_dphit_xsec, E7_dalphat_xsec;

  Double_t xsec_lst[] = {E7_x_xsec, E7_Q2_xsec, E7_W_xsec, E7_scattertheta_xsec, E7_scattermomentum_xsec, 
    E7_mesontheta_xsec, E7_mesonmomentum_xsec, E7_recoiltheta_xsec, E7_recoilmomentum_xsec, 
    E7_baryontheta_xsec, E7_baryonmomentum_xsec, E7_neutronmomentum_xsec, E7_dpt_xsec, E7_dphit_xsec, E7_dalphat_xsec};
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
  remove_neg(h_weight_E7_scattertheta);
  h_weight_E7_scattertheta->GetXaxis()->SetTitle("E7 #theta_{mu}");
  h_weight_E7_scattertheta->GetYaxis()->SetTitle("weight");
  gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
  gPad->Update();
  h_weight_E7_scattertheta->SetContour(1000);
  h_weight_E7_scattertheta->Draw("colz");
  gStyle->SetPadRightMargin(0.2);
  c1->Print("plot.pdf",h_weight_E7_scattertheta->GetTitle());
  c1->Print(Form("png/%s_%s.png", anatag, h_weight_E7_scattertheta->GetTitle()));
  */
  TString tit[]={"evtMode E7 x_{Bj}-Q^{2}", "evtMode E7 x_{Bj}-W","E7 x_{Bj}", "E7 Q^{2}(GeV^{2})", "E7 W(GeV)"};
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
    cout<<"xsec calc with : "<<var[ii%(nxsec)]<<" is "<<xsec_lst[ii]<<" cm^{2}"<<endl;
  }
  cout<<endl;
  Double_t comb_xsec{0};
  /*
  for(int ii=0; ii<(nxsec); ii++)
  {
    for(int jj=0; jj<4; jj++)
    {
      comb_xsec = comb_xsec+xsec_lst[jj*(nxsec/4)+ii];
    }
    cout<<"combind total xsec cal with "<<var[ii]<<" : "<<comb_xsec<<" cm^{2}"<<endl;
    comb_xsec = 0;
  }
  */

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

  stk_mu_theta->Draw(opt_nostack);
  stk_mu_theta->GetXaxis()->SetTitle("#theta of #mu(degree)");
  stk_mu_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_mu_theta->SetTitle("#theta of #mu of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_mu_theta->Draw();
  c1->Print("plot.pdf","#theta of #mu of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of #mu of 4 event modes"));

  stk_mu_mom->Draw(opt_nostack);
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

  stk_recoil_theta->Draw(opt_nostack);
  stk_recoil_theta->GetXaxis()->SetTitle("#theta of recoil(degree)");
  stk_recoil_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_recoil_theta->SetTitle("#theta of recoil of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_recoil_theta->Draw();
  c1->Print("plot.pdf","#theta of recoil of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of recoil of 4 event modes"));

  stk_recoil_mom->Draw(opt_nostack);
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
  
  stk_neutron_mom->Draw(opt_nostack);
  stk_neutron_mom->GetXaxis()->SetTitle("p_{n}(GeV)");
  stk_neutron_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{n}}(cm^{2}/GeV/nucleon)");
  stk_neutron_mom->SetTitle("p_{n} of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_neutron_mom->Draw();
  c1->Print("plot.pdf","p_{n} of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "p_{n} of 4 event modes"));
  
  stk_dpt->Draw(opt_nostack);
  stk_dpt->GetXaxis()->SetTitle("dpt(GeV)");
  stk_dpt->GetYaxis()->SetTitle("#frac{d#sigma}{d dpt}(cm^{2}/GeV/nucleon)");
  stk_dpt->SetTitle("dpt of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dpt->Draw();
  c1->Print("plot.pdf","dpt of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dpt of 4 event modes"));

  stk_dphit->Draw(opt_nostack);
  stk_dphit->GetXaxis()->SetTitle("dphit(degree)");
  stk_dphit->GetYaxis()->SetTitle("#frac{d#sigma}{d dphit}(cm^{2}/GeV/nucleon)");
  stk_dphit->SetTitle("dphit of 4 event modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dphit->Draw();
  c1->Print("plot.pdf","dphit of 4 event modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dphit of 4 event modes"));
  
  
  stk_dalphat->Draw(opt_nostack);
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
