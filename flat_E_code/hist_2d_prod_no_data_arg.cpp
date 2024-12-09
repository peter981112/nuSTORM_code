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
#include "TChain.h"

#include <iostream>
#include <string>

#include "HistIO.h"
#include "GeneratorUtils.h"
#include "style.h"

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
  Int_t x_bin_num;
  x_bin_num = hist->GetNbinsX();
  for (Int_t xx=1; xx<=x_bin_num; xx++){hist->SetBinContent(xx,0);}
  return;
}

void diff_xsec(TH1D* &hist, const double histnormFactor) //makes all bin contents to 0
{
  Int_t x_bin_num;
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
  //TAxis * axis = data_hist->GetXaxis();
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

void hist_2d_prod_no_data(TString fn, string root_file_path)
{

  TFile *f = new TFile(fn);
  TTree *t1 = (TTree*) f->Get("tree");
  int anaid{0};auto anatag{"outana"};string fn_string{fn};TString pdf_out_name{""};//int noff{0}

  if(fn.Contains("outAna7")){
    anaid = 2;
    //noff = 1;
    anatag = "outana7";
    cout<<"rewighting "<<anatag<<" - case with pi 0"<<endl;
  }
  else if(fn.Contains("outAna9")){
    anaid = 1;
    //noff = 2;
    anatag = "outana9";
    cout<<"rewighting "<<anatag<<" - case with no pi 0"<<endl;
  }

  cout<<"fn_string : "<<fn_string<<endl;
  size_t dot_index = fn_string.find(".");
  fn_string.erase(fn_string.begin() + dot_index, fn_string.end());
  cout<<"fn_string : "<<fn_string<<endl;
  size_t slash_index = fn_string.find("/");
  fn_string.erase(fn_string.begin(), fn_string.begin()+slash_index+1);//remove file path and extention .root from the file name
  slash_index = fn_string.find("/");
  fn_string.erase(fn_string.begin(), fn_string.begin()+slash_index+1);//remove file path and extention .root from the file name
  
  pdf_out_name = root_file_path + fn_string + "_plot.pdf";
  cout<<"fn_string : "<<fn_string<<" pdf_out_name : "<<pdf_out_name<<endl;

  TFile *fout = new TFile("tmpplot.root","recreate");
  TCanvas *c1 = new TCanvas("c","",800,600);
  TString start_pdf_out_name = pdf_out_name + "[";
  c1->Print(start_pdf_out_name);

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*13;//Scaling to per nucleon in CH 13=1+12
  
  const TString modes[]={"_all","_qe","_res","_dis","_2p2h", "_other"};

  ULong64_t minnp, maxnp;
  const bool isTopoTask = GeneratorUtils::SetTopoCut(anaid, minnp, maxnp);

  const TString npiDecomp[]={"","_0pi","_1pi","_Mpi","_other"};
  //only gointo sub-category for TopoTask
  const int nNpiCut = isTopoTask? (sizeof(npiDecomp)/sizeof(TString)) : 1;
  cout<<"n pi cut : "<<nNpiCut<<endl;

  TLegend *leg = new TLegend(0.65, 0.45, 0.85, 0.85);
  leg->SetFillStyle(0);
  //leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.023);
  gStyle-> SetOptStat(0); // turn off stat box
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.1);
  gStyle-> SetPalette(kRainbow);
  Double_t xBj, Q2, Wtrue, beamE, scattertheta, scattermomentum, mesontheta, mesonmomentum, 
    recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat, dpTT, perweight;
  Int_t prod,evtMode;
  t1->SetBranchAddress("xBj",&xBj);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("Wtrue",&Wtrue);
  t1->SetBranchAddress("beamE",&beamE);
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
  
  TString E_modes[]={"QE","RES","DIS","2p2h"};
  TString var[] = {"x", "Q2", "W", "beamE", "scattertheta", "scattermomentum", 
      "mesontheta", "mesonmomentum", "recoiltheta", "recoilmomentum", "baryontheta", "baryonmomentum",
      "neutronmomentum", "dp_{t}", "d#phi_{t}", "d#alpha_{t}"};
  //const Int_t nmode = sizeof(E_modes)/sizeof(TString);
  const Int_t n_E = sizeof(E_modes)/sizeof(TString);

  const int nh2d_temp = 2;
  const int nh1d_temp = 4;//sizeof(hh1d_temp)/sizeof(TH1D*);
  const int nh1d_mom_temp = 12;//sizeof(hh1d_mom_temp)/sizeof(TH1D*);

  TString temp_name = "";
  TString temp_tit = "";
  TH2D *hh2d[nh2d_temp*n_E]; //list of 2d hist {hxQ2_E7, hxW_E7}
  TString hh2d_name_lst[nh2d_temp*n_E];
  TString hh2d_tit_lst[nh2d_temp*n_E];
  for(int i=0; i<nh2d_temp*n_E; i++)
  {
    if(i%2 == 0)
    {
      temp_name = temp_name+"hxQ2 "+ E_modes[i/2];
      temp_tit = temp_tit+"x-Q^{2} for "+ E_modes[i/2];
    }
    else
    {
      temp_name = temp_name+"hxW "+ E_modes[i/2];
      temp_tit = temp_tit+"x-W_for "+ E_modes[i/2];
    }
    hh2d_name_lst[i] = temp_name;
    hh2d_tit_lst[i] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }

  TH1D *hh1d[nh1d_temp*n_E]; //list of 1d hist {hxQ2_E7_x, hxQ2_E7_Q2, hxQ2_E7_W, hxQ2_E7_beamE}
  TString hh1d_name_lst[nh1d_temp*n_E];
  TString hh1d_tit_lst[nh1d_temp*n_E];
  for(int i=0; i<nh1d_temp*n_E; i++)
  {
    temp_name = temp_name+"hxQ2 "+var[i%nh1d_temp]+ E_modes[i/nh1d_temp];
    temp_tit = temp_tit+var[i%nh1d_temp]+" for "+ E_modes[i/nh1d_temp];

    hh1d_name_lst[i] = temp_name;
    hh1d_tit_lst[i] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }
    
  TH1D *hh1d_mom[12*n_E];//={hxQ2_E7_scattertheta, hxQ2_E7_scattermomentum,  hxQ2_E7_mesontheta, hxQ2_E7_mesonmomentum, 
    //hxQ2_E7_recoiltheta, hxQ2_E7_recoilmomentum,  hxQ2_E7_baryontheta, hxQ2_E7_baryonmomentum, 
    //hxQ2_E7_neutronmomentum, hxQ2_E7_dpt, hxQ2_E7_dphit, hxQ2_E7_dalphat};
  TString hh1d_mom_name_lst[nh1d_mom_temp*n_E];
  TString hh1d_mom_tit_lst[nh1d_mom_temp*n_E];
  for(int i=0; i<nh1d_mom_temp*n_E; i++)
  {
    temp_name = temp_name+"hxQ2_"+var[i%nh1d_mom_temp+nh1d_temp]+ E_modes[i/nh1d_mom_temp];
    temp_tit = temp_tit+var[i%nh1d_mom_temp+nh1d_temp]+" for "+ E_modes[i/nh1d_mom_temp];

    hh1d_mom_name_lst[i] = temp_name;
    hh1d_mom_tit_lst[i] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }


  for(Int_t kk=0; kk<n_E; kk++)
  {

    Double_t hh2d_hist_x_upper_bound[] = {2,1,1,2};
    
    hh2d[kk*2] = new TH2D(hh2d_name_lst[kk*2],hh2d_tit_lst[kk*2],60,0,hh2d_hist_x_upper_bound[kk],80,0.1,10);
    hh2d[kk*2+1] = new TH2D(hh2d_name_lst[kk*2+1],hh2d_tit_lst[kk*2+1],60,0,hh2d_hist_x_upper_bound[kk],80,0.6,3);
    
    /*
    variables
    "x", "Q2", "W",
    Double_t hh1d_hist_lower_bound[] = {0,0,0.5};
    Double_t hh1d_hist_upper_bound[] = {2,4,3.5};

    "scattertheta", "scattermomentum", "mesontheta", "mesonmomentum", "recoiltheta", "recoilmomentum", "baryontheta", "baryonmomentum",
    "neutronmomentum", "dp_{t}", "d#phi_{t}", "d#alpha_{t}"};
    Double_t hh1d_mom_hist_lower_bound[] = {0,1,0,0,0,0,0,0,0,0,0,0};
    Double_t hh1d_mom_hist_upper_bound[] = {26,5.5,180,3.5,180,4.5,180,5,2,2.5,180,180};
    */

    Double_t hh1d_hist_lower_bound[] = {0,0,0.5,0};
    Double_t hh1d_hist_upper_bound[] = {2.1,2.5,3.2,8};
    for(Int_t ii=0; ii<nh1d_temp; ii++)
    {
      hh1d[kk*nh1d_temp+ii] = new TH1D(hh1d_name_lst[kk*nh1d_temp+ii],hh1d_tit_lst[kk*nh1d_temp+ii],60,hh1d_hist_lower_bound[ii],hh1d_hist_upper_bound[ii]);   
    }
    Double_t hh1d_mom_hist_lower_bound[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    Double_t hh1d_mom_hist_upper_bound[] = {180,5.5,180,3,180,3.5,180,5,1.7,1.6,180,180};
    
    for(Int_t ii=0; ii<nh1d_mom_temp; ii++)
    {
      hh1d_mom[kk*nh1d_mom_temp+ii] = new TH1D(hh1d_mom_name_lst[kk*nh1d_mom_temp+ii],hh1d_mom_tit_lst[kk*nh1d_mom_temp+ii],60,hh1d_mom_hist_lower_bound[ii],hh1d_mom_hist_upper_bound[ii]);
    }
  }

  TString tit[]={"evtMode QE x_{Bj}-Q^{2}", "evtMode QE x_{Bj}-W","QE x_{Bj}", "QE Q^{2}(GeV^{2})", "QE W(GeV)"
  , "evtMode RES x_{Bj}-Q^{2}", "evtMode RES x_{Bj}-W","RES x_{Bj}", "RES Q^{2}(GeV^{2})", "RES W(GeV)"
  , "evtMode DIS x_{Bj}-Q^{2}", "evtMode DIS x_{Bj}-W","DIS x_{Bj}", "DIS Q^{2}(GeV^{2})", "DIS W(GeV)"
  , "evtMode 2p2h x_{Bj}-Q^{2}","evtMode 2p2h x_{Bj}-W","2p2h x_{Bj}", "2p2h Q^{2}(GeV^{2})", "2p2h W(GeV)"};
  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange};
  string mode[] = {"QE", "RES", "DIS", "2p2h"};

  TH2D *h_weight_QE_scattertheta = new TH2D("h_weight_QE_scattertheta","distribution of weight and scattertheta mode : QE",80,0,86,100,0.0000001,0.001);

  THStack *stk_x = new THStack;
  TLegend *leg_x = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_x->SetFillStyle(0);
  
  THStack *stk_Q2 = new THStack;
  TLegend *leg_Q2 = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_Q2->SetFillStyle(0);
  
  THStack *stk_W = new THStack;
  TLegend *leg_W = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_W->SetFillStyle(0);

  THStack *stk_beamE = new THStack;
  TLegend *leg_beamE = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_beamE->SetFillStyle(0);
  
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
  cout<<"number of events : "<<nentries<<endl;
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);
  int npi_err{0};
  double QE_xsec{0}, RES_xsec{0}, DIS_xsec{0}, a2p2h_xsec{0}, tot_xsec{0};
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    Double_t hh1d_var_lst[] = {xBj, Q2, Wtrue, beamE};
    Double_t hh1d_mom_var_lst[] = {scattertheta, scattermomentum, mesontheta, mesonmomentum, recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat};
        
    //if(perweight<0 && i%500==0){cout<<"weight : "<<perweight<<endl;}
    if(i%1000000==0){cout<<i<<" th meson(pi) p : "<<mesonmomentum<<endl;}
    if(perweight>40){
      if(i%100==0){printf("\nAlert!!  Filling Super weight!! %d %f skiping...\n", i, perweight);}
    }
    else
    {
      if(anaid == CC1piNpID){
        if(targetZ == 1){//no hydrogen for GFS1
          cout<<"no hydrogen for GFS1"<<endl;
        }
      }
      if(anaid == 1){if(GeneratorUtils::GetNpi(npar)==0){npi_err++;}}
      else{if(!(GeneratorUtils::GetNpi(npar)==0)){npi_err++;}}

      int kk{0};
      
      tot_xsec = tot_xsec + perweight;

      if(prod==1 && evtMode==1) //fil histogram bins for quasi elastic case
      {
        kk = 0;
        hh2d[kk*2]->Fill(xBj,Q2,perweight);
        hh2d[kk*2+1]->Fill(xBj,Wtrue,perweight);
        QE_xsec = QE_xsec + perweight;
        for(Int_t ii=0; ii<4; ii++)
        {
          hh1d[kk*nh1d_temp+ii]->Fill(hh1d_var_lst[ii],perweight);
        }

        for(Int_t ii=0; ii<12; ii++)
        {
          hh1d_mom[kk*nh1d_mom_temp+ii]->Fill(hh1d_mom_var_lst[ii],perweight);;
        }
        
        h_weight_QE_scattertheta->Fill(scattertheta,perweight);
      }
      else if(prod>=2 && prod<=33 && evtMode==2) //fil histogram bins for resonance case
      {
        kk = 1;
        hh2d[kk*2]->Fill(xBj,Q2,perweight);
        hh2d[kk*2+1]->Fill(xBj,Wtrue,perweight);
        RES_xsec = RES_xsec + perweight;
        for(Int_t ii=0; ii<4; ii++)
        {
          hh1d[kk*nh1d_temp+ii]->Fill(hh1d_var_lst[ii],perweight);;
        }

        for(Int_t ii=0; ii<12; ii++)
        {
          hh1d_mom[kk*nh1d_mom_temp+ii]->Fill(hh1d_mom_var_lst[ii],perweight);;
        }
      }
      else if(prod==34 && evtMode==3) //fil histogram bins for deep inelastic case
      {
        kk = 2;
        hh2d[kk*2]->Fill(xBj,Q2,perweight);
        hh2d[kk*2+1]->Fill(xBj,Wtrue,perweight);
        DIS_xsec = DIS_xsec + perweight;
        for(Int_t ii=0; ii<4; ii++)
        {
          hh1d[kk*nh1d_temp+ii]->Fill(hh1d_var_lst[ii],perweight);;
        }

        for(Int_t ii=0; ii<12; ii++)
        {
          hh1d_mom[kk*nh1d_mom_temp+ii]->Fill(hh1d_mom_var_lst[ii],perweight);;
        }
      }
      else if(prod==35 && evtMode==4) //fil histogram bins for 2 part 2 hole case
      {
        kk = 3;
        hh2d[kk*2]->Fill(xBj,Q2,perweight);
        hh2d[kk*2+1]->Fill(xBj,Wtrue,perweight);
        a2p2h_xsec = a2p2h_xsec + perweight;
        for(Int_t ii=0; ii<4; ii++)
        {
          hh1d[kk*nh1d_temp+ii]->Fill(hh1d_var_lst[ii],perweight);;
        }

        for(Int_t ii=0; ii<12; ii++)
        {
          hh1d_mom[kk*nh1d_mom_temp+ii]->Fill(hh1d_mom_var_lst[ii],perweight);;
        }
      }
    }
  }
  cout<<" n pi err # : "<<npi_err<<endl;
  double xsec_lst_dir[] = {QE_xsec, RES_xsec, DIS_xsec, a2p2h_xsec};
  tot_xsec = tot_xsec/histnormFactor;
  cout<<"direct perweight calc total xsec : "<<tot_xsec<<" cm^{2}"<<endl;
  for(int ii=0; ii<4; ii++)
  {
    xsec_lst_dir[ii] = xsec_lst_dir[ii]/histnormFactor;
    cout<<"direct perweight calc xsec for mode : "<<mode[ii]<<" is "<<xsec_lst_dir[ii]<<" cm^{2}"<<endl;
  }

  Double_t temp_xsec{0};

  Double_t QE_x_xsec{0}, QE_Q2_xsec{0}, QE_W_xsec{0}, QE_beamE_xsec{0}, QE_scattertheta_xsec{0}, QE_scattermomentum_xsec{0}, 
    QE_mesontheta_xsec{0}, QE_mesonmomentum_xsec{0}, QE_recoiltheta_xsec{0}, QE_recoilmomentum_xsec{0}, 
    QE_baryontheta_xsec{0}, QE_baryonmomentum_xsec{0}, QE_neutronmomentum_xsec{0}, QE_dpt_xsec{0}, QE_dphit_xsec{0}, QE_dalphat_xsec{0}, 
    RES_x_xsec{0}, RES_Q2_xsec{0}, RES_W_xsec{0}, RES_beamE_xsec{0}, RES_scattertheta_xsec{0}, RES_scattermomentum_xsec{0}, 
    RES_mesontheta_xsec{0}, RES_mesonmomentum_xsec{0}, RES_recoiltheta_xsec{0}, RES_recoilmomentum_xsec{0}, 
    RES_baryontheta_xsec{0}, RES_baryonmomentum_xsec{0}, RES_neutronmomentum_xsec{0}, RES_dpt_xsec{0}, RES_dphit_xsec{0}, RES_dalphat_xsec{0}, 
    DIS_x_xsec{0}, DIS_Q2_xsec{0}, DIS_W_xsec{0}, DIS_beamE_xsec{0}, DIS_scattertheta_xsec{0}, DIS_scattermomentum_xsec{0}, 
    DIS_mesontheta_xsec{0}, DIS_mesonmomentum_xsec{0}, DIS_recoiltheta_xsec{0}, DIS_recoilmomentum_xsec{0}, 
    DIS_baryontheta_xsec{0}, DIS_baryonmomentum_xsec{0}, DIS_neutronmomentum_xsec{0}, DIS_dpt_xsec{0}, DIS_dphit_xsec{0}, DIS_dalphat_xsec{0}, 
    a2p2h_x_xsec{0}, a2p2h_Q2_xsec{0}, a2p2h_W_xsec{0}, a2p2h_beamE_xsec{0}, a2p2h_scattertheta_xsec{0}, a2p2h_scattermomentum_xsec{0}, 
    a2p2h_mesontheta_xsec{0}, a2p2h_mesonmomentum_xsec{0}, a2p2h_recoiltheta_xsec{0}, 
    a2p2h_recoilmomentum_xsec{0}, a2p2h_baryontheta_xsec{0}, a2p2h_baryonmomentum_xsec{0}, a2p2h_neutronmomentum_xsec{0}, 
    a2p2h_dpt_xsec{0}, a2p2h_dphit_xsec{0}, a2p2h_dalphat_xsec{0};

  Double_t xsec_lst[] = {QE_x_xsec, QE_Q2_xsec, QE_W_xsec, QE_beamE_xsec, QE_scattertheta_xsec, QE_scattermomentum_xsec, 
    QE_mesontheta_xsec, QE_mesonmomentum_xsec, QE_recoiltheta_xsec, QE_recoilmomentum_xsec, 
    QE_baryontheta_xsec, QE_baryonmomentum_xsec, QE_neutronmomentum_xsec, QE_dpt_xsec, QE_dphit_xsec, QE_dalphat_xsec, 
    RES_x_xsec, RES_Q2_xsec, RES_W_xsec, RES_beamE_xsec, RES_scattertheta_xsec, RES_scattermomentum_xsec, 
    RES_mesontheta_xsec, RES_mesonmomentum_xsec, RES_recoiltheta_xsec, RES_recoilmomentum_xsec, 
    RES_baryontheta_xsec, RES_baryonmomentum_xsec, RES_neutronmomentum_xsec, RES_dpt_xsec, RES_dphit_xsec, RES_dalphat_xsec, 
    DIS_x_xsec, DIS_Q2_xsec, DIS_W_xsec, DIS_beamE_xsec, DIS_scattertheta_xsec, DIS_scattermomentum_xsec, 
    DIS_mesontheta_xsec, DIS_mesonmomentum_xsec, DIS_recoiltheta_xsec, DIS_recoilmomentum_xsec, 
    DIS_baryontheta_xsec, DIS_baryonmomentum_xsec, DIS_neutronmomentum_xsec, DIS_dpt_xsec, DIS_dphit_xsec, DIS_dalphat_xsec, 
    a2p2h_x_xsec, a2p2h_Q2_xsec, a2p2h_W_xsec, a2p2h_beamE_xsec, a2p2h_scattertheta_xsec, a2p2h_scattermomentum_xsec, 
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
    c1->Print(pdf_out_name,hh2d[ii]->GetTitle());
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
    hh1d[ii]->SetLineColor(cols[ii/4]);
    hh1d[ii]->SetFillStyle(4050);
    hh1d[ii]->SetFillColorAlpha(cols[ii/4],0.35);

    temp_xsec = hh1d[ii]->Integral(0,10000,"width");
    xsec_lst[(ii/4)*15+ii%4] = temp_xsec;

    //TString temp_tit="";
    //temp_tit = tit[ii*4];
    if(ii%4==0){
      stk_x->Add(hh1d[ii]);
      leg_x->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }
    
    else if(ii%4==1){
      stk_Q2->Add(hh1d[ii]);
      leg_Q2->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }

    else if(ii%4==2){
      stk_W->Add(hh1d[ii]);
      leg_W->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }

    else{
      stk_beamE->Add(hh1d[ii]);
      leg_beamE->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");
    }
    ++aa;
  }
  //binwidth = 0;
  int count_num{12};//, stk_num{0};
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
    xsec_lst[(ii/count_num)*15+4+ii%count_num] = temp_xsec;
    switch (ii%count_num)
    {
      case 0:
        stk_mu_theta->Add(hh1d_mom[ii]);
        leg_mu_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;

      case 1:
        stk_mu_mom->Add(hh1d_mom[ii]);
        leg_mu_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;

      case 2:
        stk_pi_theta->Add(hh1d_mom[ii]);
        leg_pi_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 3:
        stk_pi_mom->Add(hh1d_mom[ii]);
        leg_pi_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 4:
        stk_recoil_theta->Add(hh1d_mom[ii]);
        leg_recoil_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 5:
        stk_recoil_mom->Add(hh1d_mom[ii]);
        leg_recoil_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 6:
        stk_baryon_theta->Add(hh1d_mom[ii]);
        leg_baryon_theta->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 7:
        stk_baryon_mom->Add(hh1d_mom[ii]);
        leg_baryon_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 8:
        stk_neutron_mom->Add(hh1d_mom[ii]);
        leg_neutron_mom->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 9:
        stk_dpt->Add(hh1d_mom[ii]);
        leg_dpt->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 10:
        stk_dphit->Add(hh1d_mom[ii]);
        leg_dphit->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
    
      case 11:
        stk_dalphat->Add(hh1d_mom[ii]);
        leg_dalphat->AddEntry(hh1d_mom[ii],hh1d_mom[ii]->GetTitle(),"fl");
        break;
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


  TString leghead = "GIBUU", leghead_ori = "GIBUU";

  const TString opt_same="same hist";
  const TString opt_err_bar="X1 E0";
  const TString opt_err_bar_same="same X1 E0";
  const TString opt_nostack="hist nostack";

  stk_x->Draw(opt_nostack);
  stk_x->GetXaxis()->SetTitle("x_{Bj}");
  stk_x->GetYaxis()->SetTitle("#frac{d#sigma}{dx_{Bj}} (cm^{2}/x_{Bj}/nucleon)");
  stk_x->SetTitle("x_{Bj} of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_x->Draw();
  c1->Print(pdf_out_name,"x_{Bj} of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "x_{Bj} of diff E modes"));

  stk_Q2->Draw(opt_nostack);
  stk_Q2->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");
  stk_Q2->GetYaxis()->SetTitle("#frac{d#sigma}{dQ^{2}}(cm^{2}/GeV^{2}/nucleon)");
  stk_Q2->SetTitle("Q^{2} of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_Q2->Draw();
  c1->Print(pdf_out_name,"Q^{2} of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "Q^{2} of diff E modes"));

  stk_W->Draw(opt_nostack);
  stk_W->GetXaxis()->SetTitle("W(GeV)");
  stk_W->GetYaxis()->SetTitle("#frac{d#sigma}{dW}(cm^{2}/GeV/nucleon)");
  stk_W->SetTitle("W of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_W->Draw();
  c1->Print(pdf_out_name,"W of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "W of diff E modes"));

  stk_beamE->Draw(opt_nostack);
  stk_beamE->GetXaxis()->SetTitle("beam E(GeV)");
  stk_beamE->GetYaxis()->SetTitle("#frac{d#sigma}{dE_{#nu}}(cm^{2}/GeV/nucleon)");
  stk_beamE->SetTitle("beam E of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_beamE->Draw();
  c1->Print(pdf_out_name,"beamE of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "beam E of diff E modes"));

  stk_mu_theta->Draw(opt_nostack);
  stk_mu_theta->GetXaxis()->SetTitle("#theta of #mu(degree)");
  stk_mu_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_mu_theta->SetTitle("#theta of #mu of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_mu_theta->Draw();
  c1->Print(pdf_out_name,"#theta of #mu of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of #mu of diff E modes"));

  stk_mu_mom->Draw(opt_nostack);
  stk_mu_mom->GetXaxis()->SetTitle("momentum of #mu(GeV)");
  stk_mu_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}}(cm^{2}/GeV/nucleon)");
  stk_mu_mom->SetTitle("momentum of #mu of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_mu_mom->Draw();
  c1->Print(pdf_out_name,"momentum of #mu of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of #mu of diff E modes"));

  stk_pi_theta->Draw(opt_nostack);
  stk_pi_theta->GetXaxis()->SetTitle("#theta of #pi(degree)");
  stk_pi_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_pi_theta->SetTitle("#theta of #pi of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_pi_theta->Draw();
  c1->Print(pdf_out_name,"#theta of #pi of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of #pi of diff E modes"));

  stk_pi_mom->Draw(opt_nostack);
  stk_pi_mom->GetXaxis()->SetTitle("momentum of #pi(GeV)");
  stk_pi_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#pi}}(cm^{2}/GeV/nucleon)");
  stk_pi_mom->SetTitle("momentum of #pi of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_pi_mom->Draw();
  c1->Print(pdf_out_name,"momentum of #pi of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of #pi of diff E modes"));

  stk_recoil_theta->Draw(opt_nostack);
  stk_recoil_theta->GetXaxis()->SetTitle("#theta of recoil(degree)");
  stk_recoil_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_recoil_theta->SetTitle("#theta of recoil of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_recoil_theta->Draw();
  c1->Print(pdf_out_name,"#theta of recoil of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of recoil of diff E modes"));

  stk_recoil_mom->Draw(opt_nostack);
  stk_recoil_mom->GetXaxis()->SetTitle("momentum of recoil(GeV)");
  stk_recoil_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{recoil}}(cm^{2}/GeV/nucleon)");
  stk_recoil_mom->SetTitle("momentum of recoil of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_recoil_mom->Draw();
  c1->Print(pdf_out_name,"momentum of recoil of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of recoil of diff E modes"));

  stk_baryon_theta->Draw(opt_nostack);
  stk_baryon_theta->GetXaxis()->SetTitle("#theta of baryon(degree)");
  stk_baryon_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)");
  stk_baryon_theta->SetTitle("#theta of baryon of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_baryon_theta->Draw();
  c1->Print(pdf_out_name,"#theta of baryon of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "#theta of baryon of diff E modes"));

  stk_baryon_mom->Draw(opt_nostack);
  stk_baryon_mom->GetXaxis()->SetTitle("momentum of baryon(GeV)");
  stk_baryon_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{baryon}}(cm^{2}/GeV/nucleon)");
  stk_baryon_mom->SetTitle("momentum of baryon of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_baryon_mom->Draw();
  c1->Print(pdf_out_name,"momentum of baryon of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "momentum of baryon of diff E modes"));
  
  stk_neutron_mom->Draw(opt_nostack);
  stk_neutron_mom->GetXaxis()->SetTitle("p_{n}(GeV)");
  stk_neutron_mom->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{n}}(cm^{2}/GeV/nucleon)");
  stk_neutron_mom->SetTitle("p_{n} of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_neutron_mom->Draw();
  c1->Print(pdf_out_name,"p_{n} of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "p_{n} of diff E modes"));
  
  stk_dpt->Draw(opt_nostack);
  stk_dpt->GetXaxis()->SetTitle("dpt(GeV)");
  stk_dpt->GetYaxis()->SetTitle("#frac{d#sigma}{d dpt}(cm^{2}/GeV/nucleon)");
  stk_dpt->SetTitle("dpt of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dpt->Draw();
  c1->Print(pdf_out_name,"dpt of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dpt of diff E modes"));

  stk_dphit->Draw(opt_nostack);
  stk_dphit->GetXaxis()->SetTitle("dphit(degree)");
  stk_dphit->GetYaxis()->SetTitle("#frac{d#sigma}{d dphit}(cm^{2}/GeV/nucleon)");
  stk_dphit->SetTitle("dphit of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dphit->Draw();
  c1->Print(pdf_out_name,"dphit of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dphit of diff E modes"));
  
  
  stk_dalphat->Draw(opt_nostack);
  stk_dalphat->GetXaxis()->SetTitle("dalphat(degree)");
  stk_dalphat->GetYaxis()->SetTitle("#frac{d#sigma}{d dalphat}(cm^{2}/GeV/nucleon)");
  stk_dalphat->SetTitle("dalphat of diff E modes");
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_dalphat->Draw();
  c1->Print(pdf_out_name,"dalphat of diff E modes");
  c1->Print(Form("png/%s_%s.png", anatag, "dalphat of diff E modes"));
  
  TString end_pdf_out_name = pdf_out_name + "]";
  c1->Print(end_pdf_out_name,"dalphat of 4 event modess");
  fout->Close();
  return;
}


int main(int argc, char* argv[])
{
  //void anaGenerator(const TString tag, const TString filelist, const int tmpana, const int nToStop=-999, const int smearBit=0)
  if(argc==3){
    hist_2d_prod_no_data(argv[1], argv[2]);
  }
  else{
    printf("wrong argc %d\n", argc); return 1;
  }
  return 0;
}