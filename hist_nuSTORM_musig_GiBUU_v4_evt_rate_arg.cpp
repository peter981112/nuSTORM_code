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
#include "TAxis.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <type_traits>

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

double integrate_const_interval(vector<array<double,2>> flux) //integrate a 1d func with constant interval with format of {x,y}
{

  double integrate_result{0};
  int len = flux.size();
  double dx = flux[1][0]-flux[0][0];
  cout<<"integration dx : "<<dx<<endl;
  double int_coeff[3] {3/8, 7/6, 23/24};
  integrate_result = integrate_result + int_coeff[0]*(flux[0][1]+flux[len-1][1])+ int_coeff[1]*(flux[1][1]+flux[len-2][1])+ int_coeff[2]*(flux[2][1]+flux[len-3][1]);
  for(int i=3; i<len-3; i++)
  {
    integrate_result = integrate_result + flux[i][1];
  }
  integrate_result = integrate_result*dx;
  return integrate_result;

}

vector<array<double,2>> read_flux_file(string flux_file) //read text file and return a vector with elements {double, double}->{nu_E, nu_flux}
{

  ifstream target_flux;
  target_flux.open(flux_file.c_str());
  string fline;
  int totlinecount=0;
  vector<array<double,2>> target_flux_vect;
  array<double, 2> nu_E_flux{ { 0.0, 0.0 } };

  if (target_flux.is_open())
  {
    while(!target_flux.eof())
    {
      getline(target_flux, fline);
      TString tmpline(fline);
      sscanf(tmpline.Data(),"%lf	%lf", &nu_E_flux[0], &nu_E_flux[1]);
      target_flux_vect.push_back(nu_E_flux);
      if(totlinecount%20==0){cout<<totlinecount<<"th nu E : "<<target_flux_vect[totlinecount][0]<<", nu flux : "<<target_flux_vect[totlinecount][1]<<endl;}
      totlinecount++;
    }
    cout << "Number of lines in the file: " << totlinecount << endl;
  }

  return target_flux_vect;
}

void hist_nuSTORM_musig_GiBUU_v4_arg(TString E3_fn, TString E4_fn, TString E5_fn, TString E6_fn, TString E7_fn, string root_file_path, 
    string E3_flux, string E4_flux, string E5_flux, string E6_flux, string E7_flux)
{
  
  TString fn_lst[] = {E3_fn, E4_fn, E5_fn, E6_fn, E7_fn};
  TString flux_lst[] = {E3_flux, E4_flux, E5_flux, E6_flux, E7_flux};

  cout<<"root files for plotting"<<endl;
  for(int i=0;i<5;i++){cout<<i<<" th root file : "<<fn_lst[i]<<endl;}

  cout<<"flux files for event rate"<<endl;
  for(int i=0;i<5;i++){cout<<i<<" th flux file : "<<flux_lst[i]<<endl;}

  vector<array<double,2>> target_flux_vect_lst[] = {read_flux_file(E3_flux), read_flux_file(E4_flux),read_flux_file(E5_flux),read_flux_file(E6_flux),read_flux_file(E7_flux)};

  double amu = 1.67e-27; //mass of nucleon in kg
  double time_interval = 5.0*365*24*3600; //time interval of 5 year in sec

  Double_t hh1d_hist_lower_bound[] = {0,0,0.5,0};
  Double_t hh1d_hist_upper_bound[] = {2.1,2.5,3.2,8};

  Double_t hh1d_mom_hist_lower_bound[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t hh1d_mom_hist_upper_bound[] = {180,5.5,180,3,180,3.5,180,5,1.7,1.6,180,180};

  int anaid{0};auto anatag{"outana"};string fn_string{fn_lst[0]};TString pdf_out_name{""};//int noff{0}

  if(fn_lst[0].Contains("outAna7")){
    anaid = 2;
    //noff = 1;
    anatag = "outana7";
    cout<<"rewighting "<<anatag<<" - case with pi 0"<<endl;

    if(fn_lst[0].Contains("Numu") || fn_lst[0].Contains("numu") ||fn_lst[0].Contains("munu")){pdf_out_name = "plot_numu_outana7.pdf";}
    else{pdf_out_name = "plot_nue_outana7.pdf";}
  }
  else if(fn_lst[0].Contains("outAna9")){
    anaid = 1;
    //noff = 2;
    anatag = "outana9";
    cout<<"rewighting "<<anatag<<" - case with no pi 0"<<endl;

    if(fn_lst[0].Contains("Numu") || fn_lst[0].Contains("numu") ||fn_lst[0].Contains("munu")){pdf_out_name = "plot_numu_outana9.pdf";}
    else{pdf_out_name = "plot_nue_outana9.pdf";}

  }
  
  pdf_out_name = root_file_path + pdf_out_name;
  cout<<" pdf_out_name : "<<pdf_out_name<<endl;
  
  const Int_t n_E = sizeof(fn_lst)/sizeof(TString);
  TCanvas *c1 = new TCanvas("c","",800,600);
  TFile *fout = new TFile("tmpplot.root","recreate");
  TString start_pdf_out_name = pdf_out_name + "[";
  c1->Print(start_pdf_out_name);

  double flux_int_lst[n_E]={0};
  for(int i=0; i<n_E; i++){flux_int_lst[i] = integrate_const_interval(target_flux_vect_lst[i]);}

  double histnormFactor_lst[n_E];
  
  string E_modes[]={"E3","E4","E5","E6","E7"};
  string var[] = {"x", "Q2", "W", "beamE", "scatter #theta", "scatter p", 
      "meson #theta", "meson p", "recoil #theta", "recoil pm", "baryon #theta", "baryon p",
      "p_{n}", "dp_{t}", "d#phi_{t}", "d#alpha_{t}"};
  string mode[] = {"QE", "RES", "DIS", "2p2h"};
  const Int_t nmode = sizeof(mode)/sizeof(string);
  const int hh2d_var_num{3}, hh1d_var_num{4}, hh1d_mom_var_num{12};

  TString temp_name = "";
  TString temp_tit = "";
  TH2D *hh2d[n_E][hh2d_var_num]; //list of 2d hist {hxQ2_E7, hxW_E7}
  TString hh2d_name_lst[n_E][hh2d_var_num];
  TString hh2d_tit_lst[n_E][hh2d_var_num];
  for(int i=0; i<hh2d_var_num*n_E; i++)
  {
    if(i%hh2d_var_num == 0)
    {
      temp_name = "hxQ2 "+ E_modes[i/hh2d_var_num];
      temp_tit = "x-Q^{2} for "+ E_modes[i/hh2d_var_num];
    }
    else if(i%hh2d_var_num == 1)
    {
      temp_name = "hxW "+ E_modes[i/hh2d_var_num];
      temp_tit = "x-W for "+ E_modes[i/hh2d_var_num];
    }
    else
    {
      temp_name = "hbeamEW "+ E_modes[i/hh2d_var_num];
      temp_tit = "beamE-W for "+ E_modes[i/hh2d_var_num];
    }
    hh2d_name_lst[i/hh2d_var_num][i%hh2d_var_num] = temp_name;
    hh2d_tit_lst[i/hh2d_var_num][i%hh2d_var_num] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }

  TH1D *hh1d[n_E][hh1d_var_num]; //list of 1d hist {hxQ2_E7_x, hxQ2_E7_Q2, hxQ2_E7_W, hxQ2_E7_beamE}
  TString hh1d_name_lst[n_E][hh1d_var_num];
  TString hh1d_tit_lst[n_E][hh1d_var_num];
  for(int i=0; i<hh1d_var_num*n_E; i++)
  {
    temp_name = "hxQ2 "+var[i%hh1d_var_num]+ E_modes[i/hh1d_var_num];
    temp_tit = var[i%hh1d_var_num]+" for "+ E_modes[i/hh1d_var_num];

    hh1d_name_lst[i/hh1d_var_num][i%hh1d_var_num] = temp_name;
    hh1d_tit_lst[i/hh1d_var_num][i%hh1d_var_num] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }
    
  TH1D *hh1d_mom[n_E][hh1d_mom_var_num];
  TString hh1d_mom_name_lst[n_E][hh1d_mom_var_num];
  TString hh1d_mom_tit_lst[n_E][hh1d_mom_var_num];
  for(int i=0; i<hh1d_mom_var_num*n_E; i++)
  {
    temp_name = "hxQ2 "+var[i%hh1d_mom_var_num+hh1d_var_num]+ E_modes[i/hh1d_mom_var_num];
    temp_tit = var[i%hh1d_mom_var_num+hh1d_var_num]+" for "+ E_modes[i/hh1d_mom_var_num];

    hh1d_mom_name_lst[i/hh1d_mom_var_num][i%hh1d_mom_var_num] = temp_name;
    hh1d_mom_tit_lst[i/hh1d_mom_var_num][i%hh1d_mom_var_num] = temp_tit;
    //cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }

  double xsec_lst[n_E][4] = {0}; //list of xsec of each E mode event mode QE, RES, DIS, 2p2h
  double event_num_lst[n_E][4] = {0}; //list of event num of each E mode event mode QE, RES, DIS, 2p2h

  const int nh1d = sizeof(hh1d)/sizeof(TH1D*);
  const int nh1d_mom = sizeof(hh1d_mom)/sizeof(TH1D*);
  cout<<"mom array len : "<<nh1d_mom<<endl;
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);

  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange, kViolet+3};

  for(Int_t kk=0; kk<n_E; kk++)
  {

    hh2d[kk][0] = new TH2D(hh2d_name_lst[kk][0],hh2d_tit_lst[kk][0],60,0,2,80,0.1,10);
    hh2d[kk][1] = new TH2D(hh2d_name_lst[kk][1],hh2d_tit_lst[kk][1],60,0,2,80,0.6,3);
    hh2d[kk][2] = new TH2D(hh2d_name_lst[kk][2],hh2d_tit_lst[kk][2],60,0,8,80,0.6,3);
    
    for(Int_t ii=0; ii<hh1d_var_num; ii++)
    {
      hh1d[kk][ii] = new TH1D(hh1d_name_lst[kk][ii],hh1d_tit_lst[kk][ii],60,hh1d_hist_lower_bound[ii],hh1d_hist_upper_bound[ii]);
    }
    
    for(Int_t ii=0; ii<hh1d_mom_var_num; ii++)
    {
      hh1d_mom[kk][ii] = new TH1D(hh1d_mom_name_lst[kk][ii],hh1d_mom_tit_lst[kk][ii],60,hh1d_mom_hist_lower_bound[ii],hh1d_mom_hist_upper_bound[ii]);
    }

    TLegend *leg = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg->SetFillStyle(0);
    //leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.023);
    gStyle-> SetOptStat(0); // turn off stat box
    gStyle->SetPadRightMargin(0.2);
    gStyle->SetPadLeftMargin(0.1);
    gStyle-> SetPalette(kRainbow);


    //fill_hist(fn_lst[kk], anaid, E_modes[kk], histnormFactor_lst[kk], hh2d[kk], hh1d[kk], hh1d_mom[kk]);
    TString fn{fn_lst[kk]};
    TFile *f = new TFile(fn);
    TTree *t1 = (TTree*) f->Get("tree");
    ULong64_t minnp, maxnp;
    const bool isTopoTask = GeneratorUtils::SetTopoCut(anaid, minnp, maxnp);

    const TString npiDecomp[]={"","_0pi","_1pi","_Mpi","_other"};
    //only gointo sub-category for TopoTask
    const int nNpiCut = isTopoTask? (sizeof(npiDecomp)/sizeof(TString)) : 1;
    cout<<"n pi cut : "<<nNpiCut<<endl;
    TString nuExp;
    int nuPDG, tarA;
    const TString gen = GeneratorUtils::SetEnv(fn_lst[kk], nuExp, nuPDG, tarA);
    //const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*13;cout<<"scaling 1/13 for CH"<<endl//Scaling to per nucleon in CH 13=1+12
    const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*12;cout<<"scaling 1/12 for Carbon target"<<endl;//Scaling to per nucleon in Carbon 12
    //const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*18;cout<<"scaling 1/18 for H2O"<<endl;//Scaling to per nucleon in H2O 18=2+16
    //const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*16;cout<<"scaling 1/16 for Oxygen target"<<endl;//Scaling to per nucleon in Oxygen 16
    histnormFactor_lst[kk] = histnormFactor;

    
    Double_t xBj, Q2, Wtrue, beamE, scattertheta, scattermomentum, mesontheta, mesonmomentum, 
      recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat, dpTT;
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

    cout<<"init hist in process"<<endl;

    //TH2D *h_weight__scattertheta = new TH2D("h_weight_E_scattertheta","distribution of weight and scattertheta mode : E",80,0,26,100,0.0000001,0.001);
  
  
    Int_t nentries = (Int_t)t1->GetEntries();
    cout<<"events for energy "<<E_modes[kk]<<" : "<<nentries<<endl;
    int npi_err{0};
    for (Int_t i=0; i<nentries; i++)
    {
      t1->GetEntry(i);
      //if(perweight<0 && i%500==0){cout<<"weight : "<<perweight<<endl;}
      if(i%1000000==0){cout<<i<<" th meson(pi) p : "<<mesonmomentum<<endl;}
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

      hh2d[kk][0]->Fill(xBj,Q2,perweight);
      hh2d[kk][1]->Fill(xBj,Wtrue,perweight);
      hh2d[kk][2]->Fill(beamE,Wtrue,perweight);
      Double_t hh1d_var_lst[] = {xBj, Q2, Wtrue, beamE};
      Double_t hh1d_mom_var_lst[] = {scattertheta, scattermomentum, mesontheta, mesonmomentum, recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat};
      for(Int_t ii=0; ii<hh1d_var_num; ii++)
      {
        hh1d[kk][ii]->Fill(hh1d_var_lst[ii],perweight);
      }

      for(Int_t ii=0; ii<hh1d_mom_var_num; ii++)
      {
        hh1d_mom[kk][ii]->Fill(hh1d_mom_var_lst[ii],perweight);
      }


      if(prod==1 && evtMode==1){xsec_lst[kk][0] = xsec_lst[kk][0]+perweight;} //add weight to get xsec for QE

      else if(prod>=2 && prod<=33 && evtMode==2){xsec_lst[kk][1] = xsec_lst[kk][1]+perweight;} //add weight to get xsec for RES

      else if(prod==34 && evtMode==3){xsec_lst[kk][2] = xsec_lst[kk][2]+perweight;} //add weight to get xsec for DIS

      else if(prod==35 && evtMode==4){xsec_lst[kk][3] = xsec_lst[kk][3]+perweight;} //add weight to get xsec for 2p2h
      
    }
    cout<<" n pi err # : "<<npi_err<<endl;

    delete t1;
    delete f;
  };

  //xsec calc

  cout<<"xsec for "<<anatag<<endl;
  double total_xsec{0};
  double total_xsec_lst[n_E] = {0};
  for(int kk=0; kk<n_E; ++kk)
  {
    for(int ii=0; ii<nmode; ++ii)
    {
      xsec_lst[kk][ii] = xsec_lst[kk][ii]/histnormFactor_lst[kk];
      total_xsec = total_xsec + xsec_lst[kk][ii];
      cout<<"xsec of "<<E_modes[kk]<<" mode : "<<mode[ii]<<" is "<<xsec_lst[kk][ii]<<" cm^{2}"<<endl;
    }
    cout<<"total xsec of "<<E_modes[kk]<<" is "<<total_xsec<<" cm^{2}"<<endl;
    total_xsec_lst[kk] = total_xsec;
    total_xsec = 0;
  }

  //event rate calc 

  double total_event_num_lst[n_E] = {0};
  for(int kk=0; kk<n_E; ++kk)
  {
    for(int ii=0; ii<nmode; ++ii)
    {
      event_num_lst[kk][ii] = xsec_lst[kk][ii]*time_interval*flux_int_lst[kk]/amu;
      cout<<"event num of "<<E_modes[kk]<<" in "<<time_interval/(365*24*3600)<<" years for mode : "
        <<mode[ii]<<" is "<<event_num_lst[kk][ii]<<" event /POT * kg"<<endl;
    }
    total_event_num_lst[kk] = total_xsec_lst[kk]*time_interval*flux_int_lst[kk]/amu;
    cout<<"total event num of "<<E_modes[kk]<<" in "<<time_interval/(365*24*3600)<<" years : "<<total_event_num_lst[kk]<<" event /POT * kg"<<endl;
  }

  cout<<"now printing cross section formatted for plotting in path : "<<root_file_path<<endl;

  cout<<"total xsec [";
  for(int ii=0; ii<n_E; ++ii)
  {
    cout<<total_xsec_lst[ii]<<", ";
  }
  cout<<"] "<<endl;

  for(int ii=0; ii<nmode; ++ii)
  {
    cout<<"xsec of "<<mode[ii]<<" [";
    for(int kk=0; kk<n_E; ++kk)
    {
      cout<<xsec_lst[kk][ii]<<", ";
    }
    cout<<"] "<<endl;
  }

  cout<<"event num units in : event /POT * kg"<<endl;
  cout<<"total event num [";
  for(int ii=0; ii<n_E; ++ii)
  {
    cout<<total_event_num_lst[ii]<<", ";
  }
  cout<<"] "<<endl;

  for(int ii=0; ii<nmode; ++ii)
  {
    cout<<"event num of "<<mode[ii]<<" [";
    for(int kk=0; kk<n_E; ++kk)
    {
      cout<<event_num_lst[kk][ii]<<", ";
    }
    cout<<"] "<<endl;
  }


  THStack stk_x;// = new THStack;
  TLegend leg_x = TLegend(0.7, 0.5, 0.9,0.9);
  leg_x.SetFillStyle(0);
    
  THStack stk_Q2;// = new THStack;
  TLegend leg_Q2 = TLegend(0.7, 0.5, 0.9,0.9);
  leg_Q2.SetFillStyle(0);
  
  THStack stk_W;// = new THStack;
  TLegend leg_W = TLegend(0.7, 0.5, 0.9,0.9);
  leg_W.SetFillStyle(0);

  THStack stk_beamE;// = new THStack;
  TLegend leg_beamE = TLegend(0.7, 0.5, 0.9,0.9);
  leg_beamE.SetFillStyle(0);
    
  THStack stk_mu_theta;// = new THStack;
  TLegend leg_mu_theta = TLegend(0.7, 0.5, 0.9,0.9);
  leg_mu_theta.SetFillStyle(0);
    
  THStack stk_mu_mom;// = new THStack;
  TLegend leg_mu_mom = TLegend(0.7, 0.5, 0.9,0.9);
  leg_mu_mom.SetFillStyle(0);

  THStack stk_pi_theta;// = new THStack;
  TLegend leg_pi_theta = TLegend(0.7, 0.5, 0.9,0.9);
  leg_pi_theta.SetFillStyle(0);
    
  THStack stk_pi_mom;// = new THStack;
  TLegend leg_pi_mom = TLegend(0.7, 0.5, 0.9,0.9);
  leg_pi_mom.SetFillStyle(0);

  THStack stk_recoil_theta;// = new THStack;
  TLegend leg_recoil_theta = TLegend(0.7, 0.5, 0.9,0.9);
  leg_recoil_theta.SetFillStyle(0);
    
  THStack stk_recoil_mom;// = new THStack;
  TLegend leg_recoil_mom = TLegend(0.7, 0.5, 0.9,0.9);
  leg_recoil_mom.SetFillStyle(0);

  THStack stk_baryon_theta;// = new THStack;
  TLegend leg_baryon_theta = TLegend(0.7, 0.5, 0.9,0.9);
  leg_baryon_theta.SetFillStyle(0);
    
  THStack stk_baryon_mom;// = new THStack;
  TLegend leg_baryon_mom = TLegend(0.7, 0.5, 0.9,0.9);
  leg_baryon_mom.SetFillStyle(0);

  THStack stk_neutron_mom;// = new THStack;
  TLegend leg_neutron_mom = TLegend(0.7, 0.5, 0.9,0.9);
  leg_neutron_mom.SetFillStyle(0);
  
  THStack stk_dpt;// = new THStack;
  TLegend leg_dpt = TLegend(0.7, 0.5, 0.9,0.9);
  leg_dpt.SetFillStyle(0);

  THStack stk_dphit;// = new THStack;
  TLegend leg_dphit = TLegend(0.7, 0.5, 0.9,0.9);
  leg_dphit.SetFillStyle(0);
  
  THStack stk_dalphat;// = new THStack;
  TLegend leg_dalphat = TLegend(0.13, 0.5, 0.33,0.9);
  leg_dalphat.SetFillStyle(0);


  int normal_choice{1}; //if 0 just remove negative cross section if 1 do column normalization
  double histnormFactor{0};
  cout<<"starting 2d hist"<<endl;
  for(int ii=0; ii<nh2d; ii++){
    histnormFactor = histnormFactor_lst[ii/2];
    if(normal_choice==0)
    {
      hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->Scale(1/histnormFactor);
      remove_neg(hh2d[ii/hh2d_var_num][ii%hh2d_var_num]);
    }
    else if(normal_choice==1)
    {
      col_normal(hh2d[ii/hh2d_var_num][ii%hh2d_var_num]);
    }
    hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetXaxis()->SetTitle("x_{Bj}");
    if(ii%hh2d_var_num==0)
    {
      hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");
      gPad->SetLogy(1); //set y axis to log-1 lin-0 scale
      gPad->Update();
    }
    else if(ii%hh2d_var_num==1)
    {
      hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetYaxis()->SetTitle("W(GeV)");
      gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
      gPad->Update();
    }
    else if(ii%hh2d_var_num==2)
    {
      hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetXaxis()->SetTitle("beam E(GeV)");
      hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetYaxis()->SetTitle("W(GeV)");
      gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
      gPad->Update();
    }
    hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->SetContour(1000);
    //hh2d[ii]->Draw("TEXT");
    hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->Draw("colz");

    gStyle->SetPadRightMargin(0.2);
    //redraw legend and update bottom margin dependent on number of entries

    // add this plot to plotbook
    c1->Print(pdf_out_name,hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetTitle());
    c1->Print(Form("png/%s_%s.png", anatag, hh2d[ii/hh2d_var_num][ii%hh2d_var_num]->GetTitle()));

    //leg->Clear();
  }
  cout<<"2d hist done"<<endl;
  /*
  remove_neg(h_weight_E7_scattertheta);
  h_weight_E7_scattertheta->GetXaxis()->SetTitle("E7 #theta_{mu}");
  h_weight_E7_scattertheta->GetYaxis()->SetTitle("weight");
  gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
  gPad->Update();
  h_weight_E7_scattertheta->SetContour(1000);
  h_weight_E7_scattertheta->Draw("colz");
  gStyle->SetPadRightMargin(0.2);
  c1->Print(pdf_out_name,h_weight_E7_scattertheta->GetTitle());
  c1->Print(Form("png/%s_%s.png", anatag, h_weight_E7_scattertheta->GetTitle()));
  */

  int aa{0};
  Double_t binwidth{0};
  cout<<"starting 1d hist"<<endl;

  for(int ii=0; ii<nh1d; ii++){
    if(aa%4==0){++aa;};
    histnormFactor = histnormFactor_lst[ii/hh1d_var_num];
    binwidth = hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->GetXaxis()->GetBinWidth(1);
    hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->Scale(1/(histnormFactor*binwidth)); //to get diff cross section
    //hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->GetXaxis()->SetTitle(tit[aa]);
    hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->GetYaxis()->SetTitle("differential cross section(cm^{2}/GeV/nucleon)");
    
    //leg->SetFillStyle(0);
    hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->SetLineStyle(kSolid);
    hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->SetLineColor(cols[ii/4]);
    hh1d[ii/hh1d_var_num][ii%hh1d_var_num]->SetFillStyle(3001);
    //hh1d[ii]->SetFillStyle(4050);//for trasparent filling
    //hh1d[ii]->SetFillColorAlpha(cols[ii/3],0.35);

    //xsec_lst[(ii/3)*14+ii%3] = temp_xsec;

    //TString temp_tit="";
    //temp_tit = tit[ii*4];
    if(ii%hh1d_var_num==0){
      stk_x.Add(hh1d[ii/4][ii%4]);
      leg_x.AddEntry(hh1d[ii/4][ii%4],hh1d[ii/4][ii%4]->GetTitle(),"fl");
    }
    
    else if(ii%hh1d_var_num==1){
      stk_Q2.Add(hh1d[ii/4][ii%4]);
      leg_Q2.AddEntry(hh1d[ii/4][ii%4],hh1d[ii/4][ii%4]->GetTitle(),"fl");
    }

    else if(ii%hh1d_var_num==2){
      stk_W.Add(hh1d[ii/4][ii%4]);
      leg_W.AddEntry(hh1d[ii/4][ii%4],hh1d[ii/4][ii%4]->GetTitle(),"fl");
    }
    else{
      stk_beamE.Add(hh1d[ii/4][ii%4]);
      leg_beamE.AddEntry(hh1d[ii/4][ii%4],hh1d[ii/4][ii%4]->GetTitle(),"fl");
    }
    ++aa;
  }
  cout<<"strting 1d kinematics hist"<<endl;
  //binwidth = 0;
    
  int count_num{hh1d_mom_var_num};
  for(int ii=0; ii<nh1d_mom; ii++)
  {
    histnormFactor = histnormFactor_lst[ii/hh1d_mom_var_num];
    diff_xsec(hh1d_mom[ii/hh1d_mom_var_num][ii%hh1d_mom_var_num], histnormFactor);
    hh1d_mom[ii/hh1d_mom_var_num][ii%hh1d_mom_var_num]->GetYaxis()->SetTitle("differential cross section(cm^{2}/GeV/nucleon)");
    
    //leg->SetFillStyle(0);
    hh1d_mom[ii/hh1d_mom_var_num][ii%hh1d_mom_var_num]->SetLineStyle(kSolid);
    hh1d_mom[ii/hh1d_mom_var_num][ii%hh1d_mom_var_num]->SetLineColor(cols[ii/count_num]);
    hh1d_mom[ii/hh1d_mom_var_num][ii%hh1d_mom_var_num]->SetFillStyle(3001);
    //hh1d_mom[ii]->SetFillStyle(4050);
    //hh1d_mom[ii]->SetFillColorAlpha(cols[ii/count_num],0.35);
    //cout<<"title : "<<hh1d_mom[ii]->GetTitle()<<endl;

    switch (ii%count_num)
    {

      case 0:
        stk_mu_theta.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_mu_theta.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;

      case 1:
        stk_mu_mom.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_mu_mom.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;

      case 2:
        stk_pi_theta.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_pi_theta.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 3:
        stk_pi_mom.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_pi_mom.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 4:
        stk_recoil_theta.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_recoil_theta.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 5:
        stk_recoil_mom.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_recoil_mom.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 6:
        stk_baryon_theta.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_baryon_theta.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 7:
        stk_baryon_mom.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_baryon_mom.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 8:
        stk_neutron_mom.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_neutron_mom.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 9:
        stk_dpt.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_dpt.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 10:
        stk_dphit.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_dphit.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;
    
      case 11:
        stk_dalphat.Add(hh1d_mom[ii/count_num][ii%count_num]);
        leg_dalphat.AddEntry(hh1d_mom[ii/count_num][ii%count_num],hh1d_mom[ii/count_num][ii%count_num]->GetTitle(),"fl");
        break;

    }
    
  }
  
  gPad->SetLogy(0);//set y axis to linear scale
  gPad->Update();


  const TString opt_same="same hist";
  const TString opt_err_bar="X1 E0";
  const TString opt_err_bar_same="same X1 E0";
  const TString opt_nostack="hist nostack";

  THStack stk_lst[] = {stk_x, stk_Q2, stk_W, stk_beamE, stk_mu_theta, stk_mu_mom, stk_pi_theta, stk_pi_mom, stk_recoil_theta, stk_recoil_mom
    , stk_baryon_theta, stk_baryon_mom, stk_neutron_mom, stk_dpt, stk_dphit, stk_dalphat};

  TLegend leg_lst[] = {leg_x, leg_Q2, leg_W, leg_beamE, leg_mu_theta, leg_mu_mom, leg_pi_theta, leg_pi_mom, leg_recoil_theta, leg_recoil_mom
    , leg_baryon_theta, leg_baryon_mom, leg_neutron_mom, leg_dpt, leg_dphit, leg_dalphat};

  const int n_stk = sizeof(stk_lst)/sizeof(THStack);
  cout<<"n_stk : "<<n_stk<<endl;

  string tit_lst[n_stk][5] = {{"x_{Bj}", "#frac{d#sigma}{dx_{Bj}} (cm^{2}/x_{Bj}/nucleon)", "x_{Bj} of diff E modes", "x_{Bj} of diff E modes", "x_{Bj} of diff E modes"},
    {"Q^{2}(GeV^{2})", "#frac{d#sigma}{dQ^{2}}(cm^{2}/GeV^{2}/nucleon)", "Q^{2} of diff E modes", "Q^{2} of diff E modes", "Q^{2} of diff E modes"},
    {"W(GeV)", "#frac{d#sigma}{dW}(cm^{2}/GeV/nucleon)", "W of diff E modes", "W of diff E modes", "W of diff E modes"},
    {"beamE(GeV)", "#frac{d#sigma}{dE_{nu}}(cm^{2}/GeV/nucleon)", "beamE of diff E modes", "beamE of diff E modes", "beamE of diff E modes"},

    {"#theta of #mu(degree)", "#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)", "#theta of #mu of diff E modes", "#theta of #mu of diff E modes", "#theta of #mu of diff E modes"},
    {"momentum of #mu(GeV)", "#frac{d#sigma}{dp_{#mu}}(cm^{2}/GeV/nucleon)", "momentum of #mu of diff E modes", "momentum of #mu of diff E modes", "momentum of #mu of diff E modes"},
    {"#theta of #pi(degree)", "#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)", "#theta of #pi of diff E modes", "#theta of #pi of diff E modes", "#theta of #pi of diff E modes"},
    {"momentum of #pi(GeV)", "#frac{d#sigma}{dp_{#pi}}(cm^{2}/GeV/nucleon)", "momentum of #pi of diff E modes", "momentum of #pi of diff E modes", "momentum of #pi of diff E modes"},
    {"#theta of recoil(degree)", "#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)", "#theta of recoil of diff E modes", "#theta of recoil of diff E modes", "#theta of recoil of diff E modes"},
    {"momentum of recoil(GeV)", "#frac{d#sigma}{dp_{recoil}}(cm^{2}/GeV/nucleon)", "momentum of recoil of diff E modes", "momentum of recoil of diff E modes", "momentum of recoil of diff E modes"},
    {"#theta of baryon(degree)", "#frac{d#sigma}{d#theta}(cm^{2}/degree/nucleon)", "#theta of baryon of diff E modes", "#theta of baryon of diff E modes", "#theta of baryon of diff E modes"},
    {"momentum of baryon(GeV)", "#frac{d#sigma}{dp_{baryon}}(cm^{2}/GeV/nucleon)", "momentum of baryon of diff E modes", "momentum of baryon of diff E modes", "momentum of baryon of diff E modes"},
    {"p_{n}(GeV)", "#frac{d#sigma}{dp_{n}}(cm^{2}/GeV/nucleon)", "p_{n} of diff E modes", "p_{n} of diff E modes", "p_{n} of diff E modes"},
    {"dp_{t}(GeV)", "#frac{d#sigma}{dp_{t}}(cm^{2}/GeV/nucleon)", "dpt of diff E modes", "dpt of diff E modes", "dpt of diff E modes"},
    {"d#phi_{t}(degree)", "#frac{d#sigma}{d#phi_{t}}(cm^{2}/GeV/nucleon)", "dphit of diff E modes", "dphit of diff E modes", "dphit of diff E modes"},
    {"d#alpha_{t}(degree)", "#frac{d#sigma}{d#alpha_{t}}(cm^{2}/GeV/nucleon)", "dalphat of diff E modes", "dalphat of diff E modes", "dalphat of diff E modes"}};

  

  for(int ii=0; ii<n_stk; ii++)
  {
    //cout<<"test for loop ii : "<<ii<<endl;
    stk_lst[ii].Draw(opt_nostack);
    stk_lst[ii].GetXaxis()->SetTitle(tit_lst[ii][0].c_str());
    stk_lst[ii].GetYaxis()->SetTitle(tit_lst[ii][1].c_str());
    stk_lst[ii].SetTitle(tit_lst[ii][2].c_str());
    gStyle->SetTitleX(0.5);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.1);
    leg_lst[ii].Draw();
    c1->Print(pdf_out_name,tit_lst[ii][3].c_str());
    c1->Print(Form("png/%s_%s.png", anatag, tit_lst[ii][4].c_str()));

  }

  TString end_pdf_out_name = pdf_out_name + "]";
  c1->Print(end_pdf_out_name,"dalphat of 4 event modess");
  fout->Close();
  cout<<"plotting done. closing file"<<endl;

  delete c1;
  for(int ii=0; ii<n_E; ii++)
  {
    for(int jj=0; jj<hh2d_var_num; jj++){delete hh2d[ii][jj];}
    for(int jj=0; jj<hh1d_var_num; jj++){delete hh1d[ii][jj];}
    for(int jj=0; jj<hh1d_mom_var_num; jj++){delete hh1d_mom[ii][jj];}
  }
  
  return;
}


int main(int argc, char* argv[])
{
  //void anaGenerator(const TString tag, const TString filelist, const int tmpana, const int nToStop=-999, const int smearBit=0)
  if(argc==12){
    hist_nuSTORM_musig_GiBUU_v4_arg(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11]);
  }
  else{
    printf("wrong argc %d\n", argc); return 1;
  }
  return 0;
}