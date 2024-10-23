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

#include "/opt/ppd/scratch/nustorm/GiBUU_2024/NuGenTKI/include/HistIO.h"
#include "/opt/ppd/scratch/nustorm/GiBUU_2024/NuGenTKI/include/GeneratorUtils.h"
#include "/opt/ppd/scratch/nustorm/GiBUU_2024/NuGenTKI/style/style.h"

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

void BinLogY(TH2D* &h)
{

   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from =  TMath::Log10(axis->GetXmin());
   Axis_t to = TMath::Log10(axis->GetXmax());
   cout<<"from : "<<from<<" to : "<<to<<endl;
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
}

void hist_simp_plot_no_data_recoil_outana7()
{
/*
  TString fn{"outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"plotting outana9 - case with no pi 0"<<endl;
   */
  TString fn{"outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root"};
  int anaid{2};
  int noff{1};
  cout<<"plotting outana7 - case with pi 0"<<endl;
  auto anatag{"outana7"};

  
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

  //xBj, Q2, Wtrue, scattertheta, scattermomentum, mesontheta, mesonmomentum, recoiltheta, recoilmomentum, baryontheta, baryonmomentum, IApN, dpt, dphit, dalphat, dpTT

  Double_t Q2, mesonmomentum, perweight, recoiltheta, recoilmomentum;//
  Int_t prod,evtMode;
  t1->SetBranchAddress("recoiltheta",&recoiltheta);
  t1->SetBranchAddress("recoilmomentum",&recoilmomentum);
  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("mesonmomentum",&mesonmomentum);

  t1->SetBranchAddress("prod",&prod);
  t1->SetBranchAddress("evtMode",&evtMode);
  t1->SetBranchAddress("perweight",&perweight);
  
  TString int_modes[]={"QE","RES","DIS","2p2h"};
  TString var = "recoiltheta";
  TString var_2d = "recoilmomentum";

  TString hh_name_lst[8];
  TString hh_tit_lst[8];
  TString temp_name{""}, temp_tit{""};
  for(int i=0; i<8; i++)
  {
    if(i%2 == 0)
    {
      temp_name = temp_name+"hxQ2 "+ int_modes[i/2];
      temp_tit = temp_tit+var+"-"+var_2d+" for "+ int_modes[i/2];
    }
    else
    {
      temp_name = temp_name+"hxW "+ int_modes[i/2];
      temp_tit = temp_tit+var+" for "+ int_modes[i/2];
    }
    hh_name_lst[i] = temp_name;
    hh_tit_lst[i] = temp_tit;
    cout<<"name : "<<temp_name<<" tit : "<<temp_tit<<endl;
    temp_name = "";
    temp_tit = "";
  }
/*
  Double_t hh_lower_bound[] = {0,0,0}; //lower bound of 1d, 2dx, 2dy
  Double_t hh_upper_bound[] = {180,180,6}; //upper bound of 1d, 2dx, 2dy
  */
  Double_t hh_lower_bound[] = {0,0,0}; //lower bound of 1d, 2dx, 2dy
  Double_t hh_upper_bound[] = {180,180,6}; //upper bound of 1d, 2dx, 2dy



  TH2D *hxQ2_QE_2d = new TH2D(hh_name_lst[0],hh_tit_lst[0],60,hh_lower_bound[1],hh_upper_bound[1],80,hh_lower_bound[2],hh_upper_bound[2]);
  TH1D *hxQ2_QE = new TH1D(hh_name_lst[1],hh_tit_lst[1],80,hh_lower_bound[0],hh_upper_bound[0]);

  TH2D *hxQ2_RES_2d = new TH2D(hh_name_lst[2],hh_tit_lst[2],60,hh_lower_bound[1],hh_upper_bound[1],80,hh_lower_bound[2],hh_upper_bound[2]);
  TH1D *hxQ2_RES = new TH1D(hh_name_lst[3],hh_tit_lst[3],80,hh_lower_bound[0],hh_upper_bound[0]);

  TH2D *hxQ2_DIS_2d = new TH2D(hh_name_lst[4],hh_tit_lst[4],60,hh_lower_bound[1],hh_upper_bound[1],80,hh_lower_bound[2],hh_upper_bound[2]);
  TH1D *hxQ2_DIS = new TH1D(hh_name_lst[5],hh_tit_lst[5],80,hh_lower_bound[0],hh_upper_bound[0]);

  TH2D *hxQ2_2p2h_2d = new TH2D(hh_name_lst[6],hh_tit_lst[6],60,hh_lower_bound[1],hh_upper_bound[1],80,hh_lower_bound[2],hh_upper_bound[2]);
  TH1D *hxQ2_2p2h = new TH1D(hh_name_lst[7],hh_tit_lst[7],80,hh_lower_bound[0],hh_upper_bound[0]);

  //BinLogY(hxQ2_QE_2d);
  //BinLogY(hxQ2_RES_2d);
  //BinLogY(hxQ2_DIS_2d);
  //BinLogY(hxQ2_2p2h_2d);

  TString h_all_2d_name{""};
  h_all_2d_name = "distribution of all mode included for "+var+"-"+var_2d;
  TH2D *h_all_2d = new TH2D("h_all_2d",h_all_2d_name,80,hh_lower_bound[1],hh_upper_bound[1],100,hh_lower_bound[2],hh_upper_bound[2]);
  //BinLogY(h_weight_var);

  TH2D *hh2d[] = {hxQ2_QE_2d, hxQ2_RES_2d, hxQ2_DIS_2d, hxQ2_2p2h_2d, h_all_2d};
  TH1D *hh1d[] = {hxQ2_QE, hxQ2_RES, hxQ2_DIS, hxQ2_2p2h};


  const Int_t cols[]={kRed-3, kBlue, kGreen+3, kOrange};
  string mode[] = {"QE", "RES", "DIS", "2p2h"};

  

  THStack *stk_1d = new THStack;
  TLegend *leg_1d = new TLegend(0.7, 0.5, 0.9,0.9);
  leg_1d->SetFillStyle(0);
 
 

  Int_t nentries = (Int_t)t1->GetEntries();
  const int nh1d = sizeof(hh1d)/sizeof(TH1D*);
  cout<<"number of events : "<<nentries<<endl;
  const int nh2d = sizeof(hh2d)/sizeof(TH2D*);
  int npi_err{0};
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    Double_t hh1d_var_lst[] = {xBj, Q2};
        
    //if(perweight<0 && i%500==0){cout<<"weight : "<<perweight<<endl;}
    if(i%1000000==0){cout<<(i/1000000)<<"000000th meson(pi) p : "<<mesonmomentum<<endl;}
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
      /*
        "recoiltheta";
        "recoilmomentum";
        scattertheta
        scattermomentum
      */
      if(prod==1 && evtMode==1) //fil histogram bins for quasi elastic case
      {
        kk = 0;
        hh2d[kk]->Fill(recoiltheta, recoilmomentum,perweight);
        hh1d[kk]->Fill(recoiltheta,perweight);        
        h_all_2d->Fill(recoiltheta, recoilmomentum,perweight);
      }
      else if(prod>=2 && prod<=33 && evtMode==2) //fil histogram bins for resonance case
      {
        kk = 1;
        hh2d[kk]->Fill(recoiltheta, recoilmomentum,perweight);
        hh1d[kk]->Fill(recoiltheta,perweight);        
        h_all_2d->Fill(recoiltheta, recoilmomentum,perweight);

      }
      else if(prod==34 && evtMode==3) //fil histogram bins for deep inelastic case
      {
        kk = 2;
        hh2d[kk]->Fill(recoiltheta, recoilmomentum,perweight);
        hh1d[kk]->Fill(recoiltheta,perweight);        
        h_all_2d->Fill(recoiltheta, recoilmomentum,perweight);
      }
      else if(prod==35 && evtMode==4) //fil histogram bins for 2 part 2 hole case
      {
        kk = 3;
        hh2d[kk]->Fill(recoiltheta, recoilmomentum,perweight);
        hh1d[kk]->Fill(recoiltheta,perweight);        
        h_all_2d->Fill(recoiltheta, recoilmomentum,perweight);
      }
    }
  }
  cout<<" n pi err # : "<<npi_err<<endl;

  Double_t temp_xsec{0};


  //const int nxsec = sizeof(xsec_lst)/sizeof(Double_t);

  int normal_choice{0}; //if 0 just remove negative cross section if 1 do column normalization

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
    hh2d[ii]->GetXaxis()->SetTitle(var);
    hh2d[ii]->GetYaxis()->SetTitle(var_2d);
    gPad->SetLogy(0); //set y axis to log-1 lin-0 scale
    gPad->Update();

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

  
  int aa{0};
  Double_t binwidth{0};
  for(int ii=0; ii<nh1d; ii++){
    if(aa%4==0){++aa;};
    binwidth = hh1d[ii]->GetXaxis()->GetBinWidth(1);
    hh1d[ii]->Scale(1/(histnormFactor*binwidth)); //to get diff cross section
    hh1d[ii]->GetXaxis()->SetTitle(var);
    hh1d[ii]->GetYaxis()->SetTitle("differential cross section(cm^{2}/GeV/nucleon)");
    
    //leg->SetFillStyle(0);
    hh1d[ii]->SetLineStyle(kSolid);
    hh1d[ii]->SetLineColor(cols[ii]);
    hh1d[ii]->SetFillStyle(4050);
    hh1d[ii]->SetFillColorAlpha(cols[ii],0.35);

    temp_xsec = hh1d[ii]->Integral(0,10000,"width");
    //xsec_lst[(ii/3)*14+ii%3] = temp_xsec;


    stk_1d->Add(hh1d[ii]);
    leg_1d->AddEntry(hh1d[ii],hh1d[ii]->GetTitle(),"fl");

    ++aa;
  }
  //binwidth = 0;

  
  gPad->SetLogy(0);//set y axis to linear scale
  gPad->Update();

  aa=0;
  /*
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
*/

  TString leghead = "GIBUU", leghead_ori = "GIBUU";

  const TString opt_same="same hist";
  const TString opt_err_bar="X1 E0";
  const TString opt_err_bar_same="same X1 E0";
  const TString opt_nostack="hist nostack";

  stk_1d->Draw(opt_nostack);
  stk_1d->GetXaxis()->SetTitle(var);
  stk_1d->GetYaxis()->SetTitle("#frac{d#sigma}{dx_{Bj}} (cm^{2}/x_{Bj}/nucleon)");
  TString stk_tit = var+"of diff E modes";
  stk_1d->SetTitle(stk_tit);
  gStyle->SetTitleX(0.5);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  leg_1d->Draw();
  c1->Print("plot.pdf","x_{Bj} of diff E modes");
  string var_str{var};
  c1->Print(Form("png/%s_%s%s.png", anatag, var_str.c_str(), " of diff E modes"));

  
  c1->Print("plot.pdf]","dalphat of 4 event modess");
  fout->Close();
}
