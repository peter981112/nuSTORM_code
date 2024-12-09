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
      //if(totlinecount%20==0){cout<<totlinecount<<"th nu E : "<<target_flux_vect[totlinecount][0]<<", nu flux : "<<target_flux_vect[totlinecount][1]<<endl;}
      totlinecount++;
    }
    cout << "Number of lines in the file: " << totlinecount << endl;
  }

  return target_flux_vect;
}

void test_event_rate()
{

  string flux_file_E7 = "E7SpectraMuSig556Numu.txt";
  vector<array<double,2>> target_flux_vect = read_flux_file(flux_file_E7);

  double flux_int{0};
  flux_int = integrate_const_interval(target_flux_vect);
  cout<<"flux integral : "<<flux_int<<endl;

  double amu = 1.67e-27; //mass of nucleon in kg
  double time_interval = 5.0*365*24*3600; //time interval of 5 year in sec


   /*
  TString fn_lst = "outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"rewighting outana9 - case with no pi 0"<<endl;

*/
  
  TString fn = "E7SpectraMuSig556Numu_numu_outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon_imitate.root";
  int anaid{2};
  int noff{1};
  auto anatag{"outana7"};
  cout<<"rewighting outana7 - case with pi 0"<<endl;

  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot_outana7.root","recreate");
  TTree *t1 = (TTree*) f->Get("tree");

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  //const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*13;cout<<"scaling 1/13 for CH"<<endl//Scaling to per nucleon in CH 13=1+12
  const double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid)*12;cout<<"scaling 1/12 for carbon"<<endl;//Scaling to per nucleon in H2O 18=2+16
  
  Double_t perweight;
  Int_t prod,evtMode;

  t1->SetBranchAddress("perweight",&perweight);
  t1->SetBranchAddress("prod",&prod);
  t1->SetBranchAddress("evtMode",&evtMode);


  Int_t nentries = (Int_t)t1->GetEntries();
  cout<<"number of events : "<<nentries<<endl;

  double QE_xsec{0}, RES_xsec{0}, DIS_xsec{0}, a2p2h_xsec{0}, tot_xsec{0};
  double QE_event_num{0}, RES_event_numc{0}, DIS_event_num{0}, a2p2h_event_num{0}, tot_event_num{0};
  int npi_err{0};
  string mode[] = {"QE", "RES", "DIS", "2p2h"};

  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i); 
    if(i%200000==0){cout<<i<<"th weight : "<<perweight<<endl;}
    if(perweight>4000){
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

      tot_xsec = tot_xsec + perweight;

      if(prod==1 && evtMode==1){ QE_xsec = QE_xsec + perweight;}//get sum of weights for quasi elastic case
      else if(prod>=2 && prod<=33 && evtMode==2){RES_xsec = RES_xsec + perweight;}//get sum of weights for resonance case
      else if(prod==34 && evtMode==3){DIS_xsec = DIS_xsec + perweight;}//get sum of weights for deep inelastic case
      else if(prod==35 && evtMode==4) {a2p2h_xsec = a2p2h_xsec + perweight;}//get sum of weights for 2 part 2 hole case
    }
  }

  double xsec_lst_dir[] = {QE_xsec, RES_xsec, DIS_xsec, a2p2h_xsec};
  tot_xsec = tot_xsec/histnormFactor;
  cout<<"direct perweight calc total xsec : "<<tot_xsec<<" cm^{2}"<<endl;
  for(int ii=0; ii<4; ii++)
  {
    xsec_lst_dir[ii] = xsec_lst_dir[ii]/histnormFactor;
    cout<<"direct perweight calc xsec for mode : "<<mode[ii]<<" is "<<xsec_lst_dir[ii]<<" cm^{2}"<<endl;
  }

  double event_num_lst[] = {QE_event_num, RES_event_numc, DIS_event_num, a2p2h_event_num};
  tot_event_num = tot_xsec*time_interval*flux_int/amu;
  cout<<"event num in "<<time_interval/(365*24*3600)<<" years : "<<tot_event_num<<" event /POT * kg"<<endl;
    for(int ii=0; ii<4; ii++)
  {
    event_num_lst[ii] = xsec_lst_dir[ii]*time_interval*flux_int/amu;
    cout<<"event num in "<<time_interval/(365*24*3600)<<" years : "<<event_num_lst[ii]<<" event /POT * kg"<<endl;
  }

  return;
}