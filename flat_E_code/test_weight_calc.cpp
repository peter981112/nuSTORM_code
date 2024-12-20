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

#include "/home/wonjong/GiBUU_2024/NuGenTKI/include/HistIO.h"
#include "/home/wonjong/GiBUU_2024/NuGenTKI/include/GeneratorUtils.h"
#include "/home/wonjong/GiBUU_2024/NuGenTKI/style/style.h"

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

double get_weight(vector<array<double,2>> weight_lst, double nu_E) //return weight at neutrino energy
{

  double dE_weight_lst = weight_lst[1][0]-weight_lst[0][0];
  double E0_weight_lst = weight_lst[0][0];
  int nu_E_index = static_cast<int>((nu_E-E0_weight_lst)/dE_weight_lst);
  double weight_nu_E{0};
  weight_nu_E = weight_lst[nu_E_index+1][1]*(nu_E-weight_lst[nu_E_index][0])+weight_lst[nu_E_index][1]*(weight_lst[nu_E_index+1][0]-nu_E);
  weight_nu_E = weight_nu_E/(weight_lst[nu_E_index+1][0]-weight_lst[nu_E_index][0]);

  cout<<"nu_E between "<<nu_E_index<<"th E in lst "<<weight_lst[nu_E_index][0]<<" and "<<(nu_E_index+1)<<"th E in lst "<<weight_lst[nu_E_index+1][0]<<endl;

  return weight_nu_E;

}

vector<array<double,2>> get_flux_weight_lst(vector<array<double,2>> target_flux_vect, vector<array<double,2>> flat_E_flux_vect)
{

  int len = target_flux_vect.size();
  vector<array<double,2>> weight_lst;
  weight_lst = target_flux_vect;
  double weight_coeff{0};
  weight_coeff = flat_E_flux_vect[flat_E_flux_vect.size()-1][0]/(integrate_const_interval(target_flux_vect));
  cout<<"weight_coeff : "<<weight_coeff<<endl;
  for(int i=0; i<len; i++)
  {
    weight_lst[i][1] = weight_coeff*target_flux_vect[i][1];
    if(i%10==0){cout<<i<<"th nu E : "<<weight_lst[i][0]<<", nu flux weight: "<<weight_lst[i][1]<<endl;}
  }
  return weight_lst;

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

void test_weight_calc()
{

  string flat_E_flux_file = "flat_E_test.txt";
  string flux_file_E7 = "E7SpectraMuSig556Numu.txt";
  vector<array<double,2>> target_flux_vect = read_flux_file(flux_file_E7);
  vector<array<double,2>> flat_E_flux_vect = read_flux_file(flat_E_flux_file);

  vector<array<double,2>> weight_lst = get_flux_weight_lst(target_flux_vect, flat_E_flux_vect);

  double test_weight_int{0};
  test_weight_int = integrate_const_interval(weight_lst);
  cout<<"weight integral : "<<test_weight_int<<endl;
  double test_nu_E{0};
  double test_weight{0};
  for(int i=0; i<20; i++)
  {
    test_nu_E = 0.1+i*0.15;
    cout<<i<<"th test nu E : "<<test_nu_E<<endl;
    test_weight = get_weight(weight_lst, test_nu_E);
    cout<<"test weight : "<<test_weight<<endl;
  }


   /*
  TString fn_lst = outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"rewighting outana9 - case with no pi 0"<<endl;

*/
  
  TString fn = "outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{2};
  int noff{1};
  auto anatag{"outana7"};
  cout<<"rewighting outana7 - case with pi 0"<<endl;

  TFile *f = new TFile(fn);
  TFile *fout = new TFile("tmpplot_outana7.root","recreate");

  TTree * theader = (TTree*) f->Get("header");
  TTree *t1 = (TTree*) f->Get("tree");
  TTree *new_t1, *new_theader;
  new_t1 = t1->CloneTree(0); 
  new_theader = theader->CloneTree(0);
  
  Double_t perweight;
  Int_t nrun;
  //Int_t prod,evtMode;

  t1->SetBranchAddress("perweight",&perweight);
  theader->SetBranchAddress("nrun",&nrun);

  Int_t nentries = (Int_t)t1->GetEntries();
  cout<<"number of events : "<<nentries<<endl;
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i); 
    if(i%200000==0){cout<<i<<"th weight : "<<perweight<<endl;}
    perweight = 1.0; // Make changes here 
    new_t1->Fill();
  }
  theader->GetEntry(0);
  new_theader->Fill();
  if (t1->GetEntries() == new_t1->GetEntries()) {
      t1->Delete();
      theader->Delete();
  };

  new_t1->Write();
  new_theader->Write();
  fout->Write();

  return;

}