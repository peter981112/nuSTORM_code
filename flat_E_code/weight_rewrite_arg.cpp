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

#include "HistIO.h"
#include "GeneratorUtils.h"
#include "style.h"

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
  if(nu_E<weight_lst[weight_lst.size()-1][0])
  {
    weight_nu_E = weight_lst[nu_E_index+1][1]*(nu_E-weight_lst[nu_E_index][0])+weight_lst[nu_E_index][1]*(weight_lst[nu_E_index+1][0]-nu_E);
    weight_nu_E = weight_nu_E/(weight_lst[nu_E_index+1][0]-weight_lst[nu_E_index][0]);
  }
  else{weight_nu_E=0;}

  //cout<<"nu_E between "<<nu_E_index<<"th E in lst "<<weight_lst[nu_E_index][0]<<" and "<<(nu_E_index+1)<<"th E in lst "<<weight_lst[nu_E_index+1][0]<<endl;

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
    //if(i%10==0){cout<<i<<"th nu E : "<<weight_lst[i][0]<<", nu flux weight: "<<weight_lst[i][1]<<endl;}
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
      //if(totlinecount%20==0){cout<<totlinecount<<"th nu E : "<<target_flux_vect[totlinecount][0]<<", nu flux : "<<target_flux_vect[totlinecount][1]<<endl;}
      totlinecount++;
    }
    cout << "Number of lines in the file: " << totlinecount << endl;
  }

  return target_flux_vect;
}

double get_hist_norm_per_nucleon(TString fn)
{
  TFile *f = new TFile(fn);
  TTree * theader = (TTree*) f->Get("header");

  double GiBUUnormFactor=-999;
  if(theader)
  {
    int tmpanr;
    theader->SetBranchAddress("nrun", &tmpanr);
    theader->GetEntry(0);
    if(tmpanr>0){
      GiBUUnormFactor=tmpanr;
    }
  }
  cout<<"=================================================getting norm factor for xsec================================================="<<endl;
  cout<<"nrun for file - "<<fn<<" is : "<<GiBUUnormFactor<<endl;
  GiBUUnormFactor = GiBUUnormFactor*1e38;
  cout<<" norm factor per nucleon : "<<GiBUUnormFactor<<endl;

  return GiBUUnormFactor;
}

void weight_rewrite(string flat_E_flux_file, string flux_file, TString fn)
{

  //string flat_E_flux_file = "flat_E_test.txt";
  //string flux_file_E7 = "E7SpectraMuSig556Numu.txt";
  vector<array<double,2>> target_flux_vect = read_flux_file(flux_file);
  vector<array<double,2>> flat_E_flux_vect = read_flux_file(flat_E_flux_file);

  vector<array<double,2>> weight_lst = get_flux_weight_lst(target_flux_vect, flat_E_flux_vect);

  double test_weight_int{0};
  test_weight_int = integrate_const_interval(weight_lst);
  cout<<"weight integral : "<<test_weight_int<<endl;


  int anaid{0};auto anatag{"outana"};string fn_string{fn};string path_string{fn};string flux_file_no_txt{flux_file};TString fout_name{""};//int noff;

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
  while (true)
  {
    if(fn_string.find("/") != fn_string.npos)
    {
      slash_index = fn_string.find("/");
      fn_string.erase(fn_string.begin(), fn_string.begin()+slash_index+1);//remove file path and extention .root from the file name
    }
    else{break;}
  }
  cout<<"fn_string : "<<fn_string<<endl;

  cout<<"flux_file_no_txt : "<<flux_file_no_txt<<endl;
  dot_index = flux_file_no_txt.find(".");
  flux_file_no_txt.erase(flux_file_no_txt.begin() + dot_index, flux_file_no_txt.end());//remove extention .txt from the file name
  cout<<"flux_file_no_txt : "<<flux_file_no_txt<<endl;
  while (true)
  {
    if(flux_file_no_txt.find("/") != flux_file_no_txt.npos)
    {
      slash_index = flux_file_no_txt.find("/");
      flux_file_no_txt.erase(flux_file_no_txt.begin(), flux_file_no_txt.begin()+slash_index+1);//remove file path and extention .txt from the file name
    }
    else{break;}
  }

  cout<<"path_string : "<<path_string<<endl;
  for (int i=0; i<3; i++)
  {
    slash_index = path_string.rfind("/");
    path_string.erase(path_string.begin() + slash_index, path_string.end());
  }



  string destination_folder{"root_target/"}; string imitate{"_imitate"};

  fout_name = path_string + "/" + destination_folder + flux_file_no_txt + "_" + fn_string + imitate + ".root";
  cout<<"fn_string : "<<fn_string<<endl<<"flux_file_no_txt : "<<flux_file_no_txt<<endl<<"path_string : "<<path_string<<endl<<"name of fout : "<<fout_name<<endl;
   /*
  TString fn_lst = outAna9_MINERvALE_GiBUU_test_GFS0PIa9nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{1};
  int noff{2};
  auto anatag{"outana9"};
  cout<<"rewighting outana9 - case with no pi 0"<<endl;


  
  TString fn = "outAna7_MINERvALE_GiBUU_test_GFSPIZEROa7nuC_Filelist_GiBUUMINERvALE_nu_T0_Carbon.root";
  int anaid{2};
  int noff{1};
  auto anatag{"outana7"};
  cout<<"rewighting outana7 - case with pi 0"<<endl;
  */

  TFile *f = new TFile(fn);
  TFile *fout = new TFile(fout_name,"recreate");

  TString nuExp;
  int nuPDG, tarA;
  const TString gen = GeneratorUtils::SetEnv(fn, nuExp, nuPDG, tarA);
  double histnormFactor = GeneratorUtils::GetHistNormPerNucleus(f, nuPDG, tarA, gen, anaid);//*13;//Scaling to per nucleon in CH 13=1+12
  cout<<"histnormFactor without nucleon num scaling"<<histnormFactor<<endl;

  TTree * theader = (TTree*) f->Get("header");
  TTree *t1 = (TTree*) f->Get("tree");
  TTree *new_t1, *new_theader;
  new_t1 = t1->CloneTree(0); 
  new_theader = theader->CloneTree(0);
  
  Double_t perweight, beamE;
  Int_t nrun;
  //Int_t prod,evtMode;

  t1->SetBranchAddress("perweight",&perweight);
  t1->SetBranchAddress("beamE",&beamE);
  theader->SetBranchAddress("nrun",&nrun);

  Int_t nentries = (Int_t)t1->GetEntries();
  cout<<"number of events : "<<nentries<<endl;
  double nu_E_weight{0};
  for (Int_t i=0; i<nentries; i++)
  {
    t1->GetEntry(i); 
    if(i%1000000==0){cout<<i<<"th weight : "<<perweight<<endl;}
    nu_E_weight = get_weight(weight_lst, beamE);
    if(nu_E_weight==0){perweight=0;}
    else
    {
      perweight = perweight*nu_E_weight; // Make changes here 
      if(i%1000000==0){cout<<i<<"th altered weight : "<<perweight<<" with nu_E weight"<<nu_E_weight<<" at nu E : "<<beamE<<endl;}
      new_t1->Fill();
    }
  }
  theader->GetEntry(0);
  new_theader->Fill();
  /*
  if (t1->GetEntries() == new_t1->GetEntries()) {
      t1->Delete();
      theader->Delete();
  };
*/
  cout<<"finished weight recalc. now writing to root file"<<endl;
  new_t1->Write();
  new_theader->Write();
  fout->Write();
  cout<<"done rewrite"<<endl;

  return;

}

int main(int argc, char* argv[])
{
  //void anaGenerator(const TString tag, const TString filelist, const int tmpana, const int nToStop=-999, const int smearBit=0)
  if(argc==4){
    weight_rewrite(argv[1], argv[2], argv[3]);
  }
  else{
    printf("wrong argc %d\n", argc); return 1;
  }
  return 0;
}
