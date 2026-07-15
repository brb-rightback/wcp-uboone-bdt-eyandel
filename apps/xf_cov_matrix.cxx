#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"
#include "WCPLEEANA/Util.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TVectorD.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./xf_cov_matrix -r[#sys 1-14]" << std::endl;
  }
  int run = 1; // run 1 ..
  int flag_spec_weights = 0;
  bool flag_osc = false;
  bool flag_v2h=false;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 's':
      flag_spec_weights = atoi(&argv[i][2]);
      break;
    case 'o':
      flag_osc = atoi(&argv[i][2]);
      break;
    case 'h':
      flag_v2h = atoi(&argv[i][2]);
      if(flag_v2h==1) std::cout<<"Will save histograms in addition to vectors"<<std::endl;
      break;
    }
  }

  TString xf_input_config_file = "./configurations/xf_input.txt";
  if (run>17) xf_input_config_file = "./configurations/rw_sys_input.txt";

  CovMatrix cov("./configurations/cov_input.txt", xf_input_config_file, "./configurations/xf_file_ch.txt", "./configurations/rw_cv_input.txt");

  if(run==3 || run==6){

    cov.add_disabled_ch_name("kdar_bdtsel_sig");
    cov.add_disabled_ch_name("kdar_bdtselns_sig");
    cov.add_disabled_ch_name("kdar_bdtselStub_sig");

    cov.add_disabled_ch_name("kdar_presel_sig");
    cov.add_disabled_ch_name("kdar_preselns_sig");
    cov.add_disabled_ch_name("kdar_preselStub_sig");

    cov.add_disabled_ch_name("kdar_antibdtsel_sig");
    cov.add_disabled_ch_name("kdar_antibdtselNp_sig");

    cov.add_disabled_ch_name("kdar_ktagged_sig");
    cov.add_disabled_ch_name("kdar_ktaggedNp_sig");

    cov.add_disabled_ch_name("kdar_lowE_bdtsel_sig");
    cov.add_disabled_ch_name("kdar_lowE_bdtselns_sig");

    cov.add_disabled_ch_name("kdar_hiE_bdtsel_sig");
    cov.add_disabled_ch_name("kdar_hiE_bdtselns_sig");
  }

  if (flag_osc) cov.add_osc_config();
  
  // special weights ...
  if (flag_spec_weights)    cov.init_spec_weights(2000,1000,0.2);
  
  if(run>17) cov.add_rw_config(run);
  cov.print_rw(cov.get_rw_info());

  // Get the file based on runno ...
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();
  // Construct the histogram ...

  // outfilename ...
  TString outfile_name;

  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::map<int, TH1F*> map_covch_hist;
  for (auto it = map_inputfile_info.begin(); it!=map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);

    
    if (period == run){
      outfile_name = out_filename;
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename, 0);
      
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	TString histoname = std::get<0>(*it1);
	Int_t nbin = std::get<1>(*it1);
	float llimit = std::get<2>(*it1);
	float hlimit = std::get<3>(*it1);
	TString var_name = std::get<4>(*it1);
	TString ch_name = std::get<5>(*it1);
	TString add_cut = std::get<6>(*it1);
	TString weight = std::get<7>(*it1);
	
	//	std::cout << std::get<0>(*it1)  << " " << std::get<1>(*it1) << " " << std::get<4>(*it1) << " " << std::get<5>(*it1) << " " << std::get<6>(*it1) << " " << std::get<7>(*it1) << std::endl;
	htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
	map_histoname_hist[histoname] = htemp;

	int covch = cov.get_covch_name(ch_name);
	if (map_covch_hist.find(covch) == map_covch_hist.end()){
	  TH1F *htemp1 = (TH1F*)htemp->Clone(Form("pred_covch_%d",covch));
	  map_covch_hist[covch] = htemp1;
	}
      }
      //  std::cout << input_filename << " " << filetype << " " << out_filename << std::endl; 
    }
  }
  std::cout << outfile_name << std::endl;
  TMatrixD* cov_add_mat = cov.get_add_cov_matrix();
  // create a covariance matrix for bootstrapping ...
  TMatrixD* cov_xf_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  TVectorD* vec_mean = new TVectorD(cov_add_mat->GetNrows());

  //std::cout << "test: " << std::endl;

  cov.gen_xf_cov_matrix(run, map_covch_hist, map_histoname_hist, vec_mean, cov_xf_mat);
  
  TMatrixD* frac_cov_xf_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  for (size_t i=0; i!= frac_cov_xf_mat->GetNrows(); i++){
    double val_1 = (*vec_mean)(i);
    for (size_t j=0; j!=frac_cov_xf_mat->GetNrows();j++){
      double val_2 = (*vec_mean)(j);
      double val = (*cov_xf_mat)(i,j);
      if (val_1 ==0 && val_2 == 0){
	(*frac_cov_xf_mat)(i,j) = 0;
      }else if (val_1 ==0 || val_2 ==0){
	if (val !=0){
	  if (i==j){
	    (*frac_cov_xf_mat)(i,j) = 0.;
	  }else{
	    (*frac_cov_xf_mat)(i,j) = 0;
	  }
	}else{
	  (*frac_cov_xf_mat)(i,j) = 0;
	}
      }else{
	(*frac_cov_xf_mat)(i,j) = val/val_1/val_2;
      }
    }
  }

  TFile *file = new TFile(outfile_name,"RECREATE");
  vec_mean->Write(Form("vec_mean_%d",run));
  cov_xf_mat->Write(Form("cov_xf_mat_%d",run));
  frac_cov_xf_mat->Write(Form("frac_cov_xf_mat_%d",run));


  if(flag_v2h){
    std::cout<<"Saving histograms in addition to vectors"<<std::endl;
    int nbin = vec_mean->GetNrows();
    TH1D* h_mean = new TH1D(Form("h_mean_%d",run), Form("h_mean_%d",run), nbin, 0, nbin);
    V2H((*vec_mean), h_mean);
    h_mean->Write();

    int nbinx = cov_xf_mat->GetNrows();
    int nbiny = cov_xf_mat->GetNcols();
    TH2D* h_cov_xf = new TH2D(Form("h_cov_xf_%d",run), Form("h_cov_xf_%d",run), nbinx, 0, nbinx, nbiny, 0, nbiny);
    M2H((*cov_xf_mat), h_cov_xf);
    h_cov_xf->Write();

    nbinx = frac_cov_xf_mat->GetNrows();
    nbiny = frac_cov_xf_mat->GetNcols();
    TH2D* h_frac_cov_xf = new TH2D(Form("h_frac_cov_xf_%d",run), Form("h_frac_cov_xf_%d",run), nbinx, 0, nbinx, nbiny, 0, nbiny);
    M2H((*frac_cov_xf_mat), h_frac_cov_xf);
    h_frac_cov_xf->Write();
   }


  // save central ... results ...
  // for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
  //  ((TH1F*)it->second)->SetDirectory(file);
  // }
  
  for (auto it = map_covch_hist.begin(); it != map_covch_hist.end(); it++){
    ((TH1F*)it->second)->SetDirectory(file);
  }
  
  file->Write();
  file->Close();
  
  return 0;
}
