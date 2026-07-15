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
    std::cout << "./det_cov_matrix -r[#sys 1-10] -g[#gp 0,1]" << std::endl;
  }
  int run = 1; // run 1 ...
  bool flag_osc = false;
  int flag_gp = 0; // gaussian process smoothing
  double seed = 0;
  bool ignore_ns=false;
  bool flag_v2h=false;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 'o':
      flag_osc = atoi(&argv[i][2]); // which run period
      break;
    case 'g':
      flag_gp = atoi(&argv[i][2]); // 1: on, 0: off.  2 = off but do GPSmoothing debugging, 3 = smoothing and debugging
      break;
     case 's':
      seed = atoi(&argv[i][2]); //Set the seed for bootsrapping
      std::cout<<"Setting bootstrapping seed to "<<seed<<std::endl;
      break;
     case 'n':
      ignore_ns = atoi(&argv[i][2]); //Set the seed for bootsrapping
      std::cout<<"Setting flag it ignore ns timing in detvar calculation to "<<ignore_ns<<std::endl;
      break;
    case 'h':
      flag_v2h = atoi(&argv[i][2]);
      if(flag_v2h==1) std::cout<<"Will save histograms in addition to vectors"<<std::endl;
      break;
    }
  }

  CovMatrix cov("./configurations/cov_input.txt", "./configurations/det_input.txt", "./configurations/det_file_ch.txt", "./configurations/rw_cv_input.txt");
  // cov.add_disabled_ch_name("nueCC_FC_nueoverlay");
  // cov.add_disabled_ch_name("nueCC_PC_nueoverlay");
  // cov.add_disabled_ch_name("LEE_FC_nueoverlay");
  // cov.add_disabled_ch_name("LEE_PC_nueoverlay");

  cov.add_xs_config();

  if (flag_osc) cov.add_osc_config();

  cov.print_rw(cov.get_rw_info());

  // Get the file based on runno ...
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();
  // Construct the histogram ...

  // outfilename ...
  TString outfile_name;

  TH1F *htemp; TH1F *htemp1;
  std::map<TString, std::tuple<TH1F*, TH1F*, int>> map_histoname_hists;
  std::map<int, std::tuple<TH1F*, TH1F*, int>> map_covch_hists;

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
        TString temp_histoname = histoname + "_CV";
	htemp = new TH1F(temp_histoname, temp_histoname, nbin, llimit, hlimit);
	temp_histoname = histoname + "_DET";
        htemp1 = new TH1F(temp_histoname, temp_histoname, nbin, llimit, hlimit);
        int htemp_n = 1+cov.is_xs_chname(ch_name);
	int covch = cov.get_covch_name(ch_name);
        map_histoname_hists[histoname] = std::make_tuple(htemp, htemp1, htemp_n);

	if (map_covch_hists.find(covch) == map_covch_hists.end()){
	  TH1F *htemp_1 = (TH1F*)htemp->Clone(Form("pred_covch_CV_%d",covch));
          TH1F *htemp1_1 = (TH1F*)htemp1->Clone(Form("pred_covch_DET_%d",covch));
          map_covch_hists[covch] = std::make_tuple(htemp_1, htemp1_1, htemp_n);
	}
      }
      //  std::cout << input_filename << " " << filetype << " " << out_filename << std::endl;
    }
  }
  std::cout << outfile_name << std::endl;

  TMatrixD* cov_add_mat = cov.get_add_cov_matrix();
  // create a covariance matrix for bootstrapping ...
  TMatrixD* cov_mat_bootstrapping = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  // create a covariance matrix for det systematics ...
  TMatrixD* cov_det_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  TVectorD* vec_mean_diff = new TVectorD(cov_add_mat->GetNrows());
  TVectorD* vec_mean = new TVectorD(cov_add_mat->GetNrows());

  cov.gen_det_cov_matrix_norm(run, map_covch_hists, map_histoname_hists, vec_mean, vec_mean_diff, cov_mat_bootstrapping, cov_det_mat, flag_gp, seed, ignore_ns);

  TMatrixD* frac_cov_det_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  for (size_t i=0; i!= frac_cov_det_mat->GetNrows(); i++){
    double val_1 = (*vec_mean)(i);
    for (size_t j=0; j!=frac_cov_det_mat->GetNrows();j++){
      double val_2 = (*vec_mean)(j);
      double val = (*cov_det_mat)(i,j);
      if (val_1 ==0 && val_2 == 0){
	(*frac_cov_det_mat)(i,j) = 0;
      }else if (val_1 ==0 || val_2 ==0){
	if (val !=0){
	  if (i==j){
	    (*frac_cov_det_mat)(i,j) = 1./16.; // 25% uncertainties ...
	  }else{
	    (*frac_cov_det_mat)(i,j) = 0;
	  }
	}else{
	  (*frac_cov_det_mat)(i,j) = 0;
	}
      }else{
	(*frac_cov_det_mat)(i,j) = val/val_1/val_2;
      }
    }
  }


  TFile *file = new TFile(outfile_name,"RECREATE");
  vec_mean->Write(Form("vec_mean_%d",run));
  vec_mean_diff->Write(Form("vec_mean_diff_%d",run));

  cov_mat_bootstrapping->Write(Form("cov_mat_boostrapping_%d",run));
  cov_det_mat->Write(Form("cov_det_mat_%d",run));
  frac_cov_det_mat->Write(Form("frac_cov_det_mat_%d",run));

  TVectorD vseed(1,1,seed,"END");
  vseed.Write("seed");

  if(flag_v2h){
    std::cout<<"Saving histograms in addition to vectors"<<std::endl;
    int nbin = vec_mean->GetNrows();
    TH1D* h_mean = new TH1D(Form("h_mean_%d",run), Form("h_mean_%d",run), nbin, 0, nbin);
    V2H((*vec_mean), h_mean);
    h_mean->Write();

    nbin = vec_mean_diff->GetNrows();
    TH1D* h_mean_diff = new TH1D(Form("h_mean_diff_%d",run), Form("h_mean_diff_%d",run), nbin, 0, nbin);
    V2H((*vec_mean_diff), h_mean_diff);
    h_mean_diff->Write();

    int nbinx = cov_det_mat->GetNrows();
    int nbiny = cov_det_mat->GetNcols();
    TH2D* h_cov_det = new TH2D(Form("h_cov_det_%d",run), Form("h_cov_det_%d",run), nbinx, 0, nbinx, nbiny, 0, nbiny);
    M2H((*cov_det_mat), h_cov_det);
    h_cov_det->Write();

    nbinx = frac_cov_det_mat->GetNrows();
    nbiny = frac_cov_det_mat->GetNcols();
    TH2D* h_frac_cov_det = new TH2D(Form("h_frac_cov_det_%d",run), Form("h_frac_cov_det_%d",run), nbinx, 0, nbinx, nbiny, 0, nbiny);
    M2H((*frac_cov_det_mat), h_frac_cov_det);
    h_frac_cov_det->Write();

    nbinx = cov_mat_bootstrapping->GetNrows();
    nbiny = cov_mat_bootstrapping->GetNcols();
    TH2D* h_cov_bootstrapping = new TH2D(Form("h_cov_bootstrapping_%d",run), Form("h_cov_bootstrapping_%d",run), nbinx, 0, nbinx, nbiny, 0, nbiny);
    M2H((*cov_mat_bootstrapping), h_cov_bootstrapping);
    h_cov_bootstrapping->Write();
   }

  for (auto it = map_covch_hists.begin(); it != map_covch_hists.end(); it++){
    auto results = it->second;
    TH1F *h1 = std::get<0>(results);
    TH1F *h2 = std::get<1>(results);
    h1->SetDirectory(file);
    h2->SetDirectory(file);
  }

  file->Write();
  file->Close();
}
