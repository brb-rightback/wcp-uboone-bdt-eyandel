#include "GPSmoothing.C"

void LEEana::CovMatrix::add_disabled_ch_name(TString name){
  disabled_ch_names.insert(name);
}
void LEEana::CovMatrix::remove_disabled_ch_name(TString name){
  disabled_ch_names.erase(name);
}

void LEEana::CovMatrix::gen_det_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TVectorD* vec_mean_diff, TMatrixD* cov_mat_bootstrapping, TMatrixD* cov_det_mat, int flag_gp=0, int seed=0){


  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ;
  std::map<int, TString> map_no_histoname;

  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);

    if (period != run) continue;
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

    for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
      int ch = map_name_ch[std::get<5>(*it1)];
      int obsch = get_obsch_name(std::get<5>(*it1));
      int covch = get_covch_name(std::get<5>(*it1));
      int flag_lee = std::get<7>(map_ch_hist[ch]);
      TString histoname = std::get<0>(*it1);
      TH1F *htemp = map_histoname_hist[histoname];
      //
      map_histoname_infos[histoname] = std::make_tuple(ncount, covch, flag_lee, input_filename);
      map_no_histoname[ncount] = histoname;
      ncount ++;

      //std::cout << histoname << obsch << " " << covch << " " << flag_lee << std::endl;
    }
  }

  // histogram and LEE flag of each histogram number, so the fill functions do not look them up by name for every entry
  std::vector<TH1F*> vec_no_hist(ncount, 0);
  std::vector<int> vec_no_lee(ncount, 0);
  for (auto it = map_no_histoname.begin(); it != map_no_histoname.end(); it++){
    auto it1 = map_histoname_hist.find(it->second);
    if (it1 != map_histoname_hist.end()) vec_no_hist.at(it->first) = it1->second;
    vec_no_lee.at(it->first) = std::get<2>(map_histoname_infos[it->second]);
  }

  // results ... filename --> re --> variable, weight, lee weight,
  std::map<TString, std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > > map_all_events;
  std::map<TString, double> map_filename_pot;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    //int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    if (period != run) continue;

    //map_all_events[input_filename];
    get_events_info(input_filename, map_all_events, map_filename_pot, map_histoname_infos);
  }

  double data_pot = 5e19; // example ...

  const int rows = cov_mat_bootstrapping->GetNcols();
  TPrincipal prin(rows, "ND");
  Double_t *x = new Double_t[rows];

  std::map<TString, TH1D*> map_filename_histo;
  // form histogram ...
  for (auto it = map_all_events.begin(); it != map_all_events.end(); it++){
    TString filename = it->first;
    int nsize = it->second.size();
    TH1D* htemp = new TH1D(filename, filename, nsize, 0.5, nsize+0.5);
    for (size_t i=0;i!=nsize;i++){
      htemp->SetBinContent(i+1, std::get<2>(it->second.at(i)) );
    }
    //std::cout << htemp->GetSum() << std::endl;
    map_filename_histo[filename] = htemp;
  }

  // working on the bootstrapping ...
  for (int qx = 0; qx != 1000; qx++){
    if (qx % 500 ==0) std::cout << qx << std::endl;

    for (int i=0;i!=rows;i++){
      x[i] = 0;
    }

    // fill the histogram with CV
    fill_det_histograms(map_filename_histo, map_all_events, vec_no_hist, vec_no_lee, map_histoname_hist);
    // merge histograms according to POTs ...
    for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
      //std::cout << it->first << std::endl;
      int covch = it->first;
      TH1F *hpred = map_covch_hist[covch];
      hpred->Reset();

      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	TH1F *htemp = (TH1F*)hpred->Clone("htemp");
	htemp->Reset();
	std::map<int, double> temp_map_mc_acc_pot;

	for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	  TString histoname = (*it2).first;
	  TString input_filename = map_histogram_inputfile[histoname];
	  auto it3 = map_inputfile_info.find(input_filename);
	  int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	  int norm_period = std::get<6>(it3->second);
	  double mc_pot = map_filename_pot[input_filename];
	  //std::cout << mc_pot << std::endl;
	  if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
	    temp_map_mc_acc_pot[norm_period] = mc_pot;
	  }else{
	    temp_map_mc_acc_pot[norm_period] += mc_pot;
	  }
	}

	for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	  TString histoname = (*it2).first;
	  TString input_filename = map_histogram_inputfile[histoname];
	  auto it3 = map_inputfile_info.find(input_filename);
	  int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	  int norm_period = std::get<6>(it3->second);
	  data_pot = std::get<5>(map_inputfile_info[input_filename]);
	  double ratio = data_pot/temp_map_mc_acc_pot[norm_period];

	  TH1F *hmc = map_histoname_hist[histoname];
	  htemp->Add(hmc, ratio);
	  //	std::cout << covch << " " << histoname << " " << ratio << std::endl;
	}

	hpred->Add(htemp);
	delete htemp;
      }

      int start_bin = map_covch_startbin[covch];
      for (int i=0;i!=hpred->GetNbinsX()+1;i++){
	x[start_bin+i] = hpred->GetBinContent(i+1) ;
	//std::cout << x[start_bin+i] << std::endl;
      }

    }

    prin.AddRow(x);

  }

  (*cov_mat_bootstrapping) = (*(TMatrixD*)prin.GetCovarianceMatrix());
  for (int i=0;i!=rows;i++){
    for (int j=0;j!=rows;j++){
      if (i<j) { (*cov_mat_bootstrapping)(i,j) = (*(TMatrixD*)prin.GetCovarianceMatrix())(j,i); }
    }
  }
  *vec_mean_diff = (*prin.GetMeanValues());

  //Do GP Smoothing Here
  GPSmoothing(vec_mean_diff, cov_mat_bootstrapping, "./configurations/gp_input.txt", flag_gp);

  // Now get the full covariance matrix ...
  TMatrixDSym DMatrix(rows);
  for (int i=0;i!=rows;i++){
    for (int j=0;j!=rows;j++){
      DMatrix(i,j) =  (*cov_mat_bootstrapping)(i,j);
    }
  }
  TMatrixDSymEigen DMatrix_eigen(DMatrix);
  TMatrixD matrix_eigenvector = DMatrix_eigen.GetEigenVectors();
  TMatrixD matrix_eigenvector_T(rows,rows);
  matrix_eigenvector_T.Transpose(matrix_eigenvector);
  TVectorD matrix_eigenvalue = DMatrix_eigen.GetEigenValues();

  TPrincipal prin_full(rows, "ND");
  TRandom3 random3(seed);
  for (int i=0;i!=16000;i++){
    TMatrixD matrix_element(rows,1);
    for (int j=0;j!=rows;j++){
      if (matrix_eigenvalue(j) >=0)
	matrix_element(j,0) = random3.Gaus(0,sqrt(matrix_eigenvalue(j)));
      else
	matrix_element(j,0) = 0;
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;
    double rel_err = random3.Gaus(0,1);
    for (int j=0;j!=rows;j++){
      //      matrix_variation(j,0) = matrix_variation(j,0)/sqrt(11.) + (*vec_mean_diff)(j); // increase MC stat by a factor of 11  x11
      matrix_variation(j,0) += (*vec_mean_diff)(j); // standard ...
      //matrix_variation(j,0) = (*vec_mean_diff)(j); // no random term
      x[j] = rel_err * matrix_variation(j,0);
      //x[j] = matrix_variation(j,0); //original no abs term
    }
    prin_full.AddRow(x);

  }
  (*cov_det_mat) = (*(TMatrixD*)prin_full.GetCovarianceMatrix());
  for (int i=0;i!=rows;i++){
    for (int j=0;j!=rows;j++){
      if (i<j) (*cov_det_mat)(i,j) = (*(TMatrixD*)prin_full.GetCovarianceMatrix())(j,i);
    }
  }


  delete[] x;

  // clean up the memory ...
   for (auto it = map_filename_histo.begin(); it != map_filename_histo.end(); it++){
     delete it->second;
   }

   // fill the histogram with CV
   fill_det_histograms(map_all_events, vec_no_hist, vec_no_lee, map_histoname_hist);
   // merge histograms according to POTs ...
   for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
     //std::cout << it->first << std::endl;
     int covch = it->first;
     TH1F *hpred = map_covch_hist[covch];
     hpred->Reset();

     for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
       TH1F *htemp = (TH1F*)hpred->Clone("htemp");
       htemp->Reset();
       std::map<int, double> temp_map_mc_acc_pot;

       for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	 TString histoname = (*it2).first;
	 TString input_filename = map_histogram_inputfile[histoname];
	 auto it3 = map_inputfile_info.find(input_filename);
	 int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	 int norm_period = std::get<6>(it3->second);
	 double mc_pot = map_filename_pot[input_filename];
	 //std::cout << mc_pot << std::endl;
	 if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
	   temp_map_mc_acc_pot[norm_period] = mc_pot;
	 }else{
	   temp_map_mc_acc_pot[norm_period] += mc_pot;
	 }
       }

       for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	 TString histoname = (*it2).first;
	 TString input_filename = map_histogram_inputfile[histoname];
	 auto it3 = map_inputfile_info.find(input_filename);
	 int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	 int norm_period = std::get<6>(it3->second);
	 data_pot = std::get<5>(map_inputfile_info[input_filename]);
	 double ratio = data_pot/temp_map_mc_acc_pot[norm_period];

	 TH1F *hmc = map_histoname_hist[histoname];
	 htemp->Add(hmc, ratio);
	 //	std::cout << covch << " " << histoname << " " << ratio << std::endl;
       }

       hpred->Add(htemp);
       delete htemp;
     }

    int start_bin = map_covch_startbin[covch];
    for (int i=0;i!=hpred->GetNbinsX()+1;i++){
      (*vec_mean)[start_bin+i] = hpred->GetBinContent(i+1);

      //std::cout << start_bin+i << " " << (*vec_mean_diff)(start_bin+i) << " " <<  hpred->GetBinContent(i+1) << std::endl;;
	//std::cout << x[start_bin+i] << std::endl;
    }

  }



}

void LEEana::CovMatrix::fill_det_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > >&map_all_events, std::vector<TH1F*>& vec_no_hist, std::vector<int>& vec_no_lee,  std::map<TString, TH1F*>& map_histoname_hist){
  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
     it->second->Reset();
   }

  // loop over files
  for (auto it = map_all_events.begin(); it!=map_all_events.end(); it++){
    TString filename = it->first;
    TH1D *hweight = map_filename_histo[filename];
    double sum = hweight->GetSum();
    for (size_t i=0;i<sum;i++){
      int global_index = hweight->FindBin(hweight->GetRandom())-1;
      double weight = std::get<2>(it->second.at(global_index));
      double weight_lee = std::get<3>(it->second.at(global_index));

      for (auto it1 = std::get<4>(it->second.at(global_index)).begin(); it1 != std::get<4>(it->second.at(global_index)).end(); it1++){
	   int no = std::get<0>(*it1);
	   double val_cv = std::get<1>(*it1);
	   bool flag_cv = std::get<2>(*it1);
	   double val_det = std::get<3>(*it1);
	   bool flag_det = std::get<4>(*it1);

       	   TH1F *htemp = vec_no_hist[no];
	   int flag_lee = vec_no_lee[no];

	   // central value ...
	   if (flag_cv){
	     if (flag_lee){
	       htemp->Fill(val_cv, -weight_lee);
	     }else{
	       htemp->Fill(val_cv, -1);
	     }
	   }
	   if (flag_det){
	     if (flag_lee){
	       htemp->Fill(val_det, weight_lee);
	     }else{
	       htemp->Fill(val_det, 1);
	     }
	   }

	   // if (no==2)
	   //std::cout << std::get<0>(it->second.at(i)) << " " << std::get<1>(it->second.at(i)) << " " << val_cv << " " << weight << std::endl;
	   // std::cout << weight << " " << weight_lee << " " << flag_lee << " " << histoname << std::endl;
	 }

    }
  }


}



 void LEEana::CovMatrix::fill_det_histograms(std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > >&map_all_events, std::vector<TH1F*>& vec_no_hist, std::vector<int>& vec_no_lee,  std::map<TString, TH1F*>& map_histoname_hist){

   for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
     it->second->Reset();
   }

   // fill central value ...

     // loop over files
     for (auto it = map_all_events.begin(); it!=map_all_events.end(); it++){
       // loop over events ...
       //std::cout << it->first << " " << it->second.size() << std::endl;

       for (size_t i=0;i!=it->second.size(); i++){
	 double weight = std::get<2>(it->second.at(i));
	 double weight_lee = std::get<3>(it->second.at(i));

	 //	 std::cout << i << " " << weight << " " << weight_lee << " " << std::get<4>(it->second.at(i)).size() << std::endl;

	 for (auto it1 = std::get<4>(it->second.at(i)).begin(); it1 != std::get<4>(it->second.at(i)).end(); it1++){
	   int no = std::get<0>(*it1);
	   double val_cv = std::get<1>(*it1);
	   bool flag_cv = std::get<2>(*it1);
	   double val_det = std::get<3>(*it1);
	   bool flag_det = std::get<4>(*it1);

       	   TH1F *htemp = vec_no_hist[no];
	   int flag_lee = vec_no_lee[no];

	   // central value ...
	   if (flag_cv){
	     if (flag_lee){
	       htemp->Fill(val_cv, weight * weight_lee);
	     }else{
	       htemp->Fill(val_cv, weight);
	     }
	   }


	   // if (flag_cv){
	   //   if (flag_lee){
	   //     htemp->Fill(val_cv, -weight*weight_lee);
	   //   }else{
	   //     htemp->Fill(val_cv, -weight);
	   //   }
	   // }
	   // if (flag_det){
	   //   if (flag_lee){
	   //     htemp->Fill(val_det, weight*weight_lee);
	   //   }else{
	   //     htemp->Fill(val_det, weight);
	   //   }
	   // }

	   // if (no==2)
	   //std::cout << std::get<0>(it->second.at(i)) << " " << std::get<1>(it->second.at(i)) << " " << val_cv << " " << weight << std::endl;
	   // std::cout << weight << " " << weight_lee << " " << flag_lee << " " << histoname << std::endl;
	 }
       }
     }


 }


 void LEEana::CovMatrix::get_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > > &map_all_events, std::map<TString, double>& map_filename_pot,  std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){

   //   std::cout << input_filename << std::endl;
  TFile *file = new TFile(input_filename);

  AnalysisTrees trees_cv = get_analysis_trees(file, "_cv");
  TTree *T_BDTvars_cv = trees_cv.T_BDTvars, *T_eval_cv = trees_cv.T_eval, *T_PFeval_cv = trees_cv.T_PFeval, *T_KINEvars_cv = trees_cv.T_KINEvars;
  TTree *T_spacepoints_cv = trees_cv.T_spacepoints, *T_pandora_cv = trees_cv.T_pandora, *T_lantern_cv = trees_cv.T_lantern, *T_glee_cv = trees_cv.T_glee;
  TTree *T_pot_cv = (TTree*)file->Get("wcpselection/T_pot_cv");

  EvalInfo eval_cv;
  POTInfo pot_cv;
  TaggerInfo tagger_cv;
  PFevalInfo pfeval_cv;
  KineInfo kine_cv;
  SpaceInfo space_cv;
  PandoraInfo pandora_cv;
  LanternInfo lantern_cv;
  GleeInfo glee_cv;

  AnalysisTrees trees_det = get_analysis_trees(file, "_det");
  TTree *T_BDTvars_det = trees_det.T_BDTvars, *T_eval_det = trees_det.T_eval, *T_PFeval_det = trees_det.T_PFeval, *T_KINEvars_det = trees_det.T_KINEvars;
  TTree *T_spacepoints_det = trees_det.T_spacepoints, *T_pandora_det = trees_det.T_pandora, *T_lantern_det = trees_det.T_lantern, *T_glee_det = trees_det.T_glee;
  TTree *T_pot_det = (TTree*)file->Get("wcpselection/T_pot_det");

  EvalInfo eval_det;
  POTInfo pot_det;
  TaggerInfo tagger_det;
  PFevalInfo pfeval_det;
  KineInfo kine_det;
  SpaceInfo space_det;
  PandoraInfo pandora_det;
  LanternInfo lantern_det;
  GleeInfo glee_det;

#include "init.txt"

  set_tree_address(T_BDTvars_cv, tagger_cv,2 );
  set_tree_address(T_eval_cv, eval_cv);
  set_tree_address(T_PFeval_cv, pfeval_cv);
  set_tree_address(T_pot_cv, pot_cv);
  set_tree_address(T_KINEvars_cv, kine_cv);
  if(T_spacepoints_cv) set_tree_address(T_spacepoints_cv, space_cv, 0);
  if(T_pandora_cv) set_tree_address(T_pandora_cv, pandora_cv);
  if(T_lantern_cv) set_tree_address(T_lantern_cv, lantern_cv);
  if(T_glee_cv) set_tree_address(T_glee_cv, glee_cv);

  set_tree_address(T_BDTvars_det, tagger_det,2 );
  set_tree_address(T_eval_det, eval_det);
  set_tree_address(T_PFeval_det, pfeval_det);
  set_tree_address(T_pot_det, pot_det);
  set_tree_address(T_KINEvars_det, kine_det);
  if(T_spacepoints_det) set_tree_address(T_spacepoints_det, space_det, 0);
  if(T_pandora_det) set_tree_address(T_pandora_det, pandora_det);
  if(T_lantern_det) set_tree_address(T_lantern_det, lantern_det);
  if(T_glee_det) set_tree_address(T_glee_det, glee_det);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot_cv->GetEntries();i++){
    T_pot_cv->GetEntry(i);
    total_pot += pot_cv.pot_tor875;
  }
  // total POT calculations ...
  map_filename_pot[input_filename] = total_pot;

  // fill histogram ...
  // switch on the branches the analysis code reads (one list for all apps, see analysis_trees.h)
  set_analysis_branch_status(trees_cv, false);
  set_analysis_branch_status(trees_det, false);

  //

  std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > vec_events;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  vec_events.resize(T_eval_cv->GetEntries());

  bool flag_first_event = true;
  for (Int_t i=0;i!=T_eval_cv->GetEntries();i++){
    T_BDTvars_cv->GetEntry(i);
    T_eval_cv->GetEntry(i);
    T_KINEvars_cv->GetEntry(i);
    T_PFeval_cv->GetEntry(i);
    if(T_spacepoints_cv) T_spacepoints_cv->GetEntry(i);
    if(T_pandora_cv) T_pandora_cv->GetEntry(i);
    if(T_lantern_cv) T_lantern_cv->GetEntry(i);
    if(T_glee_cv) T_glee_cv->GetEntry(i);

    T_BDTvars_det->GetEntry(i);
    T_eval_det->GetEntry(i);
    T_KINEvars_det->GetEntry(i);
    T_PFeval_det->GetEntry(i);
    if(T_spacepoints_det) T_spacepoints_det->GetEntry(i);
    if(T_pandora_det) T_pandora_det->GetEntry(i);
    if(T_lantern_det) T_lantern_det->GetEntry(i);
    if(T_glee_det) T_glee_det->GetEntry(i);

    if ( !(eval_cv.run == eval_det.run && eval_cv.event == eval_det.event)) std::cout <<"Wrong! " << std::endl;

    std::get<0>(vec_events.at(i)) = eval_cv.run;
    std::get<1>(vec_events.at(i)) = eval_cv.event;
    std::get<2>(vec_events.at(i)) = eval_cv.weight_cv * eval_cv.weight_spline;
    if (flag_rootino) std::get<2>(vec_events.at(i)) *= get_rootino_weight(eval_cv, glee_cv, rootino_pot_ratio); // rootino bug fix, see get_rootino_weight



    //    std::get<3>(vec_events.at(i)) = eval_cv.weight_lee;
    // hack for now ...
    std::get<3>(vec_events.at(i)) = leeweight(eval_cv.truth_nuEnergy);

    double osc_weight = 1.0;
    bool flag_updated = false;

    // event-level quantities for get_cut_pass, computed once per event
    CutEventInfo cut_info_cv;
    CutEventInfo cut_info_det;
    fill_cut_event_info(cut_info_cv, false, eval_cv, pfeval_cv, tagger_cv, kine_cv, space_cv, pandora_cv, lantern_cv);
    fill_cut_event_info(cut_info_det, false, eval_det, pfeval_det, tagger_det, kine_det, space_det, pandora_det, lantern_det);

    for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);

      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);

      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      //      TString weight = std::get<7>(*it);

      auto it3 = disabled_ch_names.find(ch_name);
      if (it3 != disabled_ch_names.end()) continue;

      bool flag_pass = get_cut_pass(ch_name, add_cut, false, cut_info_cv, eval_cv, pfeval_cv, tagger_cv, kine_cv, space_cv, pandora_cv, lantern_cv);
      bool flag_pass1 = get_cut_pass(ch_name, add_cut, false, cut_info_det, eval_det, pfeval_det, tagger_det, kine_det, space_det, pandora_det, lantern_det);

      // both variables are stored if either sample passes (always on the first event so an unknown variable name still stops the job)
      double val = 0, val1 = 0;
      if (flag_pass || flag_pass1 || flag_first_event){
        val = get_kine_var(kine_cv, eval_cv, pfeval_cv, tagger_cv, false, var_name, space_cv, pandora_cv, lantern_cv);
        val1 = get_kine_var(kine_det, eval_det, pfeval_det, tagger_det, false, var_name, space_det, pandora_det, lantern_det);
      }
      if (flag_pass || flag_pass1) {
	std::get<4>(vec_events.at(i) ).insert(std::make_tuple(no, val, flag_pass, val1, flag_pass1));
      }

      if (flag_osc && is_osc_channel(ch_name) && (!flag_updated)){
	osc_weight = get_osc_weight(eval_cv, pfeval_cv);
	flag_updated = true;
      }

    }
    flag_first_event = false;
    std::get<2>(vec_events.at(i)) *= osc_weight;
    double reweight = get_weight("add_weight", eval_cv, pfeval_cv, kine_cv, tagger_cv, glee_cv, get_rw_info(), get_time_info_allruns());//automatically 1 if reweighting is not applied
    std::get<2>(vec_events.at(i)) *= reweight;



    //std::cout << std::get<0>(vec_events.at(i)) << " " << std::get<1>(vec_events.at(i)) << " " << std::get<2>(vec_events.at(i)) << " " << std::get<3>(vec_events.at(i))  << " " << std::get<4>(vec_events.at(i) ).size() << std::endl;
  }



  map_all_events[input_filename] = vec_events;

  delete file;
}


std::pair<double, double> LEEana::CovMatrix::get_bayes_errors(double num){
  if (num > g_hlimit) return std::make_pair(sqrt(num), sqrt(num));
  else if (num < g_llimit) num = g_llimit;

  //  std::cout << num << " " << gl->Eval(num) << " " << gh->Eval(num) << std::endl;

  return std::make_pair(gl->Eval(num),gh->Eval(num));
}


std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > LEEana::CovMatrix::get_histograms(TString filename, int flag){
  if (flag == 0){
    auto it = map_inputfile_histograms.find(filename);
    if (it != map_inputfile_histograms.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }else if (flag==1){
    auto it = map_inputfile_histograms_err2.find(filename);
    if (it != map_inputfile_histograms_err2.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }else if (flag==2){
    auto it = map_inputfile_histograms_cros.find(filename);
    if (it != map_inputfile_histograms_cros.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }
}

void LEEana::CovMatrix::fill_data_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram){
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (period == run || run == 0){
      if (filetype == 5 || filetype == 15){
	std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);
	for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	  int obsch = get_obsch_name(std::get<5>(*it1));
	  TH1F *htemp = map_name_histogram[std::get<0>(*it1)].first;
	  map_obsch_histos[obsch].at(0)->Add(htemp);
	}
      }
    }
  }

}


void LEEana::CovMatrix::fill_pred_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > >& map_obsch_bayes, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > >& map_obsch_infos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram, float lee_strength, std::map<int, double> map_data_period_pot, bool flag_breakdown, std::map<int, std::vector<TH1F*> >& map_obsch_subhistos){

  for (auto it = map_pred_obsch_histos.begin(); it!=map_pred_obsch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int obsch = it->first;
    TH1F *hpred = map_obsch_histos[obsch].at(1);
    TH1F *hpred_err2 = map_obsch_histos[obsch].at(2);

    std::map<TString, std::pair<double, double> > temp_map_histo_ratios;

    std::map<TString, std::vector< std::tuple<double, double, double, int, double> > > map_histoname_values;

    // group
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      //std::cout << "sub: " << (*it1).size() << std::endl;

      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      TH1F *htemp_err2 = (TH1F*)hpred_err2->Clone("htemp_err2");
      htemp_err2->Reset();

      std::map<int, double> temp_map_data_pot;
      std::map<int, double> temp_map_mc_acc_pot;

      // get accumulated POT ...
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;

	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run && run !=0) continue; // skip ...
	auto it_data = map_data_period_pot.find(period); if (it_data == map_data_period_pot.end()) continue; // no corresponding data ...

	double data_pot = it_data->second;
	temp_map_data_pot[period] = data_pot;

	std::pair<TString, int> err2_lee = map_pred_histo_histo_err2_lee[histoname];
	TString histoname_err2 = err2_lee.first;
	int flag_lee = err2_lee.second;
	if (flag_lee == 1 && lee_strength == 0) continue; // no need to add ...

	TH1F *hmc = map_name_histogram[histoname].first;
	TH1F *hmc_err2 = map_name_histogram[histoname_err2].first;
	double mc_pot = map_name_histogram[histoname].second;
	if (temp_map_mc_acc_pot.find(period) == temp_map_mc_acc_pot.end()){
	  temp_map_mc_acc_pot[period] = mc_pot;
	}else{
	  temp_map_mc_acc_pot[period] += mc_pot;
	}
      }

      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;
	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run && run !=0) continue; // skip ...
	auto it_data = map_data_period_pot.find(period); if (it_data == map_data_period_pot.end()) continue; // no corresponding data ...
	std::pair<TString, int> err2_lee = map_pred_histo_histo_err2_lee[histoname];
	TString histoname_err2 = err2_lee.first;
	int flag_lee = err2_lee.second;
	if (flag_lee == 1 && lee_strength == 0) continue; // no need to add ...
	TH1F *hmc = map_name_histogram[histoname].first;
	TH1F *hmc_err2 = map_name_histogram[histoname_err2].first;

	double ratio = temp_map_data_pot[period]/temp_map_mc_acc_pot[period];

	temp_map_histo_ratios[histoname] = std::make_pair(ratio, lee_strength);

	if (flag_lee == 1) ratio *= lee_strength;

    /* std::string hhhhhname = hmc->GetName(); */
    /* std::istringstream sss(hhhhhname); */
    /* for(std::string line; std::getline(sss, line, '_');){ */
    /*     if(line == "CCQE") */
    /*     { */
    /*         ratio*=1.64; */
    /*     } */
    /* } */

	htemp->Add(hmc, ratio);
	htemp_err2->Add(hmc_err2, ratio*ratio);
    if(flag_breakdown){
        TH1F *hbreakdown = (TH1F*)hmc->Clone(Form("breakdown_%s", hmc->GetName()));
        hbreakdown->Scale(ratio);
        map_obsch_subhistos[obsch].push_back(hbreakdown);
        //std::cout<<"DEBUG: "<<obsch<<" "<<hmc->GetName()<<"\n";
    }

	// mean, sigma2, pot_ratio, covch, add_sys2
	std::vector< std::tuple<double, double, double, int, double> > values;
	for (int i=0;i!=hmc->GetNbinsX()+1;i++){
	  values.push_back(std::make_tuple(hmc->GetBinContent(i+1)*ratio , hmc_err2->GetBinContent(i+1)*ratio*ratio, temp_map_data_pot[period]/temp_map_mc_acc_pot[period], map_histogram_covch_add[histoname].first, pow(hmc->GetBinContent(i+1) *ratio * map_histogram_covch_add[histoname].second,2) ));
	}
	map_histoname_values[histoname] = values;

	//if (obsch==1){ //std::cout << values.at(2).first << " sep " << values.at(2).second << " " << histoname << std::endl;
	//  std::cout << histoname << " " << map_histogram_covch_add[histoname].first << " " << period << " " << temp_map_data_pot[period] << " " << temp_map_mc_acc_pot[period] << " " << hmc->GetSum() << " " << ratio << std::endl;
	// for (int i=0;i!=hmc->GetNbinsX();i++){
	//    std::cout << i << " " << hmc->GetBinContent(i+1) << std::endl;
	//  }
	//}
      }




      hpred->Add(htemp);
      hpred_err2->Add(htemp_err2);
      // delete the histograms ...
      delete htemp;
      delete htemp_err2;
    }

    for (auto it2 = map_histoname_values.begin(); it2 != map_histoname_values.end(); it2++){
      //if (obsch==1) std::cout << it2->first << " " << std::endl;
      map_obsch_infos[obsch].push_back(it2->second);
    }


    //    std::cout << map_histoname_values.size() << std::endl;
    // treat cross term ???

    for (auto it2 = map_pair_histo_histos_cros.begin(); it2 != map_pair_histo_histos_cros.end(); it2++){
      TString hist1 = it2->first.first;
      TString hist2 = it2->first.second;
      TString histogram = it2->second.first;
      int obsch1 = it2->second.second;
      if (obsch1 != obsch) continue;
      auto it3 = temp_map_histo_ratios.find(hist1);
      auto it4 = temp_map_histo_ratios.find(hist2);
      if (it3 == temp_map_histo_ratios.end() || it4 == temp_map_histo_ratios.end()) continue;
      TH1F *hmc_cros = map_name_histogram[histogram].first;
      // cross term ...
      hpred_err2->Add(hmc_cros, it3->second.first * it4->second.first * 2 * it3->second.second);

      std::vector< std::tuple<double, double, double, int, double> >& values1 = map_histoname_values[hist1];
      std::vector< std::tuple<double, double, double, int, double> >& values2 = map_histoname_values[hist2];
      for (size_t j=0;j!=values1.size();j++){
	std::get<0>(values1.at(j)) += std::get<0>(values2.at(j));
	std::get<1>(values1.at(j)) += std::get<1>(values2.at(j)) + hmc_cros->GetBinContent(j+1) * it3->second.first * it4->second.first * 2 * it3->second.second;
      }
      map_histoname_values.erase(hist2);
      //      std::cout << values1.size() << " " << values2.size() << std::endl;
      //      std::cout << hist1 << " " << hist2 << " " << histogram << " " << it3->second.first << " " << it3->second.second << " " << it4->second.first << " " << it4->second.second << std::endl;
    }
    //   std::cout <<  map_histoname_values.size() << std::endl;

    for (auto it2 = map_histoname_values.begin(); it2 != map_histoname_values.end(); it2++){
      //if (obsch==1) std::cout << it2->first << " " << std::endl;
      map_obsch_bayes[obsch].push_back(it2->second);
    }

  } // work on obs channel ...
}



int LEEana::CovMatrix::get_obsch_name(TString name){
  int ch = map_name_ch[name];
  // std::cout << ch << std::endl;
  return std::get<6>(map_ch_hist[ch]);
}

int LEEana::CovMatrix::get_covch_name(TString name){
  int ch = map_name_ch[name];
  return map_ch_covch[ch];
}

float LEEana::CovMatrix::get_ext_pot(TString filename){
  auto it = map_inputfile_info.find(filename);
  if (it != map_inputfile_info.end()){
    return std::get<3>(it->second);
  }else{
    return 0;
  }
}

void LEEana::CovMatrix::print_cvfile_info(){

  for (auto it = map_filetype_inputfiles.begin(); it!= map_filetype_inputfiles.end(); it++){
    std::cout << it->first << " " << map_filetype_name[it->first] << std::endl;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout << *it1 << " " << std::get<0>(map_inputfile_info[*it1]) << " " << std::get<1>(map_inputfile_info[*it1]) << " " << std::get<2>(map_inputfile_info[*it1]) << " " << std::get<3>(map_inputfile_info[*it1]) << " " << map_inputfile_cuts[*it1].size() << std::endl;
    }
  }
}


int LEEana::CovMatrix::get_ch(TString name){
  auto it = map_name_ch.find(name);
  if (it != map_name_ch.end()){
    return it->second;
  }else{
    return -1;
  }
}


int LEEana::CovMatrix::get_obsch(int ch){
  auto it = map_ch_obsch.find(ch);
  if (it != map_ch_obsch.end()){
    return it->second;
  }else{
    return -1;
  }
}

int LEEana::CovMatrix::get_covch(int ch){
  auto it = map_ch_covch.find(ch);
  if (it != map_ch_covch.end()){
    return it->second;
  }else{
    return -1;
  }
}

int LEEana::CovMatrix::get_obsch_fcov(int covch){
  auto it = map_covch_obsch.find(covch);
  if (it != map_covch_obsch.end()){
    return it->second;
  }else{
    return -1;
  }
}

std::pair<int, int> LEEana::CovMatrix::get_obsch_info(int obsch){
  auto it = map_obsch_startbin.find(obsch);
  if (it != map_obsch_startbin.end()){
    return std::make_pair(it->second, map_obsch_nbin[it->first]);
  }else{
    return std::make_pair(0,0);
  }
}
std::pair<int, int> LEEana::CovMatrix::get_covch_info(int covch){
  auto it = map_covch_startbin.find(covch);
  if (it != map_covch_startbin.end()){
    return std::make_pair(it->second, map_covch_nbin[it->first]);
  }else{
    return std::make_pair(0,0);
  }
}


bool LEEana::CovMatrix::get_sys_xs_flux(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<0>(it->second)==0){
      return false;
    }else{
      return true;
    }
  }else{
    return false;
  }
}

bool LEEana::CovMatrix::get_sys_det(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<1>(it->second)==0){
      return false;
    }else{
      return true;
    }
  }else{
    return false;
  }
}
std::pair<bool, float> LEEana::CovMatrix::get_sys_add(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<2>(it->second)==0){
      return std::make_pair(false,0);
    }else{
      return std::make_pair(true,std::get<2>(it->second) );
    }
  }else{
    return std::make_pair(false,0);
  }
}

int LEEana::CovMatrix::get_sys_mc_same(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    return std::get<3>(it->second);
  }else{
    return 0;
  }
}

std::vector<int> LEEana::CovMatrix::get_filetype_chs(int filetype){
  auto it = map_filetype_chs.find(filetype);
  if (it != map_filetype_chs.end()){
    return it->second;
  }else{
    std::vector<int> tmp;
    return tmp;
  }
}

int LEEana::CovMatrix::get_ch_filetype(int ch){
  auto it = map_ch_filetype.find(ch);

  if (it != map_ch_filetype.end()){
    return it->second;
  }else{
    return -1;
  }
}

std::tuple<int, float, float> LEEana::CovMatrix::get_ch_hist(int ch){
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::make_tuple(std::get<2>(it->second), std::get<3>(it->second), std::get<4>(it->second));
  }else{
    return std::make_tuple(-1,0,0);
  }
}


TString LEEana::CovMatrix::get_ch_name(int ch){
  TString result = "Null";
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::get<0>(it->second);
  }
  return result;
}

TString LEEana::CovMatrix::get_ch_var(int ch){
  TString result = "Null";
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::get<1>(it->second);
  }
  return result;
}




void LEEana::CovMatrix::print_ch_info(){
  for (auto it = map_ch_hist.begin(); it!= map_ch_hist.end(); it++){
    std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << " " << std::get<2>(it->second) << " " << std::get<3>(it->second) << " " << std::get<4>(it->second) << std::endl;
  }
}


void LEEana::CovMatrix::print_filetype_info(){
  for (auto it = map_filetype_chs.begin(); it != map_filetype_chs.end(); it++){
    std::cout << it->first << ": ";
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout << get_ch_name(*it1) << " " ;
    }
    std::cout << std::endl;
  }
}

void LEEana::CovMatrix::print_matrix(){

  for (auto it = map_ch_obsch.begin(); it != map_ch_obsch.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_ch_covch[it->first] << std::endl;
  }

  std::cout << "covch: " << std::endl;
  for (auto it = map_covch_nbin.begin(); it != map_covch_nbin.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_covch_startbin[it->first] << std::endl;
    //map_covch_nbin[it->first] = start_bin;
  }

  std::cout << "obsch: " << std::endl;
  for (auto it = map_obsch_nbin.begin(); it != map_obsch_nbin.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_obsch_startbin[it->first] << std::endl;
    //map_covch_nbin[it->first] = start_bin;
  }

  std::cout << "matrix: " << std::endl;
  for (auto it = map_covchbin_obschbin.begin(); it != map_covchbin_obschbin.end(); it++){
    std::cout << it->first << " " << it->second << std::endl;
  }
}

void LEEana::CovMatrix::print_systematics(){
  for (auto it = map_ch_systematics.begin(); it!= map_ch_systematics.end(); it++){
    std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << " " << std::get<2>(it->second) << " " << std::get<3>(it->second) << std::endl;
  }

  std::cout << "Xs+Flux: ";
  for (auto it = xfs_filetypes.begin(); it != xfs_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;

  std::cout << "Det    : ";
  for (auto it = det_filetypes.begin(); it != det_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;


  std::cout << "Add    : ";
  for (auto it = add_filetypes.begin(); it != add_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;

  std::cout << "same MC stat: " << std::endl;
  for (auto it = map_mcstat_same_chs.begin(); it != map_mcstat_same_chs.end(); it++){
    std::cout << it->first << ": ";
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout <<get_ch_name(*it1) << " ";
    }
    std::cout << std::endl;
  }

}
