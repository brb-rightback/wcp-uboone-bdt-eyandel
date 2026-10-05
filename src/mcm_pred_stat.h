void LEEana::CovMatrix::gen_pred_stat_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TMatrixD* cov_mat, int seed){

  gRandom->SetSeed(seed);

  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ;
  std::map<int, TString> map_no_histoname;


  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);

    if (filetype == 5 || filetype == 15)   continue; // if this is data ...


    if (period != run && run!=0) continue;
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
  std::map<TString, std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > > map_all_events;
  std::map<TString, double> map_filename_pot;
  std::map<int, double> map_data_period_pot;

  std::map<TString, std::set<std::pair<TString, double> > > map_filename_set_histoname_ratio;

  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);

    if (filetype == 5 || filetype == 15) {
      map_data_period_pot[period] = get_ext_pot(input_filename);
      //      std::cout << period << " " << map_data_period_pot[period] << std::endl;
      continue; // skip data ...
    }
    if (period != run && run !=0) continue;

    get_pred_events_info(input_filename, map_all_events, map_filename_pot, map_histoname_infos);
  }


  fill_pred_stat_histograms(map_all_events, vec_no_hist, vec_no_lee, map_histoname_hist);

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
    	int period = std::get<1>(it3->second);
    	if (period != run && run != 0) continue;

    	double mc_pot = map_filename_pot[input_filename];
    	//std::cout << mc_pot << std::endl;
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
    	int period = std::get<1>(it3->second);
    	if (period != run && run !=0) continue; // skip ...
    	double data_pot = map_data_period_pot[period];
    	double ratio = data_pot/temp_map_mc_acc_pot[period];

	map_filename_set_histoname_ratio[input_filename].insert(std::make_pair(histoname, ratio));

    	TH1F *hmc = map_histoname_hist[histoname];
    	htemp->Add(hmc, ratio);
      }

      hpred->Add(htemp);
      delete htemp;
    }

    int start_bin = map_covch_startbin[covch];
    for (int i=0;i!=hpred->GetNbinsX()+1;i++){
      (*vec_mean)[start_bin+i] = hpred->GetBinContent(i+1);

      // std::cout << start_bin+i << " " << (*vec_mean_diff)(start_bin+i) << " " <<  hpred->GetBinContent(i+1) << std::endl;;
      // std::cout << x[start_bin+i] << std::endl;
    }
  }


  const int rows = cov_mat->GetNcols();
  float x[rows];
  (*cov_mat).Zero();

  std::map<TString, TH1D*> map_filename_histo;
  for (auto it = map_all_events.begin(); it != map_all_events.end(); it++){
    TString filename = it->first;
    int nsize = it->second.size();
    TH1D* htemp = new TH1D(filename, filename, nsize, 0.5, nsize+0.5);
    for (size_t i=0;i!=nsize;i++){
      htemp->SetBinContent(i+1, std::get<2>(it->second.at(i)) );
    }
    map_filename_histo[filename] = htemp;
  }

  int nround = 5000;
  for (int qx = 0; qx != nround; qx++){
    if (qx % 500 ==0) std::cout << qx << std::endl;

    for (int i=0;i!=rows;i++){
      x[i] = 0;
    }

    // fill the histogram with variation ...
    fill_pred_stat_histograms(map_filename_histo, map_all_events, vec_no_hist, vec_no_lee, map_no_histoname, map_histoname_hist, map_filename_set_histoname_ratio);

     // merge histograms according to POTs ...
    for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
      //std::cout << it->first << std::endl;
      int covch = it->first;
      TH1F *hpred = map_covch_hist[covch];
      hpred->Reset();

      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	TH1F *htemp = (TH1F*)hpred->Clone("htemp");
	htemp->Reset();


	for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	  TString histoname = (*it2).first;
	  TString input_filename = map_histogram_inputfile[histoname];
	  auto it3 = map_inputfile_info.find(input_filename);
	  int period = std::get<1>(it3->second);
	  if (period != run && run !=0) continue; // skip ...

	  TH1F *hmc = map_histoname_hist[histoname];
	  htemp->Add(hmc, 1);
	}

	hpred->Add(htemp);
	delete htemp;
      }

      int start_bin = map_covch_startbin[covch];
      for (int i=0;i!=hpred->GetNbinsX()+1;i++){
	x[start_bin+i]= hpred->GetBinContent(i+1);
	//	std::cout << (*vec_mean)(start_bin+i) << " " << x[start_bin+i] << std::endl;
      }
    }


    // add covariance matrix ...
    for (size_t n = 0;n!=rows; n++){
      for (size_t m =0; m!=rows;m++){
	(*cov_mat)(n,m) += (x[n]-(*vec_mean)[n]) * (x[m] - (*vec_mean)[m]) * 1./nround;
      }
    }
  }


  // clean up the memory ...
  for (auto it = map_filename_histo.begin(); it != map_filename_histo.end(); it++){
    delete it->second;
  }



  // update the final  ones ...
   fill_pred_stat_histograms(map_all_events, vec_no_hist, vec_no_lee, map_histoname_hist);

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
    	int period = std::get<1>(it3->second);
    	if (period != run && run != 0) continue;

    	double mc_pot = map_filename_pot[input_filename];
    	//std::cout << mc_pot << std::endl;
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
    	int period = std::get<1>(it3->second);
    	if (period != run && run !=0) continue; // skip ...
    	double data_pot = map_data_period_pot[period];
    	double ratio = data_pot/temp_map_mc_acc_pot[period];

    	TH1F *hmc = map_histoname_hist[histoname];
    	htemp->Add(hmc, ratio);
      }

      hpred->Add(htemp);
      delete htemp;
    }

    //    int start_bin = map_covch_startbin[covch];
    // for (int i=0;i!=hpred->GetNbinsX()+1;i++){
    //  (*vec_mean)[start_bin+i] = hpred->GetBinContent(i+1);

      // std::cout << start_bin+i << " " << (*vec_mean_diff)(start_bin+i) << " " <<  hpred->GetBinContent(i+1) << std::endl;;
      // std::cout << x[start_bin+i] << std::endl;
    //}
  }

}

void LEEana::CovMatrix::fill_pred_stat_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::vector<TH1F*>& vec_no_hist, std::vector<int>& vec_no_lee, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist, std::map<TString, std::set<std::pair<TString, double> > >& map_filename_histoname_ratio){

  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
    it->second->Reset();
  }

  // loop over files
  for (auto it = map_all_events.begin(); it!=map_all_events.end(); it++){
    TString filename = it->first;
    TH1D *hweight = map_filename_histo[filename];

    auto it2 = map_filename_histoname_ratio.find(filename);
    if (it2 == map_filename_histoname_ratio.end()) continue;

    double max_sum = 0;
    std::map<TString, double> map_histoname_sum;

    std::map<double, double> map_ratio_sum;

    for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++){
      double ratio = it3->second;
      TString temp_histoname = it3->first;

      auto it4 = map_ratio_sum.find(ratio);
      if (it4 == map_ratio_sum.end()){
	double sum = gRandom->Poisson(hweight->GetSum() * ratio);
	if (sum > max_sum) max_sum = sum;
	map_histoname_sum[temp_histoname] = sum;
	map_ratio_sum[ratio] = sum;
      }else{
	double sum = it4->second;
	map_histoname_sum[temp_histoname] = sum;
      }
    }


    for (size_t i=0;i<max_sum;i++){
      int global_index = hweight->FindBin(hweight->GetRandom())-1;
      double weight_lee = std::get<3>(it->second.at(global_index));

       for (auto it1 = std::get<4>(it->second.at(global_index)).begin(); it1 != std::get<4>(it->second.at(global_index)).end(); it1++){
	 int no = std::get<0>(*it1);
	 double val_cv = std::get<1>(*it1);
	 bool flag_cv = std::get<2>(*it1);

	 const TString& histoname = map_no_histoname[no];
	 TH1F *htemp = vec_no_hist[no];
	 int flag_lee = vec_no_lee[no];

	 auto it4 = map_histoname_sum.find(histoname);
	 if (it4 == map_histoname_sum.end()) std::cout << "something wrong! " << std::endl;
	 double sum = it4->second;

	 if (i>=sum) continue; // go over ...

	  if (flag_cv){
	    if (flag_lee){
	      htemp->Fill(val_cv,weight_lee);
	    }else{
	      htemp->Fill(val_cv,1);
	    }
	  }
       }

    }
  }

}


void LEEana::CovMatrix::fill_pred_stat_histograms(std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::vector<TH1F*>& vec_no_hist, std::vector<int>& vec_no_lee,  std::map<TString, TH1F*>& map_histoname_hist){
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
      for (auto it1 = std::get<4>(it->second.at(i)).begin(); it1 != std::get<4>(it->second.at(i)).end(); it1++){
	int no = std::get<0>(*it1);
	double val_cv = std::get<1>(*it1);
	bool flag_cv = std::get<2>(*it1);

	TH1F *htemp = vec_no_hist[no];
	 int flag_lee = vec_no_lee[no];
	if (flag_cv){
	  if (flag_lee){
	    htemp->Fill(val_cv, weight * weight_lee);
	  }else{
	    htemp->Fill(val_cv, weight);
	  }
	}
      }
    }
  }
}

void LEEana::CovMatrix::get_pred_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){
  TFile *file = new TFile(input_filename);

  AnalysisTrees trees = get_analysis_trees(file);
  TTree *T_BDTvars = trees.T_BDTvars, *T_eval = trees.T_eval, *T_PFeval = trees.T_PFeval, *T_KINEvars = trees.T_KINEvars;
  TTree *T_spacepoints = trees.T_spacepoints, *T_pandora = trees.T_pandora, *T_lantern = trees.T_lantern, *T_glee = trees.T_glee;
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;
  SpaceInfo space;
  PandoraInfo pandora;
  LanternInfo lantern;
  GleeInfo glee;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;

  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  bool flag_data = true;
  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  set_tree_address(T_BDTvars, tagger, 2);
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  if(T_spacepoints) set_tree_address(T_spacepoints, space, 0);
  if(T_pandora) set_tree_address(T_pandora, pandora);
  if(T_lantern) set_tree_address(T_lantern, lantern);
  if(T_glee) set_tree_address(T_glee, glee);

  //  std::cout << flag_data << " " << input_filename << std::endl;

  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;

  map_filename_pot[input_filename] = total_pot;
  //std::cout << input_filename << " " << total_pot << " " << std::endl;


  // fill histogram ...
  // switch on the branches the analysis code reads (one list for all apps, see analysis_trees.h)
  set_analysis_branch_status(trees, flag_data);


  std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > vec_events;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  vec_events.resize(T_eval->GetEntries());

  bool flag_first_event = true;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    if(T_spacepoints) T_spacepoints->GetEntry(i);
    if(T_pandora) T_pandora->GetEntry(i);
    if(T_lantern) T_lantern->GetEntry(i);
    if(T_glee) T_glee->GetEntry(i);

    std::get<0>(vec_events.at(i)) = eval.run;
    std::get<1>(vec_events.at(i)) = eval.event;

    if (!flag_data){
      std::get<2>(vec_events.at(i)) = eval.weight_cv * eval.weight_spline;
      if (flag_rootino) std::get<2>(vec_events.at(i)) *= get_rootino_weight(eval, glee, rootino_pot_ratio); // rootino bug fix, see get_rootino_weight
      // hack for now ...
      std::get<3>(vec_events.at(i)) = leeweight(eval.truth_nuEnergy);
    }else{
      std::get<2>(vec_events.at(i)) = 1;
      std::get<3>(vec_events.at(i)) = 0;
    }

    double osc_weight = 1.0;
    bool flag_updated = false;

    // event-level quantities for get_cut_pass, computed once per event
    CutEventInfo cut_info;
    fill_cut_event_info(cut_info, flag_data, eval, pfeval, tagger, kine, space, pandora, lantern);

    for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);

      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);

      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);

      bool flag_pass = get_cut_pass(ch_name, add_cut, flag_data, cut_info, eval, pfeval, tagger, kine, space, pandora, lantern);
      // the variable is only stored if the event passes (always on the first event so an unknown variable name still stops the job)
      double val = 0;
      if (flag_pass || flag_first_event) val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name, space, pandora, lantern);

      if (flag_pass) std::get<4>(vec_events.at(i)).insert(std::make_tuple(no, val, flag_pass));

      if (flag_osc && is_osc_channel(ch_name) && (!flag_updated)){
	osc_weight = get_osc_weight(eval, pfeval);
	flag_updated = true;
      }
    }
    flag_first_event = false;
    if (!flag_data){
      std::get<2>(vec_events.at(i)) *= osc_weight;
      std::get<2>(vec_events.at(i)) *= get_weight("add_weight", eval, pfeval, kine, tagger, glee, get_rw_info(), get_time_info_allruns());
    }
}

  map_all_events[input_filename] = vec_events;

  delete file;

}
