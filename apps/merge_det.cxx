// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <tuple>
#include <algorithm>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMemFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "WCPLEEANA/eval.h"

#include "WCPLEEANA/tree_wrangler.h"

using namespace std;
using namespace LEEana;

//#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"
#include "WCPLEEANA/cmd_options.h"

int main( int argc, char** argv )
{
  if (argc < 4) {
    std::cout << "merge_det #input_file_cv #input_file_det #output_file " << std::endl;
    std::cout << "  -m2 -n#det_no: DetVar file without the CV events, with the flags of the CV events and subruns used (T_cv_rse_flag, T_cv_pot_flag)" << std::endl;
    std::cout << "  -m3: merge_det #input_file_cv #det_list #output_file -m3, CV file for all -m2 DetVar files, #det_list has lines \"#det_no #det_file\" ending with \"end end\"" << std::endl;
    return -1;
  }

  TString input_file_cv = argv[1];
  TString input_file_det = argv[2];
  TString out_file = argv[3];
  bool flag_config = false;
  std::string config_file_name="config.txt";
  char delimiter = ',';
  bool override_pot = false;
  int flag_mode = 1; // 1: CV and DetVar events in one file, 2: DetVar events and CV flags only, 3: CV file with the flags of the -m2 files
  int det_no = -1;
  argc = LEEana::join_option_values(argc, argv, 4); // also take "-x VALUE" for "-xVALUE"
  for (Int_t i=4;i!=argc;i++){
    switch(argv[i][1]){
    case 't':
       config_file_name = &argv[i][2];
       flag_config = true;
      break;
    case 'd':
        delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
      break;
    case 'p':
        override_pot = &argv[i][2];
        if(override_pot) std::cout<<"Overriding the pot of each subrun"<<std::endl;
      break;
    case 'm':
      flag_mode = atoi(&argv[i][2]);
      break;
    case 'n':
      det_no = atoi(&argv[i][2]);
      break;
    }
  }

  tree_wrangler wrangler_cv(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_det(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_pot_cv(flag_config, config_file_name, delimiter,true);
  tree_wrangler wrangler_pot_det(flag_config, config_file_name, delimiter,true);

  if (flag_mode == 2 && det_no < 0){
    std::cout << "ERROR! -m2 needs the DetVar number, -n#det_no" << std::endl;
    return 1;
  }

  // -m3: copy of the CV file (trees renamed _cv as in -m1), with the flags of the events and subruns each -m2 DetVar file uses:
  // T_eval_cv flag_match_rse_XX, T_pot_cv flag_match_rse_XX and pot_tor875_rse_XX (the common POT of the subrun, as T_pot_cv of -m1)
  if (flag_mode == 3){
    // list of the DetVar files, "#det_no #det_file" per line, ending with "end end"
    std::vector<int> vec_det_no;
    std::vector<TString> vec_det_file;
    std::ifstream infile(input_file_det.Data());
    if (!infile.good()){
      std::cout << "ERROR! Can not open the DetVar list " << input_file_det << std::endl;
      return 1;
    }
    TString temp_no, temp_file;
    while (infile >> temp_no >> temp_file){
      if (temp_no == "end" || temp_no == "End") break;
      if (std::find(vec_det_no.begin(), vec_det_no.end(), temp_no.Atoi()) != vec_det_no.end()){
        std::cout << "ERROR! DetVar number " << temp_no << " is listed twice in " << input_file_det << std::endl;
        return 1;
      }
      vec_det_no.push_back(temp_no.Atoi());
      vec_det_file.push_back(temp_file);
    }
    int ndet = vec_det_no.size();

    TFile *file1 = new TFile(input_file_cv);
    TTree *T_BDTvars_cv = (TTree*)file1->Get("wcpselection/T_BDTvars");
    TTree *T_eval_cv = (TTree*)file1->Get("wcpselection/T_eval");
    TTree *T_pot_cv = (TTree*)file1->Get("wcpselection/T_pot");
    TTree *T_PFeval_cv = (TTree*)file1->Get("wcpselection/T_PFeval");
    TTree *T_KINEvars_cv = (TTree*)file1->Get("wcpselection/T_KINEvars");
    TTree *T_spacepoints_cv = (TTree*)file1->Get("wcpselection/T_spacepoints");

    wrangler_cv.get_old_trees(file1);
    wrangler_pot_cv.get_old_trees(file1);

    // run, subrun, event of the CV events and run, subrun of the CV subruns, to check the flags against
    Int_t run = 0, subrun = 0, event = 0, runNo = 0, subRunNo = 0;
    T_eval_cv->SetBranchAddress("run", &run);
    T_eval_cv->SetBranchAddress("subrun", &subrun);
    T_eval_cv->SetBranchAddress("event", &event);
    T_pot_cv->SetBranchAddress("runNo", &runNo);
    T_pot_cv->SetBranchAddress("subRunNo", &subRunNo);
    std::vector<std::tuple<int, int, int> > vec_rse_cv;
    for (Long64_t i=0;i!=T_eval_cv->GetEntries();i++){
      T_eval_cv->GetBranch("run")->GetEntry(i);
      T_eval_cv->GetBranch("subrun")->GetEntry(i);
      T_eval_cv->GetBranch("event")->GetEntry(i);
      vec_rse_cv.push_back(std::make_tuple(run, subrun, event));
    }
    std::vector<std::pair<int, int> > vec_rs_pot_cv;
    for (Long64_t i=0;i!=T_pot_cv->GetEntries();i++){
      T_pot_cv->GetBranch("runNo")->GetEntry(i);
      T_pot_cv->GetBranch("subRunNo")->GetEntry(i);
      vec_rs_pot_cv.push_back(std::make_pair(runNo, subRunNo));
    }

    // flags of each DetVar file
    std::vector<std::vector<int> > vec_flag_rse(ndet);
    std::vector<std::vector<int> > vec_flag_pot(ndet);
    std::vector<std::vector<double> > vec_pot_rse(ndet);
    for (int k=0;k!=ndet;k++){
      TFile *file_det = new TFile(vec_det_file.at(k));
      TTree *T_rse_flag = (TTree*)file_det->Get("wcpselection/T_cv_rse_flag");
      TTree *T_pot_flag = (TTree*)file_det->Get("wcpselection/T_cv_pot_flag");
      if (file_det->IsZombie() || !T_rse_flag || !T_pot_flag){
        std::cout << "ERROR! No T_cv_rse_flag / T_cv_pot_flag in " << vec_det_file.at(k) << ", make it with merge_det -m2" << std::endl;
        return 1;
      }
      if (T_rse_flag->GetEntries() != T_eval_cv->GetEntries() || T_pot_flag->GetEntries() != T_pot_cv->GetEntries()){
        std::cout << "ERROR! " << vec_det_file.at(k) << " was not made with this CV file: " << T_rse_flag->GetEntries() << " " << T_eval_cv->GetEntries() << " events, "
                  << T_pot_flag->GetEntries() << " " << T_pot_cv->GetEntries() << " subruns" << std::endl;
        return 1;
      }
      Int_t f_run, f_subrun, f_event, f_flag, f_det_no;
      Double_t f_pot;
      T_rse_flag->SetBranchAddress("run", &f_run);
      T_rse_flag->SetBranchAddress("subrun", &f_subrun);
      T_rse_flag->SetBranchAddress("event", &f_event);
      T_rse_flag->SetBranchAddress("flag", &f_flag);
      T_rse_flag->SetBranchAddress("det_no", &f_det_no);
      for (Long64_t i=0;i!=T_rse_flag->GetEntries();i++){
        T_rse_flag->GetEntry(i);
        if (f_det_no != vec_det_no.at(k)){
          std::cout << "ERROR! " << vec_det_file.at(k) << " was made with -n" << f_det_no << ", listed as " << vec_det_no.at(k) << std::endl;
          return 1;
        }
        if (std::make_tuple(f_run, f_subrun, f_event) != vec_rse_cv.at(i)){
          std::cout << "ERROR! " << vec_det_file.at(k) << " was not made with this CV file: event " << i << " is " << f_run << " " << f_subrun << " " << f_event << std::endl;
          return 1;
        }
        vec_flag_rse.at(k).push_back(f_flag);
      }
      T_pot_flag->SetBranchAddress("runNo", &f_run);
      T_pot_flag->SetBranchAddress("subRunNo", &f_subrun);
      T_pot_flag->SetBranchAddress("flag", &f_flag);
      T_pot_flag->SetBranchAddress("pot_tor875", &f_pot);
      for (Long64_t i=0;i!=T_pot_flag->GetEntries();i++){
        T_pot_flag->GetEntry(i);
        if (std::make_pair(f_run, f_subrun) != vec_rs_pot_cv.at(i)){
          std::cout << "ERROR! " << vec_det_file.at(k) << " was not made with this CV file: subrun " << i << " is " << f_run << " " << f_subrun << std::endl;
          return 1;
        }
        vec_flag_pot.at(k).push_back(f_flag);
        vec_pot_rse.at(k).push_back(f_pot);
      }
      delete file_det;
    }

    TFile *file3 = new TFile(out_file,"RECREATE");

    //Setup the directories specified in the config file
    wrangler_cv.set_new_trees(file3,true,"_cv");
    wrangler_pot_cv.set_new_trees(file3,true,"_cv");

    file3->mkdir("wcpselection");
    file3->cd("wcpselection");
    TTree *t1_cv = T_eval_cv->CloneTree(0);
    t1_cv->SetObject("T_eval_cv","T_eval_cv");
    TTree *t2_cv = T_pot_cv->CloneTree(0);
    t2_cv->SetObject("T_pot_cv","T_pot_cv");
    TTree *t3_cv = T_PFeval_cv->CloneTree(0);
    t3_cv->SetObject("T_PFeval_cv","T_PFeval_cv");
    TTree *t5_cv = T_KINEvars_cv->CloneTree(0);
    t5_cv->SetObject("T_KINEvars_cv","T_KINEvars_cv");
    TTree *t4_cv = T_BDTvars_cv->CloneTree(0);
    t4_cv->SetObject("T_BDTvars_cv","T_BDTvars_cv");
    TTree *T_spacepoints_new_cv = T_spacepoints_cv->CloneTree(0);
    T_spacepoints_new_cv->SetObject("T_spacepoints_cv","T_spacepoints_cv");

    std::vector<int> flag_rse(ndet, 0);
    std::vector<int> flag_pot(ndet, 0);
    std::vector<double> pot_rse(ndet, 0);
    for (int k=0;k!=ndet;k++){
      t1_cv->Branch(Form("flag_match_rse_%d", vec_det_no.at(k)), &flag_rse.at(k), Form("flag_match_rse_%d/I", vec_det_no.at(k)));
      t2_cv->Branch(Form("flag_match_rse_%d", vec_det_no.at(k)), &flag_pot.at(k), Form("flag_match_rse_%d/I", vec_det_no.at(k)));
      t2_cv->Branch(Form("pot_tor875_rse_%d", vec_det_no.at(k)), &pot_rse.at(k), Form("pot_tor875_rse_%d/D", vec_det_no.at(k)));
    }

    int nentries = T_eval_cv->GetEntries();
    std::cout<<"Begin looping over "<<nentries<<" events"<<std::endl;
    for (int i=0;i!=nentries;i++){
      if (i%10000 == 0) std::cout << i/1000 << " k " << std::setprecision(3) << double(i)/nentries*100. << " %"<< std::endl;

      T_BDTvars_cv->GetEntry(i);
      T_eval_cv->GetEntry(i);
      T_KINEvars_cv->GetEntry(i);
      T_PFeval_cv->GetEntry(i);
      T_spacepoints_cv->GetEntry(i);
      for(auto tree_it=wrangler_cv.old_trees->begin(); tree_it!=wrangler_cv.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(i);
      }
      for (int k=0;k!=ndet;k++){
        flag_rse.at(k) = vec_flag_rse.at(k).at(i);
      }

      t1_cv->Fill();
      t3_cv->Fill();
      t4_cv->Fill();
      t5_cv->Fill();
      T_spacepoints_new_cv->Fill();
      for(auto tree_it=wrangler_cv.new_trees->begin(); tree_it!=wrangler_cv.new_trees->end(); tree_it++){
          (*tree_it)->Fill();
      }
    }

    // POT trees are copied unchanged, the common POT of each DetVar file is in pot_tor875_rse_XX
    for (int i=0;i!=T_pot_cv->GetEntries();i++){
      T_pot_cv->GetEntry(i);
      for (int k=0;k!=ndet;k++){
        flag_pot.at(k) = vec_flag_pot.at(k).at(i);
        pot_rse.at(k) = vec_pot_rse.at(k).at(i);
      }
      t2_cv->Fill();
    }
    for (size_t j=0;j!=wrangler_pot_cv.old_trees->size();j++){
      for (Long64_t i=0;i!=wrangler_pot_cv.old_trees->at(j)->GetEntries();i++){
        wrangler_pot_cv.old_trees->at(j)->GetEntry(i);
        wrangler_pot_cv.new_trees->at(j)->Fill();
      }
    }

    std::cout << out_file << std::endl;
    std::cout << "Events: " << t1_cv->GetEntries() << ", subruns: " << t2_cv->GetEntries() << std::endl;
    for (int k=0;k!=ndet;k++){
      int nevents = 0, nsubruns = 0;
      double pot = 0;
      for (size_t i=0;i!=vec_flag_rse.at(k).size();i++) nevents += vec_flag_rse.at(k).at(i);
      for (size_t i=0;i!=vec_flag_pot.at(k).size();i++){
        nsubruns += vec_flag_pot.at(k).at(i);
        if (vec_flag_pot.at(k).at(i)) pot += vec_pot_rse.at(k).at(i);
      }
      std::cout << "DetVar " << vec_det_no.at(k) << " " << vec_det_file.at(k) << ": " << nevents << " events, " << nsubruns << " subruns, POT " << pot << std::endl;
    }

    file3->Write("",TFile::kOverwrite);
    file3->Close();

    return 0;
  }

  TFile *file1 = new TFile(input_file_cv);
  TTree *T_BDTvars_cv = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_eval_cv = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot_cv = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval_cv = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars_cv = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_spacepoints_cv = (TTree*)file1->Get("wcpselection/T_spacepoints");

  TFile *file2 = new TFile(input_file_det);
  TTree *T_BDTvars_det = (TTree*)file2->Get("wcpselection/T_BDTvars");
  TTree *T_eval_det = (TTree*)file2->Get("wcpselection/T_eval");
  TTree *T_pot_det = (TTree*)file2->Get("wcpselection/T_pot");
  TTree *T_PFeval_det = (TTree*)file2->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars_det = (TTree*)file2->Get("wcpselection/T_KINEvars");
  TTree *T_spacepoints_det = (TTree*)file2->Get("wcpselection/T_spacepoints");

  //Load other trees from directories as specified by the config file
  wrangler_cv.get_old_trees(file1);
  wrangler_det.get_old_trees(file2);
  wrangler_pot_cv.get_old_trees(file1);
  wrangler_pot_det.get_old_trees(file2);

  TFile *file3 = new TFile(out_file,"RECREATE");
  // -m2: the CV trees are not saved (merge_det -m3 makes one CV file for all DetVar files), they are only made in memory
  TFile *file_cv_out = file3;
  if (flag_mode == 2) file_cv_out = new TMemFile("cv_unused","RECREATE");

  //Setup the directories specified in the config file
  wrangler_cv.set_new_trees(file_cv_out,true,"_cv");
  wrangler_det.set_new_trees(file3,true,"_det");
  wrangler_pot_cv.set_new_trees(file_cv_out,true,"_cv");
  wrangler_pot_det.set_new_trees(file3,true,"_det");

  // Build the pairs of pot trees
  wrangler_pot_cv.grow_pot_arboretum();
  wrangler_pot_det.grow_pot_arboretum();

  //Always do WC
  file3->mkdir("wcpselection");
  if (file_cv_out != file3) file_cv_out->mkdir("wcpselection");
  file_cv_out->cd("wcpselection");
  TTree *t1_cv = new TTree("T_eval_cv","T_eval_cv");
  TTree *t2_cv = new TTree("T_pot_cv","T_pot_cv");
  TTree *t3_cv = new TTree("T_PFeval_cv", "T_PFeval_cv");
  TTree *t5_cv = new TTree("T_KINEvars_cv", "T_KINEvars_cv");
  TTree *t4_cv = new TTree("T_BDTvars_cv","T_BDTvars_cv");
  TTree *T_spacepoints_new_cv = T_spacepoints_cv->CloneTree(0);
  T_spacepoints_new_cv->SetObject("T_spacepoints_cv","T_spacepoints_cv");

  file3->cd("wcpselection");
  TTree *t1_det = new TTree("T_eval_det","T_eval_det");
  TTree *t2_det = new TTree("T_pot_det","T_pot_det");
  TTree *t3_det = new TTree("T_PFeval_det", "T_PFeval_det");
  TTree *t5_det = new TTree("T_KINEvars_det", "T_KINEvars_det");
  TTree *t4_det = new TTree("T_BDTvars_det","T_BDTvars_det");
  TTree *T_spacepoints_new_det = T_spacepoints_det->CloneTree(0);
  T_spacepoints_new_det->SetObject("T_spacepoints_det","T_spacepoints_det");

  // -m2: for each CV event and each CV POT entry, whether this DetVar file uses it (input for merge_det -m3)
  Int_t flag_run = 0, flag_subrun = 0, flag_event = 0, flag_match = 0;
  Double_t flag_pot = 0;
  TTree *T_rse_flag = 0;
  TTree *T_pot_flag = 0;
  if (flag_mode == 2){
    T_rse_flag = new TTree("T_cv_rse_flag","T_cv_rse_flag");
    T_rse_flag->Branch("run", &flag_run, "run/I");
    T_rse_flag->Branch("subrun", &flag_subrun, "subrun/I");
    T_rse_flag->Branch("event", &flag_event, "event/I");
    T_rse_flag->Branch("flag", &flag_match, "flag/I");
    T_rse_flag->Branch("det_no", &det_no, "det_no/I");
    T_pot_flag = new TTree("T_cv_pot_flag","T_cv_pot_flag");
    T_pot_flag->Branch("runNo", &flag_run, "runNo/I");
    T_pot_flag->Branch("subRunNo", &flag_subrun, "subRunNo/I");
    T_pot_flag->Branch("flag", &flag_match, "flag/I");
    T_pot_flag->Branch("pot_tor875", &flag_pot, "pot_tor875/D");
  }

  EvalInfo eval_cv;
  eval_cv.file_type = new std::string();
  POTInfo pot_cv;
  TaggerInfo tagger_cv;
  PFevalInfo pfeval_cv;
  KineInfo kine_cv;
  kine_cv.kine_energy_particle = new std::vector<float>;
  kine_cv.kine_energy_info = new std::vector<int>;
  kine_cv.kine_particle_type = new std::vector<int>;
  kine_cv.kine_energy_included = new std::vector<int>;


  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_cv.ssm_kine_energy_particle = new std::vector<float>;
  tagger_cv.ssm_kine_energy_info = new std::vector<int>;
  tagger_cv.ssm_kine_particle_type = new std::vector<int>;
  tagger_cv.ssm_kine_energy_included = new std::vector<int>;
  tagger_cv.ssm_cosmict_flag_10 = new std::vector<float>;
  tagger_cv.WCPMTInfoPePred = new std::vector<double>;
  tagger_cv.WCPMTInfoPeMeas = new std::vector<double>;
  tagger_cv.WCPMTInfoPeMeasErr = new std::vector<double>;

  EvalInfo eval_det;
  eval_det.file_type = new std::string();
  POTInfo pot_det;
  TaggerInfo tagger_det;
  PFevalInfo pfeval_det;
  KineInfo kine_det;
  kine_det.kine_energy_particle = new std::vector<float>;
  kine_det.kine_energy_info = new std::vector<int>;
  kine_det.kine_particle_type = new std::vector<int>;
  kine_det.kine_energy_included = new std::vector<int>;


  tagger_det.pio_2_v_dis2 = new std::vector<float>;
  tagger_det.pio_2_v_angle2 = new std::vector<float>;
  tagger_det.pio_2_v_acc_length = new std::vector<float>;
  tagger_det.pio_2_v_flag = new std::vector<float>;
  tagger_det.sig_1_v_angle = new std::vector<float>;
  tagger_det.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_1_v_energy = new std::vector<float>;
  tagger_det.sig_1_v_energy_1 = new std::vector<float>;
  tagger_det.sig_1_v_flag = new std::vector<float>;
  tagger_det.sig_2_v_energy = new std::vector<float>;
  tagger_det.sig_2_v_shower_angle = new std::vector<float>;
  tagger_det.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_flag = new std::vector<float>;
  tagger_det.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_energy = new std::vector<float>;
  tagger_det.stw_2_v_angle = new std::vector<float>;
  tagger_det.stw_2_v_dir_length = new std::vector<float>;
  tagger_det.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_flag = new std::vector<float>;
  tagger_det.stw_3_v_angle = new std::vector<float>;
  tagger_det.stw_3_v_dir_length = new std::vector<float>;
  tagger_det.stw_3_v_energy = new std::vector<float>;
  tagger_det.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_3_v_flag = new std::vector<float>;
  tagger_det.stw_4_v_angle = new std::vector<float>;
  tagger_det.stw_4_v_dis = new std::vector<float>;
  tagger_det.stw_4_v_energy = new std::vector<float>;
  tagger_det.stw_4_v_flag = new std::vector<float>;
  tagger_det.br3_3_v_energy = new std::vector<float>;
  tagger_det.br3_3_v_angle = new std::vector<float>;
  tagger_det.br3_3_v_dir_length = new std::vector<float>;
  tagger_det.br3_3_v_length = new std::vector<float>;
  tagger_det.br3_3_v_flag = new std::vector<float>;
  tagger_det.br3_5_v_dir_length = new std::vector<float>;
  tagger_det.br3_5_v_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_det.br3_5_v_n_seg = new std::vector<float>;
  tagger_det.br3_5_v_angle = new std::vector<float>;
  tagger_det.br3_5_v_sg_length = new std::vector<float>;
  tagger_det.br3_5_v_energy = new std::vector<float>;
  tagger_det.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_det.br3_5_v_n_segs = new std::vector<float>;
  tagger_det.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_det.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag = new std::vector<float>;
  tagger_det.br3_6_v_angle = new std::vector<float>;
  tagger_det.br3_6_v_angle1 = new std::vector<float>;
  tagger_det.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.br3_6_v_direct_length = new std::vector<float>;
  tagger_det.br3_6_v_length = new std::vector<float>;
  tagger_det.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_det.br3_6_v_energy = new std::vector<float>;
  tagger_det.br3_6_v_flag = new std::vector<float>;
  tagger_det.tro_1_v_particle_type = new std::vector<float>;
  tagger_det.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_det.tro_1_v_min_dis = new std::vector<float>;
  tagger_det.tro_1_v_sg1_length = new std::vector<float>;
  tagger_det.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_1_v_tmp_length = new std::vector<float>;
  tagger_det.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_det.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_det.tro_1_v_flag = new std::vector<float>;
  tagger_det.tro_2_v_energy = new std::vector<float>;
  tagger_det.tro_2_v_stem_length = new std::vector<float>;
  tagger_det.tro_2_v_iso_angle = new std::vector<float>;
  tagger_det.tro_2_v_max_length = new std::vector<float>;
  tagger_det.tro_2_v_angle = new std::vector<float>;
  tagger_det.tro_2_v_flag = new std::vector<float>;
  tagger_det.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_det.tro_4_v_angle = new std::vector<float>;
  tagger_det.tro_4_v_angle1 = new std::vector<float>;
  tagger_det.tro_4_v_angle2 = new std::vector<float>;
  tagger_det.tro_4_v_length = new std::vector<float>;
  tagger_det.tro_4_v_length1 = new std::vector<float>;
  tagger_det.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_energy = new std::vector<float>;
  tagger_det.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.tro_4_v_flag = new std::vector<float>;
  tagger_det.tro_5_v_max_angle = new std::vector<float>;
  tagger_det.tro_5_v_min_angle = new std::vector<float>;
  tagger_det.tro_5_v_max_length = new std::vector<float>;
  tagger_det.tro_5_v_iso_angle = new std::vector<float>;
  tagger_det.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_5_v_min_count = new std::vector<float>;
  tagger_det.tro_5_v_max_count = new std::vector<float>;
  tagger_det.tro_5_v_energy = new std::vector<float>;
  tagger_det.tro_5_v_flag = new std::vector<float>;
  tagger_det.lol_1_v_energy = new std::vector<float>;
  tagger_det.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_1_v_nseg = new std::vector<float>;
  tagger_det.lol_1_v_angle = new std::vector<float>;
  tagger_det.lol_1_v_flag = new std::vector<float>;
  tagger_det.lol_2_v_length = new std::vector<float>;
  tagger_det.lol_2_v_angle = new std::vector<float>;
  tagger_det.lol_2_v_type = new std::vector<float>;
  tagger_det.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_2_v_energy = new std::vector<float>;
  tagger_det.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_det.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_det.lol_2_v_flag = new std::vector<float>;
  tagger_det.cosmict_flag_10 = new std::vector<float>;
  tagger_det.cosmict_10_flag_inside = new std::vector<float>;
  tagger_det.cosmict_10_vtx_z = new std::vector<float>;
  tagger_det.cosmict_10_flag_shower = new std::vector<float>;
  tagger_det.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_det.cosmict_10_angle_beam = new std::vector<float>;
  tagger_det.cosmict_10_length = new std::vector<float>;
  tagger_det.numu_cc_flag_1 = new std::vector<float>;
  tagger_det.numu_cc_1_particle_type = new std::vector<float>;
  tagger_det.numu_cc_1_length = new std::vector<float>;
  tagger_det.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_det.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_det.numu_cc_1_direct_length = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_det.numu_cc_flag_2 = new std::vector<float>;
  tagger_det.numu_cc_2_length = new std::vector<float>;
  tagger_det.numu_cc_2_total_length = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_det.pio_2_v_dis2 = new std::vector<float>;
  tagger_det.pio_2_v_angle2 = new std::vector<float>;
  tagger_det.pio_2_v_acc_length = new std::vector<float>;
  tagger_det.pio_2_v_flag = new std::vector<float>;
  tagger_det.sig_1_v_angle = new std::vector<float>;
  tagger_det.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_1_v_energy = new std::vector<float>;
  tagger_det.sig_1_v_energy_1 = new std::vector<float>;
  tagger_det.sig_1_v_flag = new std::vector<float>;
  tagger_det.sig_2_v_energy = new std::vector<float>;
  tagger_det.sig_2_v_shower_angle = new std::vector<float>;
  tagger_det.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_flag = new std::vector<float>;
  tagger_det.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_energy = new std::vector<float>;
  tagger_det.stw_2_v_angle = new std::vector<float>;
  tagger_det.stw_2_v_dir_length = new std::vector<float>;
  tagger_det.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_flag = new std::vector<float>;
  tagger_det.stw_3_v_angle = new std::vector<float>;
  tagger_det.stw_3_v_dir_length = new std::vector<float>;
  tagger_det.stw_3_v_energy = new std::vector<float>;
  tagger_det.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_3_v_flag = new std::vector<float>;
  tagger_det.stw_4_v_angle = new std::vector<float>;
  tagger_det.stw_4_v_dis = new std::vector<float>;
  tagger_det.stw_4_v_energy = new std::vector<float>;
  tagger_det.stw_4_v_flag = new std::vector<float>;
  tagger_det.br3_3_v_energy = new std::vector<float>;
  tagger_det.br3_3_v_angle = new std::vector<float>;
  tagger_det.br3_3_v_dir_length = new std::vector<float>;
  tagger_det.br3_3_v_length = new std::vector<float>;
  tagger_det.br3_3_v_flag = new std::vector<float>;
  tagger_det.br3_5_v_dir_length = new std::vector<float>;
  tagger_det.br3_5_v_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_det.br3_5_v_n_seg = new std::vector<float>;
  tagger_det.br3_5_v_angle = new std::vector<float>;
  tagger_det.br3_5_v_sg_length = new std::vector<float>;
  tagger_det.br3_5_v_energy = new std::vector<float>;
  tagger_det.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_det.br3_5_v_n_segs = new std::vector<float>;
  tagger_det.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_det.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag = new std::vector<float>;
  tagger_det.br3_6_v_angle = new std::vector<float>;
  tagger_det.br3_6_v_angle1 = new std::vector<float>;
  tagger_det.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.br3_6_v_direct_length = new std::vector<float>;
  tagger_det.br3_6_v_length = new std::vector<float>;
  tagger_det.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_det.br3_6_v_energy = new std::vector<float>;
  tagger_det.br3_6_v_flag = new std::vector<float>;
  tagger_det.tro_1_v_particle_type = new std::vector<float>;
  tagger_det.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_det.tro_1_v_min_dis = new std::vector<float>;
  tagger_det.tro_1_v_sg1_length = new std::vector<float>;
  tagger_det.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_1_v_tmp_length = new std::vector<float>;
  tagger_det.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_det.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_det.tro_1_v_flag = new std::vector<float>;
  tagger_det.tro_2_v_energy = new std::vector<float>;
  tagger_det.tro_2_v_stem_length = new std::vector<float>;
  tagger_det.tro_2_v_iso_angle = new std::vector<float>;
  tagger_det.tro_2_v_max_length = new std::vector<float>;
  tagger_det.tro_2_v_angle = new std::vector<float>;
  tagger_det.tro_2_v_flag = new std::vector<float>;
  tagger_det.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_det.tro_4_v_angle = new std::vector<float>;
  tagger_det.tro_4_v_angle1 = new std::vector<float>;
  tagger_det.tro_4_v_angle2 = new std::vector<float>;
  tagger_det.tro_4_v_length = new std::vector<float>;
  tagger_det.tro_4_v_length1 = new std::vector<float>;
  tagger_det.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_energy = new std::vector<float>;
  tagger_det.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.tro_4_v_flag = new std::vector<float>;
  tagger_det.tro_5_v_max_angle = new std::vector<float>;
  tagger_det.tro_5_v_min_angle = new std::vector<float>;
  tagger_det.tro_5_v_max_length = new std::vector<float>;
  tagger_det.tro_5_v_iso_angle = new std::vector<float>;
  tagger_det.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_5_v_min_count = new std::vector<float>;
  tagger_det.tro_5_v_max_count = new std::vector<float>;
  tagger_det.tro_5_v_energy = new std::vector<float>;
  tagger_det.tro_5_v_flag = new std::vector<float>;
  tagger_det.lol_1_v_energy = new std::vector<float>;
  tagger_det.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_1_v_nseg = new std::vector<float>;
  tagger_det.lol_1_v_angle = new std::vector<float>;
  tagger_det.lol_1_v_flag = new std::vector<float>;
  tagger_det.lol_2_v_length = new std::vector<float>;
  tagger_det.lol_2_v_angle = new std::vector<float>;
  tagger_det.lol_2_v_type = new std::vector<float>;
  tagger_det.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_2_v_energy = new std::vector<float>;
  tagger_det.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_det.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_det.lol_2_v_flag = new std::vector<float>;
  tagger_det.cosmict_flag_10 = new std::vector<float>;
  tagger_det.cosmict_10_flag_inside = new std::vector<float>;
  tagger_det.cosmict_10_vtx_z = new std::vector<float>;
  tagger_det.cosmict_10_flag_shower = new std::vector<float>;
  tagger_det.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_det.cosmict_10_angle_beam = new std::vector<float>;
  tagger_det.cosmict_10_length = new std::vector<float>;
  tagger_det.numu_cc_flag_1 = new std::vector<float>;
  tagger_det.numu_cc_1_particle_type = new std::vector<float>;
  tagger_det.numu_cc_1_length = new std::vector<float>;
  tagger_det.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_det.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_det.numu_cc_1_direct_length = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_det.numu_cc_flag_2 = new std::vector<float>;
  tagger_det.numu_cc_2_length = new std::vector<float>;
  tagger_det.numu_cc_2_total_length = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_det.ssm_kine_energy_particle = new std::vector<float>;
  tagger_det.ssm_kine_energy_info = new std::vector<int>;
  tagger_det.ssm_kine_particle_type = new std::vector<int>;
  tagger_det.ssm_kine_energy_included = new std::vector<int>;
  tagger_det.ssm_cosmict_flag_10 = new std::vector<float>;
  tagger_det.WCPMTInfoPePred = new std::vector<double>;
  tagger_det.WCPMTInfoPeMeas = new std::vector<double>;
  tagger_det.WCPMTInfoPeMeasErr = new std::vector<double>;


  set_tree_address(T_BDTvars_cv, tagger_cv,2 );
  set_tree_address(T_eval_cv, eval_cv);
  set_tree_address(T_PFeval_cv, pfeval_cv);
  set_tree_address(T_pot_cv, pot_cv);
  set_tree_address(T_KINEvars_cv, kine_cv);


  put_tree_address(t4_cv, tagger_cv,2);
  put_tree_address(t1_cv, eval_cv);
  put_tree_address(t3_cv, pfeval_cv);
  put_tree_address(t2_cv, pot_cv);
  put_tree_address(t5_cv, kine_cv);

  set_tree_address(T_BDTvars_det, tagger_det,2 );
  set_tree_address(T_eval_det, eval_det);
  set_tree_address(T_PFeval_det, pfeval_det);
  set_tree_address(T_pot_det, pot_det);
  set_tree_address(T_KINEvars_det, kine_det);

  put_tree_address(t4_det, tagger_det,2);
  put_tree_address(t1_det, eval_det);
  put_tree_address(t3_det, pfeval_det);
  put_tree_address(t2_det, pot_det);
  put_tree_address(t5_det, kine_det);
  // Not filled nominally, so force turn them one. 
  if(pfeval_cv.flag_ns_time_cor){
    t3_det->Branch("evtTimeNS_cor",&pfeval_det.evtTimeNS_cor);
    t3_det->Branch("cor_nu_time",&pfeval_det.cor_nu_time);
    t3_det->Branch("cor_nu_time_nospill",&pfeval_det.cor_nu_time_nospill);
    t3_det->Branch("cor_nu_time_spill",&pfeval_det.cor_nu_time_spill);
    t3_det->Branch("cor_nu_deltatime",&pfeval_det.cor_nu_deltatime);
  }

  std::map<std::pair<int, int>, int> map_re_entry_cv;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_cv;

  bool flag_presel = false;
  std::vector<std::tuple<int, int, int> > vec_rse_cv; // run, subrun, event of each CV entry (-m2)

  // Find the index of the "goofy" events for the cv, this will be common across all trees
  for (int i=0;i!=T_eval_cv->GetEntries();i++){
    T_eval_cv->GetEntry(i);
    T_BDTvars_cv->GetEntry(i);

    map_rs_re_cv[std::make_pair(eval_cv.run, eval_cv.subrun)].insert(std::make_pair(eval_cv.run, eval_cv.event));
    vec_rse_cv.push_back(std::make_tuple(eval_cv.run, eval_cv.subrun, eval_cv.event));

    int tmp_match_found = eval_cv.match_found;
    if (eval_cv.is_match_found_int) tmp_match_found = eval_cv.match_found_asInt;

    flag_presel = false;
    if (tmp_match_found != 0 && eval_cv.stm_eventtype != 0 && eval_cv.stm_lowenergy ==0 && eval_cv.stm_LM ==0 && eval_cv.stm_TGM ==0 && eval_cv.stm_STM==0 && eval_cv.stm_FullDead == 0 && eval_cv.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }

    if (tmp_match_found == -1  || (tmp_match_found == 1 && eval_cv.stm_lowenergy == -1) || (flag_presel && tagger_cv.numu_cc_flag == -1)) continue;

    map_re_entry_cv[std::make_pair(eval_cv.run, eval_cv.event)] = i;

  }

  // Find the index of the "goofy" events for the cv, this will be common across all trees
  std::map<std::pair<int, int>, int> map_re_entry_det;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_det;
  for (int i=0;i!=T_eval_det->GetEntries();i++){
    T_eval_det->GetEntry(i);
    T_BDTvars_det->GetEntry(i);

    map_rs_re_det[std::make_pair(eval_det.run, eval_det.subrun)].insert(std::make_pair(eval_det.run, eval_det.event));

    int tmp_match_found = eval_det.match_found;
    if (eval_det.is_match_found_int) tmp_match_found = eval_det.match_found_asInt;

    flag_presel = false;
    if (tmp_match_found != 0 && eval_det.stm_eventtype != 0 && eval_det.stm_lowenergy ==0 && eval_det.stm_LM ==0 && eval_det.stm_TGM ==0 && eval_det.stm_STM==0 && eval_det.stm_FullDead == 0 && eval_det.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }

    if (tmp_match_found == -1  || (tmp_match_found == 1&& eval_det.stm_lowenergy == -1)  || (flag_presel && tagger_det.numu_cc_flag == -1)) continue;

    map_re_entry_det[std::make_pair(eval_det.run, eval_det.event)] = i;

  }

  // Map out the relation between index and run-subrun for WC CV
  std::map<std::tuple<int, int>, std::pair<int, double> > map_rs_entry_pot_cv;
  std::vector<std::pair<int, int> > vec_rs_pot_cv; // run, subrun of each CV POT entry (-m2)
  for (Int_t i=0;i!=T_pot_cv->GetEntries();i++){
    T_pot_cv->GetEntry(i);
    vec_rs_pot_cv.push_back(std::make_pair(pot_cv.runNo, pot_cv.subRunNo));
    map_rs_entry_pot_cv[std::make_pair(pot_cv.runNo,pot_cv.subRunNo)] = std::make_pair(i, pot_cv.pot_tor875);
  }

  // Map out the relation between index and run-subrun for WC DetVar
  std::map<std::tuple<int, int>, std::pair<int, double> > map_rs_entry_pot_det;
  for (Int_t i=0;i!=T_pot_det->GetEntries();i++){
    T_pot_det->GetEntry(i);
    map_rs_entry_pot_det[std::make_pair(pot_det.runNo,pot_det.subRunNo)] = std::make_pair(i, pot_det.pot_tor875);
  }

  // Map out the relation between index and run-subrun for the non-WC POT CV trees
  wrangler_pot_cv.map_rs_to_entry();
  // Map out the relation between index and run-subrun for the non-WC POT DetVar trees
  wrangler_pot_det.map_rs_to_entry();
  // Check that these are one to one
  for(int arb_index=0; arb_index<wrangler_pot_cv.pot_arboretum->size(); arb_index++){
    std::string tree_name_cv = wrangler_pot_cv.pot_arboretum->at(arb_index)->old_pot_tree->GetName();
    std::string tree_name_det = wrangler_pot_det.pot_arboretum->at(arb_index)->old_pot_tree->GetName();
    if(tree_name_cv!=tree_name_det){
      std::cout<<"ERROR! Mismatch in pot trees"<<std::endl;
      return 1;
    }
    tree_name_cv = wrangler_pot_cv.pot_arboretum->at(arb_index)->new_pot_tree->GetName();
    tree_name_det = wrangler_pot_det.pot_arboretum->at(arb_index)->new_pot_tree->GetName();
    removeSubstring(tree_name_cv, "_cv");
    removeSubstring(tree_name_det, "_det");
    if(tree_name_cv!=tree_name_det){
      std::cout<<"ERROR! Mismatch in pot trees"<<std::endl;
      return 1;
    }
  }

   T_eval_cv->SetBranchStatus("*",1);
   T_PFeval_cv->SetBranchStatus("*",1);
   T_BDTvars_cv->SetBranchStatus("*",1);
   T_KINEvars_cv->SetBranchStatus("*",1);

   T_eval_det->SetBranchStatus("*",1);
   T_PFeval_det->SetBranchStatus("*",1);
   T_BDTvars_det->SetBranchStatus("*",1);
   T_KINEvars_det->SetBranchStatus("*",1);



  // fill the trees ...
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_common;
  std::map<int, int> map_cv_det_index;

  for (auto it = map_re_entry_cv.begin(); it != map_re_entry_cv.end(); it++){
    auto it1 = map_re_entry_det.find(it->first);
    if (it1 != map_re_entry_det.end()){ // common ...
      map_cv_det_index[it->second] = it1->second;
    }
  }

  int nentries = map_cv_det_index.size();
  int ientry = 0;
  std::cout<<"Begin looping over "<<nentries<<" events"<<std::endl;
  for (auto it = map_cv_det_index.begin(); it != map_cv_det_index.end(); it++){

      if (ientry%10000 == 0) std::cout << ientry/1000 << " k " << std::setprecision(3) << double(ientry)/nentries*100. << " %"<< std::endl;
      ientry++;

      // -m2 only needs the CV run/subrun/event and time correction
      if (flag_mode != 2) T_BDTvars_cv->GetEntry(it->first);
      T_eval_cv->GetEntry(it->first);
      if (flag_mode != 2) T_KINEvars_cv->GetEntry(it->first);
      T_PFeval_cv->GetEntry(it->first);

      T_BDTvars_det->GetEntry(it->second);
      T_eval_det->GetEntry(it->second);
      T_KINEvars_det->GetEntry(it->second);
      T_PFeval_det->GetEntry(it->second);

      if (flag_mode != 2) T_spacepoints_cv->GetEntry(it->first);

      T_spacepoints_det->GetEntry(it->second);

      if (flag_mode != 2){
        for(auto tree_it=wrangler_cv.old_trees->begin(); tree_it!=wrangler_cv.old_trees->end(); tree_it++){
            (*tree_it)->GetEntry(it->first);
        }
      }
      for(auto tree_it=wrangler_det.old_trees->begin(); tree_it!=wrangler_det.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(it->second);
      }


      map_rs_re_common[std::make_pair(eval_cv.run, eval_cv.subrun)].insert(std::make_pair(eval_cv.run, eval_cv.event));

      // We have to set these by hand for the detvar to make sure the beamspill time reassignment is the same. 
      // Pandora should probably have similar code?
      if(pfeval_cv.flag_ns_time_cor){
        pfeval_det.evtTimeNS_cor = pfeval_det.evtTimeNS + pfeval_cv.cor_nu_deltatime;
        pfeval_det.cor_nu_time = pfeval_cv.cor_nu_time;
        pfeval_det.cor_nu_time_nospill = pfeval_cv.cor_nu_time_nospill;
        pfeval_det.cor_nu_time_spill = pfeval_cv.cor_nu_time_spill;
        pfeval_det.cor_nu_deltatime = pfeval_cv.cor_nu_deltatime;
      }

      if (flag_mode != 2){
        t1_cv->Fill();
        t3_cv->Fill();
        t4_cv->Fill();
        t5_cv->Fill();
      }

      t1_det->Fill();
      t3_det->Fill();
      t4_det->Fill();
      t5_det->Fill();

      if (flag_mode != 2) T_spacepoints_new_cv->Fill();

      T_spacepoints_new_det->Fill();

      if (flag_mode != 2){
        for(auto tree_it=wrangler_cv.new_trees->begin(); tree_it!=wrangler_cv.new_trees->end(); tree_it++){
            (*tree_it)->Fill();
        }
      }
      for(auto tree_it=wrangler_det.new_trees->begin(); tree_it!=wrangler_det.new_trees->end(); tree_it++){
          (*tree_it)->Fill();
      }

  }

  // Add up the pot for both the cv and det files before merging
  std::vector<double> vec_cv_pot;
  std::vector<double> vec_det_pot;
  double cv_pot=0;
  double det_pot=0;
  // First for WC
  for (auto it = map_rs_entry_pot_cv.begin(); it != map_rs_entry_pot_cv.end(); it++){
    cv_pot += it->second.second;
  }
  for (auto it = map_rs_entry_pot_det.begin(); it != map_rs_entry_pot_det.end(); it++){
    det_pot += it->second.second;
  }
  vec_cv_pot.push_back(cv_pot);
  vec_det_pot.push_back(det_pot);
  // Then for the other trees
  for (auto it_map_rs_entry = wrangler_pot_cv.arboretum_map_rs_entry.begin(); it_map_rs_entry != wrangler_pot_cv.arboretum_map_rs_entry.end(); it_map_rs_entry++){
    cv_pot=0;
    for (auto it = (*it_map_rs_entry).begin(); it != (*it_map_rs_entry).end(); it++){
      cv_pot += it->second.second;
    }
    vec_cv_pot.push_back(cv_pot);
  }
  for (auto it_map_rs_entry = wrangler_pot_det.arboretum_map_rs_entry.begin(); it_map_rs_entry != wrangler_pot_det.arboretum_map_rs_entry.end(); it_map_rs_entry++){
    det_pot=0;
    for (auto it = (*it_map_rs_entry).begin(); it != (*it_map_rs_entry).end(); it++){
      det_pot += it->second.second;
    }
    vec_det_pot.push_back(cv_pot);
  }

  // Loop over each POT tree seperatly
  // Start with WireCell
  nentries =  map_rs_re_common.size();
  ientry=0;
  std::cout<<"Begin looping over WC pot tree checking "<<nentries<<" entries"<<std::endl;
  std::vector<double> vec_common_pot;
  double common_pot = 0;
  std::map<int, double> map_cv_pot_index; // CV POT entry used --> POT of the subrun saved (-m2)
  for (auto it = map_rs_re_common.begin(); it != map_rs_re_common.end(); it++){

    if (ientry%10000 == 0) std::cout << ientry/1000 << " k " << std::setprecision(3) << double(ientry)/nentries*100. << " %"<< std::endl;
    ientry++;

    auto it1 = map_rs_entry_pot_cv.find(it->first);
    auto it2 = map_rs_entry_pot_det.find(it->first);

    if (it1 != map_rs_entry_pot_cv.end() && it2 != map_rs_entry_pot_det.end()){
      T_pot_cv->GetEntry(it1->second.first);
      T_pot_det->GetEntry(it2->second.first);

      double ratio = it->second.size() * 1.0/map_rs_re_cv[it->first].size();
      pot_cv.pot_tor875 *= ratio;
      pot_cv.pot_tor875good *= ratio;

      ratio = it->second.size() * 1.0 /map_rs_re_det[it->first].size();
      pot_det.pot_tor875 *= ratio;
      pot_det.pot_tor875good *= ratio;

      common_pot += pot_cv.pot_tor875;

      if(override_pot){
        pot_cv.pot_tor875=5e18;
        pot_cv.pot_tor875good=5e18;
        pot_det.pot_tor875=5e18;
        pot_det.pot_tor875good=5e18;
      }
      map_cv_pot_index[it1->second.first] = pot_cv.pot_tor875;

      t2_cv->Fill();
      t2_det->Fill();
    }
  }
  vec_common_pot.push_back(common_pot);
  common_pot = 0;

  // Now the other trees
  for(int arb_index=0; arb_index<wrangler_pot_cv.pot_arboretum->size(); arb_index++){

    std::cout<<"Begin looping over WC pot tree checking "<<nentries<<" entries"<<std::endl;
    ientry=0;
    for (auto it_re = map_rs_re_common.begin(); it_re != map_rs_re_common.end(); it_re++){

      if (ientry%10000 == 0) std::cout << ientry/1000 << " k " << std::setprecision(3) << double(ientry)/nentries*100. << " %"<< std::endl;
      ientry++;

      auto it_cv = wrangler_pot_cv.arboretum_map_rs_entry.at(arb_index).find(it_re->first);
      auto it_det = wrangler_pot_det.arboretum_map_rs_entry.at(arb_index).find(it_re->first);

      if (it_cv != wrangler_pot_cv.arboretum_map_rs_entry.at(arb_index).end() && it_det != wrangler_pot_det.arboretum_map_rs_entry.at(arb_index).end()){
        wrangler_pot_cv.pot_arboretum->at(arb_index)->old_pot_tree->GetEntry(it_cv->second.first);      
        wrangler_pot_det.pot_arboretum->at(arb_index)->old_pot_tree->GetEntry(it_det->second.first);

        // For both the cv and the detvar, set both the double and the float, only the correct one will fill the tree
        double ratio = it_re->second.size() * 1.0/map_rs_re_cv[it_re->first].size();
        wrangler_pot_cv.pot_arboretum->at(arb_index)->fpot = wrangler_pot_cv.pot_arboretum->at(arb_index)->pot()*ratio;
        wrangler_pot_cv.pot_arboretum->at(arb_index)->dpot = wrangler_pot_cv.pot_arboretum->at(arb_index)->pot()*ratio;

        ratio = it_re->second.size() * 1.0 /map_rs_re_det[it_re->first].size();
        wrangler_pot_det.pot_arboretum->at(arb_index)->fpot = wrangler_pot_det.pot_arboretum->at(arb_index)->pot()*ratio;
        wrangler_pot_det.pot_arboretum->at(arb_index)->dpot = wrangler_pot_det.pot_arboretum->at(arb_index)->pot()*ratio;

        common_pot += wrangler_pot_cv.pot_arboretum->at(arb_index)->pot()*ratio;

      if(override_pot){
        wrangler_pot_cv.pot_arboretum->at(arb_index)->fpot = 5e18; 
        wrangler_pot_cv.pot_arboretum->at(arb_index)->dpot = 5e18;
        wrangler_pot_det.pot_arboretum->at(arb_index)->fpot = 5e18;
        wrangler_pot_det.pot_arboretum->at(arb_index)->dpot = 5e18;
      }

      wrangler_pot_cv.pot_arboretum->at(arb_index)->new_pot_tree->Fill();
      wrangler_pot_det.pot_arboretum->at(arb_index)->new_pot_tree->Fill();

    } 

  } // it_re, loop over all runs in the run subrun map 

  vec_common_pot.push_back(common_pot);
  common_pot = 0;

} // arb_index, loop over all trees in the arboretum


  int nevents_cv = 0;
  int nevents_det = 0;
  for (auto it = map_rs_re_cv.begin(); it != map_rs_re_cv.end(); it++){
    nevents_cv += it->second.size();
  }
  for (auto it = map_rs_re_det.begin(); it != map_rs_re_det.end(); it++){
    nevents_det += it->second.size();
  }
  std::cout << out_file << std::endl;
  std::cout << "Events: " << t1_det->GetEntries() << " " << nevents_cv<<"/"<<T_eval_cv->GetEntries() << " " << nevents_det << "/" << T_eval_det->GetEntries() << std::endl;
  for(int i=0; i<vec_cv_pot.size(); i++){
    std::cout << "POT:    " << vec_common_pot.at(i) << " " << vec_cv_pot.at(i) << " " << vec_det_pot.at(i) <<std::endl;
  }

  if (flag_mode == 2){
    for (size_t i=0;i!=vec_rse_cv.size();i++){
      flag_run = std::get<0>(vec_rse_cv.at(i));
      flag_subrun = std::get<1>(vec_rse_cv.at(i));
      flag_event = std::get<2>(vec_rse_cv.at(i));
      flag_match = map_cv_det_index.find(i) != map_cv_det_index.end();
      T_rse_flag->Fill();
    }
    for (size_t i=0;i!=vec_rs_pot_cv.size();i++){
      flag_run = vec_rs_pot_cv.at(i).first;
      flag_subrun = vec_rs_pot_cv.at(i).second;
      auto it = map_cv_pot_index.find(i);
      flag_match = it != map_cv_pot_index.end();
      flag_pot = flag_match ? it->second : 0;
      T_pot_flag->Fill();
    }
    std::cout << "DetVar " << det_no << ": flags of " << T_rse_flag->GetEntries() << " CV events (" << map_cv_det_index.size() << " used) and "
              << T_pot_flag->GetEntries() << " CV subruns (" << map_cv_pot_index.size() << " used)" << std::endl;
  }




  file3->Write("",TFile::kOverwrite);
  file3->Close();


  return 0;
}
