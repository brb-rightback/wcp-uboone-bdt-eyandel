#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/eval.h"
#include "WCPLEEANA/space.h"
#include "WCPLEEANA/pandora.h"
#include "WCPLEEANA/lantern.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{

  TString input_filename = argv[1];
  TString out_filename = argv[2];
  bool flag_data = true;

  bool flag_osc = false;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'o':
      flag_osc = atoi(&argv[i][2]); // run oscillation
      break;
    }
  }


  TFile *file = new TFile(input_filename,"READ");

  CovMatrix cov;
  if (flag_osc) cov.add_osc_config();

  cov.print_rw(cov.get_rw_info());

  AnalysisTrees trees = get_analysis_trees(file);
  TTree *T_BDTvars = trees.T_BDTvars, *T_eval = trees.T_eval, *T_PFeval = trees.T_PFeval, *T_KINEvars = trees.T_KINEvars;
  TTree *T_spacepoints = trees.T_spacepoints, *T_pandora = trees.T_pandora, *T_lantern = trees.T_lantern;
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");

  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;
  SpaceInfo space;
  PandoraInfo pandora;
  LanternInfo lantern;

#include "init.txt"

  set_tree_address(T_BDTvars, tagger,2 );
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);
  if(T_spacepoints) set_tree_address(T_spacepoints, space, 0);
  if(T_pandora) set_tree_address(T_pandora, pandora);
  if(T_lantern) set_tree_address(T_lantern, lantern);

  double total_pot = 0;
  //Erin
  std::map<std::pair<int, int>, bool >  already_seen;

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    auto it = already_seen.find(std::make_pair(pot.runNo,pot.subRunNo));
    if (it != already_seen.end()) continue;
    already_seen[std::make_pair(pot.runNo,pot.subRunNo)] = true;

    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;

  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;


  // prepare histograms ...
  // declare histograms ...
  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
  std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "CV:" << std::endl;
  for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);

    //    std::cout << std::get<0>( *it)  << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
  std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "Error2: " << std::endl;
  for (auto it = histo_infos_err2.begin(); it != histo_infos_err2.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
  std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
  //std::cout << "Cross: " << std::endl;
  for (auto it = histo_infos_cros.begin(); it != histo_infos_cros.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  // fill histogram ...
  // switch on the branches the analysis code reads (one list for all apps, see analysis_trees.h)
  set_analysis_branch_status(trees, flag_data);

  std::cout << "Total entries: " << T_eval->GetEntries() << std::endl;


  bool flag_first_event = true;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    if(T_spacepoints) T_spacepoints->GetEntry(i);
    if(T_pandora) T_pandora->GetEntry(i);
    if(T_lantern) T_lantern->GetEntry(i);

    if (!is_preselection(eval)) continue;

    // event-level quantities for get_cut_pass, computed once per event
    CutEventInfo cut_info;
    fill_cut_event_info(cut_info, flag_data, eval, pfeval, tagger, kine, space, pandora, lantern);

    for (auto it = all_histo_infos.begin(); it != all_histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      Int_t nbin = std::get<1>(*it);
      float llimit = std::get<2>(*it);
      float hlimit = std::get<3>(*it);
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      TString weight = std::get<7>(*it);

      htemp = map_histoname_hist[histoname];
      // get pass or not
      bool flag_pass = get_cut_pass(ch_name, add_cut, flag_data, cut_info, eval, pfeval, tagger, kine, space, pandora, lantern);
      // get kinematics variable (only needed if the event passes; always on the first event so an unknown variable name still stops the job)
      double val = 0;
      if (flag_pass || flag_first_event) val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name, space, pandora, lantern);

      double osc_weight = 1.0;



      // std::cout << weight << std::endl;
      // get weight ...
      double weight_val = get_weight(weight, eval, pfeval, kine, tagger, cov.get_rw_info(), cov.get_time_info_allruns() ,flag_data);

      if (flag_osc && cov.is_osc_channel(ch_name) && (!flag_data)){
	osc_weight = cov.get_osc_weight(eval, pfeval);
	weight_val *= osc_weight;
	if (weight == "cv_spline_cv_spline" || weight == "unity_unity" ||
	    weight == "spline_spline" )
	  weight_val *= osc_weight;
      }

      if (flag_pass)
	htemp->Fill(val,weight_val);
    }
    flag_first_event = false;
  }


  // save histograms ...
  TFile *file1 = new TFile(out_filename,"RECREATE");
  file1->cd();
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  T->Branch("pot",&total_pot,"pot/D");
  T->Fill();

  for (auto it = map_histoname_hist.begin(); it!= map_histoname_hist.end(); it++){
    //std::cout<<"DEBUG: "<<it->first<<" "<<it->second->GetName()<<" "<<it->second->GetSum()<<"\n";
    it->second->SetDirectory(file1);
  }

  file1->Write();
  file1->Close();





  return 0;
}
