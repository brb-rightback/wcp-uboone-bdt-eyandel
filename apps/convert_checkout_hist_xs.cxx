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
  TString out_filename =  argv[2];


  bool flag_data = true;

  TFile *file = new TFile(input_filename,"READ");

  CovMatrix cov;
  cov.add_xs_config();

  std::cout << "Xs mode: " << std::endl;

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
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;

  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;


  // prepare histograms ...
  // declare histograms ...
  TH1F *htemp = 0; TH1F *htemp1 = 0; TH2F *htemp2 = 0; int num = 0;
  std::map<TString, std::tuple<TH1F*, TH1F*, TH2F*, int> > map_histoname_hists;

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
    num = 1;

    bool flag_spec = cov.is_xs_chname(ch_name) ;
    if (flag_spec){
      int nbins1 = cov.get_xs_nsignals();
      TString temp_histoname = histoname + "_signal";
      htemp1 = new TH1F(temp_histoname, temp_histoname,nbins1, 0.5,nbins1+0.5);
      temp_histoname = histoname + "_R";
      htemp2 = new TH2F(temp_histoname, temp_histoname,nbin, llimit, hlimit,nbins1, 0.5,nbins1+0.5);
      num = 3;
    }
    map_histoname_hists[histoname] = std::make_tuple(htemp, htemp1, htemp2, num);
  }
  //  std::cout << std::endl;

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

    map_histoname_hists[histoname] = std::make_tuple(htemp, (TH1F*)(0), (TH2F*)(0), 1);
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
    map_histoname_hists[histoname] = std::make_tuple(htemp, (TH1F*)(0), (TH2F*)(0),1);
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

    //    if (!is_preselection(eval)) continue;

    // event-level quantities for get_cut_pass, computed once per event
    CutEventInfo cut_info;
    fill_cut_event_info(cut_info, flag_data, eval, pfeval, tagger, kine, space, pandora, lantern);
    // xs signal bin of this event, computed once (on first need) for all histograms; -2 = not computed yet
    int signal_bin_event = -2;

    for (auto it = all_histo_infos.begin(); it != all_histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      Int_t nbin = std::get<1>(*it);
      float llimit = std::get<2>(*it);
      float hlimit = std::get<3>(*it);
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      TString weight = std::get<7>(*it);

      // get pass or not
      bool flag_pass = get_cut_pass(ch_name, add_cut, flag_data, cut_info, eval, pfeval, tagger, kine, space, pandora, lantern);
      // get kinematics variable (only needed if the event passes; always on the first event so an unknown variable name still stops the job)
      double val = 0;
      if (flag_pass || flag_first_event) val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name, space, pandora, lantern);
      int signal_bin = -1;
      if (cov.is_xs_chname(ch_name)){
	if (signal_bin_event == -2) signal_bin_event = get_xs_signal_no(cov.get_cut_file(), cov.get_map_cut_xs_bin(), eval, pfeval, tagger, kine);
	signal_bin = signal_bin_event;
      }

      if (!((signal_bin != -1) || flag_pass)) continue;
      // get weight ...
      double weight_val = get_weight(weight, eval, pfeval, kine, tagger, cov.get_rw_info(), cov.get_time_info_allruns(), flag_data);

      //  if (ch_name == "numuCC_signal_Enu_FC_overlay" && weight == "cv_spline") std::cout << "Xin: " << " " << flag_pass << " " << signal_bin << " " << weight_val << " " <<eval.run << " " << eval.subrun << " " << eval.event << std::endl;


      //     htemp = map_histoname_hist[histoname];
      auto tmp_hists = map_histoname_hists[histoname];
      TH1F *h1 = std::get<0>(tmp_hists);
      TH1F *h2 = std::get<1>(tmp_hists);
      TH2F *h3 = std::get<2>(tmp_hists);
      int num = std::get<3>(tmp_hists);
      if (num==1){
       	if (flag_pass) h1->Fill(val,  weight_val);
      }else{
       	if (signal_bin != -1){
       	  if (flag_pass) h1->Fill(val, weight_val);
       	  h2->Fill(signal_bin, weight_val);
       	  if (flag_pass) h3->Fill(val, signal_bin, weight_val);
       	}else{
       	  std::cout << "[convt-hist-xs] Something wrong: cut/channel mismatch !" << std::endl;
       	}
      } // else



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

  for (auto it = map_histoname_hists.begin(); it!= map_histoname_hists.end(); it++){
  //   //std::cout<<"DEBUG: "<<it->first<<" "<<it->second->GetName()<<" "<<it->second->GetSum()<<"\n";
    TH1F *h1 = std::get<0>(it->second);
    TH1F *h2 = std::get<1>(it->second);
    TH2F *h3 = std::get<2>(it->second);
    int num = std::get<3>(it->second);
    if (num==1){
      h1->SetDirectory(file1);
    }else{
      h1->SetDirectory(file1);
      h2->SetDirectory(file1);
      h3->SetDirectory(file1);
    }
    //   it->second->SetDirectory(file1);
  }

  file1->Write();
  file1->Close();





  return 0;
}
