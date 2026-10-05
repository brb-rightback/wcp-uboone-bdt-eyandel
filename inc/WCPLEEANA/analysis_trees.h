#ifndef UBOONE_LEE_ANALYSIS_TREES
#define UBOONE_LEE_ANALYSIS_TREES

// The trees read by the histogram and covariance apps (convert_checkout_hist(_xs), xs/xf/det/stat covariance), and the
// branches the analysis code (cuts.h, the BDT bins, the weights, ...) reads from them. Kept in one place so every app
// reads the same branches: add a branch here when new code needs it.

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <vector>
#include <string>

namespace LEEana{
  struct AnalysisTrees{
    TTree *T_BDTvars;
    TTree *T_eval;
    TTree *T_PFeval;
    TTree *T_KINEvars;
    TTree *T_spacepoints;   // 0 if the file does not have it
    TTree *T_pandora;       // 0 if the file does not have it
    TTree *T_lantern;       // 0 if the file does not have it
    TTree *T_glee;          // singlephotonana/eventweight_tree, 0 if the file does not have it
  };

  // suffix: "" for normal checkout files, "_cv" or "_det" for the merged detector-variation files
  AnalysisTrees get_analysis_trees(TFile* file, TString suffix = "");

  // Switch off all branches, then switch on the ones the analysis code reads (branches a file does not have are skipped).
  // flag_data: data file, the truth branches are not switched on
  void set_analysis_branch_status(AnalysisTrees& trees, bool flag_data);

  void enable_branches(TTree* tree, const std::vector<std::string>& branches);
}

LEEana::AnalysisTrees LEEana::get_analysis_trees(TFile* file, TString suffix){
  AnalysisTrees trees;
  trees.T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars" + suffix);
  trees.T_eval = (TTree*)file->Get("wcpselection/T_eval" + suffix);
  trees.T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval" + suffix);
  trees.T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars" + suffix);
  trees.T_spacepoints = (TTree*)file->Get("wcpselection/T_spacepoints" + suffix);
  trees.T_pandora = (TTree*)file->Get("nuselection/NeutrinoSelectionFilter" + suffix);
  trees.T_lantern = (TTree*)file->Get("lantern/EventTree" + suffix);
  trees.T_glee = (TTree*)file->Get("singlephotonana/eventweight_tree" + suffix);
  return trees;
}

void LEEana::enable_branches(TTree* tree, const std::vector<std::string>& branches){
  if (!tree) return;
  for (size_t i=0; i!=branches.size(); i++){
    if (tree->GetBranch(branches.at(i).c_str())) tree->SetBranchStatus(branches.at(i).c_str(), 1);
  }
}

void LEEana::set_analysis_branch_status(AnalysisTrees& trees, bool flag_data){
  static const std::vector<std::string> bdt_branches = {
    "numu_cc_flag", "numu_score", "nue_score", "cosmict_flag", "nc_delta_score", "nc_pio_score", "all_veto_score",
    "VtxAct_bdt_score"};
  static const std::vector<std::string> eval_branches = {
    "match_isFC", "match_found", "match_found_asInt", "stm_eventtype", "stm_lowenergy", "stm_LM", "stm_TGM",
    "stm_STM", "stm_FullDead", "stm_clusterlength", "run", "event"};
  static const std::vector<std::string> eval_truth_branches = {
    "weight_spline", "weight_cv", "weight_lee", "truth_isCC", "truth_nuPdg", "truth_vtxInside", "truth_nuEnergy",
    "truth_energyInside", "truth_vtxX", "truth_vtxY", "truth_vtxZ", "truth_nuTime", "match_completeness_energy"};
  static const std::vector<std::string> pfeval_branches = {
    "run", "reco_nuvtxX", "reco_nuvtxY", "reco_nuvtxZ", "reco_muonMomentum", "reco_showerKE", "reco_larpid_pdg",
    "reco_Ntrack", "reco_id", "reco_pdg", "reco_mother", "reco_startXYZT", "reco_endXYZT", "reco_startMomentum",
    "mcs_emu_MCS", "mcs_emu_tracklen"};
  static const std::vector<std::string> pfeval_truth_branches = {
    "truth_muonMomentum", "truth_Ntrack", "truth_pdg", "truth_mother", "truth_startXYZT", "truth_startMomentum", "truth_NCDelta",
    "truth_NprimPio", "truth_nuScatType", "truth_nu_momentum", "mcflux_ntype", "mcflux_dk2gen", "mcflux_gen2vtx"};
  static const std::vector<std::string> kine_branches = {
    "kine_reco_Enu", "kine_energy_particle", "kine_particle_type", "kine_energy_info", "kine_reco_add_energy",
    "kine_pio_mass", "kine_pio_flag", "kine_pio_vtx_dis", "kine_pio_energy_1", "kine_pio_dis_1", "kine_pio_energy_2",
    "kine_pio_dis_2", "kine_pio_angle"};
  static const std::vector<std::string> space_branches = {
    "Trecchargeblob_spacepoints_x", "Trecchargeblob_spacepoints_y", "Trecchargeblob_spacepoints_z",
    "Trecchargeblob_spacepoints_q", "Trecchargeblob_spacepoints_real_cluster_id"};
  static const std::vector<std::string> pandora_branches = {
    "slice_orig_pass_id", "n_pfps", "trk_llr_pid_score_v", "pfp_generation_v", "trk_energy_proton_v", "pfpdg"};
  static const std::vector<std::string> lantern_branches = {
    "nTracks", "trackIsSecondary", "trackPID", "trackRecoE", "trackDistToVtx"};
  static const std::vector<std::string> glee_truth_branches = {
    "GTruth_ResNum"};

  TTree *all_trees[8] = {trees.T_BDTvars, trees.T_eval, trees.T_PFeval, trees.T_KINEvars, trees.T_spacepoints, trees.T_pandora, trees.T_lantern, trees.T_glee};
  for (int i=0; i!=8; i++){
    if (all_trees[i]) all_trees[i]->SetBranchStatus("*", 0);
  }

  enable_branches(trees.T_BDTvars, bdt_branches);
  enable_branches(trees.T_eval, eval_branches);
  enable_branches(trees.T_PFeval, pfeval_branches);
  enable_branches(trees.T_KINEvars, kine_branches);
  enable_branches(trees.T_spacepoints, space_branches);
  enable_branches(trees.T_pandora, pandora_branches);
  enable_branches(trees.T_lantern, lantern_branches);
  if (!flag_data){
    enable_branches(trees.T_eval, eval_truth_branches);
    enable_branches(trees.T_PFeval, pfeval_truth_branches);
    enable_branches(trees.T_glee, glee_truth_branches);
  }
}

#endif
