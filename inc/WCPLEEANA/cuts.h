#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TString.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"
#include "space.h"
#include "pandora.h"
#include "lantern.h"
#include "glee.h"
#include "particle.h"

#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

namespace LEEana{
  // this is for the real data, for fake data this should be 1 ...
  double em_charge_scale = 0.95;
  //double em_charge_scale = 1.0;

  bool wayToSort(double i, double j) { return i > j; };

  double get_mass_GeV(int pdg);
  double get_mass_MeV(int pdg);

  //truth==1 for truth, 0 for WC, else for LArPID, mother_check=-999 for all particles, mother_check=-9999 for all non-primary, which_p=0 for sum of all particles
  double get_KE(PFevalInfo& pfeval, int pdg, int truth, int mother_check, int which_p, double threshold);

  bool check_is_FC(double x, double y, double z);
  bool is_pfeval_muon(PFevalInfo& pfeval,int index, double tolerance);
  //returns (flag_FC_lepton, flag_FC_hadron); flag_FC_lepton is true if the primary muon energy is taken from range, false if from MCS
  //method=0: original choice (FC event, or muon start/end inside check_is_FC -> range)
  //method=2: two-sided, FC events with the muon inside check_is_FC -> range, others -> range if |range-MCS|/MCS<threshold, else MCS
  //method=3: one-sided, FC events with the muon inside check_is_FC -> range, others -> MCS only if (MCS-range)/MCS>threshold
  std::tuple<bool,bool> get_part_is_FC(PFevalInfo& pfeval,EvalInfo& eval, int method=0, double threshold=0.05);
  bool check_muon_range_MCS(PFevalInfo& pfeval, int method, double threshold);

  double get_muon_energy_new(PFevalInfo& pfeval, bool flag_FC_lepton, bool return_KE, bool return_MeV);
  TVector3 get_muon_momentum_new(PFevalInfo& pfeval, EvalInfo& eval);
  //drop_muon_showers: when the muon energy comes from MCS, remove particles that continue the muon past its reco end (start within drop_dist cm of the muon end, cos(angle to the muon direction)>drop_cos)
  //drop_mass: also remove the change of the masses in kine_reco_add_energy from removing those particles (false for Eavail, which removes all of kine_reco_add_energy itself)
  double get_kine_reco_Enu_new(PFevalInfo& pfeval, KineInfo& kine, SpaceInfo& space, bool flag_FC_lepton, bool flag_data, bool correct_protons, bool drop_muon_showers=false, double drop_dist=15, double drop_cos=0.9, bool drop_mass=true);
  double get_muon_continuation_energy(PFevalInfo& pfeval, KineInfo& kine, SpaceInfo& space, bool flag_data, double drop_dist, double drop_cos, double& E_mass);
  int get_reco_mother_index(PFevalInfo& pfeval, int index, const std::vector<bool>& removed);
  double get_kine_chain_end_mass(PFevalInfo& pfeval, const std::vector<int>& indices, const std::vector<bool>& removed);
  void get_principal_axis(const std::vector<TVector3>& pts, TVector3& centroid, TVector3& axis);

  int get_particle_0pNp_bdt_bin(PFevalInfo& pfeval, TaggerInfo& tagger, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern, double KE_threshold, double KE_pl_threshold, double scat_bdt_threshold, double vtxact_bdt_threshold);

  std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> get_range_proton_KE(PFevalInfo& pfeval, SpaceInfo& space, bool return_MeV);
  double get_range_proton_KE_particle(PFevalInfo& pfeval, SpaceInfo& space, int i);
  int get_reco_leading_proton(PFevalInfo& pfeval, SpaceInfo& space, double KE_threshold, double& KE_lead, int& n_protons);
  // charges of the primary muon's spacepoints, ordered from its start to its end as in create_particle (particle.h)
  std::vector<double> get_reco_muon_spacepoints_q(PFevalInfo& pfeval, SpaceInfo& space, double tolerance=0.0001);

  std::vector<double> get_pandora_proton_KE(PandoraInfo& pandora, double TRACK_SCORE_CUT, bool return_MeV);
  std::vector<double> get_lantern_KE(LanternInfo& lantern, int pdg, double vtx_cut, bool return_MeV);

  // correct reco neutrino energy and reco shower energy
  double get_reco_Enu_corr(KineInfo& kine, bool flag_data);

  double get_reco_Eproton(KineInfo& kine);

  double get_true_Eavail(EvalInfo& eval, PFevalInfo& pfeval, bool useThreshold, bool useMass);
  double get_truth_p_mu_cos(PFevalInfo& pfeval);

  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  double get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, TString var_name);

  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  // Event-level quantities of get_cut_pass (independent of the channel and add_cut), see fill_cut_event_info
  struct CutEventInfo{
    double reco_Enu;
    std::map<std::string, bool> map_cuts_flag;
    bool flag_generic;
    bool flag_numuCC;
    bool flag_numuCC_tight;
    bool flag_numuCC_1mu0p;
    bool flag_numuCC_cutbased;
    bool flag_nueCC;
    bool flag_0p;
    bool flag_cc_pi0;
    bool flag_FC;
    bool flag_FC_lepton;           // muon energy from range (true) or MCS (false), as in kine_reco_Enu_new3_5
    bool flag_FC_hadron;           // hadronic system contained (get_part_is_FC)
    int costheta_bin;
    int Enu_bin;
    bool part_bin_set = false;     // cached get_particle_0pNp_bdt_bin result for this event
    double part_bin_thresholds[4];
    int part_bin;
  };
  void fill_cut_event_info(CutEventInfo& info, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, CutEventInfo& info, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  bool get_rw_cut_pass(TString cut, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  // rootino bug fix factor (1 unless GTruth_ResNum==9), used by get_weight and the systematics
  double get_rootino_ratio(TString weight_name);
  double get_rootino_weight(EvalInfo& eval, GleeInfo& glee, double rootino_pot_ratio);
  double get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, GleeInfo& glee, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool,  std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data=false);
  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  // binning helpers of the cross-section measurements (truth: get_xs_signal_no, reco: the "*_bin" variables of get_kine_var)
  TString get_xs_bin_name(TString prefix, double value, const std::vector<double>& edges);
  TString get_xs_2d_bin_name(TString prefix, TString slice_var, double slice_value, const std::vector<double>& slice_edges, TString var, double value, const std::vector<std::vector<double>>& bin_edges);
  TString get_xs_3d_bin_name(TString prefix, TString outer_var, double outer_value, const std::vector<double>& outer_edges, TString inner_var, double inner_value, const std::vector<std::vector<double>>& inner_edges, TString var, double value, const std::vector<std::vector<std::vector<double>>>& bin_edges);
  int get_bin_index(double value, const std::vector<double>& edges);
  int get_2d_bin_index(double slice_value, const std::vector<double>& slice_edges, double value, const std::vector<std::vector<double>>& bin_edges);
  int get_3d_bin_index(double outer_value, const std::vector<double>& outer_edges, double inner_value, const std::vector<std::vector<double>>& inner_edges, double value, const std::vector<std::vector<std::vector<double>>>& bin_edges);
  double get_muon_Etot_new(PFevalInfo& pfeval, EvalInfo& eval);
  TVector3 get_reco_proton_dir(PFevalInfo& pfeval, int index);
  double get_reco_cos_mu_p(PFevalInfo& pfeval, int index);

  // Slices of the multi-differential measurements, shared by the truth and the reco binnings (inner edges, as in
  // get_bin_index). The triple-differential inner slices are given for each outer slice.
  namespace xs{
    const std::vector<double> costheta_slices = {0, 0.3, 0.5, 0.7, 0.8, 0.9};      // 12 (muon costheta), 13 (proton costheta)
    const std::vector<double> pl_slices = {0, 150, 300, 450, 625, 850, 1200};             // 16
    const std::vector<double> Eavail_slices = {75, 150, 250, 375, 550, 800};                  // 17
    const std::vector<double> Kp_slices_mup = {95, 145, 220, 345};                          // 18
    const std::vector<double> Eavail_slices_mup = {150, 250, 375, 550, 800};                  // 19
    const std::vector<double> Emu_slices_mup = {350, 550, 700, 900, 1150};                 // 20
    const std::vector<double> Eavail_outer_0p = {75, 150};                                    // 21, 22 (0p)
    const std::vector<double> Eavail_outer_Np = {150, 250, 375, 550};                         // 21, 22 (Np)
    const std::vector<std::vector<double>> costheta_inner_0p(3, {0, 0.5, 0.7, 0.9});           // 21
    const std::vector<std::vector<double>> costheta_inner_Np(5, {0, 0.3, 0.5, 0.7, 0.8, 0.9});
    const std::vector<std::vector<double>> pl_inner_0p(3, {150, 300, 450, 625, 850, 1200});     // 22
    const std::vector<std::vector<double>> pl_inner_Np(5, {0, 150, 300, 450, 625, 850, 1200});
    const std::vector<double> costheta_outer_p = {0.5, 0.7, 0.8, 0.9};                         // 23, 24
    const std::vector<std::vector<double>> Emu_inner_p = {
      {350, 450},   // costheta <= 0.5
      {350, 450, 700},   // costheta 0.5 - 0.7
      {450, 550, 900},   // costheta 0.7 - 0.8
      {450, 550, 900},   // costheta 0.8 - 0.9
      {550, 700, 900, 1150, 1600}   // costheta > 0.9
    };
    const std::vector<double> pl_outer_p = {150, 300, 450, 625, 850};                        // 25, 26
    const std::vector<std::vector<double>> pt_inner_p = {
      {150, 200, 250, 350},   // pl <= 150
      {200, 300, 450},   // pl 150 - 300
      {200, 300, 450},   // pl 300 - 450
      {250, 450},   // pl 450 - 625
      {250, 450},   // pl 625 - 850
      {250, 550}   // pl > 850
    };
    // truth bins of the multi-differential measurements in each slice (get_xs_signal_no); also the reco bins of the
    // "_bin_t" variables of get_kine_var (reco values in the truth binning)
    const std::vector<std::vector<double>> truth_Emu_bins_costheta = {
      {250, 300},   // costheta <= 0
      {250, 300},   // costheta 0 - 0.3
      {300, 350, 450},   // costheta 0.3 - 0.5
      {300, 350, 450, 700},   // costheta 0.5 - 0.7
      {300, 450, 550, 900},   // costheta 0.7 - 0.8
      {350, 450, 550, 900},   // costheta 0.8 - 0.9
      {450, 550, 700, 900, 1150, 1600}   // costheta > 0.9
    };
    const std::vector<std::vector<double>> truth_Kp_bins_costhetap = {
      {70, 95, 120, 145},   // costhetap <= 0
      {70, 95, 120, 145, 170, 220},   // costhetap 0 - 0.3
      {70, 95, 120, 145, 170, 220, 270},   // costhetap 0.3 - 0.5
      {70, 95, 120, 145, 170, 220, 270, 345},   // costhetap 0.5 - 0.7
      {70, 95, 120, 145, 170, 220, 270, 345, 470},   // costhetap 0.7 - 0.8
      {70, 95, 120, 145, 170, 220, 270, 345, 470},   // costhetap 0.8 - 0.9
      {70, 95, 120, 145, 170, 220, 270, 345}   // costhetap > 0.9
    };
    const std::vector<std::vector<double>> truth_pt_bins_pl = {
      {100, 200, 250},   // pl <= 0
      {150, 200, 250, 350},   // pl 0 - 150
      {150, 200, 300, 450},   // pl 150 - 300
      {150, 250, 350},   // pl 300 - 450
      {200, 350, 550},   // pl 450 - 625
      {200, 350, 550},   // pl 625 - 850
      {250, 450},   // pl 850 - 1200
      {350, 800}   // pl > 1200
    };
    const std::vector<std::vector<double>> truth_Emu_bins_Eavail = {
      {350, 550, 700, 900, 1150},   // Eavail <= 75
      {300, 350, 450, 550, 700, 900, 1600},   // Eavail 75 - 150
      {300, 450, 550, 900},   // Eavail 150 - 250
      {300, 450, 700, 1150},   // Eavail 250 - 375
      {250, 350, 450, 700, 1150},   // Eavail 375 - 550
      {250, 350, 550, 900},   // Eavail 550 - 800
      {300, 450, 700, 1150}   // Eavail > 800
    };
    const std::vector<std::vector<double>> truth_costhetamup_bins_Kp = {
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Kp <= 95
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Kp 95 - 145
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Kp 145 - 220
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Kp 220 - 345
      {-0.6, -0.2, 0, 0.2, 0.6}   // Kp > 345
    };
    const std::vector<std::vector<double>> truth_costhetamup_bins_Eavail = {
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Eavail <= 150
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Eavail 150 - 250
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Eavail 250 - 375
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Eavail 375 - 550
      {-0.4, 0, 0.4},   // Eavail 550 - 800
      {-0.2, 0.2, 0.6}   // Eavail > 800
    };
    const std::vector<std::vector<double>> truth_costhetamup_bins_Emu = {
      {-0.6, -0.4, -0.2, 0.2, 0.6},   // Emu <= 350
      {-0.6, -0.4, -0.2, 0, 0.2, 0.6},   // Emu 350 - 550
      {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Emu 550 - 700
      {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Emu 700 - 900
      {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // Emu 900 - 1150
      {-0.4, -0.2, 0, 0.2, 0.4, 0.6}   // Emu > 1150
    };
    const std::vector<std::vector<std::vector<double>>> truth_Emu_bins_3d_0p = {
      {
       {300},   // Eavail <= 75, costheta <= 0
       {300, 450},   // Eavail <= 75, costheta 0 - 0.5
       {350, 450, 700},   // Eavail <= 75, costheta 0.5 - 0.7
       {350, 450, 550, 700, 900, 1150},   // Eavail <= 75, costheta 0.7 - 0.9
       {450, 550, 700, 900, 1150, 1600}   // Eavail <= 75, costheta > 0.9
      },
      {
       {},   // Eavail 75 - 150, costheta <= 0
       {350},   // Eavail 75 - 150, costheta 0 - 0.5
       {550},   // Eavail 75 - 150, costheta 0.5 - 0.7
       {450, 700, 900},   // Eavail 75 - 150, costheta 0.7 - 0.9
       {700, 900, 1150}   // Eavail 75 - 150, costheta > 0.9
      },
      {
       {300},   // Eavail > 150, costheta <= 0
       {250, 350},   // Eavail > 150, costheta 0 - 0.5
       {300, 450},   // Eavail > 150, costheta 0.5 - 0.7
       {350, 450, 550, 700, 1150},   // Eavail > 150, costheta 0.7 - 0.9
       {450, 550, 700, 900, 1600}   // Eavail > 150, costheta > 0.9
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_Emu_bins_3d_Np = {
      {
       {},   // Eavail <= 150, costheta <= 0
       {},   // Eavail <= 150, costheta 0 - 0.3
       {},   // Eavail <= 150, costheta 0.3 - 0.5
       {450},   // Eavail <= 150, costheta 0.5 - 0.7
       {550, 700},   // Eavail <= 150, costheta 0.7 - 0.8
       {550, 700, 900},   // Eavail <= 150, costheta 0.8 - 0.9
       {700, 900, 1150, 1600}   // Eavail <= 150, costheta > 0.9
      },
      {
       {300},   // Eavail 150 - 250, costheta <= 0
       {},   // Eavail 150 - 250, costheta 0 - 0.3
       {450},   // Eavail 150 - 250, costheta 0.3 - 0.5
       {450, 550},   // Eavail 150 - 250, costheta 0.5 - 0.7
       {550},   // Eavail 150 - 250, costheta 0.7 - 0.8
       {550, 900},   // Eavail 150 - 250, costheta 0.8 - 0.9
       {700, 900, 1150}   // Eavail 150 - 250, costheta > 0.9
      },
      {
       {300},   // Eavail 250 - 375, costheta <= 0
       {},   // Eavail 250 - 375, costheta 0 - 0.3
       {450},   // Eavail 250 - 375, costheta 0.3 - 0.5
       {450, 700},   // Eavail 250 - 375, costheta 0.5 - 0.7
       {550, 700},   // Eavail 250 - 375, costheta 0.7 - 0.8
       {550, 700, 900},   // Eavail 250 - 375, costheta 0.8 - 0.9
       {700, 900, 1600}   // Eavail 250 - 375, costheta > 0.9
      },
      {
       {300},   // Eavail 375 - 550, costheta <= 0
       {450},   // Eavail 375 - 550, costheta 0 - 0.3
       {450},   // Eavail 375 - 550, costheta 0.3 - 0.5
       {450, 700},   // Eavail 375 - 550, costheta 0.5 - 0.7
       {700},   // Eavail 375 - 550, costheta 0.7 - 0.8
       {700},   // Eavail 375 - 550, costheta 0.8 - 0.9
       {900}   // Eavail 375 - 550, costheta > 0.9
      },
      {
       {250, 350},   // Eavail > 550, costheta <= 0
       {350},   // Eavail > 550, costheta 0 - 0.3
       {450},   // Eavail > 550, costheta 0.3 - 0.5
       {450, 700},   // Eavail > 550, costheta 0.5 - 0.7
       {550},   // Eavail > 550, costheta 0.7 - 0.8
       {700},   // Eavail > 550, costheta 0.8 - 0.9
       {900}   // Eavail > 550, costheta > 0.9
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_pt_bins_3d_0p = {
      {
       {200, 250, 300},   // Eavail <= 75, pl <= 150
       {200, 300},   // Eavail <= 75, pl 150 - 300
       {150, 250, 350},   // Eavail <= 75, pl 300 - 450
       {150, 300, 450},   // Eavail <= 75, pl 450 - 625
       {200, 350, 550},   // Eavail <= 75, pl 625 - 850
       {250, 550},   // Eavail <= 75, pl 850 - 1200
       {300}   // Eavail <= 75, pl > 1200
      },
      {
       {250},   // Eavail 75 - 150, pl <= 150
       {250},   // Eavail 75 - 150, pl 150 - 300
       {300},   // Eavail 75 - 150, pl 300 - 450
       {250, 350},   // Eavail 75 - 150, pl 450 - 625
       {300, 450},   // Eavail 75 - 150, pl 625 - 850
       {300},   // Eavail 75 - 150, pl 850 - 1200
       {450}   // Eavail 75 - 150, pl > 1200
      },
      {
       {150, 250, 350},   // Eavail > 150, pl <= 150
       {150, 200, 250, 350},   // Eavail > 150, pl 150 - 300
       {150, 200, 300, 450},   // Eavail > 150, pl 300 - 450
       {200, 300},   // Eavail > 150, pl 450 - 625
       {200, 350},   // Eavail > 150, pl 625 - 850
       {250, 450},   // Eavail > 150, pl 850 - 1200
       {350, 650}   // Eavail > 150, pl > 1200
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_pt_bins_3d_Np = {
      {
       {},   // Eavail <= 150, pl <= 0
       {300},   // Eavail <= 150, pl 0 - 150
       {250, 350},   // Eavail <= 150, pl 150 - 300
       {300},   // Eavail <= 150, pl 300 - 450
       {300, 450},   // Eavail <= 150, pl 450 - 625
       {300},   // Eavail <= 150, pl 625 - 850
       {300},   // Eavail <= 150, pl 850 - 1200
       {350}   // Eavail <= 150, pl > 1200
      },
      {
       {250},   // Eavail 150 - 250, pl <= 0
       {300},   // Eavail 150 - 250, pl 0 - 150
       {300, 450},   // Eavail 150 - 250, pl 150 - 300
       {300, 450},   // Eavail 150 - 250, pl 300 - 450
       {300},   // Eavail 150 - 250, pl 450 - 625
       {300, 450},   // Eavail 150 - 250, pl 625 - 850
       {300, 550},   // Eavail 150 - 250, pl 850 - 1200
       {350}   // Eavail 150 - 250, pl > 1200
      },
      {
       {250},   // Eavail 250 - 375, pl <= 0
       {250, 350},   // Eavail 250 - 375, pl 0 - 150
       {250, 350},   // Eavail 250 - 375, pl 150 - 300
       {300, 450},   // Eavail 250 - 375, pl 300 - 450
       {300, 450},   // Eavail 250 - 375, pl 450 - 625
       {300, 450},   // Eavail 250 - 375, pl 625 - 850
       {350, 550},   // Eavail 250 - 375, pl 850 - 1200
       {450}   // Eavail 250 - 375, pl > 1200
      },
      {
       {200, 300},   // Eavail 375 - 550, pl <= 0
       {250},   // Eavail 375 - 550, pl 0 - 150
       {250, 450},   // Eavail 375 - 550, pl 150 - 300
       {350},   // Eavail 375 - 550, pl 300 - 450
       {350, 550},   // Eavail 375 - 550, pl 450 - 625
       {450},   // Eavail 375 - 550, pl 625 - 850
       {550},   // Eavail 375 - 550, pl 850 - 1200
       {}   // Eavail 375 - 550, pl > 1200
      },
      {
       {150, 250},   // Eavail > 550, pl <= 0
       {200, 300, 450},   // Eavail > 550, pl 0 - 150
       {250, 450},   // Eavail > 550, pl 150 - 300
       {350},   // Eavail > 550, pl 300 - 450
       {450},   // Eavail > 550, pl 450 - 625
       {},   // Eavail > 550, pl 625 - 850
       {},   // Eavail > 550, pl 850 - 1200
       {650}   // Eavail > 550, pl > 1200
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_Kp_bins_costheta_Emu = {
      {
       {70, 95, 120, 145, 170, 220, 270, 345},   // costheta <= 0.5, Emu <= 350
       {95, 145, 220, 270, 345, 470},   // costheta <= 0.5, Emu 350 - 450
       {95, 145, 220, 270, 345, 470}   // costheta <= 0.5, Emu > 450
      },
      {
       {120},   // costheta 0.5 - 0.7, Emu <= 350
       {120},   // costheta 0.5 - 0.7, Emu 350 - 450
       {70, 120, 170, 220, 270, 345},   // costheta 0.5 - 0.7, Emu 450 - 700
       {120, 220, 345, 470}   // costheta 0.5 - 0.7, Emu > 700
      },
      {
       {120},   // costheta 0.7 - 0.8, Emu <= 450
       {120},   // costheta 0.7 - 0.8, Emu 450 - 550
       {95, 145, 170, 220, 270},   // costheta 0.7 - 0.8, Emu 550 - 900
       {170, 345}   // costheta 0.7 - 0.8, Emu > 900
      },
      {
       {95, 170},   // costheta 0.8 - 0.9, Emu <= 450
       {120},   // costheta 0.8 - 0.9, Emu 450 - 550
       {70, 95, 120, 145, 170, 220, 270},   // costheta 0.8 - 0.9, Emu 550 - 900
       {95, 120, 145, 220, 270, 345}   // costheta 0.8 - 0.9, Emu > 900
      },
      {
       {120},   // costheta > 0.9, Emu <= 550
       {95, 170},   // costheta > 0.9, Emu 550 - 700
       {70, 95, 120, 170},   // costheta > 0.9, Emu 700 - 900
       {70, 95, 120, 145, 170, 220},   // costheta > 0.9, Emu 900 - 1150
       {70, 95, 120, 145, 220, 270},   // costheta > 0.9, Emu 1150 - 1600
       {95, 145, 220, 345}   // costheta > 0.9, Emu > 1600
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_costhetamup_bins_costheta_Emu = {
      {
       {-0.6, -0.4, -0.2, 0.2},   // costheta <= 0.5, Emu <= 350
       {-0.6, -0.4, -0.2, 0},   // costheta <= 0.5, Emu 350 - 450
       {-0.6, -0.4, -0.2, 0}   // costheta <= 0.5, Emu > 450
      },
      {
       {0.2},   // costheta 0.5 - 0.7, Emu <= 350
       {0},   // costheta 0.5 - 0.7, Emu 350 - 450
       {-0.4, -0.2, 0, 0.2},   // costheta 0.5 - 0.7, Emu 450 - 700
       {-0.2, 0, 0.2}   // costheta 0.5 - 0.7, Emu > 700
      },
      {
       {0.2},   // costheta 0.7 - 0.8, Emu <= 450
       {0.2},   // costheta 0.7 - 0.8, Emu 450 - 550
       {-0.2, 0, 0.2, 0.4, 0.6},   // costheta 0.7 - 0.8, Emu 550 - 900
       {0, 0.2}   // costheta 0.7 - 0.8, Emu > 900
      },
      {
       {0.2, 0.6},   // costheta 0.8 - 0.9, Emu <= 450
       {0.2},   // costheta 0.8 - 0.9, Emu 450 - 550
       {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // costheta 0.8 - 0.9, Emu 550 - 900
       {-0.2, 0, 0.2, 0.4}   // costheta 0.8 - 0.9, Emu > 900
      },
      {
       {0.4},   // costheta > 0.9, Emu <= 550
       {0.2},   // costheta > 0.9, Emu 550 - 700
       {0, 0.2, 0.4, 0.6},   // costheta > 0.9, Emu 700 - 900
       {-0.2, 0.2, 0.4, 0.6},   // costheta > 0.9, Emu 900 - 1150
       {-0.2, 0, 0.2, 0.4, 0.6},   // costheta > 0.9, Emu 1150 - 1600
       {-0.2, 0.2, 0.4, 0.6}   // costheta > 0.9, Emu > 1600
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_Kp_bins_pl_pt = {
      {
       {120, 220, 345},   // pl <= 150, pt <= 150
       {120, 220},   // pl <= 150, pt 150 - 200
       {120, 220, 345},   // pl <= 150, pt 200 - 250
       {95, 145, 220, 270, 345},   // pl <= 150, pt 250 - 350
       {95, 145, 220, 270, 345, 470}   // pl <= 150, pt > 350
      },
      {
       {120},   // pl 150 - 300, pt <= 200
       {95, 145, 220},   // pl 150 - 300, pt 200 - 300
       {95, 145, 220},   // pl 150 - 300, pt 300 - 450
       {145, 270, 470}   // pl 150 - 300, pt > 450
      },
      {
       {120},   // pl 300 - 450, pt <= 200
       {95, 170},   // pl 300 - 450, pt 200 - 300
       {95, 145, 220},   // pl 300 - 450, pt 300 - 450
       {120, 220, 345}   // pl 300 - 450, pt > 450
      },
      {
       {95},   // pl 450 - 625, pt <= 250
       {70, 95, 120, 145, 170, 220},   // pl 450 - 625, pt 250 - 450
       {95, 145, 220, 345}   // pl 450 - 625, pt > 450
      },
      {
       {95},   // pl 625 - 850, pt <= 250
       {70, 95, 120, 145, 170, 220},   // pl 625 - 850, pt 250 - 450
       {95, 145, 220, 270, 345}   // pl 625 - 850, pt > 450
      },
      {
       {95},   // pl > 850, pt <= 250
       {70, 95, 120, 145, 220, 270},   // pl > 850, pt 250 - 550
       {95, 145, 220, 270, 470}   // pl > 850, pt > 550
      }
    };
    const std::vector<std::vector<std::vector<double>>> truth_costhetamup_bins_pl_pt = {
      {
       {-0.4},   // pl <= 150, pt <= 150
       {-0.6, -0.2},   // pl <= 150, pt 150 - 200
       {-0.6, -0.2, 0.2},   // pl <= 150, pt 200 - 250
       {-0.6, -0.4, -0.2, 0},   // pl <= 150, pt 250 - 350
       {-0.6, -0.4, -0.2, 0}   // pl <= 150, pt > 350
      },
      {
       {0.4},   // pl 150 - 300, pt <= 200
       {0, 0.4},   // pl 150 - 300, pt 200 - 300
       {-0.4, -0.2, 0, 0.2},   // pl 150 - 300, pt 300 - 450
       {-0.2, 0}   // pl 150 - 300, pt > 450
      },
      {
       {},   // pl 300 - 450, pt <= 200
       {0, 0.4},   // pl 300 - 450, pt 200 - 300
       {-0.2, 0, 0.2, 0.4},   // pl 300 - 450, pt 300 - 450
       {-0.2, 0, 0.2}   // pl 300 - 450, pt > 450
      },
      {
       {0.4},   // pl 450 - 625, pt <= 250
       {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // pl 450 - 625, pt 250 - 450
       {-0.2, 0, 0.2, 0.4}   // pl 450 - 625, pt > 450
      },
      {
       {0.2},   // pl 625 - 850, pt <= 250
       {-0.2, 0, 0.2, 0.4, 0.6},   // pl 625 - 850, pt 250 - 450
       {-0.2, 0, 0.2, 0.4}   // pl 625 - 850, pt > 450
      },
      {
       {0.2, 0.6},   // pl > 850, pt <= 250
       {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // pl > 850, pt 250 - 550
       {-0.2, 0, 0.2, 0.4, 0.6}   // pl > 850, pt > 550
      }
    };
  }

  // generic neutrino cuts
  // TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
  bool is_generic(EvalInfo& info);

  // preselection cuts
  // TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
  bool is_preselection(EvalInfo& info);

  // nueCC cuts
  // TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
  bool is_nueCC(TaggerInfo& tagger_info);
  bool is_loosenueCC(TaggerInfo& tagger_info);

  bool is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data);

  // numuCC cuts
  // TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
  bool is_numuCC(TaggerInfo& tagger_info);
  bool is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_numuCC_cutbased(TaggerInfo& tagger_info);

  // pio cuts (with and without vertex)
  // TCut pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_pi0(KineInfo& kine, bool flag_data);

  // must be with vertex ...
  // TCut cc_pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_cc_pi0(KineInfo& kine, bool flag_data);


  // NC cuts
  // TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
  bool is_NC(TaggerInfo& tagger_info);
  bool is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);


  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";
  bool is_FC(EvalInfo& eval);


  bool is_true_0p(PFevalInfo& pfeval,double threshold);

  const std::vector<double>& get_proton_length_bins(){
    static const std::vector<double> proton_length_bins = {7.65721e-06, 0.003103, 0.0091284, 0.0176375, 0.0283912, 0.0412807, 0.0562127, 0.0731138, 0.0919263, 0.112603, 0.135101, 0.159272, 0.185117, 0.212807, 0.24211, 0.273135, 0.305831, 0.340111, 0.376083, 0.413583, 0.452689, 0.493293, 0.535259, 0.57868, 0.62366, 0.670316, 0.718547, 0.768222, 0.819402, 0.871993, 0.926046, 0.981485, 1.03822, 1.09633, 1.15586, 1.21689, 1.27932, 1.34303, 1.40806, 1.47449, 1.54236, 1.61158, 1.68204, 1.75379, 1.82687, 1.90134, 1.97711, 2.05409, 2.13233, 2.21187, 2.29275, 2.37489, 2.45822, 2.54277, 2.62858, 2.71569, 2.80402, 2.89352, 2.98421, 3.07612, 3.16929, 3.26367, 3.35917, 3.45585, 3.55372, 3.65281, 3.75307, 3.85444, 3.95696, 4.06063, 4.1655, 4.2715, 4.37859, 4.48679, 4.59613, 4.70662, 4.81823, 4.9309, 5.04466, 5.15953, 5.27553, 5.39262, 5.51076, 5.62996, 5.75025, 5.87164, 5.9941, 6.11758, 6.24211, 6.3677, 6.49437, 6.62208, 6.75079, 6.88053, 7.0113, 7.14313, 7.27598, 7.40982, 7.54466, 7.68052, 7.81742, 7.95524, 8.09388, 8.23334, 8.37365, 8.51481, 8.65683, 8.79972, 8.94348, 9.08814, 9.23371, 9.38019, 9.52759, 9.67593, 9.82522, 9.97548, 10.1267, 10.2789, 10.4321, 10.5864, 10.7416, 10.8979, 11.0553, 11.2137, 11.3732, 11.5338, 11.6953, 11.8576, 12.0207, 12.1846, 12.3493, 12.5147, 12.681, 12.8481, 13.016, 13.1847, 13.3543, 13.5248, 13.696, 13.8682, 14.0412, 14.2152, 14.39, 14.5657, 14.7423, 14.9198, 15.0983, 15.2777, 15.4581, 15.6394, 15.8217, 16.0049, 16.1887, 16.3733, 16.5586, 16.7447, 16.9315, 17.119, 17.3073, 17.4963, 17.6862, 17.8767, 18.0681, 18.2603, 18.4532, 18.6469, 18.8415, 19.0368, 19.233, 19.43, 19.6278, 19.8265, 20.026, 20.2264, 20.4276, 20.6297, 20.8325, 21.0361, 21.2403, 21.4452, 21.6507, 21.857, 22.064, 22.2716, 22.48, 22.689, 22.8988, 23.1093, 23.3205, 23.5325, 23.7451, 23.9585, 24.1727, 24.3876, 24.6033, 24.8197, 25.0368, 25.2548, 25.4735, 25.693, 25.9133, 26.1343, 26.3559, 26.5782, 26.801, 27.0245, 27.2487, 27.4735, 27.6989, 27.925, 28.1517, 28.3791, 28.6072, 28.8359, 29.0653, 29.2953, 29.5261, 29.7575, 29.9896, 30.2223, 30.4558, 30.69, 30.9249, 31.1604, 31.3967, 31.6337, 31.8714, 32.1096, 32.3484, 32.5878, 32.8278, 33.0684, 33.3095, 33.5513, 33.7937, 34.0367, 34.2803, 34.5245, 34.7693, 35.0147, 35.2608, 35.5075, 35.7548, 36.0027, 36.2513, 36.5005, 36.7503, 37.0008, 37.2519, 37.5037, 37.7562, 38.0092, 38.2628, 38.5169, 38.7716, 39.0268, 39.2825, 39.5388, 39.7957, 40.0531, 40.3111, 40.5697, 40.8288, 41.0885, 41.3487, 41.6096, 41.871, 42.1329, 42.3955, 42.6586, 42.9224, 43.1867, 43.4516, 43.7171, 43.9832, 44.2498, 44.5171, 44.7848, 45.0531, 45.3218, 45.5911, 45.8609, 46.1312, 46.402, 46.6733, 46.9451, 47.2175, 47.4904, 47.7638, 48.0377, 48.3122, 48.5872, 48.8627, 49.1388, 49.4154, 49.6925, 49.9702, 50.2485, 50.5272, 50.8065, 51.0864, 51.3668, 51.6476, 51.9288, 52.2105, 52.4926, 52.7752, 53.0582, 53.3417, 53.6256, 53.9099, 54.1947, 54.48, 54.7657, 55.0518, 55.3384, 55.6255, 55.913, 56.201, 56.4895, 56.7784, 57.0678, 57.3576, 57.6479, 57.9387, 58.23, 58.5217, 58.8139, 59.1066, 59.3997, 59.6934, 59.9875, 60.2821, 60.5771, 60.8727, 61.1687, 61.4653, 61.7623, 62.0598, 62.3578, 62.6563, 62.9553, 63.2548, 63.5548, 63.8553, 64.1563, 64.4578, 64.7598, 65.0623, 65.3653, 65.6688, 65.9728, 66.2772, 66.5819, 66.8871, 67.1926, 67.4985, 67.8048, 68.1114, 68.4185, 68.726, 69.0338, 69.342, 69.6507, 69.9597, 70.2691, 70.5789, 70.8891, 71.1998, 71.5108, 71.8222, 72.134, 72.4462, 72.7588, 73.0719, 73.3853, 73.6991, 74.0134, 74.328, 74.6431, 74.9585, 75.2744, 75.5907, 75.9074, 76.2246, 76.5421, 76.8601, 77.1784, 77.4972, 77.8164, 78.1361, 78.4561, 78.7766, 79.0975, 79.4189, 79.7406, 80.0628, 80.3854, 80.7085, 81.032, 81.3559, 81.6802, 82.0048, 82.3298, 82.6551, 82.9808, 83.3068, 83.6331, 83.9598, 84.2868, 84.6142, 84.9419, 85.27, 85.5984, 85.9271, 86.2562, 86.5857, 86.9155, 87.2456, 87.5761, 87.907, 88.2382, 88.5697, 88.9016, 89.2339, 89.5665, 89.8995, 90.2328, 90.5665, 90.9005, 91.2349, 91.5697, 91.9048, 92.2403, 92.5761, 92.9123, 93.2489, 93.5858, 93.9231, 94.2608, 94.5988, 94.9372, 95.276, 95.6151, 95.9546, 96.2945, 96.6347, 96.9753, 97.3163, 97.6577, 97.9994, 98.3415, 98.6839, 99.0265, 99.3695, 99.7128, 100.056, 100.4, 100.744, 101.089, 101.434, 101.779, 102.124, 102.47, 102.816, 103.162, 103.509, 103.856, 104.203, 104.55, 104.898, 105.246, 105.595, 105.943, 106.292, 106.642, 106.991, 107.341, 107.692, 108.042, 108.393, 108.744, 109.096, 109.448, 109.8, 110.153, 110.505, 110.858, 111.212, 111.566, 111.92, 112.274, 112.629, 112.984, 113.339, 113.695, 114.051, 114.407, 114.764, 115.121, 115.478, 115.836, 116.194, 116.552, 116.91, 117.269, 117.628, 117.987, 118.346, 118.706, 119.066, 119.426, 119.786, 120.147, 120.508, 120.869, 121.231, 121.593, 121.955, 122.317, 122.68, 123.043, 123.406, 123.769, 124.133, 124.497, 124.861, 125.225, 125.59, 125.955, 126.321, 126.686, 127.052, 127.418, 127.785, 128.151, 128.518, 128.885, 129.253, 129.621, 129.989, 130.357, 130.726, 131.094, 131.464, 131.833, 132.203, 132.573, 132.943, 133.314, 133.684, 134.055, 134.427, 134.798, 135.17, 135.542, 135.914, 136.287, 136.659, 137.032, 137.406, 137.779, 138.153, 138.526, 138.901, 139.275, 139.649, 140.024, 140.399, 140.775, 141.15, 141.526, 141.902, 142.278, 142.655, 143.031, 143.408, 143.785, 144.163, 144.54, 144.918, 145.296, 145.675, 146.053, 146.432, 146.811, 147.191, 147.57, 147.95, 148.33, 148.71, 149.091, 149.472, 149.853, 150.234, 150.616, 150.997, 151.379, 151.762, 152.144, 152.527, 152.91, 153.293, 153.676, 154.06, 154.444, 154.828, 155.212, 155.596, 155.981, 156.366, 156.751, 157.136, 157.521, 157.907, 158.293, 158.679, 159.065, 159.452, 159.838, 160.225, 160.612, 161.0, 161.387, 161.775, 162.163, 162.551, 162.939, 163.328, 163.717, 164.106, 164.495, 164.884, 165.274, 165.664, 166.054, 166.444, 166.835, 167.225, 167.616, 168.007, 168.399, 168.79, 169.182, 169.574, 169.966, 170.359, 170.751, 171.144, 171.537, 171.93, 172.324, 172.717, 173.111, 173.505, 173.899, 174.294, 174.688, 175.083, 175.478, 175.873, 176.268, 176.664, 177.06, 177.455, 177.851, 178.248, 178.644, 179.041, 179.437, 179.834, 180.231, 180.629, 181.026, 181.424, 181.822, 182.22, 182.618, 183.017, 183.415, 183.814, 184.213, 184.612, 185.011, 185.411, 185.811, 186.211, 186.611, 187.011, 187.412, 187.812, 188.213, 188.614, 189.015, 189.417, 189.818, 190.22, 190.622, 191.024, 191.427, 191.829, 192.232, 192.635, 193.038, 193.441, 193.845, 194.248, 194.652, 195.056, 195.46, 195.864, 196.268, 196.673, 197.077, 197.482, 197.887, 198.292, 198.698, 199.103, 199.509, 199.914, 200.32, 200.726, 201.133, 201.539, 201.946, 202.352, 202.759, 203.166, 203.574, 203.981, 204.389, 204.796, 205.204, 205.612, 206.02, 206.429, 206.837, 207.246, 207.655, 208.064, 208.473, 208.882, 209.292, 209.702, 210.111, 210.521, 210.932, 211.342, 211.752, 212.163, 212.574, 212.985, 213.396, 213.807, 214.219, 214.63, 215.042, 215.454, 215.865, 216.278, 216.69, 217.102, 217.515, 217.927, 218.34, 218.753, 219.166, 219.579, 219.993, 220.406, 220.82, 221.234, 221.648, 222.062, 222.476, 222.89, 223.305, 223.72, 224.134, 224.549, 224.964, 225.38, 225.795, 226.211, 226.626, 227.042, 227.458, 227.874, 228.29, 228.707, 229.123, 229.54, 229.957, 230.374, 230.791, 231.208, 231.626, 232.043, 232.461, 232.879, 233.297, 233.715, 234.133, 234.552, 234.97, 235.389, 235.807, 236.226, 236.645, 237.064, 237.484, 237.903, 238.323, 238.742, 239.162, 239.582, 240.002, 240.422, 240.842, 241.263, 241.683, 242.104, 242.525, 242.946, 243.367, 243.788, 244.209, 244.631, 245.052, 245.474, 245.896, 246.318, 246.74, 247.162, 247.585, 248.007, 248.43, 248.852, 249.275, 249.698, 250.121, 250.545, 250.968, 251.392, 251.815, 252.239, 252.663, 253.087, 253.511, 253.935, 254.36, 254.784, 255.209, 255.634, 256.059, 256.484, 256.909, 257.334, 257.759, 258.185, 258.611, 259.036, 259.462, 259.888, 260.314, 260.74, 261.166, 261.593, 262.019, 262.446, 262.873, 263.3, 263.727, 264.154, 264.581, 265.008, 265.436, 265.863, 266.291, 266.719, 267.147, 267.575, 268.003, 268.431, 268.86, 269.288, 269.717, 270.146, 270.574, 271.003, 271.432, 271.862, 272.291, 272.72, 273.15, 273.58, 274.01, 274.44, 274.87, 275.3, 275.73, 276.16, 276.591, 277.021, 277.452, 277.883, 278.314, 278.745, 279.176, 279.607, 280.038, 280.47, 280.901, 281.333, 281.764, 282.196, 282.628, 283.06, 283.492, 283.924, 284.356, 284.789, 285.221, 285.654, 286.086, 286.519, 286.952, 287.385, 287.818, 288.251, 288.684, 289.118, 289.551, 289.985, 290.418, 290.852, 291.286, 291.72, 292.154, 292.588, 293.022, 293.457, 293.891, 294.326, 294.76, 295.195, 295.63, 296.065, 296.5, 296.935, 297.37, 297.805, 298.241, 298.676, 299.112, 299.548, 299.983, 300.419, 300.855, 301.291, 301.728, 302.164, 302.6, 303.036, 303.473, 303.91, 304.346, 304.783, 305.22, 305.657, 306.094, 306.531, 306.968, 307.406, 307.843, 308.28, 308.718, 309.156, 309.594, 310.031, 310.469, 310.907, 311.346, 311.784, 312.222, 312.661, 313.099, 313.538, 313.976, 314.415, 314.854, 315.293, 315.732, 316.171, 316.61, 317.05, 317.489, 317.929, 318.368, 328.368, 338.368, 348.368, 358.368, 368.368, 378.368, 388.368, 398.368, 408.368, 418.368, 428.368, 438.368, 448.368, 458.368, 468.368, 478.368, 488.368, 498.368, 508.368, 518.3679999999999, 528.3679999999999, 538.3679999999999, 548.3679999999999, 558.3679999999999, 568.3679999999999, 578.3679999999999, 588.3679999999999, 598.3679999999999, 608.3679999999999, 618.3679999999999, 628.3679999999999, 638.3679999999999, 648.3679999999999, 658.3679999999999, 668.3679999999999, 678.3679999999999, 688.3679999999999, 698.3679999999999, 708.3679999999999, 718.3679999999999, 728.3679999999999, 738.3679999999999, 748.3679999999999, 758.3679999999999, 768.3679999999999, 778.3679999999999, 788.3679999999999, 798.3679999999999, 808.3679999999999, 818.3679999999999, 828.3679999999999, 838.3679999999999, 848.3679999999999, 858.3679999999999, 868.3679999999999, 878.3679999999999, 888.3679999999999, 898.3679999999999, 908.3679999999999, 918.3679999999999, 928.3679999999999, 938.3679999999999, 948.3679999999999, 958.3679999999999, 968.3679999999999, 978.3679999999999, 988.3679999999999, 998.3679999999999, 1008.3679999999999, 1018.3679999999999, 1028.368, 1038.368, 1048.368, 1058.368, 1068.368, 1078.368, 1088.368, 1098.368, 1108.368};
  return proton_length_bins;
}
  const std::vector<double>& get_proton_energy_bins(){ 
    static const std::vector<double> proton_energy_bins = {0.001, 1.001, 2.001, 3.001, 4.001, 5.001, 6.001, 7.001, 8.001, 9.001, 10.001, 11.001, 12.001, 13.001, 14.001, 15.001, 16.001, 17.001, 18.001, 19.001, 20.001, 21.001, 22.001, 23.001, 24.001, 25.001, 26.001, 27.001, 28.001, 29.001, 30.001, 31.001, 32.001, 33.001, 34.001, 35.001, 36.001, 37.001, 38.001, 39.001, 40.001, 41.001, 42.001, 43.001, 44.001, 45.001, 46.001, 47.001, 48.001, 49.001, 50.001, 51.001, 52.001, 53.001, 54.001, 55.001, 56.001, 57.001, 58.001, 59.001, 60.001, 61.001, 62.001, 63.001, 64.001, 65.001, 66.001, 67.001, 68.001, 69.001, 70.001, 71.001, 72.001, 73.001, 74.001, 75.001, 76.001, 77.001, 78.001, 79.001, 80.001, 81.001, 82.001, 83.001, 84.001, 85.001, 86.001, 87.001, 88.001, 89.001, 90.001, 91.001, 92.001, 93.001, 94.001, 95.001, 96.001, 97.001, 98.001, 99.001, 100.001, 101.001, 102.001, 103.001, 104.001, 105.001, 106.001, 107.001, 108.001, 109.001, 110.001, 111.001, 112.001, 113.001, 114.001, 115.001, 116.001, 117.001, 118.001, 119.001, 120.001, 121.001, 122.001, 123.001, 124.001, 125.001, 126.001, 127.001, 128.001, 129.001, 130.001, 131.001, 132.001, 133.001, 134.001, 135.001, 136.001, 137.001, 138.001, 139.001, 140.001, 141.001, 142.001, 143.001, 144.001, 145.001, 146.001, 147.001, 148.001, 149.001, 150.001, 151.001, 152.001, 153.001, 154.001, 155.001, 156.001, 157.001, 158.001, 159.001, 160.001, 161.001, 162.001, 163.001, 164.001, 165.001, 166.001, 167.001, 168.001, 169.001, 170.001, 171.001, 172.001, 173.001, 174.001, 175.001, 176.001, 177.001, 178.001, 179.001, 180.001, 181.001, 182.001, 183.001, 184.001, 185.001, 186.001, 187.001, 188.001, 189.001, 190.001, 191.001, 192.001, 193.001, 194.001, 195.001, 196.001, 197.001, 198.001, 199.001, 200.001, 201.001, 202.001, 203.001, 204.001, 205.001, 206.001, 207.001, 208.001, 209.001, 210.001, 211.001, 212.001, 213.001, 214.001, 215.001, 216.001, 217.001, 218.001, 219.001, 220.001, 221.001, 222.001, 223.001, 224.001, 225.001, 226.001, 227.001, 228.001, 229.001, 230.001, 231.001, 232.001, 233.001, 234.001, 235.001, 236.001, 237.001, 238.001, 239.001, 240.001, 241.001, 242.001, 243.001, 244.001, 245.001, 246.001, 247.001, 248.001, 249.001, 250.001, 251.001, 252.001, 253.001, 254.001, 255.001, 256.001, 257.001, 258.001, 259.001, 260.001, 261.001, 262.001, 263.001, 264.001, 265.001, 266.001, 267.001, 268.001, 269.001, 270.001, 271.001, 272.001, 273.001, 274.001, 275.001, 276.001, 277.001, 278.001, 279.001, 280.001, 281.001, 282.001, 283.001, 284.001, 285.001, 286.001, 287.001, 288.001, 289.001, 290.001, 291.001, 292.001, 293.001, 294.001, 295.001, 296.001, 297.001, 298.001, 299.001, 300.001, 301.001, 302.001, 303.001, 304.001, 305.001, 306.001, 307.001, 308.001, 309.001, 310.001, 311.001, 312.001, 313.001, 314.001, 315.001, 316.001, 317.001, 318.001, 319.001, 320.001, 321.001, 322.001, 323.001, 324.001, 325.001, 326.001, 327.001, 328.001, 329.001, 330.001, 331.001, 332.001, 333.001, 334.001, 335.001, 336.001, 337.001, 338.001, 339.001, 340.001, 341.001, 342.001, 343.001, 344.001, 345.001, 346.001, 347.001, 348.001, 349.001, 350.001, 351.001, 352.001, 353.001, 354.001, 355.001, 356.001, 357.001, 358.001, 359.001, 360.001, 361.001, 362.001, 363.001, 364.001, 365.001, 366.001, 367.001, 368.001, 369.001, 370.001, 371.001, 372.001, 373.001, 374.001, 375.001, 376.001, 377.001, 378.001, 379.001, 380.001, 381.001, 382.001, 383.001, 384.001, 385.001, 386.001, 387.001, 388.001, 389.001, 390.001, 391.001, 392.001, 393.001, 394.001, 395.001, 396.001, 397.001, 398.001, 399.001, 400.001, 401.001, 402.001, 403.001, 404.001, 405.001, 406.001, 407.001, 408.001, 409.001, 410.001, 411.001, 412.001, 413.001, 414.001, 415.001, 416.001, 417.001, 418.001, 419.001, 420.001, 421.001, 422.001, 423.001, 424.001, 425.001, 426.001, 427.001, 428.001, 429.001, 430.001, 431.001, 432.001, 433.001, 434.001, 435.001, 436.001, 437.001, 438.001, 439.001, 440.001, 441.001, 442.001, 443.001, 444.001, 445.001, 446.001, 447.001, 448.001, 449.001, 450.001, 451.001, 452.001, 453.001, 454.001, 455.001, 456.001, 457.001, 458.001, 459.001, 460.001, 461.001, 462.001, 463.001, 464.001, 465.001, 466.001, 467.001, 468.001, 469.001, 470.001, 471.001, 472.001, 473.001, 474.001, 475.001, 476.001, 477.001, 478.001, 479.001, 480.001, 481.001, 482.001, 483.001, 484.001, 485.001, 486.001, 487.001, 488.001, 489.001, 490.001, 491.001, 492.001, 493.001, 494.001, 495.001, 496.001, 497.001, 498.001, 499.001, 500.001, 501.001, 502.001, 503.001, 504.001, 505.001, 506.001, 507.001, 508.001, 509.001, 510.001, 511.001, 512.001, 513.001, 514.001, 515.001, 516.001, 517.001, 518.001, 519.001, 520.001, 521.001, 522.001, 523.001, 524.001, 525.001, 526.001, 527.001, 528.001, 529.001, 530.001, 531.001, 532.001, 533.001, 534.001, 535.001, 536.001, 537.001, 538.001, 539.001, 540.001, 541.001, 542.001, 543.001, 544.001, 545.001, 546.001, 547.001, 548.001, 549.001, 550.001, 551.001, 552.001, 553.001, 554.001, 555.001, 556.001, 557.001, 558.001, 559.001, 560.001, 561.001, 562.001, 563.001, 564.001, 565.001, 566.001, 567.001, 568.001, 569.001, 570.001, 571.001, 572.001, 573.001, 574.001, 575.001, 576.001, 577.001, 578.001, 579.001, 580.001, 581.001, 582.001, 583.001, 584.001, 585.001, 586.001, 587.001, 588.001, 589.001, 590.001, 591.001, 592.001, 593.001, 594.001, 595.001, 596.001, 597.001, 598.001, 599.001, 600.001, 601.001, 602.001, 603.001, 604.001, 605.001, 606.001, 607.001, 608.001, 609.001, 610.001, 611.001, 612.001, 613.001, 614.001, 615.001, 616.001, 617.001, 618.001, 619.001, 620.001, 621.001, 622.001, 623.001, 624.001, 625.001, 626.001, 627.001, 628.001, 629.001, 630.001, 631.001, 632.001, 633.001, 634.001, 635.001, 636.001, 637.001, 638.001, 639.001, 640.001, 641.001, 642.001, 643.001, 644.001, 645.001, 646.001, 647.001, 648.001, 649.001, 650.001, 651.001, 652.001, 653.001, 654.001, 655.001, 656.001, 657.001, 658.001, 659.001, 660.001, 661.001, 662.001, 663.001, 664.001, 665.001, 666.001, 667.001, 668.001, 669.001, 670.001, 671.001, 672.001, 673.001, 674.001, 675.001, 676.001, 677.001, 678.001, 679.001, 680.001, 681.001, 682.001, 683.001, 684.001, 685.001, 686.001, 687.001, 688.001, 689.001, 690.001, 691.001, 692.001, 693.001, 694.001, 695.001, 696.001, 697.001, 698.001, 699.001, 700.001, 701.001, 702.001, 703.001, 704.001, 705.001, 706.001, 707.001, 708.001, 709.001, 710.001, 711.001, 712.001, 713.001, 714.001, 715.001, 716.001, 717.001, 718.001, 719.001, 720.001, 721.001, 722.001, 723.001, 724.001, 725.001, 726.001, 727.001, 728.001, 729.001, 730.001, 731.001, 732.001, 733.001, 734.001, 735.001, 736.001, 737.001, 738.001, 739.001, 740.001, 741.001, 742.001, 743.001, 744.001, 745.001, 746.001, 747.001, 748.001, 749.001, 750.001, 751.001, 752.001, 753.001, 754.001, 755.001, 756.001, 757.001, 758.001, 759.001, 760.001, 761.001, 762.001, 763.001, 764.001, 765.001, 766.001, 767.001, 768.001, 769.001, 770.001, 771.001, 772.001, 773.001, 774.001, 775.001, 776.001, 777.001, 778.001, 779.001, 780.001, 781.001, 782.001, 783.001, 784.001, 785.001, 786.001, 787.001, 788.001, 789.001, 790.001, 791.001, 792.001, 793.001, 794.001, 795.001, 796.001, 797.001, 798.001, 799.001, 800.001, 801.001, 802.001, 803.001, 804.001, 805.001, 806.001, 807.001, 808.001, 809.001, 810.001, 811.001, 812.001, 813.001, 814.001, 815.001, 816.001, 817.001, 818.001, 819.001, 820.001, 821.001, 822.001, 823.001, 824.001, 825.001, 826.001, 827.001, 828.001, 829.001, 830.001, 831.001, 832.001, 833.001, 834.001, 835.001, 836.001, 837.001, 838.001, 839.001, 840.001, 841.001, 842.001, 843.001, 844.001, 845.001, 846.001, 847.001, 848.001, 849.001, 850.001, 851.001, 852.001, 853.001, 854.001, 855.001, 856.001, 857.001, 858.001, 859.001, 860.001, 861.001, 862.001, 863.001, 864.001, 865.001, 866.001, 867.001, 868.001, 869.001, 870.001, 871.001, 872.001, 873.001, 874.001, 875.001, 876.001, 877.001, 878.001, 879.001, 880.001, 881.001, 882.001, 883.001, 884.001, 885.001, 886.001, 887.001, 888.001, 889.001, 890.001, 891.001, 892.001, 893.001, 894.001, 895.001, 896.001, 897.001, 898.001, 899.001, 900.001, 901.001, 902.001, 903.001, 904.001, 905.001, 906.001, 907.001, 908.001, 909.001, 910.001, 911.001, 912.001, 913.001, 914.001, 915.001, 916.001, 917.001, 918.001, 919.001, 920.001, 921.001, 922.001, 923.001, 924.001, 925.001, 926.001, 927.001, 928.001, 929.001, 930.001, 931.001, 932.001, 933.001, 934.001, 935.001, 936.001, 937.001, 938.001, 939.001, 940.001, 941.001, 942.001, 943.001, 944.001, 945.001, 946.001, 947.001, 948.001, 949.001, 950.001, 951.001, 952.001, 953.001, 954.001, 955.001, 956.001, 957.001, 958.001, 959.001, 960.001, 961.001, 962.001, 963.001, 964.001, 965.001, 966.001, 967.001, 968.001, 969.001, 970.001, 971.001, 972.001, 973.001, 974.001, 975.001, 976.001, 977.001, 978.001, 979.001, 980.001, 981.001, 982.001, 983.001, 984.001, 985.001, 986.001, 987.001, 988.001, 989.001, 990.001, 991.001, 992.001, 993.001, 994.001, 995.001, 996.001, 997.001, 998.001, 999.001, 1020.001, 1041.001, 1062.001, 1083.001, 1104.001, 1125.001, 1146.001, 1167.001, 1188.001, 1209.001, 1230.001, 1251.001, 1272.001, 1293.001, 1314.001, 1335.001, 1356.001, 1377.001, 1398.001, 1419.001, 1440.001, 1461.001, 1482.001, 1503.001, 1524.001, 1545.001, 1566.001, 1587.001, 1608.001, 1629.001, 1650.001, 1671.001, 1692.001, 1713.001, 1734.001, 1755.001, 1776.001, 1797.001, 1818.001, 1839.001, 1860.001, 1881.001, 1902.001, 1923.001, 1944.001, 1965.001, 1986.001, 2007.001, 2028.001, 2049.001, 2070.001, 2091.001, 2112.001, 2133.001, 2154.001, 2175.001, 2196.001, 2217.001, 2238.001, 2259.001, 2280.001, 2301.001, 2322.001, 2343.001, 2364.001, 2385.001, 2406.001, 2427.001, 2448.001, 2469.001, 2490.001, 2511.001, 2532.001, 2553.001, 2574.001, 2595.001, 2616.001, 2637.001, 2658.001};
    return proton_energy_bins;
  }


}


double LEEana::get_reco_Enu_corr(KineInfo& kine, bool flag_data){
  double reco_Enu_corr = 0;
  if (kine.kine_reco_Enu > 0){
    if (flag_data){
      for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
  	if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j) * em_charge_scale;
  	}else{
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j);
  	}
  	//	std::cout << "p: " << kine.kine_energy_particle->at(j) << " " << kine.kine_energy_info->at(j) << " " << kine.kine_particle_type->at(j) << " " << kine.kine_energy_included->at(j) << std::endl;
      }
      reco_Enu_corr += kine.kine_reco_add_energy;
      return reco_Enu_corr;
    }
  }
  return kine.kine_reco_Enu;
}


double LEEana::get_reco_Eproton(KineInfo& kine){
  double reco_Eproton=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) // proton threshold of 35 MeV
      //if(abs(pdgcode)==2212) // no proton threshold
        reco_Eproton+=kine.kine_energy_particle->at(i);

    }
  return reco_Eproton;
}

double LEEana::get_true_Eavail(EvalInfo& eval, PFevalInfo& pfeval, bool useThreshold, bool useMass){
  double Eavail=0;
  //std::cout<<pfeval.truth_Ntrack;
  for(int i=0; i<pfeval.truth_Ntrack; i++){
    int pdgcode = pfeval.truth_pdg[i];
    int mother = pfeval.truth_mother[i];
    //if (mother!=0) { break; }
    if(mother==0){
      double dx = eval.truth_vtxX-pfeval.truth_startXYZT[i][0];
      double dy = eval.truth_vtxY-pfeval.truth_startXYZT[i][1];
      double dz = eval.truth_vtxZ-pfeval.truth_startXYZT[i][2];
      double vtx_diff = sqrt( dx*dx + dy*dy + dz*dz );
      if(vtx_diff>13.25){ continue; }//catches cases with two nu
      if( (eval.truth_nuTime*1000-pfeval.truth_startXYZT[i][3])>0.01){ continue; }//catches cases with two n

      if( abs(pdgcode)==2212 && 1*(pfeval.truth_startMomentum[i][3]*1000-938.27208816)<35 && useThreshold) continue;//below threshold porton, reject if using threshold
      else if( abs(pdgcode)==211 && pfeval.truth_startMomentum[i][3]*1000-139.57039<10 && useThreshold) continue;//below threshold pion, reject if using threshold
      else if( abs(pdgcode)==2212 ) Eavail+=1*(pfeval.truth_startMomentum[i][3]*1000-938.27208816)+useMass*8.36; //proton KE only
      else if(abs(pdgcode)==2112 ) continue; //primary neutron
      else if( abs(pdgcode)==13 && abs(eval.truth_nuPdg)==14 ) continue; //primary muon from numu
      else if( abs(pdgcode)==11 && abs(eval.truth_nuPdg)==12 ) continue; //primary e from nue
      else if( abs(pdgcode)==15 && abs(eval.truth_nuPdg)==16 ) continue; //primary tau from nutau
      else if( abs(pdgcode)==12 || abs(pdgcode)==14 || abs(pdgcode)==16 ) continue; //primary nu
      else if( abs(pdgcode)>3000 && abs(pdgcode)<4000 ) continue; //strange baryon
      else if( abs(pdgcode)==211 ) Eavail+=pfeval.truth_startMomentum[i][3]*1000-139.57039+useMass*139.57039;//pion, KE only
      else if (abs(pdgcode)<10000 ) Eavail+=pfeval.truth_startMomentum[i][3]*1000;//everything else add all energy
    }
  }
  //std::cout<<"  "<<Eavail<<std::endl; 
  return Eavail;
}

double LEEana::get_truth_p_mu_cos(PFevalInfo& pfeval){
    double protonMomentum0 = -1000;
    double protonMomentum1 = -1000;
    double protonMomentum2 = -1000;
    double protonMomentum3 = -1000;
    double Ep=0;
    for(int i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] < Ep) continue;
      Ep = pfeval.truth_startMomentum[i][3];
      protonMomentum0 = pfeval.truth_startMomentum[i][0];
      protonMomentum1 = pfeval.truth_startMomentum[i][1];
      protonMomentum2 = pfeval.truth_startMomentum[i][2];
      protonMomentum3 = pfeval.truth_startMomentum[i][3];
    }

   if (protonMomentum3>0 && pfeval.truth_muonMomentum[3]>0){
      double AdotB = protonMomentum0*pfeval.truth_muonMomentum[0]+protonMomentum1*pfeval.truth_muonMomentum[1]+protonMomentum2*pfeval.truth_muonMomentum[2];
      double magA = sqrt(protonMomentum0*protonMomentum0+protonMomentum1*protonMomentum1+protonMomentum2*protonMomentum2);
      double magB = sqrt(pfeval.truth_muonMomentum[0]*pfeval.truth_muonMomentum[0]+pfeval.truth_muonMomentum[1]*pfeval.truth_muonMomentum[1]+pfeval.truth_muonMomentum[2]*pfeval.truth_muonMomentum[2]);
      double p_mu_angle = AdotB/(magA*magB);
      return p_mu_angle;
    }
    return -1000;

}

bool LEEana::is_true_0p(PFevalInfo& pfeval,double threshold=0.035){
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < threshold) continue; //Erin: CHANGE, no proton threshold
      return false;
    }
  return true;
}

//mother_check=-999 for all particles, mother_check=-9999 for all non-primary, which_p=0 for sum of all particles
double LEEana::get_KE(PFevalInfo& pfeval, int pdg, int truth, int mother_check, int which_p, double threshold){
    std::vector<double> all_KE;

    int Ntrack = pfeval.reco_Ntrack;
    if(truth==1){
      Ntrack = pfeval.truth_Ntrack;
    }
    for(size_t i=0; i<Ntrack; i++){
        int pdgcode = 0;
        double mass = 0;
        int mother = 0;
        double KE = 0;
        if(truth==1){
          pdgcode = pfeval.truth_pdg[i];
          mass = get_mass_MeV(pdgcode);
          mother = pfeval.truth_mother[i];
          KE = pfeval.truth_startMomentum[i][3]*1000-mass;
        }
        else if(truth==0){
          pdgcode = pfeval.reco_pdg[i];
          mass = get_mass_MeV(pdgcode);
          mother = pfeval.reco_mother[i];
          KE = pfeval.reco_startMomentum[i][3]*1000-mass;
        }
        else{
          pdgcode = pfeval.reco_larpid_pdg[i];
          mass = get_mass_MeV(pdgcode);
          mother = pfeval.reco_mother[i];
          KE = pfeval.reco_startMomentum[i][3]*1000-mass;
        }
        if(pdgcode!=pdg) continue;
        if(mother!=mother_check && mother_check!=-999 && mother_check!=-9999) continue;
        if(mother==0 && mother_check==-9999) continue;
        if(KE<threshold) continue;
        all_KE.push_back(KE);
     }
     
     if(all_KE.size()==0){return 0;}
     else if(all_KE.size()<which_p && which_p!=0){return 0;}
 
     std::sort(all_KE.begin(), all_KE.end(), wayToSort);

     if(which_p==0){
       double allp_energy=0;
       for(int i=0; i<all_KE.size(); i++){
         allp_energy+=all_KE.at(i);
       }
       return allp_energy;
     }
     else{ return all_KE.at(which_p-1); }

}

double LEEana::get_muon_energy_new(PFevalInfo& pfeval, bool flag_FC_lepton, bool return_KE, bool return_MeV){
  double E = pfeval.reco_muonMomentum[3];
  if (E<=0) return 0;
  if(pfeval.mcs_emu_tracklen>0 && pfeval.mcs_emu_tracklen<4 && flag_FC_lepton) E = pfeval.mcs_emu_tracklen;
  if(pfeval.mcs_emu_MCS>0 && pfeval.mcs_emu_MCS<4 && !flag_FC_lepton) E = pfeval.mcs_emu_MCS;
  if(return_KE) E = E-0.10566;
  if(return_MeV) E = E*1000;
  return E;
}

// Reco muon total energy [MeV] of kine_reco_Enu_new3_5: range or MCS with the one-sided 5% method (0 without a reco muon).
double LEEana::get_muon_Etot_new(PFevalInfo& pfeval, EvalInfo& eval){
  std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,3,0.05);
  return get_muon_energy_new(pfeval, std::get<0>(result_part_FC), false, true);
}

// Reco muon momentum vector [MeV]: magnitude from get_muon_Etot_new, direction of reco_muonMomentum. Zero vector without
// a reco muon.
TVector3 LEEana::get_muon_momentum_new(PFevalInfo& pfeval, EvalInfo& eval){
  TVector3 dir(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2]);
  if (pfeval.reco_muonMomentum[3]<=0 || dir.Mag()==0) return TVector3(0,0,0);
  double E = get_muon_Etot_new(pfeval, eval);
  double p = (E>105.66) ? sqrt(E*E-105.66*105.66) : 0;
  return p*dir.Unit();
}

// Centroid and principal axis (direction of largest spread) of a set of points
void LEEana::get_principal_axis(const std::vector<TVector3>& pts, TVector3& centroid, TVector3& axis){
  centroid.SetXYZ(0,0,0);
  for(size_t k=0; k<pts.size(); k++) centroid += pts.at(k);
  centroid *= 1./pts.size();
  TMatrixDSym cov(3);
  cov.Zero();
  for(size_t k=0; k<pts.size(); k++){
    TVector3 d = pts.at(k) - centroid;
    for(int a=0; a<3; a++) for(int b=0; b<3; b++) cov(a,b) += d[a]*d[b];
  }
  TMatrixDSymEigen eigen(cov);
  TMatrixD vectors = eigen.GetEigenVectors();   // columns, sorted by decreasing eigenvalue
  axis.SetXYZ(vectors(0,0), vectors(1,0), vectors(2,0));
}

// Kine energy of the particles that continue the primary muon past its reco end (a muon broken by the reconstruction):
// start within drop_dist cm of the muon end and cos(angle between the particle and the muon direction at its end) > drop_cos.
// Any particle type is tagged (showers and the muon-like/pion-like track pieces), except the 22/2112 pseudo-particles and
// protons (their kine energy is replaced by the range KE in get_kine_reco_Enu_new, which cannot be matched to a particle).
// Start: the reco start for showers, the endpoint closer to the muon end for tracks.
// Directions: principal axis of the muon spacepoints within 10 cm of its end (oriented start->end) and of the particle
// spacepoints (oriented away from its start); start->end if fewer than 3 spacepoints.
// Each tagged particle is matched to the kine entry of the same type with the same energy (within 0.5 MeV).
// E_mass returns the change of the masses in kine_reco_add_energy from removing the tagged particles (see
// get_kine_chain_end_mass): mostly 0, as a broken muon carries one mass however many pieces it has.
double LEEana::get_muon_continuation_energy(PFevalInfo& pfeval, KineInfo& kine, SpaceInfo& space, bool flag_data, double drop_dist, double drop_cos, double& E_mass){
  E_mass = 0;
  int mu_index = -1;
  for(int i=0; i<pfeval.reco_Ntrack; i++){
    if(is_pfeval_muon(pfeval,i,0.0001)){ mu_index = i; break; }
  }
  if(mu_index<0) return 0;
  TVector3 mu_start(pfeval.reco_startXYZT[mu_index][0], pfeval.reco_startXYZT[mu_index][1], pfeval.reco_startXYZT[mu_index][2]);
  TVector3 mu_end(pfeval.reco_endXYZT[mu_index][0], pfeval.reco_endXYZT[mu_index][1], pfeval.reco_endXYZT[mu_index][2]);

  int n_spacepoints = space.Trecchargeblob_spacepoints_real_cluster_id->size();
  std::vector<TVector3> mu_end_pts;
  for(int sp=0; sp<n_spacepoints; sp++){
    if(space.Trecchargeblob_spacepoints_real_cluster_id->at(sp)!=pfeval.reco_id[mu_index]) continue;
    TVector3 pt(space.Trecchargeblob_spacepoints_x->at(sp), space.Trecchargeblob_spacepoints_y->at(sp), space.Trecchargeblob_spacepoints_z->at(sp));
    if((pt-mu_end).Mag()<10) mu_end_pts.push_back(pt);
  }
  TVector3 mu_dir = mu_end - mu_start;
  if(mu_end_pts.size()>=3){
    TVector3 centroid;
    get_principal_axis(mu_end_pts, centroid, mu_dir);
    if(mu_dir.Dot(mu_end-mu_start)<0) mu_dir = -mu_dir;
  }
  if(mu_dir.Mag()==0) return 0;
  mu_dir = mu_dir.Unit();

  std::vector<bool> kine_used(kine.kine_energy_particle->size(), false);
  std::vector<bool> dropped(pfeval.reco_Ntrack, false);
  double E_drop = 0;
  for(int j=0; j<pfeval.reco_Ntrack; j++){
    int pdg = abs(pfeval.reco_pdg[j]);
    if(j==mu_index || pdg==22 || pdg==2112 || pdg==2212) continue;
    TVector3 p_start(pfeval.reco_startXYZT[j][0], pfeval.reco_startXYZT[j][1], pfeval.reco_startXYZT[j][2]);
    TVector3 p_end(pfeval.reco_endXYZT[j][0], pfeval.reco_endXYZT[j][1], pfeval.reco_endXYZT[j][2]);
    if(pdg!=11 && (p_end-mu_end).Mag()<(p_start-mu_end).Mag()) std::swap(p_start, p_end);
    if((p_start-mu_end).Mag()>=drop_dist) continue;
    std::vector<TVector3> p_pts;
    for(int sp=0; sp<n_spacepoints; sp++){
      if(space.Trecchargeblob_spacepoints_real_cluster_id->at(sp)!=pfeval.reco_id[j]) continue;
      p_pts.push_back(TVector3(space.Trecchargeblob_spacepoints_x->at(sp), space.Trecchargeblob_spacepoints_y->at(sp), space.Trecchargeblob_spacepoints_z->at(sp)));
    }
    if(p_pts.size()==0) continue;
    TVector3 p_dir = p_end - p_start;
    if(p_pts.size()>=3){
      TVector3 centroid;
      get_principal_axis(p_pts, centroid, p_dir);
      if(p_dir.Dot(centroid-p_start)<0) p_dir = -p_dir;
    }
    if(p_dir.Mag()==0) continue;
    if(p_dir.Unit().Dot(mu_dir)<=drop_cos) continue;
    // remove the matching kine entry (with the same EM scale as get_reco_Enu_corr for data)
    double mass = 0.000511;
    if(pdg==13) mass = 0.105658;
    else if(pdg==211) mass = 0.13957;
    double KE = (pfeval.reco_startMomentum[j][3]-mass)*1000;
    int k_match = -1;
    for(size_t k=0; k<kine.kine_energy_particle->size(); k++){
      if(kine_used.at(k) || abs(kine.kine_particle_type->at(k))!=pdg || fabs(kine.kine_energy_particle->at(k)-KE)>=0.5) continue;
      k_match = k;
      break;
    }
    if(k_match<0) continue;
    kine_used.at(k_match) = true;
    double E_kine = kine.kine_energy_particle->at(k_match);
    if(flag_data && kine.kine_energy_info->at(k_match)==2 && kine.kine_particle_type->at(k_match)==11) E_kine *= em_charge_scale;
    E_drop += E_kine;
    dropped.at(j) = true;
  }

  // masses: muon/pion particles with a kine entry, with and without the dropped particles
  std::vector<bool> kine_used_mass(kine.kine_energy_particle->size(), false);
  std::vector<int> mass_indices;
  for(int j=0; j<pfeval.reco_Ntrack; j++){
    int pdg = abs(pfeval.reco_pdg[j]);
    if(pdg!=13 && pdg!=211) continue;
    double KE = (pfeval.reco_startMomentum[j][3]-(pdg==13 ? 0.105658 : 0.13957))*1000;
    for(size_t k=0; k<kine.kine_energy_particle->size(); k++){
      if(kine_used_mass.at(k) || abs(kine.kine_particle_type->at(k))!=pdg || fabs(kine.kine_energy_particle->at(k)-KE)>=0.5) continue;
      kine_used_mass.at(k) = true;
      mass_indices.push_back(j);
      break;
    }
  }
  std::vector<bool> none_removed(pfeval.reco_Ntrack, false);
  E_mass = get_kine_chain_end_mass(pfeval, mass_indices, none_removed) - get_kine_chain_end_mass(pfeval, mass_indices, dropped);
  return E_drop;
}

// Index of the mother of particle index, skipping removed particles (their daughters count as daughters of their mother).
// -1: attached to the neutrino vertex; -2: attached through a gap (22/2112 pseudo-particle) or the mother is not found.
int LEEana::get_reco_mother_index(PFevalInfo& pfeval, int index, const std::vector<bool>& removed){
  int curr = index;
  for(int depth=0; depth<100; depth++){
    int mother_id = pfeval.reco_mother[curr];
    if(mother_id==0) return -1;
    int mother = -2;
    for(int i=0; i<pfeval.reco_Ntrack; i++){
      if(pfeval.reco_id[i]==mother_id){ mother = i; break; }
    }
    if(mother<0) return -2;
    if(abs(pfeval.reco_pdg[mother])==22 || abs(pfeval.reco_pdg[mother])==2112) return -2;
    if(!removed.at(mother)) return mother;
    curr = mother;
  }
  return -2;
}

// Muon/pion masses in kine_reco_add_energy, following WireCell (NeutrinoID::fill_kine_tree): walking the particle tree
// from the neutrino vertex, each muon/pion adds its mass, and a muon/pion continued by a muon/pion daughter gives its mass
// back. So each muon/pion chain (a broken track, pi->mu) carries one mass, of the particle at its end (one per end if it
// branches). Particles attached through a gap (22/2112 pseudo-particle in their ancestry) add no muon/pion mass.
// indices: the muon/pion particles with a kine entry; removed particles are left out (their daughters move to their mother).
double LEEana::get_kine_chain_end_mass(PFevalInfo& pfeval, const std::vector<int>& indices, const std::vector<bool>& removed){
  double E_mass = 0;
  for(size_t a=0; a<indices.size(); a++){
    int j = indices.at(a);
    if(removed.at(j)) continue;
    // attached to the neutrino vertex through the tree?
    int mother = j;
    for(int depth=0; depth<100 && mother>=0; depth++) mother = get_reco_mother_index(pfeval, mother, removed);
    if(mother!=-1) continue;
    // continued by a muon/pion daughter?
    bool flag_continued = false;
    for(size_t b=0; b<indices.size(); b++){
      int i = indices.at(b);
      if(i==j || removed.at(i)) continue;
      if(get_reco_mother_index(pfeval, i, removed)==j){ flag_continued = true; break; }
    }
    if(flag_continued) continue;
    if(abs(pfeval.reco_pdg[j])==13) E_mass += 105.658;
    else E_mass += 139.570;
  }
  return E_mass;
}

double LEEana::get_kine_reco_Enu_new(PFevalInfo& pfeval, KineInfo& kine, SpaceInfo& space, bool flag_FC_lepton, bool flag_data, bool correct_protons, bool drop_muon_showers, double drop_dist, double drop_cos, bool drop_mass){
  double E = get_reco_Enu_corr(kine, flag_data);
  if (E<0) return 0;
  double KEmuon_old = (pfeval.reco_muonMomentum[3]-0.10566)*1000;
  double KEmuon_new = get_muon_energy_new(pfeval, flag_FC_lepton, true, true);
  if(KEmuon_old>0 && KEmuon_new>0) E =  E - KEmuon_old + KEmuon_new;
  // A muon broken by the reconstruction often has its downstream pieces labelled as separate particles (showers, or
  // muon/pion-like tracks). Range only covers the labelled muon track, so those pieces complete the muon; MCS gives the full
  // muon momentum, so the pieces (and the masses added for them) would be counted twice.
  bool flag_MCS_used = !flag_FC_lepton && pfeval.mcs_emu_MCS>0 && pfeval.mcs_emu_MCS<4;
  if(drop_muon_showers && flag_MCS_used && KEmuon_old>0 && KEmuon_new>0){
    double E_mass = 0;
    E = E - get_muon_continuation_energy(pfeval, kine, space, flag_data, drop_dist, drop_cos, E_mass);
    if(drop_mass) E = E - E_mass;
  }
  if(correct_protons){
    std::vector<double> KEproton_new = std::get<1>(get_range_proton_KE(pfeval, space, true));
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
      int pdgcode = kine.kine_particle_type->at(i);
      if(abs(pdgcode)==2212) E=E-kine.kine_energy_particle->at(i);
    }
    for(size_t i=0; i<KEproton_new.size(); i++){
      E=E+KEproton_new.at(i);
    }
  }
  return E;
}

bool LEEana::check_is_FC(double x, double y, double z) {
    if (x > 3 && x < 250 &&
        y > -113 && y < 113 &&
        z > 3 && z < 1035) {
        return 1;
    } else {
        return 0;
    }
}

bool LEEana::is_pfeval_muon(PFevalInfo& pfeval,int index, double tolerance){
      bool flag_prim_mu = 0;
      if(pfeval.reco_muonMomentum[3]<0) return flag_prim_mu;
      if(pfeval.reco_pdg[index]!=13) return flag_prim_mu; 
      if(pfeval.reco_startMomentum[index][3]>pfeval.reco_muonMomentum[3]-tolerance 
         && pfeval.reco_startMomentum[index][3]<pfeval.reco_muonMomentum[3]+tolerance) flag_prim_mu = 1;
      return flag_prim_mu;
}

// For methods 2 and 3: true if the primary muon energy should come from range, false if from MCS.
// Range from an exiting muon is too low, so range well below MCS indicates a muon that leaves the detector.
bool LEEana::check_muon_range_MCS(PFevalInfo& pfeval, int method, double threshold){
  bool valid_range = pfeval.mcs_emu_tracklen>0 && pfeval.mcs_emu_tracklen<4;
  bool valid_MCS = pfeval.mcs_emu_MCS>0 && pfeval.mcs_emu_MCS<4;
  if(!valid_MCS) return true;   // no usable MCS, use range (get_muon_energy_new falls back to WC if range is not valid either)
  if(!valid_range) return false;
  double rel = (pfeval.mcs_emu_tracklen-pfeval.mcs_emu_MCS)/pfeval.mcs_emu_MCS;
  if(method==2) return fabs(rel)<threshold;   // two-sided: range only if range and MCS agree
  if(method==3) return rel>=-threshold;       // one-sided: MCS only if it is above range by more than the threshold
  return true;
}

std::tuple<bool,bool> LEEana::get_part_is_FC(PFevalInfo& pfeval,EvalInfo& eval, int method, double threshold){
  bool flag_FC_lepton = 1;
  bool flag_FC_hadron = 1;
  std::tuple<bool,bool> result = std::make_tuple(flag_FC_lepton,flag_FC_hadron);
  if(eval.match_isFC && method==0) return result;
  for(size_t i=0; i<pfeval.reco_Ntrack; i++){
    if(pfeval.reco_pdg[i]==22 || pfeval.reco_pdg[i]==2112) continue;
    double x = pfeval.reco_startXYZT[i][0];
    double y = pfeval.reco_startXYZT[i][1];
    double z = pfeval.reco_startXYZT[i][2];
    double x_end = pfeval.reco_endXYZT[i][0];
    double y_end = pfeval.reco_endXYZT[i][1];
    double z_end = pfeval.reco_endXYZT[i][2];
    bool is_FC = check_is_FC(x,y,z) * check_is_FC(x_end,y_end,z_end);
    if(is_pfeval_muon(pfeval,i,0.0001)) flag_FC_lepton = flag_FC_lepton*is_FC;
    else flag_FC_hadron = flag_FC_hadron*is_FC;
  }
  if(eval.match_isFC) flag_FC_hadron = 1;
  // Methods 2/3: only FC events whose muon passes the box keep range, all other events use the range/MCS comparison
  if((method==2 || method==3) && !(eval.match_isFC && flag_FC_lepton)) flag_FC_lepton = check_muon_range_MCS(pfeval, method, threshold);
  result = std::make_tuple(flag_FC_lepton,flag_FC_hadron);
  return result;
}


int LEEana::get_particle_0pNp_bdt_bin(PFevalInfo& pfeval, TaggerInfo& tagger, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern, double KE_threshold, double KE_pl_threshold, double scat_bdt_threshold, double vtxact_bdt_threshold){
  std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval,space,true); 
  double KE = std::get<0>(result).at(0);
  double KE_p = get_pandora_proton_KE(pandora,0.2,true).at(0);
  double KE_l = std::get<2>(result).at(0);
  double KE_ll = get_lantern_KE(lantern, 2212, 10, true).at(0);
  int bdt_bin=1;
  if(KE>=KE_threshold) bdt_bin=5;
  else if(KE_p>KE_pl_threshold || KE_l>KE_pl_threshold || KE_ll>KE_pl_threshold) bdt_bin=4;
  else if(tagger.all_veto_score>scat_bdt_threshold) bdt_bin=3;
  else if(tagger.VtxAct_bdt_score>vtxact_bdt_threshold) bdt_bin=2;
  else if(KE<1 && KE_p<1 && KE_l<1 && KE_ll<1) bdt_bin=0;
  return bdt_bin;
}

// Kinetic energy [MeV] from the range of reco particle i, assuming a proton (see get_range_proton_KE for the range).
// Returns -1 if the particle has no spacepoints.
double LEEana::get_range_proton_KE_particle(PFevalInfo& pfeval, SpaceInfo& space, int i){
  const std::vector<double>& proton_length_bins = get_proton_length_bins(); 
  const std::vector<double>& proton_energy_bins = get_proton_energy_bins();
  int n_spacepoints = space.Trecchargeblob_spacepoints_real_cluster_id->size();

  double range_proton = 0; 
  double proton_KE = 0;
  std::vector<float> temp_spacepoints_x;
  std::vector<float> temp_spacepoints_y;
  std::vector<float> temp_spacepoints_z;
  std::vector<float> temp_spacepoints_q;
  for(size_t sp=0; sp<n_spacepoints; sp++){
    if(space.Trecchargeblob_spacepoints_real_cluster_id->at(sp)==pfeval.reco_id[i]){
      temp_spacepoints_x.push_back(space.Trecchargeblob_spacepoints_x->at(sp));
      temp_spacepoints_y.push_back(space.Trecchargeblob_spacepoints_y->at(sp));
      temp_spacepoints_z.push_back(space.Trecchargeblob_spacepoints_z->at(sp));
      temp_spacepoints_q.push_back(space.Trecchargeblob_spacepoints_q->at(sp));
    }
  }
  int n_spacepoints_part = temp_spacepoints_x.size();
  if(n_spacepoints_part==0) return -1;

  // Sum the steps between consecutive spacepoints, skipping jumps >= 2 cm. Particles can be stored as
  // separate pieces (e.g. the start/end points first, then the trajectory), and summing across those
  // jumps overestimates the range.
  std::vector<bool> on_trajectory(n_spacepoints_part, false);
  for(int sp=0; sp<n_spacepoints_part-1; sp++){
    double dx = temp_spacepoints_x.at(sp) - temp_spacepoints_x.at(sp+1);
    double dy = temp_spacepoints_y.at(sp) - temp_spacepoints_y.at(sp+1);
    double dz = temp_spacepoints_z.at(sp) - temp_spacepoints_z.at(sp+1);
    double dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
    if(dist<2.0){
      range_proton+=dist;
      on_trajectory.at(sp) = true;
      on_trajectory.at(sp+1) = true;
    }
  }

  // Skipping the jumps also drops the short pieces between the track start/end and the trajectory,
  // so connect the reco start and end to the nearest trajectory point if it is within 2 cm.
  // If no steps were kept (e.g. only the start and end points are stored), fall back to the straight distance.
  bool has_trajectory = false;
  for(int sp=0; sp<n_spacepoints_part; sp++) if(on_trajectory.at(sp)) has_trajectory = true;
  if(has_trajectory){
    for(int ep=0; ep<2; ep++){
      double ex = (ep==0) ? pfeval.reco_startXYZT[i][0] : pfeval.reco_endXYZT[i][0];
      double ey = (ep==0) ? pfeval.reco_startXYZT[i][1] : pfeval.reco_endXYZT[i][1];
      double ez = (ep==0) ? pfeval.reco_startXYZT[i][2] : pfeval.reco_endXYZT[i][2];
      double dmin = 1e9;
      for(int sp=0; sp<n_spacepoints_part; sp++){
        if(!on_trajectory.at(sp)) continue;
        double d = sqrt(pow(ex-temp_spacepoints_x.at(sp),2)+pow(ey-temp_spacepoints_y.at(sp),2)+pow(ez-temp_spacepoints_z.at(sp),2));
        if(d<dmin) dmin = d;
      }
      if(dmin<2.0) range_proton+=dmin;
    }
  }else{
    range_proton = sqrt(pow(pfeval.reco_startXYZT[i][0]-pfeval.reco_endXYZT[i][0],2)
                       +pow(pfeval.reco_startXYZT[i][1]-pfeval.reco_endXYZT[i][1],2)
                       +pow(pfeval.reco_startXYZT[i][2]-pfeval.reco_endXYZT[i][2],2));
  }

  int nbb = proton_length_bins.size();
  for(int bb=0; bb<nbb-1; bb++){
    if(proton_length_bins.at(bb)<=range_proton && proton_length_bins.at(bb+1)>range_proton){
      double m = (proton_energy_bins.at(bb+1)-proton_energy_bins.at(bb))/(proton_length_bins.at(bb+1)-proton_length_bins.at(bb));
      proton_KE = proton_energy_bins.at(bb) + m * (range_proton-proton_length_bins.at(bb));
      break;
    }
  }
  return proton_KE;
}

// Leading reco proton: the primary WireCell proton (reco_pdg==2212, reco_mother==0) with the largest range KE, the one
// that sets the reco Np category in get_particle_0pNp_bdt_bin. Returns its index (-1 if none) and fills its KE [MeV]
// (0 if none) and the number of primary WireCell protons with range KE >= KE_threshold.
std::vector<double> LEEana::get_reco_muon_spacepoints_q(PFevalInfo& pfeval, SpaceInfo& space, double tolerance){
  std::vector<double> q;
  // primary muon: the reco muon matching reco_muonMomentum (flag_prim_mu of particle.h)
  int index = -1;
  for(int i=0; i<pfeval.reco_Ntrack; i++){
    if(pfeval.reco_pdg[i]==13 && pfeval.reco_startMomentum[i][3]>pfeval.reco_muonMomentum[3]-tolerance && pfeval.reco_startMomentum[i][3]<pfeval.reco_muonMomentum[3]+tolerance) index = i;
  }
  if(index<0) return q;

  // its spacepoints in stored order (as floats, as in create_particle)
  std::vector<double> x, y, z, q_stored;
  for(size_t sp=0; sp<space.Trecchargeblob_spacepoints_real_cluster_id->size(); sp++){
    if(space.Trecchargeblob_spacepoints_real_cluster_id->at(sp)!=pfeval.reco_id[index]) continue;
    x.push_back((float)space.Trecchargeblob_spacepoints_x->at(sp));
    y.push_back((float)space.Trecchargeblob_spacepoints_y->at(sp));
    z.push_back((float)space.Trecchargeblob_spacepoints_z->at(sp));
    q_stored.push_back((float)space.Trecchargeblob_spacepoints_q->at(sp));
  }
  if(x.empty()) return q;

  // order them from the start to the end: reference point the neutrino vertex for primaries, the nearer of the mother's
  // reco start/end for secondaries (see create_particle)
  double part_start[3], part_end[3];
  for(int c=0; c<3; c++){ part_start[c] = pfeval.reco_startXYZT[index][c]; part_end[c] = pfeval.reco_endXYZT[index][c]; }
  double ref_point[3] = {0, 0, 0};
  bool has_ref = false;
  if(pfeval.reco_mother[index]==0){
    ref_point[0] = pfeval.reco_nuvtxX; ref_point[1] = pfeval.reco_nuvtxY; ref_point[2] = pfeval.reco_nuvtxZ;
    has_ref = true;
  }else{
    for(int mother_part=0; mother_part<pfeval.reco_Ntrack; mother_part++){
      if(pfeval.reco_id[mother_part]!=pfeval.reco_mother[index]) continue;
      double d_mother_end = sqrt(pow(part_start[0]-pfeval.reco_endXYZT[mother_part][0],2)+pow(part_start[1]-pfeval.reco_endXYZT[mother_part][1],2)+pow(part_start[2]-pfeval.reco_endXYZT[mother_part][2],2));
      double d_mother_start = sqrt(pow(part_start[0]-pfeval.reco_startXYZT[mother_part][0],2)+pow(part_start[1]-pfeval.reco_startXYZT[mother_part][1],2)+pow(part_start[2]-pfeval.reco_startXYZT[mother_part][2],2));
      for(int c=0; c<3; c++) ref_point[c] = (d_mother_end<=d_mother_start) ? pfeval.reco_endXYZT[mother_part][c] : pfeval.reco_startXYZT[mother_part][c];
      has_ref = true;
      break;
    }
  }
  ParticleSegments segments = build_particle_segments(x, y, z, q_stored, part_start, part_end, has_ref ? ref_point : nullptr);
  for(size_t k=0; k<segments.main.size(); k++) q.push_back(q_stored.at(segments.main.at(k)));
  return q;
}

int LEEana::get_reco_leading_proton(PFevalInfo& pfeval, SpaceInfo& space, double KE_threshold, double& KE_lead, int& n_protons){
  int index = -1;
  KE_lead = 0;
  n_protons = 0;
  for(int i=0; i<pfeval.reco_Ntrack; i++){
    if(pfeval.reco_pdg[i]!=2212 || pfeval.reco_mother[i]!=0) continue;
    double KE = get_range_proton_KE_particle(pfeval, space, i);
    if(KE<0) continue;
    if(KE>=KE_threshold) n_protons++;
    if(KE>KE_lead){
      KE_lead = KE;
      index = i;
    }
  }
  return index;
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> LEEana::get_range_proton_KE(PFevalInfo& pfeval, SpaceInfo& space, bool return_MeV){

  std::vector<double> prim_proton_KEs;
  std::vector<double> proton_KEs;
  std::vector<double> prim_larpid_proton_KEs;
  std::vector<double> larpid_proton_KEs;

  for(size_t i=0; i<pfeval.reco_Ntrack; i++){

    bool flag_wc=false;
    bool flag_larpid=false;
    if (pfeval.reco_pdg[i]==2212) flag_wc=true;
    if (pfeval.reco_larpid_pdg[i]==2212) flag_larpid=true;
    if(!flag_wc && !flag_larpid) continue;

    double proton_KE = get_range_proton_KE_particle(pfeval, space, i);
    if(proton_KE<0) continue;

    if(flag_wc) proton_KEs.push_back(proton_KE);
    if(flag_larpid) larpid_proton_KEs.push_back(proton_KE);
    if(pfeval.reco_mother[i]==0){
      if(flag_wc) prim_proton_KEs.push_back(proton_KE);
      if(flag_larpid) prim_larpid_proton_KEs.push_back(proton_KE);
    }
  }

  if(proton_KEs.size()==0) proton_KEs.push_back(0);
  if(larpid_proton_KEs.size()==0) larpid_proton_KEs.push_back(0);
  if(prim_proton_KEs.size()==0) prim_proton_KEs.push_back(0);
  if(prim_larpid_proton_KEs.size()==0) prim_larpid_proton_KEs.push_back(0);
       
  std::sort(prim_proton_KEs.begin(), prim_proton_KEs.end(), wayToSort);
  std::sort(proton_KEs.begin(), proton_KEs.end(), wayToSort);
  std::sort(prim_larpid_proton_KEs.begin(), prim_larpid_proton_KEs.end(), wayToSort);
  std::sort(larpid_proton_KEs.begin(), larpid_proton_KEs.end(), wayToSort);
  if(!return_MeV){
    for (auto& val : prim_proton_KEs) { val *= 0.001;}
    for (auto& val : proton_KEs) { val *= 0.001;}
    for (auto& val : prim_larpid_proton_KEs) { val *= 0.001;}
    for (auto& val : larpid_proton_KEs) { val *= 0.001;}
  }
  std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = std::make_tuple(prim_proton_KEs,proton_KEs,prim_larpid_proton_KEs,larpid_proton_KEs);
  return result;

}

double LEEana::get_mass_GeV(int pdg){
  if (abs(pdg)==13) return 0.10566; 
  if (abs(pdg)==2212) return 0.938272; 
  if (abs(pdg)==2112) return 0.93957;
  if (abs(pdg)==211) return 0.139570;
  if (abs(pdg)==111) return 0.135;
  return 0;
}

double LEEana::get_mass_MeV(int pdg){
  return get_mass_GeV(pdg)*1000;
}

std::vector<double> LEEana::get_pandora_proton_KE(PandoraInfo& pandora, double TRACK_SCORE_CUT, bool return_MeV){

  std::vector<double> proton_KEs;
  bool pandora_inclusive_flag = 0;
  double muon_trk_llr_pid_score_v = -1;
    
  if (pandora.slice_orig_pass_id != 1){
    proton_KEs.push_back(0);
    return proton_KEs;
  }

  for(size_t part=0; part<pandora.n_pfps; part++){
        
    int generation = pandora.pfp_generation_v->at(part);
    if(generation!=2) continue;

    if(pandora.pfpdg->at(part)!=13) continue;

    //if(pandora.trk_llr_pid_score_v->at(part)>muon_trk_llr_pid_score_v) muon_trk_llr_pid_score_v = pandora.trk_llr_pid_score_v->at(part);

    if(pandora.trk_llr_pid_score_v->at(part) < TRACK_SCORE_CUT){ proton_KEs.push_back(pandora.trk_energy_proton_v->at(part)); }
  }

  if(return_MeV){
    for(size_t part=0; part<proton_KEs.size(); part++){
      proton_KEs.at(part) = proton_KEs.at(part)*1000;
    }
  }
               
  if(proton_KEs.size()==0) proton_KEs.push_back(0); 
  std::sort(proton_KEs.begin(), proton_KEs.end(), wayToSort);
  return proton_KEs;

}

std::vector<double> LEEana::get_lantern_KE(LanternInfo& lantern, int pdg, double vtx_cut, bool return_MeV){

  std::vector<double> p_KEs;

  for(size_t part=0; part<lantern.nTracks; part++){

    if(lantern.trackPID[part]!=pdg) continue;
    if(lantern.trackIsSecondary[part]!=0) continue;
    if(lantern.trackDistToVtx[part]>vtx_cut) continue;
    p_KEs.push_back(lantern.trackRecoE[part]);
  }

  if(!return_MeV){
    for(size_t part=0; part<p_KEs.size(); part++){
      p_KEs.at(part) = p_KEs.at(part)/1000;
    }
  }

  if(p_KEs.size()==0) p_KEs.push_back(0);
  std::sort(p_KEs.begin(), p_KEs.end(), wayToSort);
  return p_KEs;

}


// The rootino POT ratio is set by the weight name: cv_spline_rootino_<ratio>, <ratio> = (runs 1-5 POT)/(runs 4-5 POT).
// Also accepts the squared (err2) form cv_spline_rootino_<ratio>_cv_spline_rootino_<ratio>, and the plain cv_spline_rootino (ratio 1).
double LEEana::get_rootino_ratio(TString weight_name){
  TString key = "cv_spline_rootino";
  if (!weight_name.BeginsWith(key)) return 1.0;
  TString rest = weight_name(key.Length(), weight_name.Length()-key.Length()); // "", "_1.25", "_1.25_cv_spline_rootino_1.25", "_cv_spline_rootino"
  if (!rest.BeginsWith("_")) return 1.0;
  rest = rest(1, rest.Length()-1);
  Ssiz_t end = rest.Index("_");
  TString token = (end==kNPOS) ? rest : TString(rest(0, end));
  if (token == "cv") return 1.0; // plain cv_spline_rootino_cv_spline_rootino
  if (!token.IsFloat()){
    std::cout << "ERROR: could not read the rootino POT ratio from weight " << weight_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return token.Atof();
}

double LEEana::get_rootino_weight(EvalInfo& eval, GleeInfo& glee, double rootino_pot_ratio){
  // rootino bug fix: events with GTruth_ResNum==9 are dropped before run 4a (run 18961)
  // and the remaining ones are scaled up by rootino_pot_ratio = (runs 1-5 POT)/(runs 4-5 POT)
  if (glee.GTruth_ResNum!=9) return 1.0;
  if (eval.run < 18961) return 0.0;
  return rootino_pot_ratio;
}

double LEEana::get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, GleeInfo& glee, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data){
  double addtl_weight = 1.0;

  //Begin reweighting
  std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double> > rw_info_i;

  if(std::get<0>(rw_info) && !(flag_data)){//Are you applying any reweighting?
    for(size_t rw=0; rw<std::get<1>(rw_info).size(); rw++){
      rw_info_i = std::get<1>(rw_info)[rw];
      if(std::get<0>(rw_info_i)){//Are you reweighting this channel and cut?
        TString cut_str = std::get<1>(rw_info_i);
        TString var_str = std::get<2>(rw_info_i);
        double var = get_truth_var(kine, eval, pfeval, tagger, var_str);
        double min_var = std::get<3>(rw_info_i);
        double max_var = std::get<4>(rw_info_i);
        bool underflow = std::get<5>(rw_info_i);
        bool overflow = std::get<6>(rw_info_i);
        bool equal_binning = std::get<7>(rw_info_i);
        std::vector<double> reweight = std::get<8>(rw_info_i);

        int wbin;
        bool flag_pass = get_rw_cut_pass(cut_str, eval, pfeval, tagger, kine);
        if (flag_pass){
          if (var>max_var && overflow) addtl_weight = reweight.back();
          else if(var>max_var) addtl_weight = 1;
          else if (var<min_var && underflow) addtl_weight = reweight[0];
          else if (var>min_var){
            if(equal_binning){
              double bin_len = (max_var-min_var)/reweight.size();
              if(underflow && overflow) bin_len = (max_var-min_var)/(reweight.size()-2);
              else if (underflow || overflow) bin_len = (max_var-min_var)/(reweight.size()-1);
              wbin = floor((var-min_var)/bin_len);
            }else{
              std::vector<double> bins = std::get<9>(rw_info_i);
              for(int b=0; b<bins.size()-1; b++){
                if(var<=bins[b+1] && var>bins[b]){
                  wbin = b;
                  break;
                }
              }
            }
            if(underflow) wbin++;
            addtl_weight *= reweight[wbin];
          }
        }
      }
    }
  }

  if (weight_name == "cv_spline"){
    return addtl_weight*eval.weight_cv * eval.weight_spline;
  //rootino bug fix weights, see get_rootino_weight
  }else if (weight_name.BeginsWith("cv_spline_rootino")){
    double rootino_weight = get_rootino_weight(eval, glee, get_rootino_ratio(weight_name));
    if (weight_name.Index("_cv_spline_rootino") == kNPOS) return addtl_weight*eval.weight_cv * eval.weight_spline * rootino_weight;
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * rootino_weight,2);
  //cex bug fix weights
  }else if (weight_name == "cv_spline_cexbugfix"){
    double ratio_weight = 1.0;
    if (eval.truth_isCC == 0 && pfeval.truth_NprimPio==1)
    {
      //get the pi0 costheta and KE
      double truth_pi0_costheta = -1000.;
      double truth_pi0_KE = -1000.;
      int true_num_protons_35_MeV = 0;
      for(int jth=0; jth<pfeval.truth_Ntrack; jth++){
        int mother = pfeval.truth_mother[jth];
        if(mother != 0) continue;
        int pdgcode = pfeval.truth_pdg[jth];
        if(abs(pdgcode)==111){
          //N_th_pi0++;
          double px = pfeval.truth_startMomentum[jth][0]*1000.; // MeV
          double py = pfeval.truth_startMomentum[jth][1]*1000.; // MeV
          double pz = pfeval.truth_startMomentum[jth][2]*1000.; // MeV
          truth_pi0_costheta = pz / sqrt(px*px + py*py + pz*pz);
          truth_pi0_KE = pfeval.truth_startMomentum[jth][3]*1000. - 134.9768;
        }
        if (pdgcode==2212 && pfeval.truth_startMomentum[jth][3]*1000. - 938.272089 > 35.){
          true_num_protons_35_MeV++;
        }
      }
      //pick 0p or Np csv file
      std::string ratiofilename;
      if (true_num_protons_35_MeV>0){
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_Np1pi0.csv";
      }else{
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_0p1pi0.csv";
      }
      std::ifstream ratiofile(ratiofilename);
      if (!ratiofile.is_open()) {
          throw std::runtime_error("Could not open file");
      }
      std::string ratioline;
      std::getline(ratiofile, ratioline); // Skip header line
      while (std::getline(ratiofile, ratioline)) {
          std::istringstream ss(ratioline);
          char comma;
          double pi0_KE;
          double pi0_cos;
          double ratio;
          ss >> pi0_KE >> comma >> pi0_cos >> comma >> ratio;
          double pi0_cos_low = pi0_cos - 0.02;
          double pi0_cos_high = pi0_cos + 0.02;
          double pi0_KE_low = pi0_KE - 5.0;
          double pi0_KE_high = pi0_KE + 5.0;
          if (pi0_cos_low <= truth_pi0_costheta && truth_pi0_costheta <= pi0_cos_high && pi0_KE_low <= truth_pi0_KE && truth_pi0_KE <= pi0_KE_high){
            ratio_weight = ratio;
            break;
          }  
      }
    }
    if (ratio_weight > 10.0) ratio_weight = 10.0;
    if (ratio_weight < 0.0) ratio_weight = 0.0;
    return addtl_weight*eval.weight_cv * eval.weight_spline * ratio_weight;
  }else if (weight_name == "cv_spline_cexbugfix_cv_spline_cexbugfix"){
    double ratio_weight = 1.0;
    if (eval.truth_isCC == 0 && pfeval.truth_NprimPio==1)
    {
      //get the pi0 costheta and KE
      double truth_pi0_costheta = -1000.;
      double truth_pi0_KE = -1000.;
      int true_num_protons_35_MeV = 0;
      for(int jth=0; jth<pfeval.truth_Ntrack; jth++){
        int mother = pfeval.truth_mother[jth];
        if(mother != 0) continue;
        int pdgcode = pfeval.truth_pdg[jth];
        if(abs(pdgcode)==111){
          //N_th_pi0++;
          double px = pfeval.truth_startMomentum[jth][0]*1000.; // MeV
          double py = pfeval.truth_startMomentum[jth][1]*1000.; // MeV
          double pz = pfeval.truth_startMomentum[jth][2]*1000.; // MeV
          truth_pi0_costheta = pz / sqrt(px*px + py*py + pz*pz);
          truth_pi0_KE = pfeval.truth_startMomentum[jth][3]*1000. - 134.9768;
        }
        if (pdgcode==2212 && pfeval.truth_startMomentum[jth][3]*1000. - 938.272089 > 35.){
          true_num_protons_35_MeV++;
        }
      }
      //pick 0p or Np csv file
      std::string ratiofilename;
      if (true_num_protons_35_MeV>0){
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_Np1pi0.csv";
      }else{
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_0p1pi0.csv";
      }
      std::ifstream ratiofile(ratiofilename);
      if (!ratiofile.is_open()) {
          throw std::runtime_error("Could not open file");
      }
      std::string ratioline;
      std::getline(ratiofile, ratioline); // Skip header line
      while (std::getline(ratiofile, ratioline)) {
          std::istringstream ss(ratioline);
          char comma;
          double pi0_KE;
          double pi0_cos;
          double ratio;
          ss >> pi0_KE >> comma >> pi0_cos >> comma >> ratio;
          double pi0_cos_low = pi0_cos - 0.02;
          double pi0_cos_high = pi0_cos + 0.02;
          double pi0_KE_low = pi0_KE - 5.0;
          double pi0_KE_high = pi0_KE + 5.0;
          if (pi0_cos_low <= truth_pi0_costheta && truth_pi0_costheta <= pi0_cos_high && pi0_KE_low <= truth_pi0_KE && truth_pi0_KE <= pi0_KE_high){
            ratio_weight = ratio;
            break;
          }
      }
    }
    if (ratio_weight > 10.0) ratio_weight = 10.0;
    if (ratio_weight < 0.0) ratio_weight = 0.0;
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * ratio_weight,2);
  }else if (weight_name == "cv_spline_cv_spline"){
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "unity" || weight_name == "unity_unity"){
    return 1;
  }else if (weight_name == "lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return (eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline);
  }else if (weight_name == "lee_cv_spline_lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return pow(eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "lee_cv_spline_cv_spline" || weight_name == "cv_spline_lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return eval.weight_lee * pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "spline"){
    return eval.weight_spline;
  }else if (weight_name == "spline_spline"){
    return pow(eval.weight_spline,2);
  }else if (weight_name == "lee_spline"){
    return (eval.weight_lee * eval.weight_spline);
  }else if (weight_name == "lee_spline_lee_spline"){
    return pow(eval.weight_lee * eval.weight_spline,2);
  }else if (weight_name == "lee_spline_spline" || weight_name == "spline_lee_spline"){
    return eval.weight_lee * pow( eval.weight_spline,2);
  }else if (weight_name == "add_weight"){//for systematics
    return addtl_weight;
  }else{
    std::cout <<"Unknown weights: " << weight_name << std::endl;
  }


  return 1;
}

double LEEana::get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger , TString var_name){
  if(var_name == "truth_energyInside"){
    return eval.truth_energyInside;
  }else {std::cout<<"Unknown truth var, check configurations"<<std::endl;}
  return 0;
}


double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){
  //  if (var_name == "kine_reco_Enu"){
  //  return kine.kine_reco_Enu;
  //  }else
  if (var_name == "kine_reco_Enu"){
    return get_reco_Enu_corr(kine, flag_data);
  }else if (var_name == "kine_reco_Enu_new"){
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval);
    return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true);
  }else if (var_name == "KE_muon_new"){
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval);
    return get_muon_energy_new(pfeval, std::get<0>(result_part_FC), true, true);
  }else if (var_name == "kine_reco_Enu_new2_5" || var_name == "kine_reco_Enu_new2_10" || var_name == "kine_reco_Enu_new2_15" || var_name == "kine_reco_Enu_new2_20"
         || var_name == "kine_reco_Enu_new3_5" || var_name == "kine_reco_Enu_new3_10" || var_name == "kine_reco_Enu_new3_15" || var_name == "kine_reco_Enu_new3_20"){
    TString tag = var_name; tag.ReplaceAll("kine_reco_Enu_new","");   // e.g. "3_5": method 3, threshold 5%
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,TString(tag(0,1)).Atoi(),TString(tag(2,tag.Length()-2)).Atof()/100.);
    return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true);
  }else if (var_name == "kine_reco_Enu_new3_5_drop10_95" || var_name == "kine_reco_Enu_new3_5_drop15_90" || var_name == "kine_reco_Enu_new3_5_drop15_95" || var_name == "kine_reco_Enu_new3_5_drop30_95"){
    // new3_5, plus dropping the particles (showers and muon/pion-like pieces, with their masses) that continue the muon when its energy comes from MCS
    // "dropD_C": particle start within D cm of the muon end and cos(angle to the muon direction) > 0.C
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,3,0.05);
    if(var_name == "kine_reco_Enu_new3_5_drop10_95") return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true, true, 10, 0.95);
    if(var_name == "kine_reco_Enu_new3_5_drop15_90") return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true, true, 15, 0.90);
    if(var_name == "kine_reco_Enu_new3_5_drop15_95") return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true, true, 15, 0.95);
    return get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true, true, 30, 0.95);
  }else if (var_name == "KE_muon_new2_5" || var_name == "KE_muon_new2_10" || var_name == "KE_muon_new2_15" || var_name == "KE_muon_new2_20"
         || var_name == "KE_muon_new3_5" || var_name == "KE_muon_new3_10" || var_name == "KE_muon_new3_15" || var_name == "KE_muon_new3_20"){
    TString tag = var_name; tag.ReplaceAll("KE_muon_new","");
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,TString(tag(0,1)).Atoi(),TString(tag(2,tag.Length()-2)).Atof()/100.);
    return get_muon_energy_new(pfeval, std::get<0>(result_part_FC), true, true);
  }else if (var_name == "Ehadron_new"
         || var_name == "Ehadron_new2_5" || var_name == "Ehadron_new2_10" || var_name == "Ehadron_new2_15" || var_name == "Ehadron_new2_20"
         || var_name == "Ehadron_new3_5" || var_name == "Ehadron_new3_10" || var_name == "Ehadron_new3_15" || var_name == "Ehadron_new3_20"
         || var_name == "Ehadron_new3_5_drop10_95" || var_name == "Ehadron_new3_5_drop15_90" || var_name == "Ehadron_new3_5_drop15_95" || var_name == "Ehadron_new3_5_drop30_95"){
    // hadronic energy of the matching kine_reco_Enu_new* variable: that Enu minus the total energy (KE + mass) of the
    // muon, with the muon energy chosen the same way (range or MCS for the same method and threshold)
    if (pfeval.reco_muonMomentum[3]<=0) return -1000;
    TString enu_name = var_name; enu_name.ReplaceAll("Ehadron_new","kine_reco_Enu_new");
    TString tag = var_name; tag.ReplaceAll("Ehadron_new","");     // e.g. "", "3_5", "3_5_drop15_95"
    int method = 0;
    double threshold = 0.05;
    if (tag.Length()>0){
      method = TString(tag(0,1)).Atoi();
      TString thr = tag(2,tag.Length()-2);
      if (thr.Index("_")>=0) thr = thr(0,thr.Index("_"));
      threshold = thr.Atof()/100.;
    }
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,method,threshold);
    double E_muon = get_muon_energy_new(pfeval, std::get<0>(result_part_FC), false, true);   // total energy, MeV
    double Ehadron = get_kine_var(kine, eval, pfeval, tagger, flag_data, enu_name, space, pandora, lantern) - E_muon;
    if (Ehadron<0) return 0;   // no hadronic energy (negative only by rounding)
    return Ehadron;
  }else if (var_name == "Eavail"){
    // available energy: Enu without the added masses/binding energy (kine_reco_add_energy) and without the muon KE
    double Emu = pfeval.reco_muonMomentum[3]*1000-105.6583755;
    if (pfeval.reco_muonMomentum[3]<=0) Emu=0;
    double Enu = get_reco_Enu_corr(kine, flag_data);
    double Eadd = kine.kine_reco_add_energy;
    double Eavail = Enu - Eadd - Emu;
    if (Eavail<0) return 0;   // no available energy
    return Eavail;
  }else if (var_name == "Eavail_new"
         || var_name == "Eavail_new2_5" || var_name == "Eavail_new2_10" || var_name == "Eavail_new2_15" || var_name == "Eavail_new2_20"
         || var_name == "Eavail_new3_5" || var_name == "Eavail_new3_10" || var_name == "Eavail_new3_15" || var_name == "Eavail_new3_20"
         || var_name == "Eavail_new3_5_drop10_95" || var_name == "Eavail_new3_5_drop15_90" || var_name == "Eavail_new3_5_drop15_95" || var_name == "Eavail_new3_5_drop30_95"){
    // available energy of the matching kine_reco_Enu_new* variable: that Enu minus kine_reco_add_energy (masses and
    // binding energy) and minus the muon KE, with the muon energy chosen the same way (range or MCS for the same method and threshold)
    TString enu_name = var_name; enu_name.ReplaceAll("Eavail_new","kine_reco_Enu_new");
    TString tag = var_name; tag.ReplaceAll("Eavail_new","");      // e.g. "", "3_5", "3_5_drop15_95"
    int method = 0;
    double threshold = 0.05;
    if (tag.Length()>0){
      method = TString(tag(0,1)).Atoi();
      TString thr = tag(2,tag.Length()-2);
      if (thr.Index("_")>=0) thr = thr(0,thr.Index("_"));
      threshold = thr.Atof()/100.;
    }
    std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval,eval,method,threshold);
    double KE_muon = get_muon_energy_new(pfeval, std::get<0>(result_part_FC), true, true);   // KE, MeV (0 without a muon)
    double Enu = 0;
    if (tag.Index("drop")>=0){
      // same as kine_reco_Enu_new*_dropD_C, but keeping the masses of the dropped particles (all of kine_reco_add_energy is removed below)
      TString drop = tag(tag.Index("drop")+4, tag.Length());          // e.g. "15_95"
      double drop_dist = TString(drop(0,drop.Index("_"))).Atof();
      double drop_cos = TString(drop(drop.Index("_")+1,drop.Length())).Atof()/100.;
      Enu = get_kine_reco_Enu_new(pfeval, kine, space, std::get<0>(result_part_FC), flag_data, true, true, drop_dist, drop_cos, false);
    }else{
      Enu = get_kine_var(kine, eval, pfeval, tagger, flag_data, enu_name, space, pandora, lantern);
    }
    double Eavail = Enu - kine.kine_reco_add_energy - KE_muon;
    if (Eavail<0) return 0;   // no available energy
    return Eavail;

  }else if (var_name == "muon_pt" || var_name == "muon_pl" || var_name == "Emu_costheta_bin" || var_name == "pl_pt_bin" || var_name == "Eavail_Emu_bin"
         || var_name == "Eavail_costheta_Emu_bin" || var_name == "Eavail_pl_pt_bin"
         || var_name == "Emu_costheta_bin_t" || var_name == "pl_pt_bin_t" || var_name == "Eavail_Emu_bin_t"
         || var_name == "Eavail_costheta_Emu_bin_t" || var_name == "Eavail_pl_pt_bin_t"){
    // Reco variables of the muon (and Eavail) cross-section measurements. Muon energy / momentum as in kine_reco_Enu_new3_5
    // (get_muon_Etot_new / get_muon_momentum_new), muon direction of reco_muonMomentum, Eavail as Eavail_new3_5_drop15_95
    // (0 for negative values, in the first slice).
    //   muon_pt, muon_pl: transverse / longitudinal muon momentum w.r.t. the beam [MeV] (truth bins: cut_file 14, 15).
    //   "*_bin": the reco bin of a multi-differential measurement as the flattened bin index + 0.5, for histograms from 0
    //   to nbin. Same slices as the truth bins of get_xs_signal_no; finer bins (50 MeV), each a part of one truth bin.
    //     Emu_costheta_bin (cut_file 12): nbin = 146
    //     pl_pt_bin (16): nbin = 114
    //     Eavail_Emu_bin (17): nbin = 127
    //     Eavail_costheta_Emu_bin (21), Eavail_pl_pt_bin (22): reco 0p and reco Np (leading proton below / above 45 MeV,
    //     the reco Np category) have their own binnings, as the truth bins, so the channels need different nbin:
    //     Eavail_costheta_Emu_bin 198 (reco 0p channels) / 329 (reco Np channel), Eavail_pl_pt_bin 180 / 296.
    //   "*_bin_t": the same in the truth binning (for channels without the statistics for the finer bins), nbin =
    //     Emu_costheta_bin_t 32, pl_pt_bin_t 32, Eavail_Emu_bin_t 40, Eavail_costheta_Emu_bin_t 56 / 84 (reco 0p / Np),
    //     Eavail_pl_pt_bin_t 66 / 95.
    //   Every event the numuCC_part_bdt channels pass (reco_muonMomentum[3] >= 0) gets a bin; -1 without a reco muon.
    static const std::vector<std::vector<double>> Emu_bins_costheta = {
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 650, 850},   // costheta <= 0
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 850},   // costheta 0 - 0.3
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 950},   // costheta 0.3 - 0.5
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1250, 1450},   // costheta 0.5 - 0.7
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1050, 1150, 1300, 1550, 2100},   // costheta 0.7 - 0.8
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1400, 1500, 1600, 1750, 2000, 2400},   // costheta 0.8 - 0.9
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1600, 1700, 1800, 1900, 2050, 2200, 2400, 2750}   // costheta > 0.9
    };
    static const std::vector<std::vector<double>> pt_bins_pl = {
      {100, 150, 200, 250, 300, 350, 400, 450, 500, 600},   // pl <= 0
      {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 650},   // pl 0 - 150
      {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 650, 750},   // pl 150 - 300
      {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750},   // pl 300 - 450
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800},   // pl 450 - 625
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750, 850},   // pl 625 - 850
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 850, 950},   // pl 850 - 1200
      {200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1400, 1700}   // pl > 1200
    };
    static const std::vector<std::vector<double>> Emu_bins_Eavail = {
      {350, 450, 550, 700, 900, 1150, 1450},   // Eavail <= 75
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1200, 1300, 1400, 1600, 1850, 2300},   // Eavail 75 - 150
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 950, 1000, 1050, 1150, 1250, 1350, 1550, 1750},   // Eavail 150 - 250
      {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 850, 950, 1150, 1350, 1450, 1800},   // Eavail 250 - 375
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1150, 1350, 1600},   // Eavail 375 - 550
      {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1050, 1250, 1600},   // Eavail 550 - 800
      {150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1150, 1450}   // Eavail > 800
    };
    static const std::vector<std::vector<std::vector<double>>> Emu_bins_3d_0p = {
      {
       {250, 300, 350, 400, 500},   // Eavail <= 75, costheta <= 0
       {250, 300, 350, 400, 450, 550, 650, 800, 1000},   // Eavail <= 75, costheta 0 - 0.5
       {300, 350, 400, 450, 550, 700, 800, 950, 1150},   // Eavail <= 75, costheta 0.5 - 0.7
       {250, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1050, 1150, 1250, 1350, 1500, 1700, 2150},   // Eavail <= 75, costheta 0.7 - 0.9
       {300, 450, 550, 600, 650, 700, 750, 800, 900, 950, 1000, 1100, 1150, 1200, 1300, 1400, 1600, 1750, 2050}   // Eavail <= 75, costheta > 0.9
      },
      {
       {300, 400},   // Eavail 75 - 150, costheta <= 0
       {350, 450, 550, 700},   // Eavail 75 - 150, costheta 0 - 0.5
       {350, 450, 550, 700, 850},   // Eavail 75 - 150, costheta 0.5 - 0.7
       {300, 450, 550, 650, 700, 750, 850, 900, 950, 1050, 1150, 1300, 1650},   // Eavail 75 - 150, costheta 0.7 - 0.9
       {350, 450, 600, 700, 900, 1000, 1150, 1300, 1450, 1800}   // Eavail 75 - 150, costheta > 0.9
      },
      {
       {150, 200, 250, 300, 350, 400, 450, 550, 700},   // Eavail > 150, costheta <= 0
       {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 950, 1350},   // Eavail > 150, costheta 0 - 0.5
       {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750, 850, 1000, 1250},   // Eavail > 150, costheta 0.5 - 0.7
       {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1150, 1250, 1350, 1500, 1700, 2150},   // Eavail > 150, costheta 0.7 - 0.9
       {200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 850, 900, 950, 1000, 1050, 1150, 1250, 1350, 1600, 1750, 2050, 2550}   // Eavail > 150, costheta > 0.9
      }
    };
    static const std::vector<std::vector<std::vector<double>>> Emu_bins_3d_Np = {
      {
       {250, 300, 400},   // Eavail <= 150, costheta <= 0
       {250, 350, 450},   // Eavail <= 150, costheta 0 - 0.3
       {300, 400, 500},   // Eavail <= 150, costheta 0.3 - 0.5
       {200, 350, 450, 500, 550, 650, 800, 950},   // Eavail <= 150, costheta 0.5 - 0.7
       {300, 450, 550, 600, 700, 800, 950},   // Eavail <= 150, costheta 0.7 - 0.8
       {250, 400, 550, 650, 700, 750, 800, 850, 900, 1000, 1100, 1200, 1400},   // Eavail <= 150, costheta 0.8 - 0.9
       {250, 450, 600, 700, 800, 900, 950, 1050, 1100, 1150, 1250, 1350, 1450, 1600, 1800, 2000, 2300}   // Eavail <= 150, costheta > 0.9
      },
      {
       {200, 250, 300, 350, 400},   // Eavail 150 - 250, costheta <= 0
       {200, 300, 400, 500},   // Eavail 150 - 250, costheta 0 - 0.3
       {200, 300, 350, 450, 500, 600},   // Eavail 150 - 250, costheta 0.3 - 0.5
       {200, 300, 350, 400, 450, 550, 600, 700, 800, 950},   // Eavail 150 - 250, costheta 0.5 - 0.7
       {200, 350, 450, 550, 650, 700, 750, 850, 1000, 1250},   // Eavail 150 - 250, costheta 0.7 - 0.8
       {200, 350, 450, 550, 600, 650, 700, 750, 800, 900, 1000, 1050, 1150, 1250, 1400, 1700},   // Eavail 150 - 250, costheta 0.8 - 0.9
       {250, 400, 500, 600, 650, 700, 800, 900, 950, 1000, 1050, 1100, 1150, 1250, 1350, 1450, 1550, 1750, 1850, 2150}   // Eavail 150 - 250, costheta > 0.9
      },
      {
       {200, 250, 300, 350, 400, 450},   // Eavail 250 - 375, costheta <= 0
       {200, 300, 400, 500, 600},   // Eavail 250 - 375, costheta 0 - 0.3
       {200, 300, 450, 550, 650},   // Eavail 250 - 375, costheta 0.3 - 0.5
       {200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 850, 1000},   // Eavail 250 - 375, costheta 0.5 - 0.7
       {250, 350, 450, 550, 700, 800, 900, 1150},   // Eavail 250 - 375, costheta 0.7 - 0.8
       {200, 300, 400, 500, 550, 600, 700, 750, 850, 900, 1000, 1100, 1200, 1350},   // Eavail 250 - 375, costheta 0.8 - 0.9
       {250, 400, 500, 600, 700, 750, 800, 900, 950, 1000, 1050, 1150, 1250, 1350, 1450, 1600, 1850, 2200}   // Eavail 250 - 375, costheta > 0.9
      },
      {
       {200, 250, 300, 350, 400, 450},   // Eavail 375 - 550, costheta <= 0
       {200, 300, 450, 550},   // Eavail 375 - 550, costheta 0 - 0.3
       {200, 300, 350, 450, 600},   // Eavail 375 - 550, costheta 0.3 - 0.5
       {200, 300, 400, 450, 550, 600, 700, 800, 950},   // Eavail 375 - 550, costheta 0.5 - 0.7
       {200, 350, 450, 550, 700, 850, 1100},   // Eavail 375 - 550, costheta 0.7 - 0.8
       {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1300},   // Eavail 375 - 550, costheta 0.8 - 0.9
       {300, 400, 550, 700, 800, 900, 1000, 1150, 1300, 1550, 2100}   // Eavail 375 - 550, costheta > 0.9
      },
      {
       {150, 200, 250, 300, 350, 400, 450, 550},   // Eavail > 550, costheta <= 0
       {150, 200, 350, 450, 600},   // Eavail > 550, costheta 0 - 0.3
       {200, 300, 450, 550},   // Eavail > 550, costheta 0.3 - 0.5
       {150, 250, 350, 450, 550, 700, 850},   // Eavail > 550, costheta 0.5 - 0.7
       {250, 350, 550, 750, 950},   // Eavail > 550, costheta 0.7 - 0.8
       {200, 300, 400, 500, 600, 700, 850, 1000, 1300},   // Eavail > 550, costheta 0.8 - 0.9
       {250, 400, 550, 650, 750, 900, 1050, 1250, 1550}   // Eavail > 550, costheta > 0.9
      }
    };
    static const std::vector<std::vector<std::vector<double>>> pt_bins_3d_0p = {
      {
       {100, 150, 200, 250, 300, 350, 400, 450, 550},   // Eavail <= 75, pl <= 150
       {100, 150, 200, 250, 300, 350, 400, 450, 600},   // Eavail <= 75, pl 150 - 300
       {150, 200, 250, 300, 350, 400, 500, 650},   // Eavail <= 75, pl 300 - 450
       {150, 200, 250, 300, 350, 400, 450, 500, 550, 700},   // Eavail <= 75, pl 450 - 625
       {200, 250, 300, 350, 400, 450, 550, 650},   // Eavail <= 75, pl 625 - 850
       {250, 350, 400, 450, 500, 550, 600, 700, 850},   // Eavail <= 75, pl 850 - 1200
       {300, 350, 450, 550, 600, 700, 800, 950, 1200}   // Eavail <= 75, pl > 1200
      },
      {
       {150, 200, 250, 300, 400, 500},   // Eavail 75 - 150, pl <= 150
       {150, 250, 350, 450},   // Eavail 75 - 150, pl 150 - 300
       {150, 250, 300, 450},   // Eavail 75 - 150, pl 300 - 450
       {250, 350, 450, 550},   // Eavail 75 - 150, pl 450 - 625
       {300, 450, 550},   // Eavail 75 - 150, pl 625 - 850
       {300, 400, 500, 600},   // Eavail 75 - 150, pl 850 - 1200
       {450, 650, 900}   // Eavail 75 - 150, pl > 1200
      },
      {
       {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700},   // Eavail > 150, pl <= 150
       {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 650},   // Eavail > 150, pl 150 - 300
       {100, 150, 200, 250, 300, 350, 400, 450, 500, 600},   // Eavail > 150, pl 300 - 450
       {150, 200, 250, 300, 350, 400, 450, 500, 550, 650},   // Eavail > 150, pl 450 - 625
       {200, 300, 350, 400, 450, 500, 600, 700},   // Eavail > 150, pl 625 - 850
       {250, 300, 400, 450, 500, 550, 650, 750, 950},   // Eavail > 150, pl 850 - 1200
       {350, 450, 550, 650, 750, 850, 1050, 1400}   // Eavail > 150, pl > 1200
      }
    };
    static const std::vector<std::vector<std::vector<double>>> pt_bins_3d_Np = {
      {
       {200, 300},   // Eavail <= 150, pl <= 0
       {100, 150, 200, 250, 300, 400},   // Eavail <= 150, pl 0 - 150
       {150, 250, 350, 400, 450},   // Eavail <= 150, pl 150 - 300
       {200, 250, 300, 350, 400, 450},   // Eavail <= 150, pl 300 - 450
       {200, 300, 350, 400, 450, 550},   // Eavail <= 150, pl 450 - 625
       {200, 250, 300, 350, 400, 450, 500, 600},   // Eavail <= 150, pl 625 - 850
       {300, 350, 450, 500, 550, 650},   // Eavail <= 150, pl 850 - 1200
       {350, 450, 550, 650, 700, 850}   // Eavail <= 150, pl > 1200
      },
      {
       {150, 200, 250, 300, 400},   // Eavail 150 - 250, pl <= 0
       {100, 150, 200, 250, 300, 350, 400},   // Eavail 150 - 250, pl 0 - 150
       {100, 150, 200, 250, 300, 350, 400, 450, 500},   // Eavail 150 - 250, pl 150 - 300
       {150, 200, 250, 300, 350, 400, 450, 500, 550},   // Eavail 150 - 250, pl 300 - 450
       {200, 250, 300, 350, 400, 450, 500, 550, 650},   // Eavail 150 - 250, pl 450 - 625
       {200, 300, 350, 400, 450, 500, 550, 650},   // Eavail 150 - 250, pl 625 - 850
       {200, 250, 300, 350, 400, 450, 500, 550, 600, 750},   // Eavail 150 - 250, pl 850 - 1200
       {350, 450, 500, 550, 650, 700, 800, 1000, 1250}   // Eavail 150 - 250, pl > 1200
      },
      {
       {150, 200, 250, 300, 350},   // Eavail 250 - 375, pl <= 0
       {100, 150, 200, 250, 300, 350, 400, 500},   // Eavail 250 - 375, pl 0 - 150
       {150, 200, 250, 300, 350, 400, 450, 500, 550},   // Eavail 250 - 375, pl 150 - 300
       {150, 200, 300, 350, 400, 450, 500, 600},   // Eavail 250 - 375, pl 300 - 450
       {200, 250, 300, 350, 400, 450, 500, 550, 650},   // Eavail 250 - 375, pl 450 - 625
       {200, 250, 300, 350, 400, 450, 500, 600, 700},   // Eavail 250 - 375, pl 625 - 850
       {250, 350, 400, 450, 500, 550, 650, 800},   // Eavail 250 - 375, pl 850 - 1200
       {250, 450, 550, 650, 750, 950}   // Eavail 250 - 375, pl > 1200
      },
      {
       {100, 150, 200, 250, 300, 400},   // Eavail 375 - 550, pl <= 0
       {100, 150, 200, 250, 300, 350, 400, 500},   // Eavail 375 - 550, pl 0 - 150
       {100, 150, 200, 250, 300, 350, 450, 550},   // Eavail 375 - 550, pl 150 - 300
       {150, 250, 350, 450, 550},   // Eavail 375 - 550, pl 300 - 450
       {250, 350, 450, 550},   // Eavail 375 - 550, pl 450 - 625
       {250, 350, 450, 500, 600},   // Eavail 375 - 550, pl 625 - 850
       {300, 400, 550, 650},   // Eavail 375 - 550, pl 850 - 1200
       {450, 550, 750, 1000}   // Eavail 375 - 550, pl > 1200
      },
      {
       {100, 150, 200, 250, 300, 350, 450},   // Eavail > 550, pl <= 0
       {50, 100, 150, 200, 250, 300, 350, 450},   // Eavail > 550, pl 0 - 150
       {100, 150, 200, 250, 300, 350, 450},   // Eavail > 550, pl 150 - 300
       {150, 250, 350, 450, 550},   // Eavail > 550, pl 300 - 450
       {200, 300, 450, 550},   // Eavail > 550, pl 450 - 625
       {250, 400, 500},   // Eavail > 550, pl 625 - 850
       {300, 450, 650},   // Eavail > 550, pl 850 - 1200
       {650, 950}   // Eavail > 550, pl > 1200
      }
    };
    // "_bin_t" variants: the same reco values in the truth binning of get_xs_signal_no (LEEana::xs), nbin as the truth bins
    bool flag_truth_bins = var_name.EndsWith("_bin_t");
    TString name = flag_truth_bins ? TString(var_name(0, var_name.Length()-2)) : var_name;
    if (pfeval.reco_muonMomentum[3]<0) return (name == "muon_pt" || name == "muon_pl") ? -1000 : -1;
    TVector3 muon_p = get_muon_momentum_new(pfeval, eval);
    if (name == "muon_pt") return muon_p.Perp();
    if (name == "muon_pl") return muon_p.Z();
    if (name == "pl_pt_bin") return get_2d_bin_index(muon_p.Z(), xs::pl_slices, muon_p.Perp(), flag_truth_bins ? xs::truth_pt_bins_pl : pt_bins_pl) + 0.5;
    double Emu = get_muon_Etot_new(pfeval, eval);
    TVector3 mu_dir(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2]);
    if (name == "Emu_costheta_bin") return get_2d_bin_index(mu_dir.CosTheta(), xs::costheta_slices, Emu, flag_truth_bins ? xs::truth_Emu_bins_costheta : Emu_bins_costheta) + 0.5;
    double Eavail = get_kine_var(kine, eval, pfeval, tagger, flag_data, "Eavail_new3_5_drop15_95", space, pandora, lantern);
    if (name == "Eavail_Emu_bin") return get_2d_bin_index(Eavail, xs::Eavail_slices, Emu, flag_truth_bins ? xs::truth_Emu_bins_Eavail : Emu_bins_Eavail) + 0.5;
    double KE_lead = 0;
    int n_protons = 0;
    get_reco_leading_proton(pfeval, space, 45, KE_lead, n_protons);
    bool flag_reco_Np = KE_lead >= 45;
    if (name == "Eavail_costheta_Emu_bin"){
      if (flag_reco_Np) return get_3d_bin_index(Eavail, xs::Eavail_outer_Np, mu_dir.CosTheta(), xs::costheta_inner_Np, Emu, flag_truth_bins ? xs::truth_Emu_bins_3d_Np : Emu_bins_3d_Np) + 0.5;
      return get_3d_bin_index(Eavail, xs::Eavail_outer_0p, mu_dir.CosTheta(), xs::costheta_inner_0p, Emu, flag_truth_bins ? xs::truth_Emu_bins_3d_0p : Emu_bins_3d_0p) + 0.5;
    }
    if (flag_reco_Np) return get_3d_bin_index(Eavail, xs::Eavail_outer_Np, muon_p.Z(), xs::pl_inner_Np, muon_p.Perp(), flag_truth_bins ? xs::truth_pt_bins_3d_Np : pt_bins_3d_Np) + 0.5;
    return get_3d_bin_index(Eavail, xs::Eavail_outer_0p, muon_p.Z(), xs::pl_inner_0p, muon_p.Perp(), flag_truth_bins ? xs::truth_pt_bins_3d_0p : pt_bins_3d_0p) + 0.5;
  }else if (var_name == "reco_Np" || var_name == "reco_Kp" || var_name == "reco_costhetap" || var_name == "reco_costhetamup"
         || var_name == "Kp_costhetap_bin" || var_name == "Kp_costhetamup_bin" || var_name == "Eavail_costhetamup_bin" || var_name == "Emu_costhetamup_bin"
         || var_name == "costheta_Emu_Kp_bin" || var_name == "costheta_Emu_costhetamup_bin" || var_name == "pl_pt_Kp_bin" || var_name == "pl_pt_costhetamup_bin"
         || var_name == "Kp_costhetap_bin_t" || var_name == "Kp_costhetamup_bin_t" || var_name == "Eavail_costhetamup_bin_t" || var_name == "Emu_costhetamup_bin_t"
         || var_name == "costheta_Emu_Kp_bin_t" || var_name == "costheta_Emu_costhetamup_bin_t" || var_name == "pl_pt_Kp_bin_t" || var_name == "pl_pt_costhetamup_bin_t"){
    // Reco variables of the proton cross-section measurements. The protons are the primary WireCell protons with range KE
    // (get_reco_leading_proton); the leading one (largest range KE) sets the reco Np category of get_particle_0pNp_bdt_bin.
    // Its direction is from its end points, pointing away from the vertex (get_reco_proton_dir; better than
    // reco_startMomentum).
    //   reco_Np: number of protons with range KE >= 45 MeV (truth bins: cut_file 7, 8)
    //   reco_Kp: leading proton KE [MeV], 0 without a proton (cut_file 9)
    //   reco_costhetap, reco_costhetamup: leading proton cos(theta), cos(muon, proton) (cut_file 10, 11)
    //   "*_bin": the reco bin of a multi-differential measurement as the flattened bin index + 0.5, for histograms from 0
    //   to nbin. Same slices as the truth bins of get_xs_signal_no; finer bins (Kp: 25 MeV, cos(mu, p): 0.1), each a part
    //   of one truth bin. Muon energy / momentum as in kine_reco_Enu_new3_5, Eavail as Eavail_new3_5_drop15_95.
    //     Kp_costhetap_bin (cut_file 13): nbin = 95       Kp_costhetamup_bin (18): 98
    //     Eavail_costhetamup_bin (19): 114                   Emu_costhetamup_bin (20): 113
    //     costheta_Emu_Kp_bin (23): 209                      costheta_Emu_costhetamup_bin (24): 240
    //     pl_pt_Kp_bin (25): 228                             pl_pt_costhetamup_bin (26): 251
    //   "*_bin_t": the same in the truth binning (for channels without the statistics for the finer bins), nbin =
    //     Kp_costhetap_bin_t 58, Kp_costhetamup_bin_t 38, Eavail_costhetamup_bin_t 40, Emu_costhetamup_bin_t 42,
    //     costheta_Emu_Kp_bin_t 101, costheta_Emu_costhetamup_bin_t 83, pl_pt_Kp_bin_t 95, pl_pt_costhetamup_bin_t 88.
    //   The angles and the bins are -2 without a proton above 45 MeV (reco 0p, which goes to the channels without a reco
    //   proton).
    static const std::vector<std::vector<double>> Kp_bins_costhetap = {
      {70, 95, 120, 145, 170, 195, 220},   // costhetap <= 0
      {70, 95, 120, 145, 170, 195, 220, 245, 270},   // costhetap 0 - 0.3
      {70, 95, 120, 145, 170, 195, 220, 245, 270},   // costhetap 0.3 - 0.5
      {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370},   // costhetap 0.5 - 0.7
      {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370, 420, 470},   // costhetap 0.7 - 0.8
      {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370, 395, 420, 470},   // costhetap 0.8 - 0.9
      {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370, 395, 420, 445, 470, 520, 595}   // costhetap > 0.9
    };
    static const std::vector<std::vector<double>> costhetamup_bins_Kp = {
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Kp <= 95
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Kp 95 - 145
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Kp 145 - 220
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Kp 220 - 345
      {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}   // Kp > 345
    };
    static const std::vector<std::vector<double>> costhetamup_bins_Eavail = {
      {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // Eavail <= 150
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Eavail 150 - 250
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Eavail 250 - 375
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Eavail 375 - 550
      {-0.9, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Eavail 550 - 800
      {-0.9, -0.7, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}   // Eavail > 800
    };
    static const std::vector<std::vector<double>> costhetamup_bins_Emu = {
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Emu <= 350
      {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Emu 350 - 550
      {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Emu 550 - 700
      {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Emu 700 - 900
      {-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // Emu 900 - 1150
      {-0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}   // Emu > 1150
    };
    static const std::vector<std::vector<std::vector<double>>> Kp_bins_costheta_Emu = {
      {
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370, 395, 445, 495},   // costheta <= 0.5, Emu <= 350
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 345, 395, 470},   // costheta <= 0.5, Emu 350 - 450
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 370, 420, 470}   // costheta <= 0.5, Emu > 450
      },
      {
       {70, 95, 120, 145, 170, 195, 220, 245, 295, 370},   // costheta 0.5 - 0.7, Emu <= 350
       {70, 95, 120, 170, 220, 270},   // costheta 0.5 - 0.7, Emu 350 - 450
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 345},   // costheta 0.5 - 0.7, Emu 450 - 700
       {70, 95, 120, 145, 170, 195, 220, 245, 295, 345, 470}   // costheta 0.5 - 0.7, Emu > 700
      },
      {
       {70, 95, 120, 145, 170, 195, 220, 270},   // costheta 0.7 - 0.8, Emu <= 450
       {95, 120, 170, 220},   // costheta 0.7 - 0.8, Emu 450 - 550
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 320},   // costheta 0.7 - 0.8, Emu 550 - 900
       {70, 95, 120, 170, 220, 270, 345}   // costheta 0.7 - 0.8, Emu > 900
      },
      {
       {70, 95, 120, 145, 170, 195, 220, 270},   // costheta 0.8 - 0.9, Emu <= 450
       {70, 120, 170, 220},   // costheta 0.8 - 0.9, Emu 450 - 550
       {70, 95, 120, 145, 170, 195, 220, 245, 270},   // costheta 0.8 - 0.9, Emu 550 - 900
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 345}   // costheta 0.8 - 0.9, Emu > 900
      },
      {
       {70, 95, 120, 145, 170, 195, 220, 245, 295},   // costheta > 0.9, Emu <= 550
       {70, 95, 120, 145, 170, 220},   // costheta > 0.9, Emu 550 - 700
       {70, 95, 120, 145, 170, 220, 245},   // costheta > 0.9, Emu 700 - 900
       {70, 95, 120, 145, 170, 195, 220, 245},   // costheta > 0.9, Emu 900 - 1150
       {70, 95, 120, 145, 170, 195, 220, 270},   // costheta > 0.9, Emu 1150 - 1600
       {70, 95, 120, 145, 170, 220, 245, 345}   // costheta > 0.9, Emu > 1600
      }
    };
    static const std::vector<std::vector<std::vector<double>>> costhetamup_bins_costheta_Emu = {
      {
       {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7},   // costheta <= 0.5, Emu <= 350
       {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.3, 0.5},   // costheta <= 0.5, Emu 350 - 450
       {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.3, 0.4, 0.6}   // costheta <= 0.5, Emu > 450
      },
      {
       {-0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8},   // costheta 0.5 - 0.7, Emu <= 350
       {-0.5, -0.2, 0, 0.2, 0.4, 0.6},   // costheta 0.5 - 0.7, Emu 350 - 450
       {-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5, 0.7},   // costheta 0.5 - 0.7, Emu 450 - 700
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6}   // costheta 0.5 - 0.7, Emu > 700
      },
      {
       {-0.5, -0.3, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8},   // costheta 0.7 - 0.8, Emu <= 450
       {-0.3, 0, 0.2, 0.4},   // costheta 0.7 - 0.8, Emu 450 - 550
       {-0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.6},   // costheta 0.7 - 0.8, Emu 550 - 900
       {-0.5, -0.3, -0.1, 0, 0.1, 0.2, 0.4}   // costheta 0.7 - 0.8, Emu > 900
      },
      {
       {-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // costheta 0.8 - 0.9, Emu <= 450
       {-0.3, 0, 0.2, 0.5, 0.7},   // costheta 0.8 - 0.9, Emu 450 - 550
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // costheta 0.8 - 0.9, Emu 550 - 900
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}   // costheta 0.8 - 0.9, Emu > 900
      },
      {
       {-0.6, -0.3, -0.1, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},   // costheta > 0.9, Emu <= 550
       {-0.4, -0.2, 0.2, 0.4, 0.5, 0.7},   // costheta > 0.9, Emu 550 - 700
       {-0.5, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.8},   // costheta > 0.9, Emu 700 - 900
       {-0.6, -0.4, -0.2, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // costheta > 0.9, Emu 900 - 1150
       {-0.6, -0.4, -0.2, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // costheta > 0.9, Emu 1150 - 1600
       {-0.5, -0.3, -0.2, 0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8}   // costheta > 0.9, Emu > 1600
      }
    };
    static const std::vector<std::vector<std::vector<double>>> Kp_bins_pl_pt = {
      {
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 345, 395, 445},   // pl <= 150, pt <= 150
       {70, 95, 120, 145, 170, 195, 220, 245, 295, 345},   // pl <= 150, pt 150 - 200
       {70, 95, 120, 145, 170, 195, 220, 245, 295, 345},   // pl <= 150, pt 200 - 250
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 345, 395},   // pl <= 150, pt 250 - 350
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 345, 370, 420, 470}   // pl <= 150, pt > 350
      },
      {
       {70, 95, 120, 145, 170, 195, 220, 270},   // pl 150 - 300, pt <= 200
       {70, 95, 120, 145, 170, 195, 220, 270},   // pl 150 - 300, pt 200 - 300
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 320},   // pl 150 - 300, pt 300 - 450
       {70, 95, 120, 145, 170, 220, 270, 320, 470}   // pl 150 - 300, pt > 450
      },
      {
       {70, 95, 120, 145, 170, 195, 220},   // pl 300 - 450, pt <= 200
       {70, 95, 120, 145, 170, 195, 220},   // pl 300 - 450, pt 200 - 300
       {70, 95, 120, 145, 170, 195, 220, 245, 270},   // pl 300 - 450, pt 300 - 450
       {70, 95, 120, 145, 170, 195, 220, 245, 295, 345}   // pl 300 - 450, pt > 450
      },
      {
       {70, 95, 120, 145, 170, 220},   // pl 450 - 625, pt <= 250
       {70, 95, 120, 145, 170, 195, 220, 245, 295},   // pl 450 - 625, pt 250 - 450
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 345}   // pl 450 - 625, pt > 450
      },
      {
       {70, 95, 120, 170, 220},   // pl 625 - 850, pt <= 250
       {70, 95, 120, 145, 170, 195, 220, 245, 270},   // pl 625 - 850, pt 250 - 450
       {70, 95, 120, 145, 170, 195, 220, 270, 345}   // pl 625 - 850, pt > 450
      },
      {
       {70, 95, 120, 145, 170, 195},   // pl > 850, pt <= 250
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295},   // pl > 850, pt 250 - 550
       {70, 95, 120, 145, 170, 195, 220, 245, 270, 295, 320, 370, 470}   // pl > 850, pt > 550
      }
    };
    static const std::vector<std::vector<std::vector<double>>> costhetamup_bins_pl_pt = {
      {
       {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // pl <= 150, pt <= 150
       {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.2, 0.3, 0.6},   // pl <= 150, pt 150 - 200
       {-0.8, -0.7, -0.6, -0.5, -0.3, -0.2, 0, 0.2, 0.4},   // pl <= 150, pt 200 - 250
       {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.2, 0.4, 0.6},   // pl <= 150, pt 250 - 350
       {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.3, 0.5}   // pl <= 150, pt > 350
      },
      {
       {-0.6, -0.3, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // pl 150 - 300, pt <= 200
       {-0.6, -0.3, -0.1, 0, 0.2, 0.3, 0.4, 0.5, 0.7},   // pl 150 - 300, pt 200 - 300
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.5, 0.7},   // pl 150 - 300, pt 300 - 450
       {-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.2}   // pl 150 - 300, pt > 450
      },
      {
       {-0.4, 0, 0.3, 0.5, 0.6, 0.8},   // pl 300 - 450, pt <= 200
       {-0.4, -0.2, 0, 0.2, 0.4, 0.6},   // pl 300 - 450, pt 200 - 300
       {-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7},   // pl 300 - 450, pt 300 - 450
       {-0.5, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6}   // pl 300 - 450, pt > 450
      },
      {
       {-0.4, -0.1, 0.2, 0.4, 0.6, 0.7},   // pl 450 - 625, pt <= 250
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // pl 450 - 625, pt 250 - 450
       {-0.7, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6}   // pl 450 - 625, pt > 450
      },
      {
       {-0.2, 0.2, 0.4, 0.6},   // pl 625 - 850, pt <= 250
       {-0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // pl 625 - 850, pt 250 - 450
       {-0.6, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}   // pl 625 - 850, pt > 450
      },
      {
       {-0.2, 0.2, 0.4, 0.6, 0.8},   // pl > 850, pt <= 250
       {-0.7, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8},   // pl > 850, pt 250 - 550
       {-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}   // pl > 850, pt > 550
      }
    };
    // "_bin_t" variants: the same reco values in the truth binning of get_xs_signal_no (LEEana::xs), nbin as the truth bins
    bool flag_truth_bins = var_name.EndsWith("_bin_t");
    TString name = flag_truth_bins ? TString(var_name(0, var_name.Length()-2)) : var_name;
    double KE_lead = 0;
    int n_protons = 0;
    int index = get_reco_leading_proton(pfeval, space, 45, KE_lead, n_protons);
    if (name == "reco_Np") return n_protons;
    if (name == "reco_Kp") return KE_lead;
    if (index<0 || KE_lead<45) return -2;
    TVector3 p_dir = get_reco_proton_dir(pfeval, index);
    if (name == "reco_costhetap") return (p_dir.Mag()==0) ? -2 : p_dir.CosTheta();
    if (name == "reco_costhetamup") return (p_dir.Mag()==0 || pfeval.reco_muonMomentum[3]<=0) ? -2 : get_reco_cos_mu_p(pfeval, index);
    if (name == "Kp_costhetap_bin") return get_2d_bin_index(p_dir.CosTheta(), xs::costheta_slices, KE_lead, flag_truth_bins ? xs::truth_Kp_bins_costhetap : Kp_bins_costhetap) + 0.5;
    double cos_mu_p = get_reco_cos_mu_p(pfeval, index);
    if (name == "Kp_costhetamup_bin") return get_2d_bin_index(KE_lead, xs::Kp_slices_mup, cos_mu_p, flag_truth_bins ? xs::truth_costhetamup_bins_Kp : costhetamup_bins_Kp) + 0.5;
    if (name == "Eavail_costhetamup_bin"){
      double Eavail = get_kine_var(kine, eval, pfeval, tagger, flag_data, "Eavail_new3_5_drop15_95", space, pandora, lantern);
      return get_2d_bin_index(Eavail, xs::Eavail_slices_mup, cos_mu_p, flag_truth_bins ? xs::truth_costhetamup_bins_Eavail : costhetamup_bins_Eavail) + 0.5;
    }
    if (name == "pl_pt_Kp_bin" || name == "pl_pt_costhetamup_bin"){
      TVector3 muon_p = get_muon_momentum_new(pfeval, eval);
      if (name == "pl_pt_Kp_bin") return get_3d_bin_index(muon_p.Z(), xs::pl_outer_p, muon_p.Perp(), xs::pt_inner_p, KE_lead, flag_truth_bins ? xs::truth_Kp_bins_pl_pt : Kp_bins_pl_pt) + 0.5;
      return get_3d_bin_index(muon_p.Z(), xs::pl_outer_p, muon_p.Perp(), xs::pt_inner_p, cos_mu_p, flag_truth_bins ? xs::truth_costhetamup_bins_pl_pt : costhetamup_bins_pl_pt) + 0.5;
    }
    double Emu = get_muon_Etot_new(pfeval, eval);
    if (name == "Emu_costhetamup_bin") return get_2d_bin_index(Emu, xs::Emu_slices_mup, cos_mu_p, flag_truth_bins ? xs::truth_costhetamup_bins_Emu : costhetamup_bins_Emu) + 0.5;
    TVector3 mu_dir(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2]);
    if (name == "costheta_Emu_Kp_bin") return get_3d_bin_index(mu_dir.CosTheta(), xs::costheta_outer_p, Emu, xs::Emu_inner_p, KE_lead, flag_truth_bins ? xs::truth_Kp_bins_costheta_Emu : Kp_bins_costheta_Emu) + 0.5;
    return get_3d_bin_index(mu_dir.CosTheta(), xs::costheta_outer_p, Emu, xs::Emu_inner_p, cos_mu_p, flag_truth_bins ? xs::truth_costhetamup_bins_costheta_Emu : costhetamup_bins_costheta_Emu) + 0.5;
  }else if(var_name == "all_veto_score"){
    return tagger.all_veto_score;
  }else if(var_name == "VtxAct_bdt_score"){
  return tagger.VtxAct_bdt_score;
  }else if(var_name.BeginsWith("muon_spacepoints_q_")){
    // primary muon spacepoint charges (main sequence, ordered as in particle.h), raw charge: the particle.h BDT inputs are
    // (q+10000)*10. muon_spacepoints_q_0 ... _4: the first five, _sum5: their sum, _med: the median (dQ/dx); -999 if missing
    std::vector<double> q = get_reco_muon_spacepoints_q(pfeval, space);
    if(var_name == "muon_spacepoints_q_med"){
      if(q.empty()) return -999;
      std::sort(q.begin(), q.end());
      size_t size = q.size();
      return (size % 2 == 0) ? (q.at(size / 2 - 1) + q.at(size / 2)) / 2.0 : q.at(size / 2);
    }
    if(var_name == "muon_spacepoints_q_sum5"){
      if(q.size()<5) return -999;
      return q.at(0) + q.at(1) + q.at(2) + q.at(3) + q.at(4);
    }
    TString k_str = var_name(19, var_name.Length()-19);
    if(!k_str.IsDigit()) return -999;
    size_t k = k_str.Atoi();
    return (k<q.size()) ? q.at(k) : -999;


  }else if (var_name == "kine_reco_Eproton"){
    return get_reco_Eproton(kine);
  }else if (var_name == "kine_reco_Eproton_nothreshold"){
    double reco_Eproton=0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
        int pdgcode = kine.kine_particle_type->at(i);
        //if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) // proton threshold of 35 MeV
        if(abs(pdgcode)==2212) // no proton threshold
          reco_Eproton+=kine.kine_energy_particle->at(i);

      }
    return reco_Eproton;

  }else if (var_name == "pi0_energy"){
    double pi0_mass = 135;
    double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
    return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);

  }else if (var_name == "nue_score"){
    return (tagger.nue_score<=15.99?tagger.nue_score:15.99);
  }else if (var_name == "nc_pio_score"){
    return tagger.nc_pio_score;
  }else if (var_name == "nc_delta_score"){
    return tagger.nc_delta_score;
  }else if (var_name == "numu_score"){
    return tagger.numu_score;

  }else if (var_name == "reco_nuvtxX"){
      return pfeval.reco_nuvtxX;
  }else if (var_name == "reco_nuvtxY"){
      return pfeval.reco_nuvtxY;
  }else if (var_name == "reco_nuvtxZ"){
      return pfeval.reco_nuvtxZ;
  }else if (var_name == "reco_nuvtxU"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) - pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "reco_nuvtxV"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) + pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);

  }else if (var_name == "muon_KE"){
      return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "reco_Emuon"){
      return pfeval.reco_muonMomentum[3]*1000; // GeV --> MeV
  }else if (var_name == "muon_momentum"){
      if (pfeval.reco_muonMomentum[3] < 0) { return -1; }
      float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
      return (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));
  }else if (var_name == "muon_costheta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return TMath::Cos(muonMomentum.Theta());
      else
	return -2;
  }else if (var_name == "muon_theta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return muonMomentum.Theta()*180./TMath::Pi();
      else
	return -1000;
  }else if (var_name == "muon_phi"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return muonMomentum.Phi()/TMath::Pi()*180.;
      else
	return -1000;

  }else if (var_name == "Ehadron"){
    if (pfeval.reco_muonMomentum[3]>0)
      return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    else
      return -1000;
   
  //
  }else{
    std::cout << "No such variable: " << var_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

int get_costheta_bin (float costh) {
  int nbins = 9;
  float costheta_binning[nbins+1] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};
  if (costh == costheta_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (costh >  costheta_binning[i] && costh <= costheta_binning[i+1]) { return i; } }
  return -1;
}


int get_Enu_bin (float Enu) {
  int nbins = 4;
  float enu_binning[nbins+1] = {200, 705, 1050, 1570, 4000};
  if (Enu == enu_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (Enu >  enu_binning[i] && Enu <= enu_binning[i+1]) { return i; } }
  return -1;
}

int LEEana::get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  // every signal definition below requires a true CC interaction inside the active volume
  if (!(eval.truth_isCC==1 && eval.truth_vtxInside==1)) return -1;

  // truth kinematics
  double Emuon   = pfeval.truth_muonMomentum[3]*1000; // MeV
  double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV
  TLorentzVector muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);
  float costh = TMath::Cos(muonMomentum.Theta());

  // numuCC 0p/Np truth bins in one variable (cut_file 2-6). Signal: numu CC with the vertex in the FV, split into 0p/Np by
  // the leading primary proton KE (45 MeV), as XsnumuCCinFV0p/Np. The bins cover the whole range of each variable, so
  // every signal event gets a bin. Bin names: "numuCC.<0p|Np>.inside.<var>" + ".le.<e1>", ".le.<e(k+1)>.gt.<ek>", ".gt.<en>"
  // (inner edges ek below), e.g. "numuCC.Np.inside.Emu.le.250.gt.200". Edges follow the resolution of each variable
  // (bin width about the 68% half-width of reco - true) and the expected statistics at ~1e21 POT.
  //   2: Emu, muon total energy [MeV]
  //   3: costheta, muon cos(theta)
  //   4: nu, energy transfer Enu - Emu [MeV]; coarser for 0p, whose hadronic energy is mostly not reconstructed
  //   5: Eavail, available energy (get_true_Eavail, no thresholds, no masses) [MeV]
  //   6: Enu, neutrino energy [MeV]
  // Proton content (cut_file 7-11), counting primary protons with KE >= 45 MeV; the leading proton is the primary proton
  // with the largest KE. In 9-11 all 0p events share one bin, "numuCC.0p.inside".
  //   7: 0p / Np: "numuCC.0p.inside", "numuCC.Np.inside"
  //   8: 0p / 1p / 2p / >2p: "numuCC.0p.inside", "numuCC.1p.inside", "numuCC.2p.inside", "numuCC.gt2p.inside"
  //   9: Kp, leading proton KE [MeV] ("numuCC.Np.inside.Kp.le.70" is 45 < Kp <= 70)
  //  10: costhetap, leading proton cos(theta)
  //  11: costhetamup, cosine of the opening angle between the muon and the leading proton
  // Double differential (cut_file 12), 0p/Np as in 2-6: slices in the outer variable, each with its own bins in the inner
  // variable (more bins where the statistics and resolution allow, the same for 0p and Np). Bin names are the slice name
  // followed by the bin name, e.g. "numuCC.Np.inside.costheta.le.0.3.gt.0.Emu.le.300.gt.250".
  //  12: costheta slices x Emu bins
  //  13: costhetap slices x Kp bins of the leading proton (Np only; all 0p events in one bin, "numuCC.0p.inside")
  // Muon transverse / longitudinal momentum (w.r.t. the beam, z) [MeV], 0p/Np as in 2-6:
  //  14: pt = p sin(theta)
  //  15: pl = p cos(theta)
  //  16: pl slices x pt bins
  //  17: Eavail slices x Emu bins, 0p/Np
  //  18: Kp slices x costhetamup bins (Np only; all 0p events in one bin, "numuCC.0p.inside")
  //  19: Eavail slices x costhetamup bins (Np only; all 0p events in one bin)
  //  20: Emu slices x costhetamup bins (Np only; all 0p events in one bin)
  // Triple differential: outer slices, inner slices, bins (helpers get_xs_3d_bin_name / get_3d_bin_index). 0p and Np
  // have different binnings: true 0p events are concentrated at low Eavail (~2/3 at Eavail <= 150 MeV), Np events extend
  // to high Eavail, so common binnings would starve one of them. Inner slices are the same in all outer slices.
  //  21: Eavail (outer) x costheta (inner) x Emu
  //  22: Eavail (outer) x pl (inner) x pt
  // Proton triple differential (Np only; all 0p events in one bin, "numuCC.0p.inside"). The inner slices follow the
  // statistics in each outer slice (e.g. more energetic muons in the forward slices), the same for 23/24 and for 25/26.
  //  23: costheta (outer) x Emu (inner) x Kp
  //  24: costheta (outer) x Emu (inner) x costhetamup
  //  25: pl (outer) x pt (inner) x Kp
  //  26: pl (outer) x pt (inner) x costhetamup
  // The edges of all binnings come from one grid per variable (the 1D edges; pt also 50 MeV steps), so the slices and bins
  // of different measurements line up; the slices shared with the reco binnings are in LEEana::xs.
  if (cut_file>=2 && cut_file<=26){
    if (eval.truth_nuPdg!=14) return -1;
    // 1D
    static const std::vector<double> Emu_edges = {200, 250, 300, 350, 450, 550, 700, 900, 1150, 1600};
    static const std::vector<double> costheta_edges = {-0.3, 0, 0.15, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
    static const std::vector<double> nu_edges_0p = {100, 200, 500, 950};
    static const std::vector<double> nu_edges_Np = {100, 200, 325, 500, 700, 950, 1300, 1700};
    static const std::vector<double> Eavail_edges = {75, 150, 250, 375, 550, 800, 1100};
    static const std::vector<double> Enu_edges = {400, 500, 650, 800, 1000, 1300, 1700, 2300};
    static const std::vector<double> Kp_edges = {70, 95, 120, 145, 170, 220, 270, 345, 470};
    static const std::vector<double> costhetap_edges = {-0.5, -0.2, 0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    static const std::vector<double> costhetamup_edges = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
    static const std::vector<double> pt_edges = {100, 200, 300, 450, 650};
    static const std::vector<double> pl_edges = {-150, 0, 75, 150, 300, 450, 625, 850, 1200, 1700};
    // multi-differential: bins in each slice (bins and slices in LEEana::xs, shared with the "_bin_t" reco variables)
    const std::vector<std::vector<double>>& Emu_bins_costheta = xs::truth_Emu_bins_costheta;
    const std::vector<std::vector<double>>& Kp_bins_costhetap = xs::truth_Kp_bins_costhetap;
    const std::vector<std::vector<double>>& pt_bins_pl = xs::truth_pt_bins_pl;
    const std::vector<std::vector<double>>& Emu_bins_Eavail = xs::truth_Emu_bins_Eavail;
    const std::vector<std::vector<double>>& costhetamup_bins_Kp = xs::truth_costhetamup_bins_Kp;
    const std::vector<std::vector<double>>& costhetamup_bins_Eavail = xs::truth_costhetamup_bins_Eavail;
    const std::vector<std::vector<double>>& costhetamup_bins_Emu = xs::truth_costhetamup_bins_Emu;
    const std::vector<std::vector<std::vector<double>>>& Emu_bins_3d_0p = xs::truth_Emu_bins_3d_0p;
    const std::vector<std::vector<std::vector<double>>>& Emu_bins_3d_Np = xs::truth_Emu_bins_3d_Np;
    const std::vector<std::vector<std::vector<double>>>& pt_bins_3d_0p = xs::truth_pt_bins_3d_0p;
    const std::vector<std::vector<std::vector<double>>>& pt_bins_3d_Np = xs::truth_pt_bins_3d_Np;
    const std::vector<std::vector<std::vector<double>>>& Kp_bins_costheta_Emu = xs::truth_Kp_bins_costheta_Emu;
    const std::vector<std::vector<std::vector<double>>>& costhetamup_bins_costheta_Emu = xs::truth_costhetamup_bins_costheta_Emu;
    const std::vector<std::vector<std::vector<double>>>& Kp_bins_pl_pt = xs::truth_Kp_bins_pl_pt;
    const std::vector<std::vector<std::vector<double>>>& costhetamup_bins_pl_pt = xs::truth_costhetamup_bins_pl_pt;
    // leading primary proton (as get_KE(pfeval, 2212, 1, 0, 1, 0)) and the number of primary protons above 45 MeV
    int lead_p = -1;
    double Kp = 0;
    int n_p = 0;
    for (int i=0; i<pfeval.truth_Ntrack; i++){
      if (pfeval.truth_mother[i]!=0 || pfeval.truth_pdg[i]!=2212) continue;
      double KE = pfeval.truth_startMomentum[i][3]*1000 - get_mass_MeV(2212);
      if (KE>=45) n_p++;
      if (KE>Kp){
        Kp = KE;
        lead_p = i;
      }
    }
    bool flag_0p = Kp < 45;
    TString prefix = flag_0p ? "numuCC.0p.inside." : "numuCC.Np.inside.";
    TVector3 muon_p(pfeval.truth_muonMomentum[0]*1000, pfeval.truth_muonMomentum[1]*1000, pfeval.truth_muonMomentum[2]*1000);
    double Eavail_true = (cut_file==5 || cut_file==17 || cut_file==19 || cut_file==21 || cut_file==22) ? get_true_Eavail(eval, pfeval, false, false) : 0;
    TString bin_name;
    // 0p and Np bins
    if (cut_file==2) bin_name = get_xs_bin_name(prefix+"Emu", Emuon, Emu_edges);
    else if (cut_file==3) bin_name = get_xs_bin_name(prefix+"costheta", costh, costheta_edges);
    else if (cut_file==4) bin_name = get_xs_bin_name(prefix+"nu", Ehadron, flag_0p ? nu_edges_0p : nu_edges_Np);
    else if (cut_file==5) bin_name = get_xs_bin_name(prefix+"Eavail", Eavail_true, Eavail_edges);
    else if (cut_file==6) bin_name = get_xs_bin_name(prefix+"Enu", eval.truth_nuEnergy, Enu_edges);
    else if (cut_file==7) bin_name = flag_0p ? "numuCC.0p.inside" : "numuCC.Np.inside";
    else if (cut_file==8) bin_name = (n_p==0) ? "numuCC.0p.inside" : (n_p==1) ? "numuCC.1p.inside" : (n_p==2) ? "numuCC.2p.inside" : "numuCC.gt2p.inside";
    else if (cut_file==12) bin_name = get_xs_2d_bin_name(prefix, "costheta", costh, xs::costheta_slices, "Emu", Emuon, Emu_bins_costheta);
    else if (cut_file==14) bin_name = get_xs_bin_name(prefix+"pt", muon_p.Perp(), pt_edges);
    else if (cut_file==15) bin_name = get_xs_bin_name(prefix+"pl", muon_p.Z(), pl_edges);
    else if (cut_file==16) bin_name = get_xs_2d_bin_name(prefix, "pl", muon_p.Z(), xs::pl_slices, "pt", muon_p.Perp(), pt_bins_pl);
    else if (cut_file==17) bin_name = get_xs_2d_bin_name(prefix, "Eavail", Eavail_true, xs::Eavail_slices, "Emu", Emuon, Emu_bins_Eavail);
    else if (cut_file==21) bin_name = flag_0p ? get_xs_3d_bin_name(prefix, "Eavail", Eavail_true, xs::Eavail_outer_0p, "costheta", costh, xs::costheta_inner_0p, "Emu", Emuon, Emu_bins_3d_0p)
                                              : get_xs_3d_bin_name(prefix, "Eavail", Eavail_true, xs::Eavail_outer_Np, "costheta", costh, xs::costheta_inner_Np, "Emu", Emuon, Emu_bins_3d_Np);
    else if (cut_file==22) bin_name = flag_0p ? get_xs_3d_bin_name(prefix, "Eavail", Eavail_true, xs::Eavail_outer_0p, "pl", muon_p.Z(), xs::pl_inner_0p, "pt", muon_p.Perp(), pt_bins_3d_0p)
                                              : get_xs_3d_bin_name(prefix, "Eavail", Eavail_true, xs::Eavail_outer_Np, "pl", muon_p.Z(), xs::pl_inner_Np, "pt", muon_p.Perp(), pt_bins_3d_Np);
    // one bin for all 0p events, Np bins of the leading proton
    else if (flag_0p) bin_name = "numuCC.0p.inside";
    else{
      TVector3 p_dir(pfeval.truth_startMomentum[lead_p][0], pfeval.truth_startMomentum[lead_p][1], pfeval.truth_startMomentum[lead_p][2]);
      double costhp = p_dir.CosTheta();
      double cos_mu_p = p_dir.Unit().Dot(muon_p.Unit());
      if (cut_file==9) bin_name = get_xs_bin_name("numuCC.Np.inside.Kp", Kp, Kp_edges);
      else if (cut_file==10) bin_name = get_xs_bin_name("numuCC.Np.inside.costhetap", costhp, costhetap_edges);
      else if (cut_file==11) bin_name = get_xs_bin_name("numuCC.Np.inside.costhetamup", cos_mu_p, costhetamup_edges);
      else if (cut_file==13) bin_name = get_xs_2d_bin_name("numuCC.Np.inside.", "costhetap", costhp, xs::costheta_slices, "Kp", Kp, Kp_bins_costhetap);
      else if (cut_file==18) bin_name = get_xs_2d_bin_name("numuCC.Np.inside.", "Kp", Kp, xs::Kp_slices_mup, "costhetamup", cos_mu_p, costhetamup_bins_Kp);
      else if (cut_file==19) bin_name = get_xs_2d_bin_name("numuCC.Np.inside.", "Eavail", Eavail_true, xs::Eavail_slices_mup, "costhetamup", cos_mu_p, costhetamup_bins_Eavail);
      else if (cut_file==20) bin_name = get_xs_2d_bin_name("numuCC.Np.inside.", "Emu", Emuon, xs::Emu_slices_mup, "costhetamup", cos_mu_p, costhetamup_bins_Emu);
      else if (cut_file==23) bin_name = get_xs_3d_bin_name("numuCC.Np.inside.", "costheta", costh, xs::costheta_outer_p, "Emu", Emuon, xs::Emu_inner_p, "Kp", Kp, Kp_bins_costheta_Emu);
      else if (cut_file==24) bin_name = get_xs_3d_bin_name("numuCC.Np.inside.", "costheta", costh, xs::costheta_outer_p, "Emu", Emuon, xs::Emu_inner_p, "costhetamup", cos_mu_p, costhetamup_bins_costheta_Emu);
      else if (cut_file==25) bin_name = get_xs_3d_bin_name("numuCC.Np.inside.", "pl", muon_p.Z(), xs::pl_outer_p, "pt", muon_p.Perp(), xs::pt_inner_p, "Kp", Kp, Kp_bins_pl_pt);
      else bin_name = get_xs_3d_bin_name("numuCC.Np.inside.", "pl", muon_p.Z(), xs::pl_outer_p, "pt", muon_p.Perp(), xs::pt_inner_p, "costhetamup", cos_mu_p, costhetamup_bins_pl_pt);   // 26
    }
    auto it = map_cut_xs_bin.find(bin_name);
    if (it != map_cut_xs_bin.end()) return it->second;
    std::cout << "get_xs_signal_no: no bin " << bin_name << " in the xs bin file!" << std::endl;
    return -1;
  }

  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    const TString& cut_name = it->first;
    int number = it->second;

    if (cut_file == 1){
      if (cut_name == "numuCC.inside.Enu.le.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400.gt.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 && eval.truth_nuEnergy>300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500.gt.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 && eval.truth_nuEnergy>400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.600.gt.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=600 && eval.truth_nuEnergy>500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.700.gt.600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=700 && eval.truth_nuEnergy>600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.800.gt.700"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=800 && eval.truth_nuEnergy>700) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.900.gt.800"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=900 && eval.truth_nuEnergy>800) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1000.gt.900"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1000 && eval.truth_nuEnergy>900) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1100.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1100 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1100) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1500.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1500 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2100.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2100 && eval.truth_nuEnergy>1500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1400.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1400 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1600.gt.1400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1600 && eval.truth_nuEnergy>1400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2000.gt.1600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2000 && eval.truth_nuEnergy>1600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2500.gt.2000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2500 && eval.truth_nuEnergy>2000) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2500) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2100) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>1500) return number;
      }else{
	         std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
  }

  return -1;
}

// Index of the bin of value for the inner bin edges e1 < ... < en: 0 for value <= e1, k for ek < value <= e(k+1),
// n for value > en (the same convention as get_xs_bin_name).
int LEEana::get_bin_index(double value, const std::vector<double>& edges){
  int index = 0;
  while (index < (int)edges.size() && value > edges.at(index)) index++;
  return index;
}

// Flattened index of a double-differential bin: the bins of all lower slices, then the bin of value within its slice
// (slices from the inner edges slice_edges, bins in slice k from the inner edges bin_edges[k]).
int LEEana::get_2d_bin_index(double slice_value, const std::vector<double>& slice_edges, double value, const std::vector<std::vector<double>>& bin_edges){
  int slice = get_bin_index(slice_value, slice_edges);
  int index = 0;
  for (int k=0; k<slice; k++) index += bin_edges.at(k).size() + 1;
  return index + get_bin_index(value, bin_edges.at(slice));
}

// Flattened index of a triple-differential bin: outer slices (outer_edges), inner slices within outer slice k
// (inner_edges[k]), bins within each inner slice (bin_edges[outer][inner]).
int LEEana::get_3d_bin_index(double outer_value, const std::vector<double>& outer_edges, double inner_value, const std::vector<std::vector<double>>& inner_edges, double value, const std::vector<std::vector<std::vector<double>>>& bin_edges){
  int outer = get_bin_index(outer_value, outer_edges);
  int index = 0;
  for (int k=0; k<outer; k++)
    for (size_t l=0; l<bin_edges.at(k).size(); l++) index += bin_edges.at(k).at(l).size() + 1;
  return index + get_2d_bin_index(inner_value, inner_edges.at(outer), value, bin_edges.at(outer));
}

// Name of the bin of value for the inner bin edges e1 < ... < en: prefix + ".le.e1" (value <= e1),
// ".le.e(k+1).gt.ek" (ek < value <= e(k+1)) or ".gt.en" (value > en). Just prefix without edges (a single bin).
TString LEEana::get_xs_bin_name(TString prefix, double value, const std::vector<double>& edges){
  if (edges.empty()) return prefix;
  if (value <= edges.front()) return prefix + Form(".le.%g", edges.front());
  for (size_t k=1; k<edges.size(); k++){
    if (value <= edges.at(k)) return prefix + Form(".le.%g.gt.%g", edges.at(k), edges.at(k-1));
  }
  return prefix + Form(".gt.%g", edges.back());
}

// Name of a double-differential bin: the slice of slice_value (inner slice edges slice_edges) followed by the bin of value
// within that slice (inner edges bin_edges[slice index]), e.g. prefix + "costheta.le.0.3.gt.0" + ".Emu.le.300.gt.250".
TString LEEana::get_xs_2d_bin_name(TString prefix, TString slice_var, double slice_value, const std::vector<double>& slice_edges, TString var, double value, const std::vector<std::vector<double>>& bin_edges){
  TString slice_name = get_xs_bin_name(prefix + slice_var, slice_value, slice_edges);
  return get_xs_bin_name(slice_name + "." + var, value, bin_edges.at(get_bin_index(slice_value, slice_edges)));
}

// Name of a triple-differential bin: outer slice, inner slice (inner_edges[outer]), bin, e.g.
// prefix + "Eavail.le.75" + ".costheta.le.0.5.gt.0" + ".Emu.le.450.gt.350".
TString LEEana::get_xs_3d_bin_name(TString prefix, TString outer_var, double outer_value, const std::vector<double>& outer_edges, TString inner_var, double inner_value, const std::vector<std::vector<double>>& inner_edges, TString var, double value, const std::vector<std::vector<std::vector<double>>>& bin_edges){
  int outer = get_bin_index(outer_value, outer_edges);
  TString outer_name = get_xs_bin_name(prefix + outer_var, outer_value, outer_edges);
  return get_xs_2d_bin_name(outer_name + ".", inner_var, inner_value, inner_edges.at(outer), var, value, bin_edges.at(outer));
}

// Direction of reco particle index, from its end points, pointing away from the reco neutrino vertex.
TVector3 LEEana::get_reco_proton_dir(PFevalInfo& pfeval, int index){
  TVector3 p_start(pfeval.reco_startXYZT[index][0], pfeval.reco_startXYZT[index][1], pfeval.reco_startXYZT[index][2]);
  TVector3 p_end(pfeval.reco_endXYZT[index][0], pfeval.reco_endXYZT[index][1], pfeval.reco_endXYZT[index][2]);
  TVector3 vtx(pfeval.reco_nuvtxX, pfeval.reco_nuvtxY, pfeval.reco_nuvtxZ);
  return ((p_end-vtx).Mag() >= (p_start-vtx).Mag()) ? p_end-p_start : p_start-p_end;
}

// Cosine of the opening angle between the reco muon (reco_muonMomentum) and reco particle index (get_reco_proton_dir).
double LEEana::get_reco_cos_mu_p(PFevalInfo& pfeval, int index){
  TVector3 mu_dir(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2]);
  return get_reco_proton_dir(pfeval, index).Unit().Dot(mu_dir.Unit());
}

// Event-level quantities used by get_cut_pass that do not depend on the channel or add_cut.
// Fill once per event (after reading the entry) and pass to get_cut_pass for every histogram of that event.
void LEEana::fill_cut_event_info(CutEventInfo& info, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){
  double& reco_Enu = info.reco_Enu;
  double KE_muon;
  double Pmuon;
  double Emuon;
  double Ehadron;
  TLorentzVector truth_muonMomentum;
  bool flag_truth_inside;
  double true_proton_KE;
  std::map<std::string, bool>& map_cuts_flag = info.map_cuts_flag;
  bool& flag_generic = info.flag_generic;
  bool& flag_numuCC = info.flag_numuCC;
  bool& flag_numuCC_tight = info.flag_numuCC_tight;
  bool& flag_numuCC_1mu0p = info.flag_numuCC_1mu0p;
  bool& flag_numuCC_cutbased = info.flag_numuCC_cutbased;
  bool& flag_nueCC = info.flag_nueCC;
  bool& flag_0p = info.flag_0p;
  bool& flag_cc_pi0 = info.flag_cc_pi0;
  bool& flag_FC = info.flag_FC;
  bool& flag_FC_lepton = info.flag_FC_lepton;
  bool& flag_FC_hadron = info.flag_FC_hadron;
  TLorentzVector muonMomentum;
  int& costheta_bin = info.costheta_bin;
  int& Enu_bin = info.Enu_bin;

  reco_Enu = get_reco_Enu_corr(kine, flag_data);

  KE_muon = pfeval.truth_muonMomentum[3]*1000.-105.66; // MeV
  Pmuon = (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));

  Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
  Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

  truth_muonMomentum = TLorentzVector(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);

  flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  true_proton_KE = get_KE(pfeval, 2212, 1, 0, 1, 0);
  // definition of additional cuts
  map_cuts_flag.clear();
  if(is_far_sideband(kine, tagger, flag_data)) map_cuts_flag["farsideband"] = true;
  else map_cuts_flag["farsideband"] = false;

  if(is_near_sideband(kine, tagger, flag_data)) map_cuts_flag["nearsideband"] = true;
  else map_cuts_flag["nearsideband"] = false;

  if(is_nueCC(tagger)) map_cuts_flag["nueCC"] = true;
  else map_cuts_flag["nueCC"] = false;

  if(is_loosenueCC(tagger)) map_cuts_flag["loosenueCC"] = true;
  else map_cuts_flag["loosenueCC"] = false;

  if(is_generic(eval)) map_cuts_flag["generic"] = true;
  else map_cuts_flag["generic"] = false;


  if(eval.truth_nuEnergy<=400) map_cuts_flag["LowEnu"] = true;
  else map_cuts_flag["LowEnu"] = false;

  if(!(eval.truth_nuEnergy<=400)) map_cuts_flag["antiLowEnu"] = true;
  else map_cuts_flag["antiLowEnu"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["badmatch"] = true;
  else map_cuts_flag["badmatch"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["numuCCinFV"] = true;
  else map_cuts_flag["numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["RnumuCCinFV"] = true;
  else map_cuts_flag["RnumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["AnumuCCinFV"] = true;
  else map_cuts_flag["AnumuCCinFV"] = false;

  // Xs related cuts ...

  map_cuts_flag["XsnumuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1;


  map_cuts_flag["XsnumuCCinFV0p"] = map_cuts_flag["XsnumuCCinFV"] && true_proton_KE<45;
  map_cuts_flag["XsnumuCCinFVNp"] = map_cuts_flag["XsnumuCCinFV"] && true_proton_KE>=45;

  map_cuts_flag["Xs_Enu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  map_cuts_flag["Xs_Enu_mu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon > 105.7 && Emuon<=2506) map_cuts_flag["Xs_Emu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Emu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon > 0 && Pmuon<=2500) map_cuts_flag["Xs_Pmu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Pmu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron > 30 && Ehadron <=2500) map_cuts_flag["Xs_Ehad_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Ehad_numuCCinFV"] = false;

  // xs breakdown mode
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCinFV"] = true;
  else map_cuts_flag["XsecNumuCCinFV"] = false;

  //if(true_proton_KE<15 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCC0pinFV"] = true;
  //else map_cuts_flag["XsecNumuCC0pinFV"] = false;
  if(true_proton_KE<45 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCC0pinFV"] = true;
  else map_cuts_flag["XsecNumuCC0pinFV"] = false;
  //if(true_proton_KE>15 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCNpinFV"] = true;
  //else map_cuts_flag["XsecNumuCCNpinFV"] = false;
  if(true_proton_KE>45 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCNpinFV"] = true;
  else map_cuts_flag["XsecNumuCCNpinFV"] = false;

  if(true_proton_KE<15 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCStub0pinFV"] = true;
  else map_cuts_flag["XsecNumuCCStub0pinFV"] = false;
  if(true_proton_KE>15 && true_proton_KE<45 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCStubinFV"] = true;
  else map_cuts_flag["XsecNumuCCStubinFV"] = false;
  if(true_proton_KE>45 && eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCStubNpinFV"] = true;
  else map_cuts_flag["XsecNumuCCStubNpinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0) map_cuts_flag["XsecNC"] = true;
  else map_cuts_flag["XsecNC"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["XsecCosmic"] = true;
  else map_cuts_flag["XsecCosmic"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==1 && !(eval.truth_nuPdg==14 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200)) map_cuts_flag["XsecBkgCC"] = true;
  else map_cuts_flag["XsecBkgCC"] = false;

  // finish Xs related cuts ...


  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["NCinFV"] = true;
  else map_cuts_flag["NCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["outFV"] = true;
  else map_cuts_flag["outFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==1) map_cuts_flag["CCpi0inFV"] = true;
  else map_cuts_flag["CCpi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==1) map_cuts_flag["NCpi0inFV"] = true;
  else map_cuts_flag["NCpi0inFV"] = false;

  // breakdown categories for NC Delta analysis
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1) map_cuts_flag["NCDeltainFV"] = true;
  else map_cuts_flag["NCDeltainFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["NC1Pi0inFV"] = true;
  else map_cuts_flag["NC1Pi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCC1Pi0inFV"] = true;
  else map_cuts_flag["numuCC1Pi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCCotherinFV"] = true;
  else map_cuts_flag["numuCCotherinFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["NCotherinFV"] = true;
  else map_cuts_flag["NCotherinFV"] = false;
  // done with NC Delta breakdown categories

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  pfeval.reco_muonMomentum[3] > 0) map_cuts_flag["muon"] = true;
  else map_cuts_flag["muon"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  !(pfeval.reco_muonMomentum[3] > 0)) map_cuts_flag["nomuon"] = true;
  else map_cuts_flag["nomuon"] = false;

  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCMEC"] = true;
  else map_cuts_flag["CCMEC"] = false;

  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCMEC"] = true;
  else map_cuts_flag["NCMEC"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCQE"] = true;
  else map_cuts_flag["CCQE"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCQE"] = true;
  else map_cuts_flag["NCQE"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCRES"] = true;
  else map_cuts_flag["CCRES"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCRES"] = true;
  else map_cuts_flag["NCRES"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCDIS"] = true;
  else map_cuts_flag["CCDIS"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCDIS"] = true;
  else map_cuts_flag["NCDIS"] = false;

  if(pfeval.truth_nuScatType!=10 && pfeval.truth_nuScatType!=1 && pfeval.truth_nuScatType!=3 && pfeval.truth_nuScatType!=4 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["OTHER"] = true;
  else map_cuts_flag["OTHER"] = false;

  map_cuts_flag["none"] = false;
  map_cuts_flag["LEE"] = true;


  flag_generic = is_generic(eval);
  flag_numuCC = is_numuCC(tagger);
  //bool flag_numuCC = is_numuCC(tagger) && (is_far_sideband(kine, tagger, flag_data) || is_near_sideband(kine, tagger, flag_data) );
  flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  flag_nueCC = is_nueCC(tagger);

  flag_0p = is_0p(tagger, kine, pfeval);

  flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  flag_FC = is_FC(eval);
  // lepton FC/PC = muon energy from range/MCS with the one-sided 5% method (kine_reco_Enu_new3_5); hadron FC/PC from the
  // containment of the other particles (always FC in FC events)
  std::tuple<bool,bool> result_part_FC = get_part_is_FC(pfeval, eval, 3, 0.05);
  flag_FC_lepton = std::get<0>(result_part_FC);
  flag_FC_hadron = std::get<1>(result_part_FC);

  muonMomentum = TLorentzVector(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);

  costheta_bin = get_costheta_bin(TMath::Cos(muonMomentum.Theta()));
  Enu_bin = get_Enu_bin(reco_Enu);

  info.part_bin_set = false;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){
  CutEventInfo info;
  fill_cut_event_info(info, flag_data, eval, pfeval, tagger, kine, space, pandora, lantern);
  return get_cut_pass(ch_name, add_cut, flag_data, info, eval, pfeval, tagger, kine, space, pandora, lantern);
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, CutEventInfo& info, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){

  // event-level quantities, computed once per event in fill_cut_event_info
  double reco_Enu = info.reco_Enu;
  std::map<std::string, bool>& map_cuts_flag = info.map_cuts_flag;   // only existing keys are read below

  // figure out additional cuts and flag_data ...
  bool flag_add = true;
  if(add_cut == "all") flag_add = true;
  else if( (flag_data && (add_cut=="none" || add_cut=="farsideband" || add_cut=="nearsideband" || add_cut=="nueCC" || add_cut=="generic" || add_cut=="loosenueCC")) || !flag_data ){
      std::istringstream sss(add_cut.Data());
      for(std::string line; std::getline(sss, line, '_');){
          if(map_cuts_flag.find(line)!=map_cuts_flag.end()){
              flag_add *= map_cuts_flag[line];
          }
          else{
              std::cout<<"ERROR: add_cut "<<line<<" not defined!\n";
              exit(EXIT_FAILURE);
          }
      }
  }
  else{
    std::cout<<"ERROR: add_cut "<<add_cut<<" of channel "<< ch_name <<" is not assigned to sample "<<flag_data<<" [1: data; 0: mc]\n";
    std::cout<<"Please modify inc/WCPLEEANA/cuts.h\n";
    exit(EXIT_FAILURE);
  }

  if (!flag_add) return false;

  bool flag_generic = info.flag_generic;
  bool flag_numuCC = info.flag_numuCC;
  bool flag_numuCC_tight = info.flag_numuCC_tight;
  bool flag_numuCC_1mu0p = info.flag_numuCC_1mu0p;
  bool flag_numuCC_cutbased = info.flag_numuCC_cutbased;
  bool flag_nueCC = info.flag_nueCC;
  bool flag_0p = info.flag_0p;
  bool flag_cc_pi0 = info.flag_cc_pi0;
  bool flag_FC = info.flag_FC;
  int costheta_bin = info.costheta_bin;
  int Enu_bin = info.Enu_bin;

  std::string ch_name_string(ch_name.Data());
  std::string sequence_to_find = "numuCC_part_bdt";  
  size_t pos = ch_name_string.find(sequence_to_find);
  // optional FC/PC split, given as the last "_" token of the channel name (removed before matching the name):
  //   _FC / _PC                              event FC/PC (match_isFC)
  //   _LFC / _LPC                            lepton: muon energy from range / MCS (flag_FC_lepton)
  //   _HFC / _HPC                            hadronic system contained or not (flag_FC_hadron)
  //   _LFCHFC / _LFCHPC / _LPCHFC / _LPCHPC  lepton and hadronic system
  //   _LPCorHPC                              lepton or hadronic system PC (everything but _LFCHFC)
  // and before it an optional copy index "_<n>" (the very last token), e.g. numuCC_part_bdt_sig_LFCHFC_1: the same selection
  // in more than one channel of a cov_input.txt (channel names must be unique), e.g. with different variables.
  int require_FC = -1;       // -1: no requirement, 1: FC, 0: PC
  int require_lepton = -1;
  int require_hadron = -1;
  bool require_any_PC = false;
  if(pos != std::string::npos){
    size_t posCopy = ch_name_string.rfind('_');
    if(posCopy != std::string::npos && posCopy+1 < ch_name_string.size()
       && ch_name_string.find_first_not_of("0123456789", posCopy+1) == std::string::npos) ch_name_string.erase(posCopy);
    size_t posSuffix = ch_name_string.rfind('_');
    std::string suffix = (posSuffix != std::string::npos) ? ch_name_string.substr(posSuffix+1) : "";
    bool flag_suffix = true;
    if(suffix == "FC") require_FC = 1;
    else if(suffix == "PC") require_FC = 0;
    else if(suffix == "LFC") require_lepton = 1;
    else if(suffix == "LPC") require_lepton = 0;
    else if(suffix == "HFC") require_hadron = 1;
    else if(suffix == "HPC") require_hadron = 0;
    else if(suffix == "LPCorHPC") require_any_PC = true;
    else if(suffix.size() == 6 && (suffix.substr(0,3) == "LFC" || suffix.substr(0,3) == "LPC") && (suffix.substr(3) == "HFC" || suffix.substr(3) == "HPC")){
      require_lepton = (suffix[1] == 'F');
      require_hadron = (suffix[4] == 'F');
    }
    else flag_suffix = false;
    if(flag_suffix) ch_name_string.erase(posSuffix);
  }

  if(ch_name_string == "numuCC_part_bdt_Np_sig"     || ch_name_string == "numuCC_part_bdt_Np_bck"     || ch_name_string == "numuCC_part_bdt_Np_ext"     || ch_name_string == "numuCC_part_bdt_Np_dirt"     || ch_name_string == "numuCC_part_bdt_Np"
   ||ch_name_string == "numuCC_part_bdt_PorLNp_sig" || ch_name_string == "numuCC_part_bdt_PorLNp_bck" || ch_name_string == "numuCC_part_bdt_PorLNp_ext" || ch_name_string == "numuCC_part_bdt_PorLNp_dirt" || ch_name_string == "numuCC_part_bdt_PorLNp"
   ||ch_name_string == "numuCC_part_bdt_Scat0p_sig" || ch_name_string == "numuCC_part_bdt_Scat0p_bck" || ch_name_string == "numuCC_part_bdt_Scat0p_ext" || ch_name_string == "numuCC_part_bdt_Scat0p_dirt" || ch_name_string == "numuCC_part_bdt_Scat0p"
   ||ch_name_string == "numuCC_part_bdt_Act0p_sig"  || ch_name_string == "numuCC_part_bdt_Act0p_bck"  || ch_name_string == "numuCC_part_bdt_Act0p_ext"  || ch_name_string == "numuCC_part_bdt_Act0p_dirt"  || ch_name_string == "numuCC_part_bdt_Act0p"
   ||ch_name_string == "numuCC_part_bdt_Trk0p_sig"  || ch_name_string == "numuCC_part_bdt_Trk0p_bck"  || ch_name_string == "numuCC_part_bdt_Trk0p_ext"  || ch_name_string == "numuCC_part_bdt_Trk0p_dirt"  || ch_name_string == "numuCC_part_bdt_Trk0p"
   ||ch_name_string == "numuCC_part_bdt_G0p_sig"    || ch_name_string == "numuCC_part_bdt_G0p_bck"    || ch_name_string == "numuCC_part_bdt_G0p_ext"    || ch_name_string == "numuCC_part_bdt_G0p_dirt"    || ch_name_string == "numuCC_part_bdt_G0p"
   ||ch_name_string == "numuCC_part_bdt_sig"        || ch_name_string == "numuCC_part_bdt_bck"        || ch_name_string == "numuCC_part_bdt_ext"        || ch_name_string == "numuCC_part_bdt_dirt"        || ch_name_string == "numuCC_part_bdt"
   ||ch_name_string == "numuCC_part_bdt_0p_sig"     || ch_name_string == "numuCC_part_bdt_0p_bck"     || ch_name_string == "numuCC_part_bdt_0p_ext"     || ch_name_string == "numuCC_part_bdt_0p_dirt"     || ch_name_string == "numuCC_part_bdt_0p"
   || ch_name_string == "numuCC_part_bdt_Np_open" || ch_name_string == "numuCC_part_bdt_PorLNp_open" || ch_name_string == "numuCC_part_bdt_Scat0p_open" || ch_name_string == "numuCC_part_bdt_Act0p_open" || ch_name_string == "numuCC_part_bdt_Trk0p_open" || ch_name_string == "numuCC_part_bdt_G0p_open" || ch_name_string == "numuCC_part_bdt_0p_open"){

    if(pfeval.run<20700 && (ch_name_string == "numuCC_part_bdt_Np_open" || ch_name_string == "numuCC_part_bdt_PorLNp_open" || ch_name_string == "numuCC_part_bdt_Scat0p_open" || ch_name_string == "numuCC_part_bdt_Act0p_open" || ch_name_string == "numuCC_part_bdt_Trk0p_open" || ch_name_string == "numuCC_part_bdt_G0p_open" || ch_name_string == "numuCC_part_bdt_0p_open")) return false;

    if( (ch_name_string == "numuCC_part_bdt_Np_sig"    || ch_name_string == "numuCC_part_bdt_PorLNp_sig" || ch_name_string == "numuCC_part_bdt_Scat0p_sig" 
      || ch_name_string == "numuCC_part_bdt_Act0p_sig" || ch_name_string == "numuCC_part_bdt_Trk0p_sig"  || ch_name_string == "numuCC_part_bdt_G0p_sig" || ch_name_string == "numuCC_part_bdt_sig" || ch_name_string == "numuCC_part_bdt_0p_sig") 
      && map_cuts_flag["XsnumuCCinFV"]==false) return false;

    if( (ch_name_string == "numuCC_part_bdt_Np_bck" || ch_name_string == "numuCC_part_bdt_Np_dirt" || ch_name_string == "numuCC_part_bdt_PorLNp_bck" || ch_name_string == "numuCC_part_bdt_PorLNp_dirt"
       ||ch_name_string == "numuCC_part_bdt_Scat0p_bck" || ch_name_string == "numuCC_part_bdt_Scat0p_dirt" || ch_name_string == "numuCC_part_bdt_Act0p_bck" || ch_name_string == "numuCC_part_bdt_Act0p_dirt"
       ||ch_name_string == "numuCC_part_bdt_Trk0p_bck" || ch_name_string == "numuCC_part_bdt_Trk0p_dirt" || ch_name_string == "numuCC_part_bdt_G0p_bck" || ch_name_string == "numuCC_part_bdt_G0p_dirt"
       || ch_name_string == "numuCC_part_bdt_bck" || ch_name_string == "numuCC_part_bdt_dirt" || ch_name_string == "numuCC_part_bdt_0p_bck" || ch_name_string == "numuCC_part_bdt_0p_dirt") 
        && map_cuts_flag["XsnumuCCinFV"]==true) return false; 

    if(require_FC >= 0 && eval.match_isFC != require_FC) return false;
    if(require_lepton >= 0 && info.flag_FC_lepton != (require_lepton == 1)) return false;
    if(require_hadron >= 0 && info.flag_FC_hadron != (require_hadron == 1)) return false;
    if(require_any_PC && info.flag_FC_lepton && info.flag_FC_hadron) return false;

    if(tagger.numu_score<0.9 || pfeval.reco_muonMomentum[3]<0) return false; 

    if((ch_name_string == "numuCC_part_bdt_sig" || ch_name_string == "numuCC_part_bdt_bck" || ch_name_string == "numuCC_part_bdt_ext"
     || ch_name_string == "numuCC_part_bdt_dirt" || ch_name_string == "numuCC_part_bdt")) return true;

    //int part_bin = get_particle_0pNp_bdt_bin(pfeval, tagger, space, pandora, lantern, 45, 45, -0.65, 1.60); 
    // computed once per event and cached in info (recomputed if the thresholds change)
    double part_bin_thresholds[4] = {45, 45, 0.65, 1.60};
    if(!info.part_bin_set || !std::equal(part_bin_thresholds, part_bin_thresholds+4, info.part_bin_thresholds)){
      info.part_bin = get_particle_0pNp_bdt_bin(pfeval, tagger, space, pandora, lantern, part_bin_thresholds[0], part_bin_thresholds[1], part_bin_thresholds[2], part_bin_thresholds[3]);
      std::copy(part_bin_thresholds, part_bin_thresholds+4, info.part_bin_thresholds);
      info.part_bin_set = true;
    }
    int part_bin = info.part_bin;

    if((ch_name_string == "numuCC_part_bdt_Np_sig" || ch_name_string == "numuCC_part_bdt_Np_bck" || ch_name_string == "numuCC_part_bdt_Np_ext" 
     || ch_name_string == "numuCC_part_bdt_Np_dirt" || ch_name_string == "numuCC_part_bdt_Np" || ch_name_string == "numuCC_part_bdt_Np_open") && part_bin==5) return true;

    if((ch_name_string == "numuCC_part_bdt_PorLNp_sig" || ch_name_string == "numuCC_part_bdt_PorLNp_bck" || ch_name_string == "numuCC_part_bdt_PorLNp_ext" 
     || ch_name_string == "numuCC_part_bdt_PorLNp_dirt" || ch_name_string == "numuCC_part_bdt_PorLNp" || ch_name_string == "numuCC_part_bdt_PorLNp_open") && part_bin==4) return true;

    if((ch_name_string == "numuCC_part_bdt_Scat0p_sig" || ch_name_string == "numuCC_part_bdt_Scat0p_bck" || ch_name_string == "numuCC_part_bdt_Scat0p_ext" 
     || ch_name_string == "numuCC_part_bdt_Scat0p_dirt" || ch_name_string == "numuCC_part_bdt_Scat0p" || ch_name_string == "numuCC_part_bdt_Scat0p_open") && part_bin==3) return true;

    if((ch_name_string == "numuCC_part_bdt_Act0p_sig" || ch_name_string == "numuCC_part_bdt_Act0p_bck" || ch_name_string == "numuCC_part_bdt_Act0p_ext" 
     || ch_name_string == "numuCC_part_bdt_Act0p_dirt" || ch_name_string == "numuCC_part_bdt_Act0p" || ch_name_string == "numuCC_part_bdt_Act0p_open") && part_bin==2) return true;

    if((ch_name_string == "numuCC_part_bdt_Trk0p_sig" || ch_name_string == "numuCC_part_bdt_Trk0p_bck" || ch_name_string == "numuCC_part_bdt_Trk0p_ext" 
     || ch_name_string == "numuCC_part_bdt_Trk0p_dirt" || ch_name_string == "numuCC_part_bdt_Trk0p" || ch_name_string == "numuCC_part_bdt_Trk0p_open") && part_bin==1) return true;

    if((ch_name_string == "numuCC_part_bdt_G0p_sig" || ch_name_string == "numuCC_part_bdt_G0p_bck" || ch_name_string == "numuCC_part_bdt_G0p_ext" 
     || ch_name_string == "numuCC_part_bdt_G0p_dirt" || ch_name_string == "numuCC_part_bdt_G0p" || ch_name_string == "numuCC_part_bdt_G0p_open") && part_bin==0) return true;
     
    if((ch_name_string == "numuCC_part_bdt_0p_sig" || ch_name_string == "numuCC_part_bdt_0p_bck" || ch_name_string == "numuCC_part_bdt_0p_ext"
     || ch_name_string == "numuCC_part_bdt_0p_dirt" || ch_name_string == "numuCC_part_bdt_0p" || ch_name_string == "numuCC_part_bdt_0p_open") && part_bin<4) return true;

    return false;
  }else if (ch_name == "numuCC_theta_all_0p_FC_overlay" || ch_name == "BG_numuCC_theta_all_0p_FC_ext" || ch_name =="BG_numuCC_theta_all_0p_FC_dirt" || ch_name == "numuCC_theta_all_0p_FC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_Np_FC_overlay" || ch_name == "BG_numuCC_theta_all_Np_FC_ext" || ch_name =="BG_numuCC_theta_all_Np_FC_dirt" || ch_name == "numuCC_theta_all_Np_FC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all_FC_overlay" || ch_name == "BG_numuCC_all_FC_ext" || ch_name =="BG_numuCC_all_FC_dirt" || ch_name == "numuCC_all_FC_bnb") {
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all2_FC_overlay" || ch_name == "BG_numuCC_all2_FC_ext" || ch_name =="BG_numuCC_all2_FC_dirt" || ch_name == "numuCC_all2_FC_bnb") {
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_0p_FC_overlay" || ch_name == "BG_numuCC_0p_FC_ext" || ch_name =="BG_numuCC_0p_FC_dirt" || ch_name == "numuCC_0p_FC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_FC_overlay" || ch_name == "BG_numuCC_Np_FC_ext" || ch_name =="BG_numuCC_Np_FC_dirt" || ch_name == "numuCC_Np_FC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_PC_overlay" || ch_name == "BG_numuCC_theta_all_PC_ext" || ch_name =="BG_numuCC_theta_all_PC_dirt" || ch_name == "numuCC_theta_all_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_0p_PC_overlay" || ch_name == "BG_numuCC_theta_all_0p_PC_ext" || ch_name =="BG_numuCC_theta_all_0p_PC_dirt" || ch_name == "numuCC_theta_all_0p_PC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_Np_PC_overlay" || ch_name == "BG_numuCC_theta_all_Np_PC_ext" || ch_name =="BG_numuCC_theta_all_Np_PC_dirt" || ch_name == "numuCC_theta_all_Np_PC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all_PC_overlay" || ch_name == "BG_numuCC_all_PC_ext" || ch_name =="BG_numuCC_all_PC_dirt" || ch_name == "numuCC_all_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all2_PC_overlay" || ch_name == "BG_numuCC_all2_PC_ext" || ch_name =="BG_numuCC_all2_PC_dirt" || ch_name == "numuCC_all2_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_0p_PC_overlay" || ch_name == "BG_numuCC_0p_PC_ext" || ch_name =="BG_numuCC_0p_PC_dirt" || ch_name == "numuCC_0p_PC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_PC_overlay" || ch_name == "BG_numuCC_Np_PC_ext" || ch_name =="BG_numuCC_Np_PC_dirt" || ch_name == "numuCC_Np_PC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_both_overlay" || ch_name == "BG_numuCC_Np_both_ext" || ch_name =="BG_numuCC_Np_both_dirt" || ch_name == "numuCC_Np_both_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb" || ch_name == "numuCC_nopi0_nonueCC_FC_numu2nueoverlay"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb" || ch_name == "numuCC_nopi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay" || ch_name =="BG_CCpi0_nonueCC_FC_ext" || ch_name == "BG_CCpi0_nonueCC_FC_dirt" || ch_name == "CCpi0_nonueCC_FC_bnb" || ch_name == "CCpi0_nonueCC_FC_numu2nueoverlay"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay" || ch_name == "BG_CCpi0_nonueCC_PC_ext" || ch_name == "BG_CCpi0_nonueCC_PC_dirt" || ch_name == "CCpi0_nonueCC_PC_bnb" || ch_name == "CCpi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;


 // generic selection nu PC+FC 1 obs channel
}else if (ch_name == "generic_nu_overlay" || ch_name == "BG_generic_nu_ext" || ch_name =="BG_generic_nu_dirt" || ch_name == "generic_nu_bnb" ||
          ch_name == "generic_nu_overlay_2" || ch_name == "BG_generic_nu_ext_2" || ch_name =="BG_generic_nu_dirt_2" || ch_name == "generic_nu_bnb_2" ||
          ch_name == "generic_nu_overlay_3" || ch_name == "BG_generic_nu_ext_3" || ch_name =="BG_generic_nu_dirt_3" || ch_name == "generic_nu_bnb_3" ||
          ch_name == "generic_nu_overlay_4" || ch_name == "BG_generic_nu_ext_4" || ch_name =="BG_generic_nu_dirt_4" || ch_name == "generic_nu_bnb_4"){
    if (flag_generic) return true;
    else return false;
 // numuCC selection PC+FC 1 obs channel
  }else if (ch_name == "numuCC_overlay" || ch_name == "BG_numuCC_ext" || ch_name =="BG_numuCC_dirt" || ch_name == "numuCC_bnb"){
    if (flag_numuCC) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc" || ch_name == "BG_numuCC_ext_fc" || ch_name =="BG_numuCC_dirt_fc" || ch_name == "numuCC_bnb_fc"){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc_np" || ch_name == "BG_numuCC_ext_fc_np" || ch_name =="BG_numuCC_dirt_fc_np" || ch_name == "numuCC_bnb_fc_np"){
    if (flag_numuCC && flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc_0p" || ch_name == "BG_numuCC_ext_fc_0p" || ch_name =="BG_numuCC_dirt_fc_0p" || ch_name == "numuCC_bnb_fc_0p"){
    if (flag_numuCC && flag_FC && (flag_0p)) return true;
    else return false;
 // cutbased numuCC selection PC+FC 1 obs channel
  }else if (ch_name == "numuCC_cutbased_overlay" || ch_name == "BG_numuCC_cutbased_ext" || ch_name =="BG_numuCC_cutbased_dirt" || ch_name == "numuCC_cutbased_bnb"){
    if (flag_numuCC_cutbased) return true;
    else return false;
    // add some cuts for Xs related cases ...
  }else if (ch_name == "numuCC_FC_bnb" || ch_name == "BG_numuCC_FC_ext" || ch_name == "BG_numuCC_FC_dirt"
	    || ch_name == "numuCC1_FC_bnb" || ch_name == "BG_numuCC1_FC_ext" || ch_name == "BG_numuCC1_FC_dirt"
	    || ch_name == "numuCC2_FC_bnb" || ch_name == "BG_numuCC2_FC_ext" || ch_name == "BG_numuCC2_FC_dirt"
	    ){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb" || ch_name == "BG_numuCC_PC_ext" || ch_name == "BG_numuCC_PC_dirt"
	    || ch_name == "numuCC1_PC_bnb" || ch_name == "BG_numuCC1_PC_ext" || ch_name == "BG_numuCC1_PC_dirt"
	    || ch_name == "numuCC2_PC_bnb" || ch_name == "BG_numuCC2_PC_ext" || ch_name == "BG_numuCC2_PC_dirt"
	    ){
    if (flag_numuCC && (!flag_FC)) return true;
    else return false;
  // 1D Enu channel with the same inclusive signal definition as the 3D selection
  }else if (ch_name == "numuCC_Enu_mu_FC_bnb" || ch_name == "BG_numuCC_Enu_mu_FC_ext" || ch_name == "BG_numuCC_Enu_mu_FC_dirt"){
    if (flag_numuCC &&   flag_FC  && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;
  }else if (ch_name == "numuCC_Enu_mu_PC_bnb" || ch_name == "BG_numuCC_Enu_mu_PC_ext" || ch_name == "BG_numuCC_Enu_mu_PC_dirt"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;

  }else if (ch_name == "numuCC_FC_bnb_L800MeV" || ch_name == "BG_numuCC_FC_ext_L800MeV" || ch_name == "BG_numuCC_FC_dirt_L800MeV"
	    || ch_name == "numuCC1_FC_bnb_L800MeV" || ch_name == "BG_numuCC1_FC_ext_L800MeV" || ch_name == "BG_numuCC1_FC_dirt_L800MeV"
	    || ch_name == "numuCC2_FC_bnb_L800MeV" || ch_name == "BG_numuCC2_FC_ext_L800MeV" || ch_name == "BG_numuCC2_FC_dirt_L800MeV"
	    ){
    if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb_L800MeV" || ch_name == "BG_numuCC_PC_ext_L800MeV" || ch_name == "BG_numuCC_PC_dirt_L800MeV"
	    || ch_name == "numuCC1_PC_bnb_L800MeV" || ch_name == "BG_numuCC1_PC_ext_L800MeV" || ch_name == "BG_numuCC1_PC_dirt_L800MeV"
	    || ch_name == "numuCC2_PC_bnb_L800MeV" || ch_name == "BG_numuCC2_PC_ext_L800MeV" || ch_name == "BG_numuCC2_PC_dirt_L800MeV"
	    ){
    if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC_PC_overlay_L800MeV"
	    || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV"
	    || ch_name == "numuCC2_FC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV"   ){
    if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC2_FC_overlay_L800MeV"){
      if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    }else if (ch_name == "numuCC_PC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV" ){
      if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    }
    return false;
  }else if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC_PC_overlay"
	    || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC1_PC_overlay"
	    || ch_name == "numuCC2_FC_overlay" || ch_name == "numuCC2_PC_overlay"   ){
    if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC2_FC_overlay"){
      if (flag_numuCC && flag_FC ) return true;
    }else if (ch_name == "numuCC_PC_overlay" || ch_name == "numuCC1_PC_overlay" || ch_name == "numuCC2_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) ) return true;
    }
    return false;
  }else if (ch_name == "generic_nu_ext" || ch_name == "generic_nu_dirt" ||ch_name == "generic_nu_bnb_LEE"){
    if (flag_generic) return true;
    else return false;

  }else{
    std::cout << "Not sure what cut: " << ch_name << std::endl;
  }

  return false;
}

bool LEEana::get_rw_cut_pass(TString cut, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  if(cut == "NCPi0"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_Np"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && !(is_true_0p(pfeval)) && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_0p"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && is_true_0p(pfeval) && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCDeltaNp_scale"){
    if(is_NCdelta_sel(tagger, pfeval) && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && !(is_0p(tagger, kine, pfeval)) ) return true;
    return false;
  }else if(cut == "NCDelta0p_scale"){
    if(is_NCdelta_sel(tagger, pfeval) && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && is_0p(tagger, kine, pfeval)) return true;
    return false;
  }else if(cut == "NCPi0_NCDelta_Np"){
    if(eval.truth_isCC==0 && !(is_true_0p(pfeval)) && (pfeval.truth_NprimPio==1 || pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_NCDelta_0p"){
    if(eval.truth_isCC==0 && is_true_0p(pfeval) && (pfeval.truth_NprimPio==1 || pfeval.truth_NCDelta==1)) return true;
    return false;
  }else{
    std::cout<<"No matching reweighting cut, check reweight configuration file"<<std::endl;
  }
return false;
}

bool LEEana::is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;

  bool flag_numuCC = is_numuCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);

  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  if ((reco_Enu>=800 && tagger.nue_score >=0) ||
      (tagger.nue_score<=0 && (flag_numuCC || (flag_pi0 && flag_NC) ))) flag = true;
  return flag;
}
bool LEEana::is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  if (reco_Enu < 800 && tagger.nue_score>0 && (reco_Enu>=600 || tagger.nue_score<=7)) flag = true;

  return flag ;
}

bool LEEana::is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if (reco_Enu < 600 && tagger.nue_score>7) flag = true;
  return flag;
}


bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){
    return true;
  }else{
    return false;
  }
}

bool LEEana::is_cc_pi0(KineInfo& kine, bool flag_data){

  bool flag = false;

  if (flag_data){
    if (kine.kine_pio_mass>0){
      //     TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      // TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      //pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;

      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }

  return flag;
}


bool LEEana::is_pi0(KineInfo& kine, bool flag_data){
  bool flag = false;

  if (flag_data){
    if (kine.kine_pio_mass>0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      // pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;

      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }

  return flag;
}


bool LEEana::is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){ // includes all cuts except FC
  bool flag = false;
  if (tagger_info.nc_delta_score > 2.61 && tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}


//


bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;

  return flag;
}


bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;

  return flag;
}

bool LEEana::is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0)
    flag = true;

  return flag;
}

bool LEEana::is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          //if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>0) Nproton++; //Erin: CHANGE, "actual 0p" aka no 35 MeV threshold
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Nproton==0) flag = true;
  }

  return flag;
}


bool LEEana::is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Npion==0) flag = true;
  }

  return flag;
}

bool LEEana::is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Nproton==0) flag = true;
  }

  return flag;
}


bool LEEana::is_numuCC_cutbased(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag==1 && tagger_info.cosmict_flag==0)
    flag = true;

  return flag;
}


bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    //  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score <= 7.0 && tagger_info.nue_score > 0)
    flag = true;

  return flag;
}

bool LEEana::is_loosenueCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 4.0)
    flag = true;

  return flag;
}

bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){ // == T_BDTvars.numu_cc_flag >= 0
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;


  return flag;
}


#endif
