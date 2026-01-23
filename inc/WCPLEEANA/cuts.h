#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"

#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TRandom3.h"

namespace LEEana{
  // this is for the real data, for fake data this should be 1 ...
  double em_charge_scale = 0.95;
  //double em_charge_scale = 1.0;

  double get_angle(double momentum_0, double momentum_1, double momentum_2, int option, bool to_numi);
  double get_angle_to_absorber(double momentum_0, double momentum_1, double momentum_2);
  double get_true_nu_angle_to_absorber(PFevalInfo& pfeval);
  double get_true_nu_angle(PFevalInfo& pfeval, int option, bool to_numi);

  bool is_true_kdar_all(EvalInfo& eval, PFevalInfo& pfeval);
  bool is_true_kdar(EvalInfo& eval, PFevalInfo& pfeval);
  bool is_true_kdar_outFV(EvalInfo& eval, PFevalInfo& pfeval);

  double get_ssmE(TaggerInfo& tagger_info, bool flag_2track=false);
  bool is_kdar_presel(TaggerInfo& tagger_info,EvalInfo& eval);
  bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=0.1523, double hiE_cut=1.409);
  bool is_lowE_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=0.1523);
  bool is_hiE_kdar_bdtsel(TaggerInfo& tagger_info, double hiE_cut=1.409);

  //bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=-0.077, double hiE_cut=1.409);
  //bool is_lowE_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=-0.077);
  //bool is_hiE_kdar_bdtsel(TaggerInfo& tagger_info, double hiE_cut=1.409);




  //bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=0.69, double hiE_cut=1.57);
  //bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=-0.3, double hiE_cut=1.409);
  //bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=-0.3, double hiE_cut=1.02);
  //bool is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut=-0.3, double hiE_cut=1.77);
  double get_oldflux_weight(EvalInfo& eval,PFevalInfo& pfeval);

  double recover_nsbeamtiming(PFevalInfo& pfeval, bool flag_data, bool flag_kdar=false);
  double get_dE_dx_range(double R, int pdg);
  std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > get_extrapolated_times(PFevalInfo& pfeval, int particle, double mother_time);

  // correct reco neutrino energy and reco shower energy
  double get_reco_Enu_corr(KineInfo& kine, bool flag_data);
  double get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data);

  double get_reco_Eproton(KineInfo& kine);
  double get_reco_Epion(KineInfo& kine);

  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name="kine_reco_Enu", bool flag_kdar_file=false);
  double get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, TString var_name);

  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  bool get_rw_cut_pass(TString cut, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  double get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool,  std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data=false);
  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);

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
  bool is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data);
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
  bool is_NCpio_sel(TaggerInfo& tagger_info, KineInfo& kine);
  bool is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);

  //Erin
  bool is_singlephoton_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singleshower_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singleshower_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_numu_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_other_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_ncpi0_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue_sel_allshw(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_nsbeam(PFevalInfo& pfeval, EvalInfo& eval);
  bool is_nsbeam_photon(PFevalInfo& pfeval, EvalInfo& eval);
  bool is_singlephoton_pre(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_other(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_other(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_oneshw(TaggerInfo& tagger_info, PFevalInfo& pfeval);


  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";
  bool is_FC(EvalInfo& eval);


  // TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
  bool is_truth_nueCC_inside(EvalInfo& eval);
  bool is_truth_numuCC_inside(EvalInfo& eval);

  bool is_true_0p(PFevalInfo& pfeval);
  double get_true_Kp(PFevalInfo& pfeval,bool MeV=true,bool total=false);
  double get_true_p_angle_absorber(PFevalInfo& pfeval);

  int mcc8_pmuon_costheta_bin(float pmuon, float costh);
  int alt_var_index(std::string var1, float val1, std::string var2, float val2, std::string config="./configurations/alt_var_xbins.txt");
  std::map<std::string, TH1F> map_var_hist; // variable name and binning
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

double LEEana::get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data){
    if (flag_data){
        return pfeval.reco_showerKE * em_charge_scale;
    } else {
        return pfeval.reco_showerKE;
    }
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

double LEEana::get_reco_Epion(KineInfo& kine){
  double reco_Epion=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      //if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10)  // KE threshold: 10 keV
      if(abs(pdgcode)==211)  // no threshold
        reco_Epion+=kine.kine_energy_particle->at(i);
    }
  return reco_Epion;
}


bool LEEana::is_true_0p(PFevalInfo& pfeval){
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue; //Erin: CHANGE, no proton threshold
      return false;
    }
  return true;
}

double LEEana::get_true_Kp(PFevalInfo& pfeval,bool MeV,bool total){
    double Kp=0;
    double units=1;
    if(MeV) units=1000;
    double mass=0.938272*units;
    if(total) mass=0;
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] >Kp) Kp=pfeval.truth_startMomentum[i][3]; 
    }
  if (Kp==0.938272) return 0;
  return Kp*units-mass;
}

double LEEana::get_true_p_angle_absorber(PFevalInfo& pfeval){
    double Kp=0;
    double p0=0;
    double p1=0;
    double p2=0;
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] >Kp){
        Kp=pfeval.truth_startMomentum[i][3];
	p0 = pfeval.truth_startMomentum[i][0];
	p1 = pfeval.truth_startMomentum[i][1];
	p2 = pfeval.truth_startMomentum[i][2];
      }
    }
  if (Kp==0) return -999;
  return get_angle_to_absorber(p0,p1,p2);
}


double LEEana::get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data){
  double addtl_weight = 1.0;
  // CV correction from numuCC cross section data
  // if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1){
  //   if (eval.truth_nuEnergy>200 && eval.truth_nuEnergy<=540) addtl_weight = 1.28043;
  //   else if (eval.truth_nuEnergy>540 && eval.truth_nuEnergy<=705) addtl_weight = 1.21158;
  //   else if (eval.truth_nuEnergy>705 && eval.truth_nuEnergy<=805) addtl_weight = 1.19091;
  //   else if (eval.truth_nuEnergy>805 && eval.truth_nuEnergy<=920) addtl_weight = 1.17733;
  //   else if (eval.truth_nuEnergy>920 && eval.truth_nuEnergy<=1050) addtl_weight = 1.13983;
  //   else if (eval.truth_nuEnergy>1050 && eval.truth_nuEnergy<=1200) addtl_weight = 1.07864;
  //   else if (eval.truth_nuEnergy>1200 && eval.truth_nuEnergy<=1375) addtl_weight = 1.00722;
  //   else if (eval.truth_nuEnergy>1375 && eval.truth_nuEnergy<=1570) addtl_weight = 0.93857;
  //   else if (eval.truth_nuEnergy>1570 && eval.truth_nuEnergy<=2050) addtl_weight = 0.886241;
  //   else if (eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) addtl_weight = 0.858724;
  //   else if (eval.truth_nuEnergy>4000) addtl_weight = 0.858724;
  // }
  // std::cout << "energy: " << eval.truth_nuEnergy << " addtl_weight: " << addtl_weight << std::endl;
  // end of data correction

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
  }else if (weight_name == "cv_spline_oldflux"){
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    return eval.weight_cv * eval.weight_spline * oldflux_weight;
  }else if (weight_name == "cv_spline_oldflux_cv_spline_oldflux"){
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    return pow(eval.weight_cv * eval.weight_spline * oldflux_weight,2);
  }else if (weight_name == "cv_spline_cv_spline_oldflux" || weight_name == "cv_spline_oldflux_cv_spline"){
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    return oldflux_weight * pow(eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "add_weight"){//for systematics
    return addtl_weight;

    //Erin - ns beam time scaling
  }else if (weight_name == "cv_spline_nsbeam"){
    float beam_scale = 0.86;

    bool has_muon = false; //set to true to turn off
    if (pfeval.reco_muonMomentum[3] > 0){has_muon = true;}

    if (pfeval.run >= 13697){ 
      beam_scale=std::get<0>(time_info[3]);//0.845;}//0.931503;}//0.913671 - 0.0812331; }
      if (!has_muon){ beam_scale = 0.82;}
    }
    else if (pfeval.run >= 8321){ 
      beam_scale=std::get<0>( time_info[2]);
      if (!has_muon){ beam_scale = 0.86;}
    }//0.88;}//0.919618;}//0.900644 - 0.044328;}
    else if (pfeval.run > 0 ){ 
      beam_scale=std::get<0>( time_info[1]);//0.912832;}//0.885887 - 0.0315298;}
      if (!has_muon){ beam_scale = 0.78;}
    }
    //beam_scale = beam_scale - 0.03;
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}
    //0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}
    //0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}
    //0.527815;}//0.471911 + 0.0315298;}
    if(eval.match_completeness_energy<=0.1*eval.truth_energyInside){beam_scale = 1.0-ext_rej;}
    return addtl_weight*eval.weight_cv * eval.weight_spline * beam_scale;
  }else if (weight_name == "dirt_nsbeam"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    return addtl_weight*eval.weight_cv * eval.weight_spline * (1.0-ext_rej);
  }else if (weight_name == "nsbeam_ext"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}
    //0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}
    //0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}
    //0.527815;}//0.471911 + 0.0315298;}
    //ext_rej = ext_rej + 0.03;
    float ext_scale = 1.0 - ext_rej;
    return ext_scale;
  }else if (weight_name == "cv_spline_nsbeam_cv_spline_nsbeam"){
    float beam_scale = 0.86;

    bool has_muon = false; //set to true to turn off
    if (pfeval.reco_muonMomentum[3] > 0){has_muon = true;}

    if (pfeval.run >= 13697){ 
      beam_scale=std::get<0>(time_info[3]);//0.845;}//0.931503;}//0.913671 - 0.0812331; }
      if (!has_muon){ beam_scale = 0.82;}
    }
    else if (pfeval.run >= 8321){ 
      beam_scale=std::get<0>( time_info[2]);//0.88;}//0.919618;}//0.900644 - 0.044328;}
      if (!has_muon){ beam_scale = 0.86;}
    }
    else if (pfeval.run > 0 ){ 
      beam_scale=std::get<0>( time_info[1]);//0.912832;}//0.885887 - 0.0315298;}
      if (!has_muon){ beam_scale = 0.78;}
    }
    //beam_scale = beam_scale - 0.03;
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    if(eval.match_completeness_energy<=0.1*eval.truth_energyInside){beam_scale = 1.0-ext_rej;}
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * beam_scale,2);
  }else if (weight_name == "dirt_nsbeam_dirt_nsbeam"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * (1.0-ext_rej),2);
  }else if (weight_name == "nsbeam_ext_nsbeam_ext"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    //ext_rej = ext_rej + 0.03;
    float ext_scale = 1.0 - ext_rej;
    return pow(ext_scale,2);
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
  }else if (weight_name == "cv_spline_oldflux"){
    std::cout<<"Entered"<<std::endl;
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    std::cout<<oldflux_weight<<std::endl;
    return eval.weight_cv * eval.weight_spline * oldflux_weight;
  }else if (weight_name == "cv_spline_oldflux_cv_spline_oldflux"){
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    return pow(eval.weight_cv * eval.weight_spline * oldflux_weight,2);
  }else if (weight_name == "cv_spline_cv_spline_oldflux" || weight_name == "cv_spline_oldflux_cv_spline"){
    double oldflux_weight = get_oldflux_weight(eval,pfeval);
    return oldflux_weight * pow(eval.weight_cv * eval.weight_spline,2);
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


double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data , TString var_name, bool flag_kdar_file){
  //  if (var_name == "kine_reco_Enu"){
  //  return kine.kine_reco_Enu;
  //  }else
  if (var_name == "kine_reco_Enu"){
    return get_reco_Enu_corr(kine, flag_data);
  }else if (var_name == "ssm_kine_reco_Enu"){
    return tagger.ssm_kine_reco_Enu;
  }else if (var_name == "ssm_kine_energy"){
    return tagger.ssm_kine_energy;
  }else if (var_name == "ssm_kine_energy_force_min"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_kine_energy<30) return 30.1;
    return tagger.ssm_kine_energy;
  }else if (var_name == "ssm_angle_to_absorber"){
    if(std::isnan(tagger.ssm_angle_to_absorber)) return -999;
    return tagger.ssm_angle_to_absorber;
  }else if (var_name == "ssm_angle_to_absorber_deg"){
    if(std::isnan(tagger.ssm_angle_to_absorber)) return -999;
    if (tagger.ssm_kine_energy<0) return -999;
    return tagger.ssm_angle_to_absorber*180/3.14159;
  }else if (var_name == "ssm_pl" || var_name == "ssm_pt" || var_name == "ssm_q" || var_name == "ssm_Q2" || var_name == "ssm_sqrtQ2"){
    if(std::isnan(tagger.ssm_angle_to_absorber)) return -999;
    if (tagger.ssm_kine_energy<0) return -999;
    double k = sqrt(tagger.ssm_kine_energy*tagger.ssm_kine_energy+2*tagger.ssm_kine_energy*105.7);
    double var = -999;
    if(var_name == "ssm_pl") var = k*TMath::Cos(tagger.ssm_angle_to_absorber);
    if(var_name == "ssm_pt") var = k*TMath::Sin(tagger.ssm_angle_to_absorber);
    if(var_name == "ssm_q") var = sqrt(235.5*235.5+k*k-2*235.5*k*cos(tagger.ssm_angle_to_absorber));
    double Q2 = 2*235.5*(105.7+tagger.ssm_kine_energy-k*cos(tagger.ssm_angle_to_absorber))-105.7*105.7;		   
    if(var_name == "ssm_Q2") var = Q2;
    if(var_name == "ssm_sqrtQ2") var = sqrt(Q2);
    return var;
  }else if (var_name == "ssm_prim_track1_kine_energy_range"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_prim_track1_kine_energy_range<0) return -0.01;
    return tagger.ssm_prim_track1_kine_energy_range;
  }else if (var_name == "ssm_prim_track1_kine_energy_range_vtx"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_prim_track1_kine_energy_range<0 && tagger.ssm_vtx_activity) return -0.01;
    if (tagger.ssm_prim_track1_kine_energy_range<0) return -0.01-20;
    return tagger.ssm_prim_track1_kine_energy_range;
  }else if (var_name == "ssm_prim_track1_kine_energy_range_vtx50"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_prim_track1_kine_energy_range<0 && tagger.ssm_vtx_activity) return -0.01;
    if (tagger.ssm_prim_track1_kine_energy_range<0) return -0.01-50;
    if (tagger.ssm_prim_track1_kine_energy_range<50) return 10;
    return tagger.ssm_prim_track1_kine_energy_range;
  }else if (var_name == "ssm_prim_track1_kine_energy_range_vtx60"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_prim_track1_kine_energy_range<0 && tagger.ssm_vtx_activity) return -0.01;
    if (tagger.ssm_prim_track1_kine_energy_range<0) return -0.01-60;
    if (tagger.ssm_prim_track1_kine_energy_range<60) return 10;
    return tagger.ssm_prim_track1_kine_energy_range;
  }else if (var_name == "ssm_prim_track1_kine_energy_range_vtx40"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (tagger.ssm_prim_track1_kine_energy_range<0 && tagger.ssm_vtx_activity) return -0.01;
    if (tagger.ssm_prim_track1_kine_energy_range<0) return -0.01-40;
    return tagger.ssm_prim_track1_kine_energy_range;

  }else if (var_name == "ssm_E" || var_name == "ssm_KE" || var_name == "ssm_KE_2track"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (var_name == "ssm_KE_2track") return get_ssmE(tagger,true)-105.7;
    double E = get_ssmE(tagger);
    if(var_name == "ssm_E") return E;
    return E-105.7;

  }else if (var_name == "ssm_cosP" || var_name == "ssm_angle_P" || var_name == "ssm_angle_P_deg" || var_name == "ssm_angle_P_deg_vtx"){
    if (tagger.ssm_kine_energy<0) return -999;
    if (var_name == "ssm_angle_P_deg_vtx" && tagger.ssm_prim_track1_kine_energy_range<0 && tagger.ssm_vtx_activity) return -0.01;
    if (var_name == "ssm_angle_P_deg_vtx" && tagger.ssm_prim_track1_kine_energy_range<0) return -50-0.1;
    if (var_name != "ssm_angle_P_deg_vtx" && tagger.ssm_prim_track1_kine_energy_range<0) return -999;

    double theta = get_angle_to_absorber(tagger.ssm_prim_track1_x_dir, tagger.ssm_prim_track1_y_dir, tagger.ssm_prim_track1_z_dir);

    if (std::isnan(theta) && tagger.ssm_vtx_activity && var_name == "ssm_angle_P_deg_vtx") return -0.1;
    if (std::isnan(theta) && var_name == "ssm_angle_P_deg_vtx") return -50-0.1;
    if (std::isnan(theta)) return -999;

    if (var_name == "ssm_cosP") return TMath::Cos(theta);
    if (var_name == "ssm_angle_P") return theta;
    if (var_name == "ssm_angle_P_deg" || var_name == "ssm_angle_P_deg_vtx") return theta*180/3.14159;

  }else if (var_name == "evtTimeNS_recover"){
    return recover_nsbeamtiming(pfeval, flag_data); 
  }else if (var_name == "evtTimeNS"){
    if (flag_data) return pfeval.evtTimeNS;
    else return pfeval.evtTimeNS_cor;
  }else if(var_name=="merge_time_random"){
    gRandom = new TRandom3(eval.run*100000000000+eval.subrun*1000000+eval.event);
    double merge_time = gRandom->Uniform(-9.42, 9.42);
    //std::cout<<merge_time<<std::endl;
    return merge_time;
  }else if(var_name=="merge_time_recover_numi"){
    bool temp_flag_kdar_file = flag_kdar_file;
    if(!flag_data){if (flag_kdar_file==false && is_true_kdar(eval,pfeval)){temp_flag_kdar_file=true;}}// This is a dirty hack that works under the assumption that all the KDARs come from the KDAR overlay.
    double evtTimeNS = recover_nsbeamtiming(pfeval, flag_data, temp_flag_kdar_file);
    double TThelp=0;
    double gap = 18.8305;
    double fShift=0.9378;
    if(flag_data){
      // run4a
      if                        (eval.run<=19672) fShift=0.9378-3;
      else if (eval.run>19672 && eval.run<=19683) fShift=0.9378-7;
      else if (eval.run>19683 && eval.run<=19700) fShift=0.9378-3;
      else if (eval.run>19700 && eval.run<=19730) fShift=0.9378-6.5;
      else if (eval.run>19730 && eval.run<=19752) fShift=0.9378-5.5;
      // run4b
      else if (eval.run>19752 && eval.run<=19855) fShift=0.9378-7.2836;
      else if (eval.run>20800 && eval.run<=20855) fShift=0.9378-1.0278;
      else if (eval.run>20855 && eval.run<=20870) fShift=0.9378-10.1;
      else if (eval.run>20870 && eval.run<=20890) fShift=0.9378-2.48;
      else if (eval.run>20890 && eval.run<=20915) fShift=0.9378-5.02;
      else if (eval.run>20915 && eval.run<=21000) fShift=0.9378+3.3475;
      else if (eval.run>21000 && eval.run<=21286) fShift=0.9378-6.3796;
      //run4c
      else if (eval.run>21286 && eval.run<=22100) fShift=0.9378-6.3796;
      else if (eval.run>22100 && eval.run<=22125) fShift=0.9378-7.7;
      else if (eval.run>22125 && eval.run<=22140) fShift=0.9378-9.5;
      else if (eval.run>22140 && eval.run<=22160) fShift=0.9378-6;
      else if (eval.run>22160 && eval.run<=22190) fShift=0.9378-7.7;
      else if (eval.run>22190 && eval.run<=22240) fShift=0.9378-10.4;
      else if (eval.run>22240 && eval.run<=22269) fShift=0.9378-9.5;
      //run4d
      else if (eval.run>22269 && eval.run<=22290) fShift=0.9378-9.5;
      else if (eval.run>22290 && eval.run<=22300) fShift=0.9378-2.5;
      else if (eval.run>22300 && eval.run<=22320) fShift=0.9378-5.2;
      else if (eval.run>22320 && eval.run<=22350) fShift=0.9378-4;
      else if (eval.run>22350 && eval.run<=23150) fShift=0.9378-7.7;
      else if (eval.run>23150 && eval.run<=24319) fShift=0.9378-3.9;
      //run5
      else if (eval.run>24319) fShift=0.9378+4.7;
    }else if (!temp_flag_kdar_file){fShift=7.1;}
    else{fShift=0;}
    TThelp = evtTimeNS-fShift+gap*0.5;
    double TT_merged = -9999.;
    if (evtTimeNS>0) TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
//std::cout<<" "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<" "<<TT_merged<<std::endl;
    return TT_merged;
  }else if (var_name == "kdar_score_lowE"){
    return tagger.ssm_kdar_score_lowE;
  }else if (var_name == "kdar_score_hiE"){
    return tagger.ssm_kdar_score_hiE;

  }else if (var_name == "ssm_nu_angle_target"){
    if(std::isnan(tagger.ssm_nu_angle_target)) return -999;
    return tagger.ssm_nu_angle_target;
  }else if (var_name == "ssm_nu_angle_target_deg"){
    if(std::isnan(tagger.ssm_nu_angle_target)) return -999;
    if (tagger.ssm_nu_angle_target<0) return -999;
    return tagger.ssm_nu_angle_target*180/3.14159;

  }else if (var_name == "ssm_nu_angle_absorber"){
    if(std::isnan(tagger.ssm_nu_angle_absorber)) return -999;
    return tagger.ssm_nu_angle_absorber;
  }else if (var_name == "ssm_nu_angle_absorber_deg"){
    if(std::isnan(tagger.ssm_nu_angle_absorber)) return -999;
    if (tagger.ssm_nu_angle_absorber<0) return -999;
    return tagger.ssm_nu_angle_absorber*180/3.14159;

  }else if (var_name=="ssm_vtxX"){
    return tagger.ssm_vtxX;
  }else if (var_name=="ssm_vtxY"){
    return tagger.ssm_vtxY;
  }else if (var_name=="ssm_vtxZ"){
    return tagger.ssm_vtxZ;

  }else if (var_name == "reco_showerKE"){
    return get_reco_showerKE_corr(pfeval, flag_data) * 1000.;
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
  }else if (var_name == "kine_reco_Epion"){
    return get_reco_Epion(kine);
  }else if (var_name == "kine_pio_energy_1"){
    return kine.kine_pio_energy_1;
  }else if (var_name == "kine_pio_energy_2"){
    return kine.kine_pio_energy_2;
  }else if (var_name == "kine_pio_energy_max"){
    if(flag_data)
      return std::max(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::max(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }else if (var_name == "kine_pio_energy_min"){
    if(flag_data)
      return std::min(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::min(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }else if (var_name == "kine_pio_angle" || var_name == "kine_pio_costheta"){
    if (var_name == "kine_pio_angle")
      return kine.kine_pio_angle;
    else
      return TMath::Cos(kine.kine_pio_angle/180.*TMath::Pi());
  }else if (var_name == "match_energy"){
    return eval.match_energy;
  }else if (var_name == "pi0_energy"){
    double pi0_mass = 135;
    double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
    return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
  }else if (var_name == "pi0_mass"){
    if (kine.kine_pio_mass >0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      if (flag_data) {
	return kine.kine_pio_mass * em_charge_scale;
      }else{
	return kine.kine_pio_mass;
      }
    }else{
      return kine.kine_pio_mass;
    }
    //  }else if (var_name == "pi0_mass"){
    // return kine.kine_pio_mass;
  }else if (var_name == "nue_score"){
    return (tagger.nue_score<=15.99?tagger.nue_score:15.99);
  }else if (var_name == "nc_pio_score"){
    return tagger.nc_pio_score;
  }else if (var_name == "nc_delta_score"){
    return tagger.nc_delta_score;
  }else if (var_name == "numu_score"){
    return tagger.numu_score;
  }else if (var_name == "shower_energy"){
    if(flag_data)
      return tagger.mip_energy*em_charge_scale;
    else
      return tagger.mip_energy;
  }else if (var_name == "electron_energy"){
    if (flag_data){
      return pfeval.reco_showerMomentum[3] * em_charge_scale *1000;
    }else{
      return pfeval.reco_showerMomentum[3] * 1000;
    }
  }else if (var_name == "electron_polar_angle"){
    return pfeval.reco_showerMomentum[2]/pfeval.reco_showerMomentum[3];
  }else if (var_name == "shower_angle_beam"){
    return tagger.mip_angle_beam;
  }else if (var_name == "shower_angle_vertical"){
    return tagger.spt_angle_vertical;
  }else if (var_name == "shwvtx_nuvtx_dis"){
    return sqrt(pow(pfeval.reco_nuvtxX-pfeval.reco_showervtxX,2)+pow(pfeval.reco_nuvtxY-pfeval.reco_showervtxY,2)+pow(pfeval.reco_nuvtxZ-pfeval.reco_showervtxZ,2));
  }else if (var_name == "median_dQdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    return vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
  }else if (var_name == "median_dEdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    return median_dedx; // MeV/cm
  }else if (var_name == "reco_showervtxX"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxX;
    }else{
      Float_t x = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          x = pfeval.reco_startXYZT[i_p][0];
        }
      }
      return x;
    }
  }else if (var_name == "reco_showervtxY"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxY;
    }else{
      Float_t y = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          y = pfeval.reco_startXYZT[i_p][1];
        }
      }
      return y;
    }
  }else if (var_name == "reco_showervtxZ"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxZ;
    }else{
      Float_t z = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          z = pfeval.reco_startXYZT[i_p][2];
        }
      }
      return z;
    }
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
  }else if (var_name == "mip_quality_n_tracks"){
      return tagger.mip_quality_n_tracks;
  }else if (var_name == "mip_quality_n_showers"){
      return tagger.mip_quality_n_showers;
  }else if (var_name == "gap_n_bad"){
      return tagger.gap_n_bad;
  }else if (var_name == "muon_KE"){
      return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "reco_Emuon"){
      return pfeval.reco_muonMomentum[3]*1000; // GeV --> MeV
  /*}else if(var_name == "reco_Emuon_hybrid"){
      if (eval.match_isFC) {
        return 1000.0*kine.vlne_v4_numu_full_primaryE;
      }
      else {
        return 1000.0*pfeval.reco_muonMomentum[3];
      }
  }else if(var_name == "reco_Emuon_dlnew"){
      // std::cout << "vlne_v4_numu_full_primaryE: " << kine.vlne_v4_numu_full_primaryE << std::endl;
      return 1000.0*kine.vlne_v4_numu_full_primaryE;*/
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
   }else if (var_name == "reco_Eqe_muon" || var_name == "reco_Eqe_muon_Enu_diff" || var_name == "reco_Eqe_electron" || var_name == "reco_Eqe_electron_Enu_diff"){
      // everything is in MeV
      float neutron_mass = 939.57;
      float binding_energy = 30.0;
      float muon_mass = 105.66;
      float electron_mass = 0.511;
      float proton_mass = 938.27;

      float muon_costheta;
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
        muon_costheta = TMath::Cos(muonMomentum.Theta());

      float shower_costheta;
      TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);
      if (pfeval.reco_showerMomentum[3]>0)
        shower_costheta = TMath::Cos(showerMomentum.Theta());

      shower_costheta = TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi());

      float reco_Emuon =  pfeval.reco_muonMomentum[3]*1000.; // GeV --> MeV
      float reco_Eqe_muon = 0.5 * (2*(neutron_mass-binding_energy)*reco_Emuon - (pow(neutron_mass-binding_energy,2) + pow(muon_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Emuon + sqrt(pow(reco_Emuon,2)-pow(muon_mass,2))*muon_costheta);

      float reco_Eelectron= pfeval.reco_showerMomentum[3]*1000.;

      reco_Eelectron = get_reco_showerKE_corr(pfeval, flag_data) * 1000.;


      float reco_Eqe_electron = 0.5 * (2*(neutron_mass-binding_energy)*reco_Eelectron - (pow(neutron_mass-binding_energy,2) + pow(electron_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Eelectron + sqrt(pow(reco_Eelectron,2)-pow(electron_mass,2))*shower_costheta);

      if (shower_costheta > 0.999 || shower_costheta < -0.999) {
          reco_Eqe_electron = -1;
      }

      if(var_name=="reco_Eqe_muon") return reco_Eqe_muon;
      else if(var_name=="reco_Eqe_electron") return reco_Eqe_electron;
      else if(var_name=="reco_Eqe_muon_Enu_diff") return reco_Eqe_muon - get_reco_Enu_corr(kine, flag_data);
      else if(var_name=="reco_Eqe_electron_Enu_diff") return reco_Eqe_electron - get_reco_Enu_corr(kine, flag_data);
  }else if (var_name == "proton_KE"){
      return pfeval.reco_protonMomentum[3]*1000.-938.27; // GeV--> MeV
  }else if (var_name == "proton_theta"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Theta()/TMath::Pi()*180.;
  }else if (var_name == "proton_phi"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Phi()/TMath::Pi()*180.;
   }else if (var_name == "shower_theta" || var_name == "shower_costheta" || var_name == "shower_costheta" || var_name == "shower_phi"){
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(var_name == "shower_theta")
      return tagger.mip_angle_beam;
      //return showerMomentum.Theta()/TMath::Pi()*180.;

    if(var_name == "shower_costheta"){
      return TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi());
    }

    if(var_name == "shower_phi"){
      if (pfeval.reco_showerMomentum[3]>0)
        return showerMomentum.Phi()/TMath::Pi()*180.;
      else
        return -1000;
    }
   }else if (var_name=="shower_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(pfeval.reco_showerMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return showerMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else
      return -1000;
  }
  else if (var_name=="muon_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);

    if(pfeval.reco_muonMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return muonMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else
      return -1000;
  }else if (var_name == "Ehadron"){
    if (pfeval.reco_muonMomentum[3]>0)
      return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    else
      return -1000;
    //  }else if (var_name == "Ehadron"){
      /* Float_t Ehadron = kine.kine_reco_Enu; */
      /* for(size_t i=0; i<kine.kine_energy_particle->size(); i++) */
      /* { */
      /*     int pdgcode = kine.kine_particle_type->at(i); */
      /*     if(abs(pdgcode)==13) Ehadron = Ehadron - kine.kine_energy_particle->at(i) - 105.658; */
      /*     //if(abs(pdgcode)==11) Ehadron = Ehadron - kine.kine_energy_particle->at(i); */
      /* } */
    // return kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
  /*}else if (var_name == "Ehadron_hybrid"){
    if (eval.match_isFC)
      return 1000.0*(kine.vlne_v4_numu_full_totalE - kine.vlne_v4_numu_full_primaryE);
    else {
      if (pfeval.reco_muonMomentum[3]>0)
        return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
      else
        return -1000;
    }
  }else if (var_name == "Ehadron_dlnew"){
      return 1000.0*(kine.vlne_v4_numu_full_totalE - kine.vlne_v4_numu_full_primaryE);*/
  }else if (var_name == "Q2"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
    //  }else if (var_name == "Q2"){
    // Float_t Enu = kine.kine_reco_Enu;
    //Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    //Float_t Ehadron = Enu - Emu;
    //Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    //TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    //Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    //return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
  }else if (var_name == "x_Bjorken"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
    //  }else if (var_name == "x_Bjorken"){
    // Float_t Enu = kine.kine_reco_Enu;
    // Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    // Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    // TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    // Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    // return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
  }else if (var_name == "N_tracks"){
      int N_tracks = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==11) continue;
          if(kine.kine_energy_particle->at(i)<10) continue;
          if(abs(pdgcode)==13 || abs(pdgcode)==211){
            N_tracks += 1;
          }
          else if(kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_tracks += 1;
          }
      }
      return N_tracks;
   }else if (var_name == "N_other_tracks"){
        int Nothertracks = 0;
        for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
                int pdgcode = kine.kine_particle_type->at(i);
                if((abs(pdgcode)==211 || abs(pdgcode)==13) && kine.kine_energy_particle->at(i)>10) Nothertracks++; // KE threshold: 10 MeV
        }
        return Nothertracks;
   }else if (var_name == "N_showers"){
      int N_showers = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)!=11) continue;
          if(kine.kine_energy_particle->at(i)>10) N_showers += 1;
      }
      return N_showers;
  }else if (var_name == "N_protons"){
      int N_protons = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)== 2212 && kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_protons += 1;
          }
      }
      return N_protons;
  }else if (var_name == "N_true_protons"){
    if (flag_data) return -1;
    int np = 0;
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
      np++;
    }
    return np;
  }else if (var_name == "EhadShwrFrac"){
      double EhadShwr=0, EhadTot=0;

      if (pfeval.reco_muonMomentum[3]>0) {
        EhadTot = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
        for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
          if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
            EhadShwr +=  kine.kine_energy_particle->at(j);
          }

        }
        return EhadShwr/ EhadTot;

      }
      else return -1;
  }else if (var_name == "reco_mcc8_pmuoncosth_Enu"){
    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon momentum
    float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66;
    float pmuon = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66) / 1000.0;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());
    // flattened Enu
    int indx = mcc8_pmuon_costheta_bin(pmuon, costh); // index of pmuon-costh
    double reco_Enu = get_reco_Enu_corr(kine, flag_data) / 100.0;
    if (reco_Enu<0) return -10000;
    else if (reco_Enu>25.0) return 10000; // overflow bin
    else return (indx-1)*25.0 + reco_Enu;
  }
  else if (var_name == "reco_concatenated_Pmuon"){

    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon momentum
    float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66;
    float pmuon = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66) / 1000.0;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    float pmuon_MeV = pmuon * 1000.0;
    if (pmuon_MeV > 1500.0) return -10000;
    if (costh>=-1 and costh<-0.5) { return pmuon_MeV; }
    else if (costh>=-0.5 and costh<0){ return pmuon_MeV + 1500*1; }
    else if (costh>=0 and costh<0.27){ return pmuon_MeV + 1500*2; }
    else if (costh>=0.27 and costh<0.45){ return pmuon_MeV + 1500*3; }
    else if (costh>=0.45 and costh<0.62){ return pmuon_MeV + 1500*4; }
    else if (costh>=0.62 and costh<0.76){ return pmuon_MeV + 1500*5; }
    else if (costh>=0.76 and costh<0.86){ return pmuon_MeV + 1500*6; }
    else if (costh>=0.86 and costh<0.94){ return pmuon_MeV + 1500*7; }
    else if (costh>=0.94 and costh<=1.00){ return pmuon_MeV + 1500*8; }

    return -10000;
  }
  else if (var_name == "muon_momentum_costheta"){
    float muon_momentum = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_momentum");
    float costheta = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_costheta");
    int bin = alt_var_index("muon_momentum",muon_momentum, "costheta", costheta);
    return bin;
  }
  else if (var_name == "Ehad_muon_costheta"){
    float Ehadron = -1000;
    if (pfeval.reco_muonMomentum[3]>0)
      Ehadron = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    float costheta = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_costheta");
    int bin = alt_var_index("Ehadron",Ehadron, "costheta", costheta);
    return bin;
  }
  else if (var_name == "reco_concatenated_Ehad"){

    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    float Ehad = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.0;

    if (Ehad > 1500.0) return -10000;
    if (costh>=-1 and costh<-0.5) { return Ehad; }
    else if (costh>=-0.5 and costh<0){ return Ehad + 1500*1; }
    else if (costh>=0 and costh<0.27){ return Ehad + 1500*2; }
    else if (costh>=0.27 and costh<0.45){ return Ehad + 1500*3; }
    else if (costh>=0.45 and costh<0.62){ return Ehad + 1500*4; }
    else if (costh>=0.62 and costh<0.76){ return Ehad + 1500*5; }
    else if (costh>=0.76 and costh<0.86){ return Ehad + 1500*6; }
    else if (costh>=0.86 and costh<0.94){ return Ehad + 1500*7; }
    else if (costh>=0.94 and costh<=1.00){ return Ehad + 1500*8; }

    return -10000;
  }else if (var_name == "proton_pi0_total_momentum" || var_name == "proton_pi0_invariant_mass") {
        bool debug_pf_info = 0;
        TLorentzVector max_energy_proton_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_1_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_2_momentum(-1., -1., -1., -1.);
        float max_proton_energy = 0.;
        if (debug_pf_info) {
                std::cout << "********************************* starting event  ************************\n";
                std::cout << "(100 hard coded) looping over " << pfeval.reco_Ntrack << " reco particles (" << pfeval.truth_Ntrack << " true particles)\n";
        }
        for(size_t i=0; i<100; i++) // looping over all reconstructed particles
        {
                int pdgcode = pfeval.reco_pdg[i];
                if (debug_pf_info) std::cout << "investigating particle" << i << "//" << pfeval.reco_Ntrack << " : " << pdgcode << "\n";
                if(abs(pdgcode)==2212 && pfeval.reco_startMomentum[i][3] > max_proton_energy){ // new max energy proton
                        if (debug_pf_info) std::cout << "new max energy proton\n";
                        max_energy_proton_momentum = pfeval.reco_startMomentum[i];
                        max_proton_energy = pfeval.reco_startMomentum[i][3];
                }
                if(abs(pdgcode)==22 || abs(pdgcode)==11) { // reconstructed shower
                        float shower_energy = 1000. * pfeval.reco_startMomentum[i][3];
                        if (abs(shower_energy - kine.kine_pio_energy_1) / kine.kine_pio_energy_1 < 0.01) { // very close to gamma 1 energy
                                if (debug_pf_info) std::cout << "gamma 1 matched";
                                if (flag_data) {
                                        gamma_1_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][1],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][2],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][3]);
                                } else {
                                        gamma_1_momentum = pfeval.reco_startMomentum[i];
                                }
                        }
                        if (abs(shower_energy - kine.kine_pio_energy_2) / kine.kine_pio_energy_2 < 0.01) { // very close to gamma 2 energy
                                if (debug_pf_info) std::cout << "gamma 2 matched";
                                if (flag_data) {
                                        gamma_2_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][1],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][2],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][3]);
                                } else {
                                        gamma_2_momentum = pfeval.reco_startMomentum[i];
                                }
                        }

                }
        }

        float proton_pi0_invariant_mass = -1.;
        float proton_pi0_total_momentum = -1.;

        float pi0_mass = 134.9768;
        float proton_mass = 938.272;

        if (max_energy_proton_momentum[3]>0 && gamma_1_momentum[3]>0 && gamma_2_momentum[3]>0) { // found proton and pi0 in PF tree
                if (debug_pf_info) std::cout << "matched proton and pi0!" << std::endl;
                TLorentzVector pi0_momentum = gamma_1_momentum + gamma_2_momentum;
                proton_pi0_invariant_mass = sqrt(pi0_mass * pi0_mass + proton_mass * proton_mass
                                + 2. * 1000. * 1000. * (max_energy_proton_momentum[3] * pi0_momentum[3]
                                                      - max_energy_proton_momentum[0] * pi0_momentum[0]
                                                      - max_energy_proton_momentum[1] * pi0_momentum[1]
                                                      - max_energy_proton_momentum[2] * pi0_momentum[2]));
                proton_pi0_total_momentum = 1000. * sqrt((max_energy_proton_momentum[0] + pi0_momentum[0]) * (max_energy_proton_momentum[0] + pi0_momentum[0])
                                                       + (max_energy_proton_momentum[1] + pi0_momentum[1]) * (max_energy_proton_momentum[1] + pi0_momentum[1])
                                                       + (max_energy_proton_momentum[2] + pi0_momentum[2]) * (max_energy_proton_momentum[2] + pi0_momentum[2]));
        }

        if (debug_pf_info) std::cout << "******************* ending event *************\n";

        if (var_name == "proton_pi0_total_momentum") {
                return proton_pi0_total_momentum;
        } else if (var_name == "proton_pi0_invariant_mass") {
                return proton_pi0_invariant_mass;
        } else {
                std::cout << "No such proton-pi0 variable: " << var_name << std::endl;
        }
  //Erin
  }else if (var_name == "single_photon_numu_score"){
    return tagger.single_photon_numu_score;
  }else if (var_name == "single_photon_other_score"){
    return tagger.single_photon_other_score;
  }else if (var_name == "single_photon_ncpi0_score"){
    return tagger.single_photon_ncpi0_score;
  }else if (var_name == "single_photon_nue_score"){
    return tagger.single_photon_nue_score;
  }else if (var_name == "shower_energy_sp"){
    if(flag_data)
      return tagger.shw_sp_energy*em_charge_scale;
    else
      return tagger.shw_sp_energy;
  }else if (var_name == "shower_angle_beam_sp"){
    return tagger.shw_sp_angle_beam;
  }else if (var_name == "cos_shower_angle_beam_sp"){
    return TMath::Cos(tagger.shw_sp_angle_beam/180.*TMath::Pi());
  }else if (var_name == "num_shower_sp"){
    return tagger.shw_sp_n_20br1_showers;
  }else if (var_name == "median_dEdx_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_median_dedx;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_median_dedx;
  }else if (var_name == "median_dEdx_sp_15"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_2);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_3);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_4);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_5);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_6);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_7);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    if(flag_data)
      return median_dedx;//*em_charge_scale;
    else
      return median_dedx; // MeV/cm
  }else if (var_name == "dQdx_0_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_dQ_dx_0;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_dQ_dx_0;
  }else if (var_name == "dQdx_1_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_dQ_dx_1;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_dQ_dx_1;
  }else if (var_name == "shw_vtx_dis_sp"){
    return tagger.shw_sp_shw_vtx_dis;
  }else if (var_name == "max_shw_dis_sp"){
    return tagger.shw_sp_max_shw_dis;
  }else if (var_name == "shw_sp_br3_1_n_shower_segments"){
    return tagger.shw_sp_br3_1_n_shower_segments;
  }else if (var_name == "shw_sp_E_indirect_max_energy"){
    return tagger.shw_sp_E_indirect_max_energy;
  }else if (var_name == "shw_sp_lem_n_3seg"){
    return tagger.shw_sp_lem_n_3seg;
  }else if (var_name == "shw_sp_lem_n_3seg"){
    return tagger.shw_sp_lem_n_3seg;
  }else if (var_name == "kine_pio_phi_2"){
    return kine.kine_pio_phi_2;
  }else if (var_name == "shw_sp_pio_flag_pio"){
    return tagger.shw_sp_pio_flag_pio;
  }else if (var_name == "shw_sp_length_total"){
    return tagger.shw_sp_length_total;
  }else if (var_name == "shw_sp_n_vertex"){
    return tagger.shw_sp_n_vertex;
  }else if (var_name == "shw_backwards_projected_dist"){
    float backwards_projected_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_backwards_projected_dist = 1e9;

        //projecting to x walls
        if (shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
              min_backwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2];
        }
        if (isinf(min_backwards_projected_dist)) min_backwards_projected_dist = -99999.0;

        backwards_projected_dist = min_backwards_projected_dist;
      }

      return backwards_projected_dist;
  }else if (var_name == "shw_forwards_projected_dist"){
    float forwards_projected_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_forwards_projected_dist = 1e9;

        //projecting to x walls
        if (inv_shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
              min_forwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (inv_shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / inv_shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / inv_shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (inv_shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2];
        }
        if (isinf(min_forwards_projected_dist)) min_forwards_projected_dist = -99999.0;

        forwards_projected_dist = min_forwards_projected_dist;
      }

      return forwards_projected_dist;
    }else if (var_name == "shw_min_dist"){
    float minimum_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_dist = 1e9;

        //x walls
        if (sqrt(pow(pfeval.reco_showervtxX,2) - pow(-1.0,2)) < min_dist)
            min_dist =  sqrt(pow(pfeval.reco_showervtxX,2) - pow(-1.0,2));
        if (sqrt(pow(pfeval.reco_showervtxX,2) - pow(254.3,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxX,2) - pow(254.3,2));
      
        //y walls
        if (sqrt(pow(pfeval.reco_showervtxY,2) - pow(-115.0,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxY,2) - pow(-115.0,2));
        if (sqrt(pow(pfeval.reco_showervtxY,2) - pow(117.0,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxY,2) - pow(117.0,2));
        
        //projecting to z walls
        if (sqrt(pow(pfeval.reco_showervtxZ,2) - pow(0.6,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxZ,2) - pow(0.6,2));
        if (sqrt(pow(pfeval.reco_showervtxZ,2) - pow(1036.4,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxZ,2) - pow(1036.4,2));
        
        if (isinf(min_dist)) min_dist = -99999.0;

        minimum_dist = min_dist;
      }

      return minimum_dist;
    }else if (var_name == "shw_projected_dist"){
      float projected_dist = -99999.0;
      float forwards_projected_dist = -99999.0;
      float backwards_projected_dist = -99999.0;

      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);

        //forwards projected dist
        float min_forwards_projected_dist = 1e9;

        //projecting to x walls
        if (inv_shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
              min_forwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (inv_shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / inv_shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / inv_shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (inv_shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2];
        }
        if (isinf(min_forwards_projected_dist)) min_forwards_projected_dist = -99999.0;

        forwards_projected_dist = min_forwards_projected_dist;

        //backwards projected dist
        float min_backwards_projected_dist = 1e9;

        //projecting to x walls
        if (shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
              min_backwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2];
        }
        if (isinf(min_backwards_projected_dist)) min_backwards_projected_dist = -99999.0;

        backwards_projected_dist = min_backwards_projected_dist;

        if (forwards_projected_dist < backwards_projected_dist) projected_dist = forwards_projected_dist;
        else projected_dist = backwards_projected_dist;
      }

      return projected_dist;
  }else if (var_name == "shw_inwardness_3d"){
    }else if (var_name == "shw_backwards_projected_dist"){
      float inwardness = -99999.0;
    
      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);
        inwardness = inwardness_3d;
      }
      return inwardness;
  }else if (var_name == "shw_inwardness_2d"){
    }else if (var_name == "shw_backwards_projected_dist"){
      float inwardness = -99999.0;
    
      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);
        inwardness = inwardness_2d;
      }
      return inwardness;
  }else if (var_name == "ns_beam_time"){
    if(flag_data){
      double delta_time_calc = -9999.;
      //Merge Peaks
      double gap=18.936;
      double Shift=0;
      double TThelp=0;
      if (pfeval.run >= 17380){ Shift=2916.0; }
      else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
      else if (pfeval.run >= 10812){ Shift = 3568.5; }
      else if (pfeval.run >= 8321){ Shift = 3610.7;}
      else if (pfeval.run >= 5800){ Shift = 3164.4;}
      else if (pfeval.run >= 0){ Shift = 3168.9;}
      //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
      /*if (pfeval.run >= 13697){ Shift = 3166.9;}
      else if(pfeval.run>=10812){ Shift = 3568.5; }
      else if (pfeval.run >= 8321){ Shift = 3610.7;}
      else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}*/
      //9.43;}
      //if(run>8000 && run<10812){Shift=3610.7; }
      //if(run>=10812 && run <12500){Shift=3568.5; }
      TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
      double TT_merged = -9999.;

      //merge peaks
      if(TThelp>=0 && TThelp<gap*81.0){
        TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
      }

      delta_time_calc = TT_merged;

      return delta_time_calc;
    }else{
      return -9999.;
    }
  }else if (var_name == "flag_0p"){
    if (is_0p(tagger, kine, pfeval)){
      return 1;
    }else{
      return 0;
    }
  }else if (var_name == "flag_FC"){
    if (is_FC(eval)){
      return 1;
    }else{
      return 0;
    }
  }else if (var_name == "run_period"){
    if (pfeval.run >= 13697){ 
      return 3;
    }
    else if (pfeval.run >= 8321){ 
      return 2;
    }
    else if (pfeval.run > 0 ){ 
      return 1;
    }
  }else if (var_name == "run_period_100_300"){
    float shwen = -1;
    if(flag_data)
      shwen = tagger.shw_sp_energy*em_charge_scale;
    else
      shwen = tagger.shw_sp_energy;
    
    if (shwen >= 100.0 && shwen < 300){
      if (pfeval.run >= 13697){ 
        return 3;
      }
      else if (pfeval.run >= 8321){ 
        return 2;
      }
      else if (pfeval.run > 0 ){ 
        return 1;
      }
    }
    return -1;
    
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

int get_Pmuon_bin (float Pmuon) {
  int nbins = 6;
  float pmuon_binning[nbins+1] = {0, 180, 300, 450, 770, 1280, 2500};
  if (Pmuon == pmuon_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (Pmuon >  pmuon_binning[i] && Pmuon <= pmuon_binning[i+1]) { return i; } }
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
  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    TString cut_name = it->first;
    int number = it->second;

    double KE_muon = pfeval.truth_muonMomentum[3]*1000.-105.66; // MeV
    double pmuon   = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66); // MeV
    double Pmuon   = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66); // MeV
    double Emuon   = pfeval.truth_muonMomentum[3]*1000; // MeV
    double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

    float pmuon_binning[7] = {0, 180, 300, 450, 770, 1280, 2500};

    float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};			//fine binning
    //float costheta_binning[5]  = {-1,         .27,      .62,      .86,      1};		//coarse binning
    //float costheta_binning[3]    = {-1,                   .62,                1};		//very coarse binning
    TLorentzVector muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    double theta_absorber = get_angle_to_absorber(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2]);
    float costh_absorber = TMath::Cos(theta_absorber);
    float sinth_absorber = TMath::Sin(theta_absorber);
    double pl_absorber = pmuon*costh_absorber;
    double pt_absorber = pmuon*sinth_absorber;
    double q_absorber = sqrt(235.5*235.2+pmuon*pmuon-2*235.5*pmuon*costh_absorber);
    double Q2_absorber = sqrt(2*235.5*(Emuon-pmuon*costh_absorber)-105.7*105.7 );

    double Kp = get_true_Kp(pfeval);
    double p_angle_absorber_deg = get_true_p_angle_absorber(pfeval)*180/3.14159;

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
    else if (cut_file == 2) {
      if (cut_name == "numuCC.inside.Emuon.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=100 && Emuon>0) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=200 && Emuon>100) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=300 && Emuon>200) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=400 && Emuon>300) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=500 && Emuon>400) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=600 && Emuon>500) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=700 && Emuon>600) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=800 && Emuon>700) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=900 && Emuon>800) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1000 && Emuon>900) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1200.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1200 && Emuon>1000) return number;
      }
      else if (cut_name == "numuCC.inside.Emuon.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1200) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 3) {
      if (cut_name == "numuCC.inside.Ehadron.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=200 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=300 && Ehadron>200) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=400 && Ehadron>300) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=500 && Ehadron>400) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=600 && Ehadron>500) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=700 && Ehadron>600) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=800 && Ehadron>700) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=900 && Ehadron>800) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1000 && Ehadron>900) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1000) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 4){
      bool pre_cut = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && muonMomentum[3]>0 && Pmuon > 0 && Pmuon <= 2500;
      if (cut_name == "numuCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
	if (pre_cut && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){
	if (pre_cut && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.805.gt.705"){
	if (pre_cut && eval.truth_nuEnergy<=805 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.920.gt.805"){
	if (pre_cut && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>805) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1050.gt.920"){
	if (pre_cut && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1050"){
	if (pre_cut && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1375.gt.1200"){
	if (pre_cut && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1375"){
	if (pre_cut && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){
	if (pre_cut && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
	if (pre_cut && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 5) {
      if (cut_name == "numuCC.inside.Emuon.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 6) {
      if (cut_name == "numuCC.inside.Ehadron.le.100.gt.30"){ // 30 MeV - 100 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100 && Ehadron>30) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.150.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=150 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.225.gt.150"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=225 && Ehadron>150) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.275.gt.225"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=275 && Ehadron>225) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.336.gt.275"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=336 && Ehadron>275) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.411.gt.336"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=411 && Ehadron>336) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.502.gt.411"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=502 && Ehadron>411) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.614.gt.502"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=614 && Ehadron>502) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.750.gt.614"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=750 && Ehadron>614) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1120.gt.750"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1120 && Ehadron>750) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.2500.gt.1120"){ // 1120 - 2500 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1120 && Ehadron<= 2500) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 7){
      if (cut_name == "numuCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.920.gt.705"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.920"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 8){
      if (cut_name == "numuCC.inside.Enu.le.4000.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 9){
      if (cut_name == "numuCC.inside.Enu.le.1200.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>1200) return number;
      }
      else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 10) {
      if (cut_name == "numuCC.inside.Ehadron.le.150.gt.30"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=150 && Ehadron>30) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.275.gt.150"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=275 && Ehadron>150) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.411.gt.275"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=411 && Ehadron>275) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.502.gt.411"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=502 && Ehadron>411) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.614.gt.502"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=614 && Ehadron>502) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.750.gt.614"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=750 && Ehadron>614) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1120.gt.750"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1120 && Ehadron>750) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.2500.gt.1120"){ // 1120 - 2500 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1120 && Ehadron<= 2500) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 11){
      if (cut_name == "nueCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.705.gt.540"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.920.gt.705"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.1200.gt.920"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.1570.gt.1200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.2050.gt.1570"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 12){
      if (cut_name == "nueCC.inside.Enu.le.1200.gt.200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.4000.gt.1200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>1200) return number;
      }
      else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }

    //MCC8 binning -> 40 bins, muon mometum
    else if (cut_file == 13){
      if  (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;

      }else{
	std::cout << "get_xs_signal_no: no cut found!   cut_name = " << cut_name << std::endl;
      }
    }

    //MCC8 binning -> 36 bins, muon mometum
    else if (cut_file==14) {

      if       (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;

      }else{
	std::cout << "get_xs_signal_no: no cut found!   cut_name = " << cut_name << std::endl;
      }

    }

    //very coarse angle binning
    else if (cut_file == 15){
      if (number==-1) { std::cout << "cut_name, number = " << cut_name << ", " << number << std::endl; }
      if       (cut_name == "numuCC.inside.Emuon.theta0.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }
    }

    else if (cut_file == 16){ // pmuon, costh

      if (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<2500.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=0.00 && costh<0.27 ) return number;
        }
        else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<1280.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=1280.00 && pmuon<2500.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<1280.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=1280.00 && pmuon<2500.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else{
      std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }

    else if (cut_file == 17){ // Enu, costheta, Pmuon

      int Enu_bin      = get_Enu_bin(eval.truth_nuEnergy);
      int costheta_bin = get_costheta_bin(costh);
      int Pmuon_bin    = get_Pmuon_bin(Pmuon);
      bool pre_cut     = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && (muonMomentum[3]>0);

      if      (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.450.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.450.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else { std::cout << "get_xs_signal_no: no cut found!" << std::endl; }

    }

    //1D Enu truth binning using the same inclusive selection as the 3D binning
    else if (cut_file == 18){

      int Enu_bin      = get_Enu_bin(eval.truth_nuEnergy);
      int costheta_bin = get_costheta_bin(costh);
      int Pmuon_bin    = get_Pmuon_bin(Pmuon);
      bool pre_cut     = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && (muonMomentum[3]>0) && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8 && Pmuon_bin>=0 && Pmuon_bin<=5;

      if      (cut_name == "numuCC.inside.Enu.le.540.gt.200"){   if (pre_cut && eval.truth_nuEnergy<=540  && eval.truth_nuEnergy>200)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){   if (pre_cut && eval.truth_nuEnergy<=705  && eval.truth_nuEnergy>540)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.805.gt.705"){   if (pre_cut && eval.truth_nuEnergy<=805  && eval.truth_nuEnergy>705)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.920.gt.805"){   if (pre_cut && eval.truth_nuEnergy<=920  && eval.truth_nuEnergy>805)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1050.gt.920"){  if (pre_cut && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1050"){ if (pre_cut && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1375.gt.1200"){ if (pre_cut && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1375"){ if (pre_cut && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){ if (pre_cut && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ if (pre_cut && eval.truth_nuEnergy>2050  && eval.truth_nuEnergy<=4000) { return number; } }
      else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl; }
    }

    else if (cut_file==769){
      int bin_width=10;
      int min = 30;
      int max = 120;
      int nbins = int(max/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.nuwro.Emu.le.";
      if(!is_true_kdar(eval,pfeval)  || eval.event%10<5 ){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && KE_muon>max)   { return number; }
        cut_string = base_cut_string+std::to_string(min+bin_width)+".gt."+std::to_string(min);
        if(cut_name == cut_string && KE_muon<min)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==770){
      int bin_width=10;
      int min = 30;
      int max = 120;
      int nbins = int(max/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Emu.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && KE_muon>max)   { return number; }
        cut_string = base_cut_string+std::to_string(min+bin_width)+".gt."+std::to_string(min);
        if(cut_name == cut_string && KE_muon<min)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==771){
      int bin_width=20;
      int min = 0;
      int max = 120;
      int nbins = int(max/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Emu.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && KE_muon>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }


    else if (cut_file==772){
      int bin_width=10;
      int min = 0;
      int max = 130;
      int nbins = int(max/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.nuwro.Emu.le.";
      if(!is_true_kdar(eval,pfeval)  || eval.event%10<5 ){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && KE_muon>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }


    else if (cut_file==773){
      int bin_width=10;
      int min = 0;
      int max = 130;
      int nbins = int(max/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Emu.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && KE_muon>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==774){
      bool found_cut = false;
      std::string cut_string = "kdar.gibuu.tot";
      if(!is_true_kdar(eval,pfeval) || eval.event%10<5){found_cut=true;}
      else{
        if(cut_name == cut_string){
          found_cut = true;
          return number;
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }
    else if (cut_file==775){
      bool found_cut = false;
      std::string cut_string = "kdar.nuwro.tot";
      if(!is_true_kdar(eval,pfeval) || eval.event%10<5){found_cut=true;}
      else{
        if(cut_name == cut_string){
          found_cut = true;
          return number;
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }
    else if (cut_file==776){
      bool found_cut = false;
      std::string cut_string = "kdar.tot";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{         
        if(cut_name == cut_string){
          found_cut = true;
          return number; 
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }
    else if (cut_file==777){
      int bin_width=5;
      int min = 0;
      int max = 130; 
      int nbins = int(max/bin_width);  
      bool found_cut = false;   
      std::string base_cut_string = "kdar.Emu.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
	// Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(bin*bin_width+bin_width)+".gt."+std::to_string(bin*bin_width);
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
	  if(cut_name == cut_string){
	    found_cut = true;
	    if(KE_muon<=bin*bin_width+bin_width  && KE_muon>bin*bin_width)   { return number; } 
          }
  	}
	// Double check the overflow
	std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
	if(cut_name == cut_string && KE_muon>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }


    else if (cut_file==778){
      int bin_width=10;
      int min = -200;
      int max = 200;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.pl.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(pl_absorber<=bin*bin_width+bin_width+min  && pl_absorber>bin*bin_width+min)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && pl_absorber>max)   { return number; }
        // Double check the underflow
        cut_string = base_cut_string+std::to_string(min+bin_width)+".gt."+std::to_string(min);
        if(cut_name == cut_string && pl_absorber<min)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==779){
      int bin_width=8;
      int min = 0;
      int max = 200;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.pt.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(pt_absorber<=bin*bin_width+bin_width+min  && pt_absorber>bin*bin_width+min)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && pt_absorber>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==780){
      int bin_width=12;
      int min = 0;
      int max = 480;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.q.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(q_absorber<=bin*bin_width+bin_width+min  && q_absorber>bin*bin_width+min)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && q_absorber>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==781){
      int bin_width=12;
      int min = 0;
      int max = 420;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Q2.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(Q2_absorber<=bin*bin_width+bin_width+min  && Q2_absorber>bin*bin_width+min)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && Q2_absorber>max)   { return number; }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==782){
      int bin_width=5;
      int min = 0;
      int max = 130;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Kp.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(Kp<bin*bin_width+bin_width+min && Kp>=bin*bin_width+min)   { return number; }
          }
        }
        // Double check the overflow
        std::string cut_string = base_cut_string+std::to_string(max)+".gt."+std::to_string(max-bin_width);
        if(cut_name == cut_string && Kp>max)   { return number; }
        // Double check the underflow
        cut_string = base_cut_string+std::to_string(min+bin_width)+".gt."+std::to_string(min);
        if(cut_name == cut_string && Kp<min)   { return number; }        
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }
    else if (cut_file==783){
      int slice_width=12;
      int slice_min = 0;
      int slice_max = 120;
      int slice_nbins = int( (slice_max-slice_min)/slice_width);
      int bin_width=8;
      int min = 0;
      int max = 128;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Kp.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
	// Check all slices
	for(int slice_bin=0; slice_bin<slice_nbins; slice_bin++){
          std::string cut_string_slice = base_cut_string+std::to_string(abs(slice_bin*slice_width+slice_width+slice_min))+".gt."+std::to_string(abs(slice_bin*slice_width+slice_min))+".Kmu.le.";
          // Check all the bins
          for(int bin=0; bin<nbins; bin++){
            std::string cut_string = cut_string_slice+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
            //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
            if(cut_name == cut_string){
              found_cut = true;
              if(Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min)   { return number; }
	      //check overflow and underflow for the slice while looping bins
              if(slice_bin==slice_nbins-1 && Kp>slice_max && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
              if(slice_bin==0 && Kp<slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
	    }
          }
          // Double check the overflow
          std::string cut_string = cut_string_slice+std::to_string(max)+".gt."+std::to_string(max-bin_width);
          if(cut_name == cut_string && KE_muon>max && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          if(cut_name == cut_string && KE_muon>max && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          if(cut_name == cut_string && KE_muon>max && slice_bin==0 && Kp<slice_min){ return number; }
          // Double check the underflow
          cut_string = cut_string_slice+std::to_string(min+bin_width)+".gt."+std::to_string(min);
          if(cut_name == cut_string && KE_muon<min && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          if(cut_name == cut_string && KE_muon<min && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          if(cut_name == cut_string && KE_muon<min && slice_bin==0 && Kp<slice_min){ return number; }
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==784){
      int i_slice_width=125;
      int i_slice_min = -1000;
      int i_slice_max = 1000;
      int slice_nbins = int( (i_slice_max-i_slice_min)/i_slice_width);
      double slice_width = double(i_slice_width)/1000;
      double slice_min = double(i_slice_min)/1000;
      double slice_max = double(i_slice_max)/1000;
      int bin_width=8;
      int min = 0;
      int max = 128;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.CosMu.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all slices
        for(int slice_bin=0; slice_bin<slice_nbins; slice_bin++){
          std::string cut_string_slice = base_cut_string+std::to_string(abs(slice_bin*i_slice_width+i_slice_width+i_slice_min))+".gt."+std::to_string(abs(slice_bin*i_slice_width+i_slice_min))+".Kmu.le.";
          // Check all the bins
          for(int bin=0; bin<nbins; bin++){
            std::string cut_string = cut_string_slice+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
            //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
            if(cut_name == cut_string){
              found_cut = true;
              if(costh_absorber<=slice_bin*slice_width+slice_width+slice_min && costh_absorber>=slice_bin*slice_width+slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min)   { return number; }
              //check overflow and underflow for the slice while looping bins
              if(slice_bin==slice_nbins-1 && costh_absorber>slice_max && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
              if(slice_bin==0 && costh_absorber<slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
            }
          }
          // Double check the overflow
          std::string cut_string = cut_string_slice+std::to_string(max)+".gt."+std::to_string(max-bin_width);
          if(cut_name == cut_string && KE_muon>max && costh_absorber<slice_bin*slice_width+slice_width+slice_min && costh_absorber>=slice_bin*slice_width+slice_min)   { return number; }
          // Double check the underflow
          cut_string = cut_string_slice+std::to_string(min+bin_width)+".gt."+std::to_string(min);
          if(cut_name == cut_string && KE_muon<min && costh_absorber<slice_bin*slice_width+slice_width+slice_min && costh_absorber>=slice_bin*slice_width+slice_min)   { return number; }
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }
    else if (cut_file==785){
      int slice_width=8;
      int slice_min = 0;
      int slice_max = 128;
      int slice_nbins = int( (slice_max-slice_min)/slice_width);
      int bin_width=8;
      int min = 0;
      int max = 128;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Kp.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all slices
        for(int slice_bin=0; slice_bin<slice_nbins; slice_bin++){
          std::string cut_string_slice = base_cut_string+std::to_string(abs(slice_bin*slice_width+slice_width+slice_min))+".gt."+std::to_string(abs(slice_bin*slice_width+slice_min))+".Kmu.le.";
          // Check all the bins
          for(int bin=0; bin<nbins; bin++){
            std::string cut_string = cut_string_slice+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
            //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
            if(cut_name == cut_string){
              found_cut = true;
              if(Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min)   { return number; }
              //check overflow and underflow for the slice while looping bins
              if(slice_bin==slice_nbins-1 && Kp>slice_max && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
              if(slice_bin==0 && Kp<slice_min && KE_muon<bin*bin_width+bin_width+min && KE_muon>=bin*bin_width+min){ return number; }
            }
          }
          // Double check the overflow
          std::string cut_string = cut_string_slice+std::to_string(max)+".gt."+std::to_string(max-bin_width);
          if(cut_name == cut_string && KE_muon>max && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          if(cut_name == cut_string && KE_muon>max && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          if(cut_name == cut_string && KE_muon>max && slice_bin==0 && Kp<slice_min){ return number; }
          // Double check the underflow
          cut_string = cut_string_slice+std::to_string(min+bin_width)+".gt."+std::to_string(min);
          if(cut_name == cut_string && KE_muon<min && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          if(cut_name == cut_string && KE_muon<min && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          if(cut_name == cut_string && KE_muon<min && slice_bin==0 && Kp<slice_min){ return number; }
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==786){
      int bin_width=5;
      int min = 0;
      int max = 180;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.ThetaP.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all the bins
        for(int bin=0; bin<nbins; bin++){
          std::string cut_string = base_cut_string+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
          //std::cout<<cut_string<<"  "<<cut_name<<" "<<p_angle_absorber_deg<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
          if(cut_name == cut_string){
            found_cut = true;
            if(p_angle_absorber_deg<=bin*bin_width+bin_width+min && p_angle_absorber_deg>bin*bin_width+min)   { return number; }
          }
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }

    else if (cut_file==787){
      int slice_width=8;
      int slice_min = 0;
      int slice_max = 128;
      int slice_nbins = int( (slice_max-slice_min)/slice_width);
      int bin_width=10;
      int min = 0;
      int max = 180;
      int nbins = int( (max-min)/bin_width);
      bool found_cut = false;
      std::string base_cut_string = "kdar.Kp.le.";
      if(!is_true_kdar(eval,pfeval)){found_cut=true;}
      else{
        // Check all slices
        for(int slice_bin=0; slice_bin<slice_nbins; slice_bin++){
          std::string cut_string_slice = base_cut_string+std::to_string(abs(slice_bin*slice_width+slice_width+slice_min))+".gt."+std::to_string(abs(slice_bin*slice_width+slice_min))+".ThetaP.le.";
          // Check all the bins
          for(int bin=0; bin<nbins; bin++){
            std::string cut_string = cut_string_slice+std::to_string(abs(bin*bin_width+bin_width+min))+".gt."+std::to_string(abs(bin*bin_width+min));
            //std::cout<<cut_string<<"  "<<cut_name<<" "<<KE_muon<<" "<<is_true_kdar(eval,pfeval)<<std::endl;
            if(cut_name == cut_string){
              found_cut = true;
              if(Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min && p_angle_absorber_deg<bin*bin_width+bin_width+min && p_angle_absorber_deg>=bin*bin_width+min)   { return number; }
              //check overflow and underflow for the slice while looping bins
              if(slice_bin==slice_nbins-1 && Kp>slice_max && p_angle_absorber_deg<bin*bin_width+bin_width+min && p_angle_absorber_deg>=bin*bin_width+min){ return number; }
              if(slice_bin==0 && Kp<slice_min && p_angle_absorber_deg<bin*bin_width+bin_width+min && p_angle_absorber_deg>=bin*bin_width+min){ return number; }
            }
          }
          // Double check the overflow
          std::string cut_string = cut_string_slice+std::to_string(max)+".gt."+std::to_string(max-bin_width);
          //if(cut_name == cut_string && KE_muon>max && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          //if(cut_name == cut_string && KE_muon>max && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          //if(cut_name == cut_string && KE_muon>max && slice_bin==0 && Kp<slice_min){ return number; }
          // Double check the underflow
          cut_string = cut_string_slice+std::to_string(min+bin_width)+".gt."+std::to_string(min);
          //if(cut_name == cut_string && p_angle_absorber_deg<min && Kp<slice_bin*slice_width+slice_width+slice_min && Kp>=slice_bin*slice_width+slice_min)   { return number; }
          //if(cut_name == cut_string && KE_muon<min && slice_bin==slice_nbins-1 && Kp>slice_max){ return number; }
          if(cut_name == cut_string && p_angle_absorber_deg<min && slice_bin==0 && Kp<slice_min){ return number; }
        }
      }
      if(!found_cut) std::cout << "get_xs_signal_no: no cut found! " << cut_name <<std::endl;
    }


  }

  return -1;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){


  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  double KE_muon = pfeval.truth_muonMomentum[3]*1000.-105.66; // MeV
  double Pmuon   = (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));

  double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
  double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

  TLorentzVector truth_muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);

  bool flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  // definition of additional cuts
  std::map<std::string, bool> map_cuts_flag;
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

  if(eval.truth_nuEnergy <=400) map_cuts_flag["LowEintnueCC"] = true;
  else map_cuts_flag["LowEintnueCC"] = false;

  if (!(eval.truth_nuEnergy <=400)) map_cuts_flag["antiLowEintnueCC"] = true;
  else map_cuts_flag["antiLowEintnueCC"] = false;

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

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0) map_cuts_flag["XsecNC"] = true;
  else map_cuts_flag["XsecNC"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["XsecCosmic"] = true;
  else map_cuts_flag["XsecCosmic"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==1 && !(eval.truth_nuPdg==14 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200)) map_cuts_flag["XsecBkgCC"] = true;
  else map_cuts_flag["XsecBkgCC"] = false;

  // finish Xs related cuts ...
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200) map_cuts_flag["Xs_Enu_nueCCinFV"] = true;
  else map_cuts_flag["Xs_Enu_nueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["nueCCinFV"] = true;
  else map_cuts_flag["nueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["RnueCCinFV"] = true;
  else map_cuts_flag["RnueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["AnueCCinFV"] = true;
  else map_cuts_flag["AnueCCinFV"] = false;

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

  //Erin
  // breakdown categories for single photon analysis
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && eval.truth_vtxInside==1) map_cuts_flag["SPNCDeltaSig"] = true;
  else map_cuts_flag["SPNCDeltaSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && (pfeval.truth_showerMother==111) && eval.truth_vtxInside==1) map_cuts_flag["SPNCPi0Sig"] = true;
  else map_cuts_flag["SPNCPi0Sig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_showerMother!=111 && pfeval.truth_NCDelta==0 && eval.truth_vtxInside==1) map_cuts_flag["SPNCOtherSig"] = true;
  else map_cuts_flag["SPNCOtherSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1 && eval.truth_vtxInside==1) map_cuts_flag["SPNumuCCSig"] = true;
  else map_cuts_flag["SPNumuCCSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && (eval.truth_isCC==0 || (eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1)) && eval.truth_vtxInside==0) map_cuts_flag["SPOutFVSig"] = true;
  else map_cuts_flag["SPOutFVSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  pfeval.reco_muonMomentum[3] > 0) map_cuts_flag["muon"] = true;
  else map_cuts_flag["muon"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  !(pfeval.reco_muonMomentum[3] > 0)) map_cuts_flag["nomuon"] = true;
  else map_cuts_flag["nomuon"] = false;

  //for testing, no FV check
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1) map_cuts_flag["SPNCDeltaSigNoFV"] = true;
  else map_cuts_flag["SPNCDeltaSigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && (pfeval.truth_showerMother==111)) map_cuts_flag["SPNCPi0SigNoFV"] = true;
  else map_cuts_flag["SPNCPi0SigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_showerMother!=111 && pfeval.truth_NCDelta==0) map_cuts_flag["SPNCOtherSigNoFV"] = true;
  else map_cuts_flag["SPNCOtherSigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1) map_cuts_flag["SPNumuCCSigNoFV"] = true;
  else map_cuts_flag["SPNumuCCSigNoFV"] = false;
  //

  map_cuts_flag["SPdirtBkg"] = false;
  map_cuts_flag["SPoutFVBkg"] = false;
  map_cuts_flag["SPnumuCCBkg"] = false;
  map_cuts_flag["SPnumuCCpi0Bkg"] = false;
  map_cuts_flag["SPNCBkg"] = false;
  map_cuts_flag["SPNCpi0Bkg"] = false;

  if(!(map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] || map_cuts_flag["SPNumuCCSig"] || map_cuts_flag["SPOutFVSig"])){
      map_cuts_flag["SPdirtBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["SPoutFVBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_Npi0==0) map_cuts_flag["SPnumuCCBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_Npi0>0) map_cuts_flag["SPnumuCCpi0Bkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_Npi0==0) map_cuts_flag["SPNCBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_Npi0>0) map_cuts_flag["SPNCpi0Bkg"] = true;
  }

  if((map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] || map_cuts_flag["SPNumuCCSig"] || map_cuts_flag["SPOutFVSig"])){
      if (is_true_0p(pfeval)==1) {
        map_cuts_flag["SP0p"] = true;
        map_cuts_flag["SPNp"] = false;
      }
      else {
        map_cuts_flag["SPNp"] = true;
        map_cuts_flag["SP0p"] = false;
      }
  }else{
      map_cuts_flag["SPNp"] = false;
      map_cuts_flag["SP0p"] = false;
  }
  // done with single photon breakdown categories
  //


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

  if(is_true_kdar(eval,pfeval)) map_cuts_flag["kdar"] = true;
  else map_cuts_flag["kdar"] = false;
  if(is_true_kdar_outFV(eval,pfeval)) map_cuts_flag["kdarOutFV"] = true;
  else map_cuts_flag["kdaroutFV"] = false;
  if(is_true_kdar(eval,pfeval) && get_true_Kp(pfeval)>0) map_cuts_flag["kdarNp"] = true;
  else map_cuts_flag["kdarNp"] = false;
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

  bool flag_generic = is_generic(eval);
  bool flag_numuCC = is_numuCC(tagger);
  //bool flag_numuCC = is_numuCC(tagger) && (is_far_sideband(kine, tagger, flag_data) || is_near_sideband(kine, tagger, flag_data) );
  bool flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  bool flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  bool flag_numuCC_lowEhad = is_numuCC_lowEhad(tagger, kine, pfeval, flag_data);
  bool flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  bool flag_nueCC = is_nueCC(tagger);
  bool flag_nueCC_loose = is_loosenueCC(tagger);

  bool flag_0p = is_0p(tagger, kine, pfeval);
  bool flag_1p = is_1p(tagger, kine, pfeval);
  bool flag_0pi = is_0pi(tagger, kine, pfeval);

  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);
  bool flag_FC = is_FC(eval);

  //bool flag_ncpio_sel = is_NCpio_bdt(tagger) && (!flag_0p);
  bool flag_ncpio_sel = is_NCpio_sel(tagger, kine);
  bool flag_ncdelta_sel = is_NCdelta_sel(tagger, pfeval);

  //Erin
  bool flag_singlephoton_sel = is_singlephoton_sel(tagger, pfeval);
  bool flag_singlephoton_eff_sel = is_singlephoton_eff_sel(tagger, pfeval);
  bool flag_singleshower_sel = is_singleshower_sel(tagger, pfeval);
  bool flag_singleshower_eff_sel = is_singleshower_eff_sel(tagger, pfeval);
  bool flag_singlephoton_numu_sel = is_singlephoton_numu_sel(tagger, pfeval);
  bool flag_singlephoton_other_sel = is_singlephoton_other_sel(tagger, pfeval);
  bool flag_singlephoton_ncpi0_sel = is_singlephoton_ncpi0_sel(tagger, pfeval);
  bool flag_singlephoton_nue_sel = is_singlephoton_nue_sel(tagger, pfeval);
  bool flag_singlephoton_nue_sel_allshw = is_singlephoton_nue_sel_allshw(tagger, pfeval);
  bool flag_nsbeam = is_nsbeam_photon(pfeval, eval); //set all cuts to shifted
  bool flag_nsbeam_photon = is_nsbeam_photon(pfeval, eval);
  bool flag_singlephoton_pre = is_singlephoton_pre(tagger, pfeval);
  bool flag_singlephoton_numu = is_singlephoton_numu(tagger, pfeval);
  bool flag_singlephoton_other = is_singlephoton_other(tagger, pfeval);
  bool flag_singlephoton_ncpi0 = is_singlephoton_ncpi0(tagger, pfeval);
  bool flag_singlephoton_nue = is_singlephoton_nue(tagger, pfeval);
  bool flag_singlephoton_eff_numu = is_singlephoton_eff_numu(tagger, pfeval);
  bool flag_singlephoton_eff_other = is_singlephoton_eff_other(tagger, pfeval);
  bool flag_singlephoton_eff_ncpi0 = is_singlephoton_eff_ncpi0(tagger, pfeval);
  bool flag_singlephoton_eff_nue = is_singlephoton_eff_nue(tagger, pfeval);
  bool flag_singlephoton_oneshw = is_singlephoton_oneshw(tagger, pfeval);

  bool flag_kdar_presel = is_kdar_presel(tagger,eval);
  bool flag_kdar_bdtsel = is_kdar_bdtsel(tagger);
  bool flag_kdar_lowE_bdtsel = is_lowE_kdar_bdtsel(tagger);
  bool flag_kdar_hiE_bdtsel = is_hiE_kdar_bdtsel(tagger);
  //

  float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};		// PeLEE binning
  //float costheta_binning[7]  = {-1,         .27,      .62, .76, .86, .94, 1};		// coarse binning
  //float costheta_binning[3]    = {-1,                   .62,                1};	//very coarse binning
  TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
  float reco_pmuon = TMath::Sqrt(pow(pfeval.reco_muonMomentum[0],2)+pow(pfeval.reco_muonMomentum[1],2)+pow(pfeval.reco_muonMomentum[2],2))*1000;

  int costheta_bin = get_costheta_bin(TMath::Cos(muonMomentum.Theta()));
  int Pmu_bin      = get_Pmuon_bin(reco_pmuon);
  int Enu_bin      = get_Enu_bin(reco_Enu);

  bool flag_kdar_file = false;
  if(ch_name == "kdar_bdtsel_sig" || ch_name == "kdar_bdtsel_outFV_sig" || ch_name == "kdar_bdtsel_sigNp" || ch_name == "kdar_bdtselns_fake_bck"
  || ch_name == "kdar_bdtselns_sig" || ch_name == "kdar_bdtselns_outFV_sig" || ch_name == "kdar_bdtselns_sigNp" 
  || ch_name == "kdar_bdtsel_nuwro_train" || ch_name == "kdar_bdtsel_gibuu_train" || ch_name == "kdar_bdtsel_nuwro_trainNp" || ch_name == "kdar_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_bdtselns_nuwro_train" || ch_name == "kdar_bdtselns_gibuu_train" || ch_name == "kdar_bdtselns_nuwro_trainNp" || ch_name == "kdar_bdtselns_gibuu_trainNp"

  || ch_name == "kdar_lowE_bdtsel_sig" || ch_name == "kdar_lowE_bdtsel_outFV_sig" || ch_name == "kdar_lowE_bdtsel_sigNp"
  || ch_name == "kdar_lowE_bdtselns_sig" || ch_name == "kdar_lowE_bdtselns_outFV_sig" || ch_name == "kdar_lowE_bdtselns_sigNp"
  || ch_name == "kdar_lowE_bdtsel_nuwro_train" || ch_name == "kdar_lowE_bdtsel_gibuu_train" || ch_name == "kdar_lowE_bdtsel_nuwro_trainNp" || ch_name == "kdar_lowE_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_lowE_bdtselns_nuwro_train" || ch_name == "kdar_lowE_bdtselns_gibuu_train" || ch_name == "kdar_lowE_bdtselns_nuwro_trainNp" || ch_name == "kdar_lowE_bdtselns_gibuu_trainNp"

  || ch_name == "kdar_hiE_bdtsel_sig" || ch_name == "kdar_hiE_bdtsel_outFV_sig" || ch_name == "kdar_hiE_bdtsel_sigNp"
  || ch_name == "kdar_hiE_bdtselns_sig" || ch_name == "kdar_hiE_bdtselns_outFV_sig" || ch_name == "kdar_hiE_bdtselns_sigNp"
  || ch_name == "kdar_hiE_bdtsel_nuwro_train" || ch_name == "kdar_hiE_bdtsel_gibuu_train" || ch_name == "kdar_hiE_bdtsel_nuwro_trainNp" || ch_name == "kdar_hiE_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_hiE_bdtselns_nuwro_train" || ch_name == "kdar_hiE_bdtselns_gibuu_train" || ch_name == "kdar_hiE_bdtselns_nuwro_trainNp" || ch_name == "kdar_hiE_bdtselns_gibuu_trainNp"

  || ch_name == "kdar_presel_sig" || ch_name == "kdar_presel_outFV_sig" || ch_name == "kdar_presel_sigNp"
  || ch_name == "kdar_preselns_sig" || ch_name == "kdar_preselns_outFV_sig" || ch_name == "kdar_preselns_sigNp" 
  || ch_name == "kdar_presel_nuwro_train" || ch_name == "kdar_presel_gibuu_train" || ch_name == "kdar_presel_nuwro_trainNp" || ch_name == "kdar_presel_gibuu_trainNp"
  || ch_name == "kdar_preselns_nuwro_train" || ch_name == "kdar_preselns_gibuu_train" || ch_name == "kdar_preselns_nuwro_trainNp" || ch_name == "kdar_preselns_gibuu_trainNp"

  || ch_name == "kdar_ktagged_sig" || ch_name == "kdar_ktagged_outFV_sig" || ch_name == "kdar_ktagged_sigNp"
  || ch_name == "kdar_ktaggedns_sig" || ch_name == "kdar_ktaggedns_outFV_sig" || ch_name == "kdar_ktaggedns_sigNp"
  || ch_name == "kdar_ktagged_nuwro_train" || ch_name == "kdar_ktagged_gibuu_train" || ch_name == "kdar_ktagged_nuwro_trainNp" || ch_name == "kdar_ktagged_gibuu_trainNp"
  || ch_name == "kdar_ktaggedns_nuwro_train" || ch_name == "kdar_ktaggedns_gibuu_train" || ch_name == "kdar_ktaggedns_nuwro_trainNp" || ch_name == "kdar_ktaggedns_gibuu_trainNp"

  || ch_name == "kdar_bdtselNp_sig" || ch_name == "kdar_bdtselNpns_sig" || ch_name == "kdar_preselNp_sig" || ch_name == "kdar_preselNpns_sig" || ch_name == "kdar_ktaggedNp_sig" || ch_name == "kdar_ktaggedNpns_sig")
  {flag_kdar_file=true;}
  double merge_time = 0;//
  if(!(ch_name=="kdar_sideband_overlap1_bck" || ch_name=="kdar_sideband_overlap2_bck" || ch_name=="kdar_sideband_overlap3_bck" || ch_name=="kdar_sideband_overlap4_bck" || ch_name == "kdar_bdtselns_nons_bck")) merge_time = get_kine_var(kine, eval, pfeval, tagger, flag_data, "merge_time_recover_numi",flag_kdar_file);
  if(ch_name == "kdar_bdtselnsrand_dirt" || ch_name == "kdar_hiE_bdtselnsrand_dirt" || ch_name == "kdar_lowE_bdtselnsrand_dirt" || ch_name == "kdar_preselnsrand_dirt" || ch_name == "kdar_ktaggednsrand_dirt" || ch_name == "kdar_bdtselNpnsrand_dirt" || ch_name == "kdar_preselNprand_dirt" || ch_name == "kdar_ktaggedNpnsrand_dirt") merge_time = get_kine_var(kine, eval, pfeval, tagger, flag_data, "merge_time_random");
  //std::cout<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<" "<<merge_time<<std::endl;








  if(ch_name == "kdar_bdtsel_bck"  || ch_name == "kdar_bdtsel_dirt" || ch_name == "kdar_bdtsel_bckNp"
  || ch_name == "kdar_bdtsel_sig" || ch_name == "kdar_bdtsel_outFV_sig" || ch_name == "kdar_bdtsel_sigNp"
  || ch_name == "kdar_bdtsel_nuwro_train" || ch_name == "kdar_bdtsel_gibuu_train" || ch_name == "kdar_bdtsel_nuwro_trainNp" || ch_name == "kdar_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_bdtsel" || ch_name == "kdar_bdtsel_ext" 

  || ch_name == "kdar_bdtselns_bck"  || ch_name == "kdar_bdtselns_dirt" || ch_name == "kdar_bdtselnsrand_dirt" || ch_name == "kdar_bdtselns_bckNp" || ch_name == "kdar_bdtselns_nons_bck" || ch_name == "kdar_bdtselns_fake_bck"
  || ch_name == "kdar_bdtselns_sig" || ch_name == "kdar_bdtselns_outFV_sig" || ch_name == "kdar_bdtselns_sigNp"
  || ch_name == "kdar_bdtselns_nuwro_train" || ch_name == "kdar_bdtselns_gibuu_train" || ch_name == "kdar_bdtselns_nuwro_trainNp" || ch_name == "kdar_bdtselns_gibuu_trainNp"
  || ch_name == "kdar_bdtselns" || ch_name == "kdar_bdtselns_ext"

  || ch_name == "kdar_bdtselNp_sig" || ch_name == "kdar_bdtselNp_bck" || ch_name == "kdar_bdtselNp_dirt" || ch_name == "kdar_bdtselNp_ext" || ch_name == "kdar_bdtselNp" 
  || ch_name == "kdar_bdtselNpns_sig" || ch_name == "kdar_bdtselNpns_bck" || ch_name == "kdar_bdtselNpns_dirt" || ch_name == "kdar_bdtselNpns_ext" || ch_name == "kdar_bdtselNpns" || ch_name == "kdar_bdtselNpnsrand_dirt"){
    bool flag_pass = flag_kdar_bdtsel && flag_kdar_presel;

    if((ch_name == "kdar_bdtselns_bck"  || ch_name == "kdar_bdtselns_dirt" || ch_name == "kdar_bdtselnsrand_dirt" || ch_name == "kdar_bdtselns_bckNp" || ch_name == "kdar_bdtselns_nons_bck"
    || ch_name == "kdar_bdtselns_sig" || ch_name == "kdar_bdtselns_outFV_sig" || ch_name == "kdar_bdtselns_sigNp"
    || ch_name == "kdar_bdtselns_nuwro_train" || ch_name == "kdar_bdtselns_gibuu_train" || ch_name == "kdar_bdtselns_nuwro_trainNp" || ch_name == "kdar_bdtselns_gibuu_trainNp"
    || ch_name == "kdar_bdtselNpns_dirt" || ch_name == "kdar_bdtselNpnsrand_dirt"
    //|| ch_name == "kdar_bdtselns" || ch_name == "kdar_bdtselns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    //|| ch_name == "kdar_bdtselns" || ch_name == "kdar_bdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}
    || ch_name == "kdar_bdtselns" || ch_name == "kdar_bdtselns_ext" || ch_name == "kdar_bdtselns_fake_bck") && ((merge_time>-4.5 && merge_time<4.5) )){flag_pass=false;}

    if(ch_name == "kdar_bdtsel_sig" || ch_name == "kdar_bdtsel_nuwro_train" || ch_name == "kdar_bdtsel_gibuu_train" || ch_name == "kdar_bdtselns_sig" || ch_name == "kdar_bdtselns_nuwro_train" || ch_name == "kdar_bdtselns_gibuu_train" || ch_name == "kdar_bdtselNp_sig" || ch_name == "kdar_bdtselNpns_sig")  flag_pass = flag_pass && map_cuts_flag["kdar"];

    if(ch_name == "kdar_bdtsel_outFV_sig" || ch_name == "kdar_bdtselns_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];

    if(ch_name == "kdar_bdtsel_sigNp" || ch_name == "kdar_bdtsel_nuwro_trainNp" || ch_name == "kdar_bdtsel_gibuu_trainNp" || ch_name == "kdar_bdtselns_sigNp" || ch_name == "kdar_bdtselns_nuwro_trainNp" || ch_name == "kdar_bdtselns_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];

    if(ch_name == "kdar_bdtsel_bck" || ch_name == "kdar_bdtselns_bck" || ch_name == "kdar_bdtselns_nons_bck" || ch_name == "kdar_bdtselNp_bck" || ch_name == "kdar_bdtselNpns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }

    if(ch_name == "kdar_bdtsel_bckNp" || ch_name == "kdar_bdtselns_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }

    if(ch_name == "kdar_bdtsel_nuwro_train" || ch_name == "kdar_bdtsel_gibuu_train"|| ch_name == "kdar_bdtsel_nuwro_trainNp" || ch_name == "kdar_bdtsel_gibuu_trainNp" || ch_name == "kdar_bdtselns_nuwro_train" || ch_name == "kdar_bdtselns_gibuu_train"|| ch_name == "kdar_bdtselns_nuwro_trainNp" || ch_name == "kdar_bdtselns_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }

    if(ch_name == "kdar_bdtselNp_sig" || ch_name == "kdar_bdtselNp_bck" || ch_name == "kdar_bdtselNp_dirt" || ch_name == "kdar_bdtselNp_ext" || ch_name == "kdar_bdtselNp"        
  || ch_name == "kdar_bdtselNpns_sig" || ch_name == "kdar_bdtselNpns_bck" || ch_name == "kdar_bdtselNpns_dirt" || ch_name == "kdar_bdtselNpns_ext" || ch_name == "kdar_bdtselNpns" || ch_name == "kdar_bdtselNpnsrand_dirt"){if(tagger.ssm_prim_track1_kine_energy_range<0) flag_pass=false;}

    return flag_pass;



  }else if(ch_name == "kdar_lowE_bdtsel_bck"  || ch_name == "kdar_lowE_bdtsel_dirt" || ch_name == "kdar_lowE_bdtsel_bckNp"
  || ch_name == "kdar_lowE_bdtsel_sig" || ch_name == "kdar_lowE_bdtsel_outFV_sig" || ch_name == "kdar_lowE_bdtsel_sigNp"
  || ch_name == "kdar_lowE_bdtsel_nuwro_train" || ch_name == "kdar_lowE_bdtsel_gibuu_train" || ch_name == "kdar_lowE_bdtsel_nuwro_trainNp" || ch_name == "kdar_lowE_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_lowE_bdtsel" || ch_name == "kdar_lowE_bdtsel_ext"

  || ch_name == "kdar_lowE_bdtselns_bck"  || ch_name == "kdar_lowE_bdtselns_dirt" || ch_name == "kdar_lowE_bdtselns_bckNp" || ch_name == "kdar_lowE_bdtselnsrand_dirt"
  || ch_name == "kdar_lowE_bdtselns_sig" || ch_name == "kdar_lowE_bdtselns_outFV_sig" || ch_name == "kdar_lowE_bdtselns_sigNp"
  || ch_name == "kdar_lowE_bdtselns_nuwro_train" || ch_name == "kdar_lowE_bdtselns_gibuu_train" || ch_name == "kdar_lowE_bdtselns_nuwro_trainNp" || ch_name == "kdar_lowE_bdtselns_gibuu_trainNp"
  || ch_name == "kdar_lowE_bdtselns" || ch_name == "kdar_lowE_bdtselns_ext"){

    bool flag_pass = flag_kdar_lowE_bdtsel && flag_kdar_presel;

    if((ch_name == "kdar_lowE_bdtselns_bck"  || ch_name == "kdar_lowE_bdtselns_dirt" || ch_name == "kdar_lowE_bdtselnsrand_dirt" || ch_name == "kdar_lowE_bdtselns_bckNp"
    || ch_name == "kdar_lowE_bdtselns_sig" || ch_name == "kdar_lowE_bdtselns_outFV_sig" || ch_name == "kdar_lowE_bdtselns_sigNp"
    || ch_name == "kdar_lowE_bdtselns_nuwro_train" || ch_name == "kdar_lowE_bdtselns_gibuu_train" || ch_name == "kdar_lowE_bdtselns_nuwro_trainNp" || ch_name == "kdar_lowE_bdtselns_gibuu_trainNp"
    //|| ch_name == "kdar_lowE_bdtselns" || ch_name == "kdar_lowE_bdtselns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    //|| ch_name == "kdar_lowE_bdtselns" || ch_name == "kdar_lowE_bdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}
    || ch_name == "kdar_lowE_bdtselns" || ch_name == "kdar_lowE_bdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5)) ){flag_pass=false;}

    if(ch_name == "kdar_lowE_bdtsel_sig" || ch_name == "kdar_lowE_bdtsel_nuwro_train" || ch_name == "kdar_lowE_bdtsel_gibuu_train" || ch_name == "kdar_lowE_bdtselns_sig" || ch_name == "kdar_lowE_bdtselns_nuwro_train" || ch_name == "kdar_lowE_bdtselns_gibuu_train")  flag_pass = flag_pass && map_cuts_flag["kdar"];
    if(ch_name == "kdar_lowE_bdtsel_outFV_sig" || ch_name == "kdar_lowE_bdtselns_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];
    if(ch_name == "kdar_lowE_bdtsel_sigNp" || ch_name == "kdar_lowE_bdtsel_nuwro_trainNp" || ch_name == "kdar_lowE_bdtsel_gibuu_trainNp" || ch_name == "kdar_lowE_bdtselns_sigNp" || ch_name == "kdar_lowE_bdtselns_nuwro_trainNp" || ch_name == "kdar_lowE_bdtselns_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];
    if(ch_name == "kdar_lowE_bdtsel_bck" || ch_name == "kdar_lowE_bdtselns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }
    if(ch_name == "kdar_lowE_bdtsel_bckNp" || ch_name == "kdar_lowE_bdtselns_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }
    if(ch_name == "kdar_lowE_bdtsel_nuwro_train" || ch_name == "kdar_lowE_bdtsel_gibuu_train"|| ch_name == "kdar_lowE_bdtsel_nuwro_trainNp" || ch_name == "kdar_lowE_bdtsel_gibuu_trainNp" || ch_name == "kdar_lowE_bdtselns_nuwro_train" || ch_name == "kdar_lowE_bdtselns_gibuu_train"|| ch_name == "kdar_lowE_bdtselns_nuwro_trainNp" || ch_name == "kdar_lowE_bdtselns_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }
    return flag_pass;




  }else if(ch_name == "kdar_hiE_bdtsel_bck"  || ch_name == "kdar_hiE_bdtsel_dirt" || ch_name == "kdar_hiE_bdtsel_bckNp"
  || ch_name == "kdar_hiE_bdtsel_sig" || ch_name == "kdar_hiE_bdtsel_outFV_sig" || ch_name == "kdar_hiE_bdtsel_sigNp"
  || ch_name == "kdar_hiE_bdtsel_nuwro_train" || ch_name == "kdar_hiE_bdtsel_gibuu_train" || ch_name == "kdar_hiE_bdtsel_nuwro_trainNp" || ch_name == "kdar_hiE_bdtsel_gibuu_trainNp"
  || ch_name == "kdar_hiE_bdtsel" || ch_name == "kdar_hiE_bdtsel_ext"

  || ch_name == "kdar_hiE_bdtselns_bck"  || ch_name == "kdar_hiE_bdtselns_dirt" || ch_name == "kdar_hiE_bdtselns_bckNp" || ch_name == "kdar_hiE_bdtselnsrand_dirt"
  || ch_name == "kdar_hiE_bdtselns_sig" || ch_name == "kdar_hiE_bdtselns_outFV_sig" || ch_name == "kdar_hiE_bdtselns_sigNp"
  || ch_name == "kdar_hiE_bdtselns_nuwro_train" || ch_name == "kdar_hiE_bdtselns_gibuu_train" || ch_name == "kdar_hiE_bdtselns_nuwro_trainNp" || ch_name == "kdar_hiE_bdtselns_gibuu_trainNp"
  || ch_name == "kdar_hiE_bdtselns" || ch_name == "kdar_hiE_bdtselns_ext"){

    bool flag_pass = flag_kdar_hiE_bdtsel && flag_kdar_presel;

    if((ch_name == "kdar_hiE_bdtselns_bck"  || ch_name == "kdar_hiE_bdtselns_dirt" || ch_name == "kdar_hiE_bdtselns_bckNp" || ch_name == "kdar_hiE_bdtselnsrand_dirt"
    || ch_name == "kdar_hiE_bdtselns_sig" || ch_name == "kdar_hiE_bdtselns_outFV_sig" || ch_name == "kdar_hiE_bdtselns_sigNp"
    || ch_name == "kdar_hiE_bdtselns_nuwro_train" || ch_name == "kdar_hiE_bdtselns_gibuu_train" || ch_name == "kdar_hiE_bdtselns_nuwro_trainNp" || ch_name == "kdar_hiE_bdtselns_gibuu_trainNp"
    //|| ch_name == "kdar_hiE_bdtselns" || ch_name == "kdar_hiE_bdtselns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    //|| ch_name == "kdar_hiE_bdtselns" || ch_name == "kdar_hiE_bdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}
    || ch_name == "kdar_hiE_bdtselns" || ch_name == "kdar_hiE_bdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5) )){flag_pass=false;}

    if(ch_name == "kdar_hiE_bdtsel_sig" || ch_name == "kdar_hiE_bdtsel_nuwro_train" || ch_name == "kdar_hiE_bdtsel_gibuu_train" || ch_name == "kdar_hiE_bdtselns_sig" || ch_name == "kdar_hiE_bdtselns_nuwro_train" || ch_name == "kdar_hiE_bdtselns_gibuu_train")  flag_pass = flag_pass && map_cuts_flag["kdar"];
    if(ch_name == "kdar_hiE_bdtsel_outFV_sig" || ch_name == "kdar_hiE_bdtselns_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];
    if(ch_name == "kdar_hiE_bdtsel_sigNp" || ch_name == "kdar_hiE_bdtsel_nuwro_trainNp" || ch_name == "kdar_hiE_bdtsel_gibuu_trainNp" || ch_name == "kdar_hiE_bdtselns_sigNp" || ch_name == "kdar_hiE_bdtselns_nuwro_trainNp" || ch_name == "kdar_hiE_bdtselns_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];
    if(ch_name == "kdar_hiE_bdtsel_bck" || ch_name == "kdar_hiE_bdtselns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }
    if(ch_name == "kdar_hiE_bdtsel_bckNp" || ch_name == "kdar_hiE_bdtselns_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }
    if(ch_name == "kdar_hiE_bdtsel_nuwro_train" || ch_name == "kdar_hiE_bdtsel_gibuu_train"|| ch_name == "kdar_hiE_bdtsel_nuwro_trainNp" || ch_name == "kdar_hiE_bdtsel_gibuu_trainNp" || ch_name == "kdar_hiE_bdtselns_nuwro_train" || ch_name == "kdar_hiE_bdtselns_gibuu_train"|| ch_name == "kdar_hiE_bdtselns_nuwro_trainNp" || ch_name == "kdar_hiE_bdtselns_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }
    return flag_pass;


  }else if(ch_name == "kdar_presel_bck"  || ch_name == "kdar_presel_dirt" || ch_name == "kdar_presel_bckNp"
  || ch_name == "kdar_presel_sig" || ch_name == "kdar_presel_outFV_sig" || ch_name == "kdar_presel_sigNp"
  || ch_name == "kdar_presel_nuwro_train" || ch_name == "kdar_presel_gibuu_train" || ch_name == "kdar_presel_nuwro_trainNp" || ch_name == "kdar_presel_gibuu_trainNp"
  || ch_name == "kdar_presel" || ch_name == "kdar_presel_ext"

  || ch_name == "kdar_preselns_bck"  || ch_name == "kdar_preselns_dirt" || ch_name == "kdar_preselns_bckNp" || ch_name == "kdar_preselnsrand_dirt"
  || ch_name == "kdar_preselns_sig" || ch_name == "kdar_preselns_outFV_sig" || ch_name == "kdar_preselns_sigNp"
  || ch_name == "kdar_preselns_nuwro_train" || ch_name == "kdar_preselns_gibuu_train" || ch_name == "kdar_preselns_nuwro_trainNp" || ch_name == "kdar_preselns_gibuu_trainNp"
  || ch_name == "kdar_preselns" || ch_name == "kdar_preselns_ext"

  || ch_name == "kdar_preselNp_sig" || ch_name == "kdar_preselNp_bck" || ch_name == "kdar_preselNp_dirt" || ch_name == "kdar_preselNp_ext" || ch_name == "kdar_preselNp"
  || ch_name == "kdar_preselNpns_sig" || ch_name == "kdar_preselNpns_bck" || ch_name == "kdar_preselNpns_dirt" || ch_name == "kdar_preselNpns_ext" || ch_name == "kdar_preselNpns" || ch_name == "kdar_preselNpnsrand_dirt"){

    bool flag_pass = flag_kdar_presel;

    if((ch_name == "kdar_preselns_bck"  || ch_name == "kdar_preselns_dirt" || ch_name == "kdar_preselns_bckNp" || ch_name == "kdar_preselnsrand_dirt"
    || ch_name == "kdar_preselns_sig" || ch_name == "kdar_preselns_outFV_sig" || ch_name == "kdar_preselns_sigNp"
    || ch_name == "kdar_preselns_nuwro_train" || ch_name == "kdar_preselns_gibuu_train" || ch_name == "kdar_preselns_nuwro_trainNp" || ch_name == "kdar_preselns_gibuu_trainNp"
      || ch_name == "kdar_preselNpns_sig" || ch_name == "kdar_preselNpns_bck" || ch_name == "kdar_preselNpns_dirt" || ch_name == "kdar_preselNpns_ext" || ch_name == "kdar_preselNpns" || ch_name == "kdar_preselNpnsrand_dirt"
    //|| ch_name == "kdar_preselns" || ch_name == "kdar_preselns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    || ch_name == "kdar_preselns" || ch_name == "kdar_preselns_ext") && ((merge_time>-4.5 && merge_time<4.5) ) ){flag_pass=false;}
    //|| ch_name == "kdar_preselns" || ch_name == "kdar_preselns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}

    if(ch_name == "kdar_presel_sig" || ch_name == "kdar_presel_nuwro_train" || ch_name == "kdar_presel_gibuu_train" || ch_name == "kdar_preselNp_sig" || ch_name == "kdar_preselNpns_sig")  flag_pass = flag_pass && map_cuts_flag["kdar"];
    if(ch_name == "kdar_presel_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];
    if(ch_name == "kdar_presel_sigNp" || ch_name == "kdar_presel_nuwro_trainNp" || ch_name == "kdar_presel_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];
    if(ch_name == "kdar_presel_bck" || ch_name == "kdar_preselNp_bck" || ch_name == "kdar_preselNpns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }
    if(ch_name == "kdar_presel_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }
    if(ch_name == "kdar_presel_nuwro_train" || ch_name == "kdar_presel_gibuu_train"|| ch_name == "kdar_presel_nuwro_trainNp" || ch_name == "kdar_presel_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }
    return flag_pass;


  }else if(ch_name == "kdar_sideband_bck"  || ch_name == "kdar_sideband_dirt" || ch_name == "kdar_sideband_sig" || ch_name == "kdar_sideband" || ch_name == "kdar_sideband_ext"
     || ch_name == "kdar_sideband1_bck"  || ch_name == "kdar_sideband1_dirt" || ch_name == "kdar_sideband1_sig" || ch_name == "kdar_sideband1" || ch_name == "kdar_sideband1_ext"
     || ch_name == "kdar_sideband2_bck"  || ch_name == "kdar_sideband2_dirt" || ch_name == "kdar_sideband2_sig" || ch_name == "kdar_sideband2" || ch_name == "kdar_sideband2_ext"
     || ch_name == "kdar_sideband3_bck"  || ch_name == "kdar_sideband3_dirt" || ch_name == "kdar_sideband3_sig" || ch_name == "kdar_sideband3" || ch_name == "kdar_sideband3_ext"
     || ch_name == "kdar_sideband4_bck"  || ch_name == "kdar_sideband4_dirt" || ch_name == "kdar_sideband4_sig" || ch_name == "kdar_sideband4" || ch_name == "kdar_sideband4_ext"){
    bool flag_pass = tagger.numu_score>0.9 && eval.match_isFC==1 && !(flag_kdar_bdtsel && flag_kdar_presel);
    if( (ch_name == "kdar_sideband1_bck"  || ch_name == "kdar_sideband1_dirt" || ch_name == "kdar_sideband1_sig" || ch_name == "kdar_sideband1" || ch_name == "kdar_sideband1_ext") && reco_Enu>400) flag_pass = false;
    if( (ch_name == "kdar_sideband2_bck"  || ch_name == "kdar_sideband2_dirt" || ch_name == "kdar_sideband2_sig" || ch_name == "kdar_sideband2" || ch_name == "kdar_sideband2_ext") && (reco_Enu>700 || reco_Enu<400)) flag_pass = false;
    if( (ch_name == "kdar_sideband3_bck"  || ch_name == "kdar_sideband3_dirt" || ch_name == "kdar_sideband3_sig" || ch_name == "kdar_sideband3" || ch_name == "kdar_sideband3_ext") && (reco_Enu>1000 || reco_Enu<700)) flag_pass = false;
    if( (ch_name == "kdar_sideband4_bck"  || ch_name == "kdar_sideband4_dirt" || ch_name == "kdar_sideband4_sig" || ch_name == "kdar_sideband4" || ch_name == "kdar_sideband4_ext") && reco_Enu<1000) flag_pass = false;
    return flag_pass;


  }else if(ch_name == "kdar_sideband_beam_bck"  || ch_name == "kdar_sideband_beam_dirt" || ch_name == "kdar_sideband_beam_sig" || ch_name == "kdar_sideband_beam" || ch_name == "kdar_sideband_beam_ext"
     || ch_name == "kdar_sideband_beam1_bck"  || ch_name == "kdar_sideband_beam1_dirt" || ch_name == "kdar_sideband_beam1_sig" || ch_name == "kdar_sideband_beam1" || ch_name == "kdar_sideband_beam1_ext"
     || ch_name == "kdar_sideband_beam2_bck"  || ch_name == "kdar_sideband_beam2_dirt" || ch_name == "kdar_sideband_beam2_sig" || ch_name == "kdar_sideband_beam2" || ch_name == "kdar_sideband_beam2_ext"
     || ch_name == "kdar_sideband_beam3_bck"  || ch_name == "kdar_sideband_beam3_dirt" || ch_name == "kdar_sideband_beam3_sig" || ch_name == "kdar_sideband_beam3" || ch_name == "kdar_sideband_beam3_ext"
     || ch_name == "kdar_sideband_beam4_bck"  || ch_name == "kdar_sideband_beam4_dirt" || ch_name == "kdar_sideband_beam4_sig" || ch_name == "kdar_sideband_beam4" || ch_name == "kdar_sideband_beam4_ext"){
    bool flag_pass = tagger.numu_score>0.9 && eval.match_isFC==1 && !(flag_kdar_bdtsel && flag_kdar_presel);
    if(!(merge_time>-4.5 && merge_time<4.5)){flag_pass=false;}
    if( (ch_name == "kdar_sideband_beam1_bck"  || ch_name == "kdar_sideband_beam1_dirt" || ch_name == "kdar_sideband_beam1_sig" || ch_name == "kdar_sideband_beam1" || ch_name == "kdar_sideband_beam1_ext") && reco_Enu>400) flag_pass = false;
    if( (ch_name == "kdar_sideband_beam2_bck"  || ch_name == "kdar_sideband_beam2_dirt" || ch_name == "kdar_sideband_beam2_sig" || ch_name == "kdar_sideband_beam2" || ch_name == "kdar_sideband_beam2_ext") && (reco_Enu>700 || reco_Enu<400)) flag_pass = false;
    if( (ch_name == "kdar_sideband_beam3_bck"  || ch_name == "kdar_sideband_beam3_dirt" || ch_name == "kdar_sideband_beam3_sig" || ch_name == "kdar_sideband_beam3" || ch_name == "kdar_sideband_beam3_ext") && (reco_Enu>1000 || reco_Enu<700)) flag_pass = false;
    if( (ch_name == "kdar_sideband_beam4_bck"  || ch_name == "kdar_sideband_beam4_dirt" || ch_name == "kdar_sideband_beam4_sig" || ch_name == "kdar_sideband_beam4" || ch_name == "kdar_sideband_beam4_ext") && reco_Enu<1000) flag_pass = false;
    return flag_pass;

  }else if(ch_name == "kdar_sideband_decay_bck"  || ch_name == "kdar_sideband_decay_dirt" || ch_name == "kdar_sideband_decay_sig" || ch_name == "kdar_sideband_decay" || ch_name == "kdar_sideband_decay_ext"
     || ch_name == "kdar_sideband_decay1_bck"  || ch_name == "kdar_sideband_decay1_dirt" || ch_name == "kdar_sideband_decay1_sig" || ch_name == "kdar_sideband_decay1" || ch_name == "kdar_sideband_decay1_ext"
     || ch_name == "kdar_sideband_decay2_bck"  || ch_name == "kdar_sideband_decay2_dirt" || ch_name == "kdar_sideband_decay2_sig" || ch_name == "kdar_sideband_decay2" || ch_name == "kdar_sideband_decay2_ext"
     || ch_name == "kdar_sideband_decay3_bck"  || ch_name == "kdar_sideband_decay3_dirt" || ch_name == "kdar_sideband_decay3_sig" || ch_name == "kdar_sideband_decay3" || ch_name == "kdar_sideband_decay3_ext"
     || ch_name == "kdar_sideband_decay4_bck"  || ch_name == "kdar_sideband_decay4_dirt" || ch_name == "kdar_sideband_decay4_sig" || ch_name == "kdar_sideband_decay4" || ch_name == "kdar_sideband_decay4_ext"){
    bool flag_pass = tagger.numu_score>0.9 && eval.match_isFC==1 && !(flag_kdar_bdtsel && flag_kdar_presel);
    if(merge_time>-4.5 && merge_time<4.5){flag_pass=false;}
    if( (ch_name == "kdar_sideband_decay1_bck"  || ch_name == "kdar_sideband_decay1_dirt" || ch_name == "kdar_sideband_decay1_sig" || ch_name == "kdar_sideband_decay1" || ch_name == "kdar_sideband_decay1_ext") && reco_Enu>400) flag_pass = false;
    if( (ch_name == "kdar_sideband_decay2_bck"  || ch_name == "kdar_sideband_decay2_dirt" || ch_name == "kdar_sideband_decay2_sig" || ch_name == "kdar_sideband_decay2" || ch_name == "kdar_sideband_decay2_ext") && (reco_Enu>700 || reco_Enu<400)) flag_pass = false;
    if( (ch_name == "kdar_sideband_decay3_bck"  || ch_name == "kdar_sideband_decay3_dirt" || ch_name == "kdar_sideband_decay3_sig" || ch_name == "kdar_sideband_decay3" || ch_name == "kdar_sideband_decay3_ext") && (reco_Enu>1000 || reco_Enu<700)) flag_pass = false;
    if( (ch_name == "kdar_sideband_decay4_bck"  || ch_name == "kdar_sideband_decay4_dirt" || ch_name == "kdar_sideband_decay4_sig" || ch_name == "kdar_sideband_decay4" || ch_name == "kdar_sideband_decay4_ext") && reco_Enu<1000) flag_pass = false;
    return flag_pass;



  }else if(ch_name == "kdar_sideband_overlap_bck"  || ch_name == "kdar_sideband_overlap_dirt" || ch_name == "kdar_sideband_overlap_sig" || ch_name == "kdar_sideband_overlap" || ch_name == "kdar_sideband_overlap_ext"
     || ch_name == "kdar_sideband_overlap1_bck"  || ch_name == "kdar_sideband_overlap1_dirt" || ch_name == "kdar_sideband_overlap1_sig" || ch_name == "kdar_sideband_overlap1" || ch_name == "kdar_sideband_overlap1_ext"
     || ch_name == "kdar_sideband_overlap2_bck"  || ch_name == "kdar_sideband_overlap2_dirt" || ch_name == "kdar_sideband_overlap2_sig" || ch_name == "kdar_sideband_overlap2" || ch_name == "kdar_sideband_overlap2_ext"
     || ch_name == "kdar_sideband_overlap3_bck"  || ch_name == "kdar_sideband_overlap3_dirt" || ch_name == "kdar_sideband_overlap3_sig" || ch_name == "kdar_sideband_overlap3" || ch_name == "kdar_sideband_overlap3_ext"
     || ch_name == "kdar_sideband_overlap4_bck"  || ch_name == "kdar_sideband_overlap4_dirt" || ch_name == "kdar_sideband_overlap4_sig" || ch_name == "kdar_sideband_overlap4" || ch_name == "kdar_sideband_overlap4_ext"){
    bool flag_pass = tagger.numu_score>0.9 && eval.match_isFC==1;
    if( (ch_name == "kdar_sideband_overlap1_bck"  || ch_name == "kdar_sideband_overlap1_dirt" || ch_name == "kdar_sideband_overlap1_sig" || ch_name == "kdar_sideband_overlap1" || ch_name == "kdar_sideband_overlap1_ext") && reco_Enu>400) flag_pass = false;
    if( (ch_name == "kdar_sideband_overlap2_bck"  || ch_name == "kdar_sideband_overlap2_dirt" || ch_name == "kdar_sideband_overlap2_sig" || ch_name == "kdar_sideband_overlap2" || ch_name == "kdar_sideband_overlap2_ext") && (reco_Enu>700 || reco_Enu<400)) flag_pass = false;
    if( (ch_name == "kdar_sideband_overlap3_bck"  || ch_name == "kdar_sideband_overlap3_dirt" || ch_name == "kdar_sideband_overlap3_sig" || ch_name == "kdar_sideband_overlap3" || ch_name == "kdar_sideband_overlap3_ext") && (reco_Enu>1000 || reco_Enu<700)) flag_pass = false;
    if( (ch_name == "kdar_sideband_overlap4_bck"  || ch_name == "kdar_sideband_overlap4_dirt" || ch_name == "kdar_sideband_overlap4_sig" || ch_name == "kdar_sideband_overlap4" || ch_name == "kdar_sideband_overlap4_ext") && reco_Enu<1000) flag_pass = false;
    return flag_pass;




  }else if(ch_name == "kdar_ktagged_bck"  || ch_name == "kdar_ktagged_dirt" || ch_name == "kdar_ktagged_bckNp"
  || ch_name == "kdar_ktagged_sig" || ch_name == "kdar_ktagged_outFV_sig" || ch_name == "kdar_ktagged_sigNp"
  || ch_name == "kdar_ktagged_nuwro_train" || ch_name == "kdar_ktagged_gibuu_train" || ch_name == "kdar_ktagged_nuwro_trainNp" || ch_name == "kdar_ktagged_gibuu_trainNp"
  || ch_name == "kdar_ktagged" || ch_name == "kdar_ktagged_ext"

  || ch_name == "kdar_ktaggedns_bck"  || ch_name == "kdar_ktaggedns_dirt" || ch_name == "kdar_ktaggedns_bckNp" || ch_name == "kdar_ktaggednsrand_dirt"
  || ch_name == "kdar_ktaggedns_sig" || ch_name == "kdar_ktaggedns_outFV_sig" || ch_name == "kdar_ktaggedns_sigNp"
  || ch_name == "kdar_ktaggedns_nuwro_train" || ch_name == "kdar_ktaggedns_gibuu_train" || ch_name == "kdar_ktaggedns_nuwro_trainNp" || ch_name == "kdar_ktaggedns_gibuu_trainNp"
  || ch_name == "kdar_ktaggedns" || ch_name == "kdar_ktaggedns_ext"

  || ch_name == "kdar_ktaggedNp_sig" || ch_name == "kdar_ktaggedNp_bck" || ch_name == "kdar_ktaggedNp_dirt" || ch_name == "kdar_ktaggedNp_ext" || ch_name == "kdar_ktaggedNp"
  || ch_name == "kdar_ktaggedNpns_sig" || ch_name == "kdar_ktaggedNpns_bck" || ch_name == "kdar_ktaggedNpns_dirt" || ch_name == "kdar_ktaggedNpns_ext" || ch_name == "kdar_ktaggedNpns" || ch_name == "kdar_ktaggedNpnsrand_dirt"){

    bool flag_pass = false;
    if(tagger.ssm_kine_energy>0) flag_pass = true;

    if((ch_name == "kdar_ktaggedns_bck"  || ch_name == "kdar_ktaggedns_dirt" || ch_name == "kdar_ktaggedns_bckNp" || ch_name == "kdar_ktaggednsrand_dirt"
    || ch_name == "kdar_ktaggedns_sig" || ch_name == "kdar_ktaggedns_outFV_sig" || ch_name == "kdar_ktaggedns_sigNp"
    || ch_name == "kdar_ktaggedns_nuwro_train" || ch_name == "kdar_ktaggedns_gibuu_train" || ch_name == "kdar_ktaggedns_nuwro_trainNp" || ch_name == "kdar_ktaggedns_gibuu_trainNp"
      || ch_name == "kdar_ktaggedNpns_sig" || ch_name == "kdar_ktaggedNpns_bck" || ch_name == "kdar_ktaggedNpns_dirt" || ch_name == "kdar_ktaggedNpns_ext" || ch_name == "kdar_ktaggedNpns" || ch_name == "kdar_ktaggedNpnsrand_dirt"
    //|| ch_name == "kdar_ktaggedns" || ch_name == "kdar_ktaggedns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    || ch_name == "kdar_ktaggedns" || ch_name == "kdar_ktaggedns_ext") && ((merge_time>-4.5 && merge_time<4.5) ) ){flag_pass=false;}
    //|| ch_name == "kdar_ktaggedns" || ch_name == "kdar_ktaggedns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}

    if(ch_name == "kdar_ktagged_sig" || ch_name == "kdar_ktagged_nuwro_train" || ch_name == "kdar_ktagged_gibuu_train" || ch_name == "kdar_ktaggedNp_sig" || ch_name == "kdar_ktaggedNpns_sig")  flag_pass = flag_pass && map_cuts_flag["kdar"];
    if(ch_name == "kdar_ktagged_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];
    if(ch_name == "kdar_ktagged_sigNp" || ch_name == "kdar_ktagged_nuwro_trainNp" || ch_name == "kdar_ktagged_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];
    if(ch_name == "kdar_ktagged_bck" || ch_name == "kdar_ktaggedNp_bck" || ch_name == "kdar_ktaggedNpns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }
    if(ch_name == "kdar_ktagged_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }
    if(ch_name == "kdar_ktagged_nuwro_train" || ch_name == "kdar_ktagged_gibuu_train"|| ch_name == "kdar_ktagged_nuwro_trainNp" || ch_name == "kdar_ktagged_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }
    return flag_pass;



  }else if(ch_name == "kdar_antibdtsel_bck"  || ch_name == "kdar_antibdtsel_dirt" || ch_name == "kdar_antibdtsel_bckNp"
  || ch_name == "kdar_antibdtsel_sig" || ch_name == "kdar_antibdtsel_outFV_sig" || ch_name == "kdar_antibdtsel_sigNp"
  || ch_name == "kdar_antibdtsel_nuwro_train" || ch_name == "kdar_antibdtsel_gibuu_train" || ch_name == "kdar_antibdtsel_nuwro_trainNp" || ch_name == "kdar_antibdtsel_gibuu_trainNp"
  || ch_name == "kdar_antibdtsel" || ch_name == "kdar_antibdtsel_ext"

  || ch_name == "kdar_antibdtselns_bck"  || ch_name == "kdar_antibdtselns_dirt" || ch_name == "kdar_antibdtselnsrand_dirt" || ch_name == "kdar_antibdtselns_bckNp" || ch_name == "kdar_antibdtselns_nons_bck" || ch_name == "kdar_antibdtselns_fake_bck"
  || ch_name == "kdar_antibdtselns_sig" || ch_name == "kdar_antibdtselns_outFV_sig" || ch_name == "kdar_antibdtselns_sigNp"
  || ch_name == "kdar_antibdtselns_nuwro_train" || ch_name == "kdar_antibdtselns_gibuu_train" || ch_name == "kdar_antibdtselns_nuwro_trainNp" || ch_name == "kdar_antibdtselns_gibuu_trainNp"
  || ch_name == "kdar_antibdtselns" || ch_name == "kdar_antibdtselns_ext"

  || ch_name == "kdar_antibdtselNp_sig" || ch_name == "kdar_antibdtselNp_bck" || ch_name == "kdar_antibdtselNp_dirt" || ch_name == "kdar_antibdtselNp_ext" || ch_name == "kdar_antibdtselNp"
  || ch_name == "kdar_antibdtselNpns_sig" || ch_name == "kdar_antibdtselNpns_bck" || ch_name == "kdar_antibdtselNpns_dirt" || ch_name == "kdar_antibdtselNpns_ext" || ch_name == "kdar_antibdtselNpns" || ch_name == "kdar_antibdtselNpnsrand_dirt"){
    bool flag_pass = !(flag_kdar_bdtsel) && flag_kdar_presel;

    if((ch_name == "kdar_antibdtselns_bck"  || ch_name == "kdar_antibdtselns_dirt" || ch_name == "kdar_antibdtselnsrand_dirt" || ch_name == "kdar_antibdtselns_bckNp" || ch_name == "kdar_antibdtselns_nons_bck"
    || ch_name == "kdar_antibdtselns_sig" || ch_name == "kdar_antibdtselns_outFV_sig" || ch_name == "kdar_antibdtselns_sigNp"
    || ch_name == "kdar_antibdtselns_nuwro_train" || ch_name == "kdar_antibdtselns_gibuu_train" || ch_name == "kdar_antibdtselns_nuwro_trainNp" || ch_name == "kdar_antibdtselns_gibuu_trainNp"
    || ch_name == "kdar_antibdtselNpns_dirt" || ch_name == "kdar_antibdtselNpnsrand_dirt"
    //|| ch_name == "kdar_antibdtselns" || ch_name == "kdar_antibdtselns_ext") && merge_time>-3.14 && merge_time<3.14 ){flag_pass=false;}
    //|| ch_name == "kdar_antibdtselns" || ch_name == "kdar_antibdtselns_ext") && ((merge_time>-4.5 && merge_time<4.5) || merge_time<-10) ){flag_pass=false;}
    || ch_name == "kdar_antibdtselns" || ch_name == "kdar_antibdtselns_ext" || ch_name == "kdar_antibdtselns_fake_bck") && ((merge_time>-4.5 && merge_time<4.5) )){flag_pass=false;}

    if(ch_name == "kdar_antibdtsel_sig" || ch_name == "kdar_antibdtsel_nuwro_train" || ch_name == "kdar_antibdtsel_gibuu_train" || ch_name == "kdar_antibdtselns_sig" || ch_name == "kdar_antibdtselns_nuwro_train" || ch_name == "kdar_antibdtselns_gibuu_train" || ch_name == "kdar_antibdtselNp_sig" || ch_name == "kdar_antibdtselNpns_sig")  flag_pass = flag_pass && map_cuts_flag["kdar"];

    if(ch_name == "kdar_antibdtsel_outFV_sig" || ch_name == "kdar_antibdtselns_outFV_sig")  flag_pass = flag_pass && map_cuts_flag["kdarOutFV"];

    if(ch_name == "kdar_antibdtsel_sigNp" || ch_name == "kdar_antibdtsel_nuwro_trainNp" || ch_name == "kdar_antibdtsel_gibuu_trainNp" || ch_name == "kdar_antibdtselns_sigNp" || ch_name == "kdar_antibdtselns_nuwro_trainNp" || ch_name == "kdar_antibdtselns_gibuu_trainNp")  flag_pass = flag_pass && map_cuts_flag["kdarNp"];

    if(ch_name == "kdar_antibdtsel_bck" || ch_name == "kdar_antibdtselns_bck" || ch_name == "kdar_antibdtselns_nons_bck" || ch_name == "kdar_antibdtselNp_bck" || ch_name == "kdar_antibdtselNpns_bck")  { if(flag_pass && map_cuts_flag["kdar"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdar"]; }

    if(ch_name == "kdar_antibdtsel_bckNp" || ch_name == "kdar_antibdtselns_bckNp")  { if(flag_pass && map_cuts_flag["kdarNp"]){std::cout<<"cutting KDAR "<<eval.run<<" "<<eval.subrun<<" "<<eval.event<<std::endl;} flag_pass = flag_pass && !map_cuts_flag["kdarNp"]; }

    if(ch_name == "kdar_antibdtsel_nuwro_train" || ch_name == "kdar_antibdtsel_gibuu_train"|| ch_name == "kdar_antibdtsel_nuwro_trainNp" || ch_name == "kdar_antibdtsel_gibuu_trainNp" || ch_name == "kdar_antibdtselns_nuwro_train" || ch_name == "kdar_antibdtselns_gibuu_train"|| ch_name == "kdar_antibdtselns_nuwro_trainNp" || ch_name == "kdar_antibdtselns_gibuu_trainNp")  {if(eval.event%10<5){flag_pass=false;} }

    if(ch_name == "kdar_antibdtselNp_sig" || ch_name == "kdar_antibdtselNp_bck" || ch_name == "kdar_antibdtselNp_dirt" || ch_name == "kdar_antibdtselNp_ext" || ch_name == "kdar_antibdtselNp"
  || ch_name == "kdar_antibdtselNpns_sig" || ch_name == "kdar_antibdtselNpns_bck" || ch_name == "kdar_antibdtselNpns_dirt" || ch_name == "kdar_antibdtselNpns_ext" || ch_name == "kdar_antibdtselNpns" || ch_name == "kdar_antibdtselNpnsrand_dirt"){if(tagger.ssm_prim_track1_kine_energy_range<0) flag_pass=false;}

    return flag_pass;



  }else if (ch_name == "LEE_FC_nueoverlay"  || ch_name == "nueCC_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "nueCC_FC_nueoverlay_numi"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if ( ch_name == "nueCC_FC_numu2nueoverlay" ){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if ( ch_name == "nueCC_FC_numu2nueoverlay_numi" ){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext" || ch_name == "BG_nueCC_FC_dirt" || ch_name =="nueCC_FC_bnb"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext_numi" || ch_name == "BG_nueCC_FC_dirt_numi" || ch_name =="nueCC_FC_numi"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay_numi"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay" ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "nueCC_PC_nueoverlay_numi" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if ( ch_name == "nueCC_PC_numu2nueoverlay" ){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if ( ch_name == "nueCC_PC_numu2nueoverlay_numi" ){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext" || ch_name == "BG_nueCC_PC_dirt" || ch_name == "nueCC_PC_bnb"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext_numi" || ch_name == "BG_nueCC_PC_dirt_numi" || ch_name == "nueCC_PC_numi"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay_numi"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta0_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta0_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta0_FC_overlay"     || ch_name == "numuCC_signal_nu_theta0_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta0_FC_overlay" || ch_name == "numuCC_background_Emu_theta0_FC_overlay" || ch_name == "numuCC_background_Pmu_theta0_FC_overlay" || ch_name == "numuCC_background_nu_theta0_FC_overlay"
         || ch_name == "BG_numuCC_theta0_FC_ext"                 || ch_name =="BG_numuCC_theta0_FC_dirt"                 || ch_name == "numuCC_theta0_FC_bnb"
         || ch_name == "BG_numuCC_theta0_FC_ext_v2"              || ch_name =="BG_numuCC_theta0_FC_dirt_v2"              || ch_name == "numuCC_theta0_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) {
      if      (ch_name == "numuCC_signal_Enu_theta0_FC_overlay" || ch_name == "numuCC_background_Enu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta0_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta0_FC_overlay" || ch_name == "numuCC_background_Emu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta0_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta0_FC_overlay" || ch_name == "numuCC_background_Pmu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta0_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta0_FC_overlay" || ch_name ==  "numuCC_background_nu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta0_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta1_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta1_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta1_FC_overlay"     || ch_name == "numuCC_signal_nu_theta1_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta1_FC_overlay" || ch_name == "numuCC_background_Emu_theta1_FC_overlay" || ch_name == "numuCC_background_Pmu_theta1_FC_overlay" || ch_name == "numuCC_background_nu_theta1_FC_overlay"
         || ch_name == "BG_numuCC_theta1_FC_ext"                 || ch_name =="BG_numuCC_theta1_FC_dirt"                 || ch_name == "numuCC_theta1_FC_bnb"
         || ch_name == "BG_numuCC_theta1_FC_ext_v2"              || ch_name =="BG_numuCC_theta1_FC_dirt_v2"              || ch_name == "numuCC_theta1_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) {
      if      (ch_name == "numuCC_signal_Enu_theta1_FC_overlay" || ch_name == "numuCC_background_Enu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta1_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta1_FC_overlay" || ch_name == "numuCC_background_Emu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta1_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta1_FC_overlay" || ch_name == "numuCC_background_Pmu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta1_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta1_FC_overlay" || ch_name ==  "numuCC_background_nu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta1_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta2_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta2_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta2_FC_overlay"     || ch_name == "numuCC_signal_nu_theta2_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta2_FC_overlay" || ch_name == "numuCC_background_Emu_theta2_FC_overlay" || ch_name == "numuCC_background_Pmu_theta2_FC_overlay" || ch_name == "numuCC_background_nu_theta2_FC_overlay"
         || ch_name == "BG_numuCC_theta2_FC_ext"                 || ch_name =="BG_numuCC_theta2_FC_dirt"                 || ch_name == "numuCC_theta2_FC_bnb"
         || ch_name == "BG_numuCC_theta2_FC_ext_v2"              || ch_name =="BG_numuCC_theta2_FC_dirt_v2"              || ch_name == "numuCC_theta2_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) {
      if      (ch_name == "numuCC_signal_Enu_theta2_FC_overlay" || ch_name == "numuCC_background_Enu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta2_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta2_FC_overlay" || ch_name == "numuCC_background_Emu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta2_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta2_FC_overlay" || ch_name == "numuCC_background_Pmu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta2_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta2_FC_overlay" || ch_name ==  "numuCC_background_nu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta2_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta3_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta3_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta3_FC_overlay"     || ch_name == "numuCC_signal_nu_theta3_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta3_FC_overlay" || ch_name == "numuCC_background_Emu_theta3_FC_overlay" || ch_name == "numuCC_background_Pmu_theta3_FC_overlay" || ch_name == "numuCC_background_nu_theta3_FC_overlay"
         || ch_name == "BG_numuCC_theta3_FC_ext"                 || ch_name =="BG_numuCC_theta3_FC_dirt"                 || ch_name == "numuCC_theta3_FC_bnb"
         || ch_name == "BG_numuCC_theta3_FC_ext_v2"              || ch_name =="BG_numuCC_theta3_FC_dirt_v2"              || ch_name == "numuCC_theta3_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) {
      if      (ch_name == "numuCC_signal_Enu_theta3_FC_overlay" || ch_name == "numuCC_background_Enu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta3_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta3_FC_overlay" || ch_name == "numuCC_background_Emu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta3_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta3_FC_overlay" || ch_name == "numuCC_background_Pmu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta3_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta3_FC_overlay" || ch_name ==  "numuCC_background_nu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta3_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta4_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta4_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta4_FC_overlay"     || ch_name == "numuCC_signal_nu_theta4_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta4_FC_overlay" || ch_name == "numuCC_background_Emu_theta4_FC_overlay" || ch_name == "numuCC_background_Pmu_theta4_FC_overlay" || ch_name == "numuCC_background_nu_theta4_FC_overlay"
         || ch_name == "BG_numuCC_theta4_FC_ext"                 || ch_name =="BG_numuCC_theta4_FC_dirt"                 || ch_name == "numuCC_theta4_FC_bnb"
         || ch_name == "BG_numuCC_theta4_FC_ext_v2"              || ch_name =="BG_numuCC_theta4_FC_dirt_v2"              || ch_name == "numuCC_theta4_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) {
      if      (ch_name == "numuCC_signal_Enu_theta4_FC_overlay" || ch_name == "numuCC_background_Enu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta4_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta4_FC_overlay" || ch_name == "numuCC_background_Emu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta4_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta4_FC_overlay" || ch_name == "numuCC_background_Pmu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta4_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta4_FC_overlay" || ch_name ==  "numuCC_background_nu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta4_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta5_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta5_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta5_FC_overlay"     || ch_name == "numuCC_signal_nu_theta5_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta5_FC_overlay" || ch_name == "numuCC_background_Emu_theta5_FC_overlay" || ch_name == "numuCC_background_Pmu_theta5_FC_overlay" || ch_name == "numuCC_background_nu_theta5_FC_overlay"
         || ch_name == "BG_numuCC_theta5_FC_ext"                 || ch_name =="BG_numuCC_theta5_FC_dirt"                 || ch_name == "numuCC_theta5_FC_bnb"
         || ch_name == "BG_numuCC_theta5_FC_ext_v2"              || ch_name =="BG_numuCC_theta5_FC_dirt_v2"              || ch_name == "numuCC_theta5_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) {
      if      (ch_name == "numuCC_signal_Enu_theta5_FC_overlay" || ch_name == "numuCC_background_Enu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta5_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta5_FC_overlay" || ch_name == "numuCC_background_Emu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta5_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta5_FC_overlay" || ch_name == "numuCC_background_Pmu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta5_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta5_FC_overlay" || ch_name ==  "numuCC_background_nu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta5_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta6_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta6_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta6_FC_overlay"     || ch_name == "numuCC_signal_nu_theta6_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta6_FC_overlay" || ch_name == "numuCC_background_Emu_theta6_FC_overlay" || ch_name == "numuCC_background_Pmu_theta6_FC_overlay" || ch_name == "numuCC_background_nu_theta6_FC_overlay"
         || ch_name == "BG_numuCC_theta6_FC_ext"                 || ch_name =="BG_numuCC_theta6_FC_dirt"                 || ch_name == "numuCC_theta6_FC_bnb"
         || ch_name == "BG_numuCC_theta6_FC_ext_v2"              || ch_name =="BG_numuCC_theta6_FC_dirt_v2"              || ch_name == "numuCC_theta6_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) {
      if      (ch_name == "numuCC_signal_Enu_theta6_FC_overlay" || ch_name == "numuCC_background_Enu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta6_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta6_FC_overlay" || ch_name == "numuCC_background_Emu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta6_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta6_FC_overlay" || ch_name == "numuCC_background_Pmu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta6_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta6_FC_overlay" || ch_name ==  "numuCC_background_nu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta6_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta7_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta7_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta7_FC_overlay"     || ch_name == "numuCC_signal_nu_theta7_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta7_FC_overlay" || ch_name == "numuCC_background_Emu_theta7_FC_overlay" || ch_name == "numuCC_background_Pmu_theta7_FC_overlay" || ch_name == "numuCC_background_nu_theta7_FC_overlay"
         || ch_name == "BG_numuCC_theta7_FC_ext"                 || ch_name =="BG_numuCC_theta7_FC_dirt"                 || ch_name == "numuCC_theta7_FC_bnb"
         || ch_name == "BG_numuCC_theta7_FC_ext_v2"              || ch_name =="BG_numuCC_theta7_FC_dirt_v2"              || ch_name == "numuCC_theta7_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) {
      if      (ch_name == "numuCC_signal_Enu_theta7_FC_overlay" || ch_name == "numuCC_background_Enu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta7_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta7_FC_overlay" || ch_name == "numuCC_background_Emu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta7_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta7_FC_overlay" || ch_name == "numuCC_background_Pmu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta7_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta7_FC_overlay" || ch_name ==  "numuCC_background_nu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta7_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta8_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta8_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta8_FC_overlay"     || ch_name == "numuCC_signal_nu_theta8_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta8_FC_overlay" || ch_name == "numuCC_background_Emu_theta8_FC_overlay" || ch_name == "numuCC_background_Pmu_theta8_FC_overlay" || ch_name == "numuCC_background_nu_theta8_FC_overlay"
         || ch_name == "BG_numuCC_theta8_FC_ext"                 || ch_name =="BG_numuCC_theta8_FC_dirt"                 || ch_name == "numuCC_theta8_FC_bnb"
         || ch_name == "BG_numuCC_theta8_FC_ext_v2"              || ch_name =="BG_numuCC_theta8_FC_dirt_v2"              || ch_name == "numuCC_theta8_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) {
      if      (ch_name == "numuCC_signal_Enu_theta8_FC_overlay" || ch_name == "numuCC_background_Enu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta8_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta8_FC_overlay" || ch_name == "numuCC_background_Emu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta8_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta8_FC_overlay" || ch_name == "numuCC_background_Pmu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta8_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta8_FC_overlay" || ch_name ==  "numuCC_background_nu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta8_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu0_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta0_Pmu_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu0_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu0_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu0_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu0_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu0_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu0_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu0_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu0_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu1_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu1_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu1_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu1_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu1_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu1_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu1_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu1_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu1_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu2_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu2_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu2_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu2_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu2_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu2_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu2_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu2_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu2_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu3_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu3_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu3_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu3_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu3_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu3_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;

  }else if (ch_name == "numuCC_signal_Enu3_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu3_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu3_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu3_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Pmu0_FC_overlay" || ch_name == "numuCC_background_Pmu0_FC_overlay"
         || ch_name == "BG_numuCC_Pmu0_FC_ext"         || ch_name =="BG_numuCC_Pmu0_FC_dirt" || ch_name == "numuCC_Pmu0_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==0) {
      if      (ch_name == "numuCC_signal_Pmu0_FC_overlay" || ch_name == "numuCC_background_Pmu0_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu0_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu1_FC_overlay" || ch_name == "numuCC_background_Pmu1_FC_overlay"
         || ch_name == "BG_numuCC_Pmu1_FC_ext"         || ch_name =="BG_numuCC_Pmu1_FC_dirt" || ch_name == "numuCC_Pmu1_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==1) {
      if      (ch_name == "numuCC_signal_Pmu1_FC_overlay" || ch_name == "numuCC_background_Pmu1_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu1_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu2_FC_overlay" || ch_name == "numuCC_background_Pmu2_FC_overlay"
         || ch_name == "BG_numuCC_Pmu2_FC_ext"         || ch_name =="BG_numuCC_Pmu2_FC_dirt" || ch_name == "numuCC_Pmu2_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==2) {
      if      (ch_name == "numuCC_signal_Pmu2_FC_overlay" || ch_name == "numuCC_background_Pmu2_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu2_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu3_FC_overlay" || ch_name == "numuCC_background_Pmu3_FC_overlay"
         || ch_name == "BG_numuCC_Pmu3_FC_ext"         || ch_name =="BG_numuCC_Pmu3_FC_dirt" || ch_name == "numuCC_Pmu3_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==3) {
      if      (ch_name == "numuCC_signal_Pmu3_FC_overlay" || ch_name == "numuCC_background_Pmu3_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu3_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu4_FC_overlay" || ch_name == "numuCC_background_Pmu4_FC_overlay"
         || ch_name == "BG_numuCC_Pmu4_FC_ext"         || ch_name =="BG_numuCC_Pmu4_FC_dirt" || ch_name == "numuCC_Pmu4_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==4) {
      if      (ch_name == "numuCC_signal_Pmu4_FC_overlay" || ch_name == "numuCC_background_Pmu4_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu4_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu5_FC_overlay" || ch_name == "numuCC_background_Pmu5_FC_overlay"
         || ch_name == "BG_numuCC_Pmu5_FC_ext"         || ch_name =="BG_numuCC_Pmu5_FC_dirt" || ch_name == "numuCC_Pmu5_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==5) {
      if      (ch_name == "numuCC_signal_Pmu5_FC_overlay" || ch_name == "numuCC_background_Pmu5_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu5_FC_overlay")); }
      else return true;
    }  else return false;
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
  }else if (ch_name == "numuCC_signal_Enu_theta0_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta0_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta0_PC_overlay"     || ch_name == "numuCC_signal_nu_theta0_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta0_PC_overlay" || ch_name == "numuCC_background_Emu_theta0_PC_overlay" || ch_name == "numuCC_background_Pmu_theta0_PC_overlay" || ch_name == "numuCC_background_nu_theta0_PC_overlay"
         || ch_name == "BG_numuCC_theta0_PC_ext"                 || ch_name =="BG_numuCC_theta0_PC_dirt"                 || ch_name == "numuCC_theta0_PC_bnb"
         || ch_name == "BG_numuCC_theta0_PC_ext_v2"              || ch_name =="BG_numuCC_theta0_PC_dirt_v2"              || ch_name == "numuCC_theta0_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) {
      if      (ch_name == "numuCC_signal_Enu_theta0_PC_overlay" || ch_name == "numuCC_background_Enu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta0_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta0_PC_overlay" || ch_name == "numuCC_background_Emu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta0_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta0_PC_overlay" || ch_name == "numuCC_background_Pmu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta0_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta0_PC_overlay" || ch_name ==  "numuCC_background_nu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta0_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta1_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta1_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta1_PC_overlay"     || ch_name == "numuCC_signal_nu_theta1_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta1_PC_overlay" || ch_name == "numuCC_background_Emu_theta1_PC_overlay" || ch_name == "numuCC_background_Pmu_theta1_PC_overlay" || ch_name == "numuCC_background_nu_theta1_PC_overlay"
         || ch_name == "BG_numuCC_theta1_PC_ext"                 || ch_name =="BG_numuCC_theta1_PC_dirt"                 || ch_name == "numuCC_theta1_PC_bnb"
         || ch_name == "BG_numuCC_theta1_PC_ext_v2"              || ch_name =="BG_numuCC_theta1_PC_dirt_v2"              || ch_name == "numuCC_theta1_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) {
      if      (ch_name == "numuCC_signal_Enu_theta1_PC_overlay" || ch_name == "numuCC_background_Enu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta1_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta1_PC_overlay" || ch_name == "numuCC_background_Emu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta1_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta1_PC_overlay" || ch_name == "numuCC_background_Pmu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta1_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta1_PC_overlay" || ch_name ==  "numuCC_background_nu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta1_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta2_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta2_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta2_PC_overlay"     || ch_name == "numuCC_signal_nu_theta2_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta2_PC_overlay" || ch_name == "numuCC_background_Emu_theta2_PC_overlay" || ch_name == "numuCC_background_Pmu_theta2_PC_overlay" || ch_name == "numuCC_background_nu_theta2_PC_overlay"
         || ch_name == "BG_numuCC_theta2_PC_ext"                 || ch_name =="BG_numuCC_theta2_PC_dirt"                 || ch_name == "numuCC_theta2_PC_bnb"
         || ch_name == "BG_numuCC_theta2_PC_ext_v2"              || ch_name =="BG_numuCC_theta2_PC_dirt_v2"              || ch_name == "numuCC_theta2_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) {
      if      (ch_name == "numuCC_signal_Enu_theta2_PC_overlay" || ch_name == "numuCC_background_Enu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta2_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta2_PC_overlay" || ch_name == "numuCC_background_Emu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta2_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta2_PC_overlay" || ch_name == "numuCC_background_Pmu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta2_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta2_PC_overlay" || ch_name ==  "numuCC_background_nu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta2_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta3_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta3_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta3_PC_overlay"     || ch_name == "numuCC_signal_nu_theta3_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta3_PC_overlay" || ch_name == "numuCC_background_Emu_theta3_PC_overlay" || ch_name == "numuCC_background_Pmu_theta3_PC_overlay" || ch_name == "numuCC_background_nu_theta3_PC_overlay"
         || ch_name == "BG_numuCC_theta3_PC_ext"                 || ch_name =="BG_numuCC_theta3_PC_dirt"                 || ch_name == "numuCC_theta3_PC_bnb"
         || ch_name == "BG_numuCC_theta3_PC_ext_v2"              || ch_name =="BG_numuCC_theta3_PC_dirt_v2"              || ch_name == "numuCC_theta3_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) {
      if      (ch_name == "numuCC_signal_Enu_theta3_PC_overlay" || ch_name == "numuCC_background_Enu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta3_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta3_PC_overlay" || ch_name == "numuCC_background_Emu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta3_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta3_PC_overlay" || ch_name == "numuCC_background_Pmu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta3_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta3_PC_overlay" || ch_name ==  "numuCC_background_nu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta3_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta4_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta4_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta4_PC_overlay"     || ch_name == "numuCC_signal_nu_theta4_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta4_PC_overlay" || ch_name == "numuCC_background_Emu_theta4_PC_overlay" || ch_name == "numuCC_background_Pmu_theta4_PC_overlay" || ch_name == "numuCC_background_nu_theta4_PC_overlay"
         || ch_name == "BG_numuCC_theta4_PC_ext"                 || ch_name =="BG_numuCC_theta4_PC_dirt"                 || ch_name == "numuCC_theta4_PC_bnb"
         || ch_name == "BG_numuCC_theta4_PC_ext_v2"              || ch_name =="BG_numuCC_theta4_PC_dirt_v2"              || ch_name == "numuCC_theta4_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) {
      if      (ch_name == "numuCC_signal_Enu_theta4_PC_overlay" || ch_name == "numuCC_background_Enu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta4_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta4_PC_overlay" || ch_name == "numuCC_background_Emu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta4_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta4_PC_overlay" || ch_name == "numuCC_background_Pmu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta4_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta4_PC_overlay" || ch_name ==  "numuCC_background_nu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta4_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta5_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta5_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta5_PC_overlay"     || ch_name == "numuCC_signal_nu_theta5_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta5_PC_overlay" || ch_name == "numuCC_background_Emu_theta5_PC_overlay" || ch_name == "numuCC_background_Pmu_theta5_PC_overlay" || ch_name == "numuCC_background_nu_theta5_PC_overlay"
         || ch_name == "BG_numuCC_theta5_PC_ext"                 || ch_name =="BG_numuCC_theta5_PC_dirt"                 || ch_name == "numuCC_theta5_PC_bnb"
         || ch_name == "BG_numuCC_theta5_PC_ext_v2"              || ch_name =="BG_numuCC_theta5_PC_dirt_v2"              || ch_name == "numuCC_theta5_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) {
      if      (ch_name == "numuCC_signal_Enu_theta5_PC_overlay" || ch_name == "numuCC_background_Enu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta5_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta5_PC_overlay" || ch_name == "numuCC_background_Emu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta5_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta5_PC_overlay" || ch_name == "numuCC_background_Pmu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta5_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta5_PC_overlay" || ch_name ==  "numuCC_background_nu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta5_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta6_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta6_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta6_PC_overlay"     || ch_name == "numuCC_signal_nu_theta6_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta6_PC_overlay" || ch_name == "numuCC_background_Emu_theta6_PC_overlay" || ch_name == "numuCC_background_Pmu_theta6_PC_overlay" || ch_name == "numuCC_background_nu_theta6_PC_overlay"
         || ch_name == "BG_numuCC_theta6_PC_ext"                 || ch_name =="BG_numuCC_theta6_PC_dirt"                 || ch_name == "numuCC_theta6_PC_bnb"
         || ch_name == "BG_numuCC_theta6_PC_ext_v2"              || ch_name =="BG_numuCC_theta6_PC_dirt_v2"              || ch_name == "numuCC_theta6_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) {
      if      (ch_name == "numuCC_signal_Enu_theta6_PC_overlay" || ch_name == "numuCC_background_Enu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta6_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta6_PC_overlay" || ch_name == "numuCC_background_Emu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta6_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta6_PC_overlay" || ch_name == "numuCC_background_Pmu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta6_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta6_PC_overlay" || ch_name ==  "numuCC_background_nu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta6_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta7_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta7_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta7_PC_overlay"     || ch_name == "numuCC_signal_nu_theta7_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta7_PC_overlay" || ch_name == "numuCC_background_Emu_theta7_PC_overlay" || ch_name == "numuCC_background_Pmu_theta7_PC_overlay" || ch_name == "numuCC_background_nu_theta7_PC_overlay"
         || ch_name == "BG_numuCC_theta7_PC_ext"                 || ch_name =="BG_numuCC_theta7_PC_dirt"                 || ch_name == "numuCC_theta7_PC_bnb"
         || ch_name == "BG_numuCC_theta7_PC_ext_v2"              || ch_name =="BG_numuCC_theta7_PC_dirt_v2"              || ch_name == "numuCC_theta7_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) {
      if      (ch_name == "numuCC_signal_Enu_theta7_PC_overlay" || ch_name == "numuCC_background_Enu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta7_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta7_PC_overlay" || ch_name == "numuCC_background_Emu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta7_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta7_PC_overlay" || ch_name == "numuCC_background_Pmu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta7_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta7_PC_overlay" || ch_name ==  "numuCC_background_nu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta7_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta8_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta8_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta8_PC_overlay"     || ch_name == "numuCC_signal_nu_theta8_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta8_PC_overlay" || ch_name == "numuCC_background_Emu_theta8_PC_overlay" || ch_name == "numuCC_background_Pmu_theta8_PC_overlay" || ch_name == "numuCC_background_nu_theta8_PC_overlay"
         || ch_name == "BG_numuCC_theta8_PC_ext"                 || ch_name =="BG_numuCC_theta8_PC_dirt"                 || ch_name == "numuCC_theta8_PC_bnb"
         || ch_name == "BG_numuCC_theta8_PC_ext_v2"              || ch_name =="BG_numuCC_theta8_PC_dirt_v2"              || ch_name == "numuCC_theta8_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) {
      if      (ch_name == "numuCC_signal_Enu_theta8_PC_overlay" || ch_name == "numuCC_background_Enu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta8_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta8_PC_overlay" || ch_name == "numuCC_background_Emu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta8_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta8_PC_overlay" || ch_name == "numuCC_background_Pmu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta8_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta8_PC_overlay" || ch_name ==  "numuCC_background_nu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta8_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu0_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu0_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu0_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu0_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu0_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu0_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu0_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu0_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu0_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu1_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu1_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu1_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu1_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu1_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu1_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu1_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu1_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu1_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu2_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu2_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu2_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu2_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu2_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu2_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu2_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu2_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu2_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu3_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu3_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu3_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu3_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu3_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu3_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu3_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu3_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu3_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Pmu0_PC_overlay" || ch_name == "numuCC_background_Pmu0_PC_overlay"
         || ch_name == "BG_numuCC_Pmu0_PC_ext"         || ch_name =="BG_numuCC_Pmu0_PC_dirt" || ch_name == "numuCC_Pmu0_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==0) {
      if      (ch_name == "numuCC_signal_Pmu0_PC_overlay" || ch_name == "numuCC_background_Pmu0_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu0_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu1_PC_overlay" || ch_name == "numuCC_background_Pmu1_PC_overlay"
         || ch_name == "BG_numuCC_Pmu1_PC_ext"         || ch_name =="BG_numuCC_Pmu1_PC_dirt" || ch_name == "numuCC_Pmu1_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==1) {
      if      (ch_name == "numuCC_signal_Pmu1_PC_overlay" || ch_name == "numuCC_background_Pmu1_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu1_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu2_PC_overlay" || ch_name == "numuCC_background_Pmu2_PC_overlay"
         || ch_name == "BG_numuCC_Pmu2_PC_ext"         || ch_name =="BG_numuCC_Pmu2_PC_dirt" || ch_name == "numuCC_Pmu2_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==2) {
      if      (ch_name == "numuCC_signal_Pmu2_PC_overlay" || ch_name == "numuCC_background_Pmu2_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu2_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu3_PC_overlay" || ch_name == "numuCC_background_Pmu3_PC_overlay"
         || ch_name == "BG_numuCC_Pmu3_PC_ext"         || ch_name =="BG_numuCC_Pmu3_PC_dirt" || ch_name == "numuCC_Pmu3_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==3) {
      if      (ch_name == "numuCC_signal_Pmu3_PC_overlay" || ch_name == "numuCC_background_Pmu3_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu3_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu4_PC_overlay" || ch_name == "numuCC_background_Pmu4_PC_overlay"
         || ch_name == "BG_numuCC_Pmu4_PC_ext"         || ch_name =="BG_numuCC_Pmu4_PC_dirt" || ch_name == "numuCC_Pmu4_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==4) {
      if      (ch_name == "numuCC_signal_Pmu4_PC_overlay" || ch_name == "numuCC_background_Pmu4_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu4_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu5_PC_overlay" || ch_name == "numuCC_background_Pmu5_PC_overlay"
         || ch_name == "BG_numuCC_Pmu5_PC_ext"         || ch_name =="BG_numuCC_Pmu5_PC_dirt" || ch_name == "numuCC_Pmu5_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==5) {
      if      (ch_name == "numuCC_signal_Pmu5_PC_overlay" || ch_name == "numuCC_background_Pmu5_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu5_PC_overlay")); }
      else return true;
    }  else return false;
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
  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay_numi" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext_numi" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt_numi" || ch_name == "numuCC_nopi0_nonueCC_FC_numi" || ch_name == "numuCC_nopi0_nonueCC_FC_numu2nueoverlay_numi"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb" || ch_name == "numuCC_nopi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay_numi" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext_numi" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt_numi" || ch_name == "numuCC_nopi0_nonueCC_PC_numi" || ch_name == "numuCC_nopi0_nonueCC_PC_numu2nueoverlay_numi"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay" || ch_name =="BG_CCpi0_nonueCC_FC_ext" || ch_name == "BG_CCpi0_nonueCC_FC_dirt" || ch_name == "CCpi0_nonueCC_FC_bnb" || ch_name == "CCpi0_nonueCC_FC_numu2nueoverlay"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay_numi" || ch_name =="BG_CCpi0_nonueCC_FC_ext_numi" || ch_name == "BG_CCpi0_nonueCC_FC_dirt_numi" || ch_name == "CCpi0_nonueCC_FC_numi" || ch_name == "CCpi0_nonueCC_FC_numu2nueoverlay_numi"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay" || ch_name == "BG_CCpi0_nonueCC_PC_ext" || ch_name == "BG_CCpi0_nonueCC_PC_dirt" || ch_name == "CCpi0_nonueCC_PC_bnb" || ch_name == "CCpi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay_numi" || ch_name == "BG_CCpi0_nonueCC_PC_ext_numi" || ch_name == "BG_CCpi0_nonueCC_PC_dirt_numi" || ch_name == "CCpi0_nonueCC_PC_numi" || ch_name == "CCpi0_nonueCC_PC_numu2nueoverlay_numi"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay" || ch_name == "BG_NCpi0_nonueCC_ext" || ch_name == "BG_NCpi0_nonueCC_dirt" || ch_name == "NCpi0_nonueCC_bnb" || ch_name == "NCpi0_nonueCC_numu2nueoverlay"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    // if (flag_NC && flag_pi0 && (!flag_nueCC) && flag_FC && (!flag_0p) ) return true; // a test ...
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay_numi" || ch_name == "BG_NCpi0_nonueCC_ext_numi" || ch_name == "BG_NCpi0_nonueCC_dirt_numi" || ch_name == "NCpi0_nonueCC_numi" || ch_name == "NCpi0_nonueCC_numu2nueoverlay_numi"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "nueCC_bnb" || ch_name == "nueCC_nueoverlay"){   // side band ...
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay" || ch_name == "nueCC_bnb") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb" || ch_name == "all_but_nueCC_overlay" || ch_name == "all_but_nueCC_ext" || ch_name == "all_but_nueCC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay" || ch_name != "all_but_nueCC_overlay") return true;
    else return false;
  }else if (ch_name == "nueCC_bnb1" || ch_name == "nueCC_nueoverlay1"){
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay1" || ch_name == "nueCC_bnb1") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb1" || ch_name == "all_but_nueCC_overlay1" || ch_name == "all_but_nueCC_ext1" || ch_name == "all_but_nueCC_dirt1"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay1" || ch_name != "all_but_nueCC_overlay1") return true;
    else return false;
  }else if (ch_name == "testA_bnb" || ch_name == "testA_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testA_nueoverlay" || ch_name == "testA_bnb") return true;
    else return false;
  }else if (ch_name == "testA_overlay" || ch_name == "testA_ext" || ch_name == "testA_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testA_overlay" || ch_name != "testA_overlay") return true;
    else return false;
  }else if (ch_name == "testB_bnb" || ch_name == "testB_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testB_nueoverlay" || ch_name == "testB_bnb") return true;
    else return false;
  }else if (ch_name == "testB_overlay" || ch_name == "testB_ext" || ch_name == "testB_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testB_overlay" || ch_name != "testB_overlay") return true;
    else return false;
  }else if (ch_name == "testC_bnb" || ch_name == "testC_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testC_nueoverlay" || ch_name == "testC_bnb") return true;
    else return false;
  }else if (ch_name == "testC_overlay" || ch_name == "testC_ext" || ch_name == "testC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testC_overlay" || ch_name != "testC_overlay") return true;
    else return false;
  }else if (ch_name == "testD_bnb" || ch_name == "testD_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testD_nueoverlay" || ch_name == "testD_bnb") return true;
    else return false;
  }else if (ch_name == "testD_overlay" || ch_name == "testD_ext" || ch_name == "testD_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testD_overlay" || ch_name != "testD_overlay") return true;
    else return false;
 // Janet's requests: <600 MeV numuCC PC, FC for three variables = 6 obs channels
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;

  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;

  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // Mike Shaevitz >800 MeV nueCC PC+FC 1 obs channel
  }else if (ch_name == "nueCC_extra_nueoverlay"){
    if (flag_nueCC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext" || ch_name == "BG_nueCC_extra_dirt" || ch_name =="nueCC_extra_bnb"){
    if (flag_nueCC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;

  }else if (ch_name == "nueCC_extra_nueoverlay_fc"){
    if (flag_nueCC && flag_truth_inside && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc" || ch_name == "BG_nueCC_extra_dirt_fc" || ch_name =="nueCC_extra_bnb_fc"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_FC) return true;
    else return false;

    // FC and Np
  }else if (ch_name == "nueCC_extra_nueoverlay_fc_np"){
    if (flag_nueCC && flag_truth_inside && flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc_np" || ch_name == "BG_nueCC_extra_dirt_fc_np" || ch_name =="nueCC_extra_bnb_fc_np"){
    if (flag_nueCC&& flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc_np"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)&& flag_FC && (!flag_0p)) return true;
    else return false;
    // FC and 0p
  }else if (ch_name == "nueCC_extra_nueoverlay_fc_0p"){
    if (flag_nueCC && flag_truth_inside && flag_FC && (flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc_0p" || ch_name == "BG_nueCC_extra_dirt_fc_0p" || ch_name =="nueCC_extra_bnb_fc_0p"){
    if (flag_nueCC&& flag_FC && (flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc_0p"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)&& flag_FC && (flag_0p)) return true;
    else return false;

 // cut-based numuCC FC/PC 2 obs channels
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_cutbased && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_cutbased && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
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
 // nueCC 3 variables: n_trakcs, n_showers, gap_n_bad, FC/PC x3 = 6 channels; 4 additional channels
  }else if (ch_name == "nueCC2_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_ext" || ch_name == "BG_nueCC2_FC_dirt" || ch_name =="nueCC2_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC2_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_ext" || ch_name == "BG_nueCC2_PC_dirt" || ch_name == "nueCC2_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_ext" || ch_name == "BG_nueCC3_FC_dirt" || ch_name =="nueCC3_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_ext" || ch_name == "BG_nueCC3_PC_dirt" || ch_name == "nueCC3_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
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
  }else if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC_background_PC_overlay"
  	    || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC1_background_PC_overlay"
  	    || ch_name == "numuCC2_signal_FC_overlay" || ch_name == "numuCC2_signal_PC_overlay" || ch_name == "numuCC2_background_FC_overlay" || ch_name == "numuCC2_background_PC_overlay"
  	    ){
    if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC2_signal_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC2_signal_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC2_background_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_PC_overlay" || ch_name == "numuCC1_background_PC_overlay" || ch_name == "numuCC2_background_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }
    return false;
  } else if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC_background_Enu_PC_overlay"
	     || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay"
	     || ch_name == "numuCC2_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"
       || ch_name == "numuCC_signal_Enu_overlay" || ch_name == "numuCC_background_Enu_overlay"
	    ){
    if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_signal_Enu_overlay"){
      if (flag_numuCC && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Enu_overlay"){
      if (flag_numuCC && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }
    return false;
  // 1D Enu channel with the same inclusive signal definition as the 3D selection
  }else if (ch_name == "numuCC_Enu_mu_FC_bnb" || ch_name == "BG_numuCC_Enu_mu_FC_ext" || ch_name == "BG_numuCC_Enu_mu_FC_dirt"){
    if (flag_numuCC &&   flag_FC  && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;
  }else if (ch_name == "numuCC_Enu_mu_PC_bnb" || ch_name == "BG_numuCC_Enu_mu_PC_ext" || ch_name == "BG_numuCC_Enu_mu_PC_dirt"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;
  } else if (ch_name == "numuCC_signal_Enu_mu_FC_overlay"    || ch_name == "numuCC_signal_Enu_mu_PC_overlay"
         || ch_name == "numuCC_background_Enu_mu_FC_overlay" || ch_name == "numuCC_background_Enu_mu_PC_overlay" ){
    bool pre_cut = flag_numuCC && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8;
    if      (ch_name == "numuCC_signal_Enu_mu_FC_overlay"     && pre_cut &&   flag_FC  &&   map_cuts_flag["Xs_Enu_mu_numuCCinFV"])  { return true; }
    else if (ch_name == "numuCC_signal_Enu_mu_PC_overlay"     && pre_cut && (!flag_FC) &&   map_cuts_flag["Xs_Enu_mu_numuCCinFV"])  { return true; }
    else if (ch_name == "numuCC_background_Enu_mu_FC_overlay" && pre_cut &&   flag_FC  && (!map_cuts_flag["Xs_Enu_mu_numuCCinFV"])) { return true; }
    else if (ch_name == "numuCC_background_Enu_mu_PC_overlay" && pre_cut && (!flag_FC) && (!map_cuts_flag["Xs_Enu_mu_numuCCinFV"])) { return true; }
    return false;
  // ------
  }else if (ch_name == "nueCC_signal_Enu_FC_overlay" || ch_name == "nueCC_signal_Enu_PC_overlay" || ch_name == "nueCC_background_Enu_FC_overlay" || ch_name == "nueCC_background_Enu_PC_overlay"
       || ch_name == "nueCC_signal_Enu_overlay" || ch_name == "nueCC_background_Enu_overlay"
      ){
    if (ch_name == "nueCC_signal_Enu_FC_overlay"){
      if (flag_nueCC && flag_FC && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_signal_Enu_PC_overlay"){
      if (flag_nueCC && (!flag_FC) && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_background_Enu_FC_overlay"){
      if (flag_nueCC && flag_FC && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }else if (ch_name == "nueCC_background_Enu_PC_overlay" ){
      if (flag_nueCC && (!flag_FC) && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }else if (ch_name == "nueCC_signal_Enu_overlay"){
      if (flag_nueCC && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_background_Enu_overlay"){
      if (flag_nueCC && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }
    return false;

  } else if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC_background_Emu_PC_overlay"
	     || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay"
	     || ch_name == "numuCC2_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"
	     ){
    if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }
    return false;

  } else if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC_background_Ehad_PC_overlay"
	     || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay"
	     || ch_name == "numuCC2_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"
	    ){
    if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }
    return false;

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
  }else if (ch_name == "nc_pio_energy_FC" || ch_name == "nc_pio_score_FC"
	    || ch_name == "nc_pio_energy_FC_ncpio_overlay" || ch_name == "nc_pio_score_FC_ncpio_overlay"
	    || ch_name == "nc_pio_energy_FC_ncdelta_overlay" || ch_name == "nc_pio_score_FC_ncdelta_overlay"
	    || ch_name == "nc_pio_energy_FC_overlay" || ch_name == "nc_pio_score_FC_overlay"
	    || ch_name == "nc_pio_energy_FC_ext" || ch_name == "nc_pio_score_FC_ext"
	    || ch_name == "nc_pio_energy_FC_dirt" || ch_name == "nc_pio_score_FC_dirt"
	    ){
    if (ch_name == "nc_pio_energy_FC"
	|| ch_name == "nc_pio_energy_FC_ext"
	|| ch_name == "nc_pio_energy_FC_dirt" ){
      if (flag_ncpio_sel && flag_FC) return true;
    }else if (ch_name == "nc_pio_energy_FC_ncpio_overlay" ){
      if (flag_ncpio_sel && flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_pio_energy_FC_ncdelta_overlay" ){
      if (flag_ncpio_sel && flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_energy_FC_overlay" ){
      if (flag_ncpio_sel && flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //if (flag_ncpio_sel && flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC"
	|| ch_name == "nc_pio_score_FC_ext"
	|| ch_name == "nc_pio_score_FC_dirt"){
      if (flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC_ncpio_overlay"){
      if (flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_pio_score_FC_ncdelta_overlay"){
      if (flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_score_FC_overlay"){
      if (flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
			&& (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //  if (flag_FC) return true;
    }

    return false;
  // start cuts from Lee's cuts.h
  // NC delta (1g) channels:
  }else if (ch_name == "nc_delta_0p_01" || ch_name == "nc_delta_0p_02" || ch_name == "nc_delta_0p_03" || ch_name == "nc_delta_0p_04"
                    || ch_name == "nc_delta_0p_05" || ch_name == "nc_delta_0p_06" || ch_name == "nc_delta_0p_07" || ch_name == "nc_delta_0p_08"
                    || ch_name == "nc_delta_0p_09" || ch_name == "nc_delta_0p_10" || ch_name == "nc_delta_0p_11" || ch_name == "nc_delta_0p_12"
                    || ch_name == "nc_delta_0p_13" || ch_name == "nc_delta_0p_14" || ch_name == "nc_delta_0p_15" || ch_name == "nc_delta_0p_16"
                    || ch_name == "nc_delta_0p_17" || ch_name == "nc_delta_0p_18" || ch_name == "nc_delta_0p_19" || ch_name == "nc_delta_0p_20"){
            if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01" || ch_name == "nc_delta_Np_02" || ch_name == "nc_delta_Np_03" || ch_name == "nc_delta_Np_04"
            || ch_name == "nc_delta_Np_05" || ch_name == "nc_delta_Np_06" || ch_name == "nc_delta_Np_07" || ch_name == "nc_delta_Np_08"
            || ch_name == "nc_delta_Np_09" || ch_name == "nc_delta_Np_10" || ch_name == "nc_delta_Np_11" || ch_name == "nc_delta_Np_12"
            || ch_name == "nc_delta_Np_13" || ch_name == "nc_delta_Np_14" || ch_name == "nc_delta_Np_15" || ch_name == "nc_delta_Np_16"
            || ch_name == "nc_delta_Np_17" || ch_name == "nc_delta_Np_18" || ch_name == "nc_delta_Np_19" || ch_name == "nc_delta_Np_20"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_ext" || ch_name == "nc_delta_0p_02_ext" || ch_name == "nc_delta_0p_03_ext" || ch_name == "nc_delta_0p_04_ext"
            || ch_name == "nc_delta_0p_05_ext" || ch_name == "nc_delta_0p_06_ext" || ch_name == "nc_delta_0p_07_ext" || ch_name == "nc_delta_0p_08_ext"
            || ch_name == "nc_delta_0p_09_ext" || ch_name == "nc_delta_0p_10_ext" || ch_name == "nc_delta_0p_11_ext" || ch_name == "nc_delta_0p_12_ext"
            || ch_name == "nc_delta_0p_13_ext" || ch_name == "nc_delta_0p_14_ext" || ch_name == "nc_delta_0p_15_ext" || ch_name == "nc_delta_0p_16_ext"
            || ch_name == "nc_delta_0p_17_ext" || ch_name == "nc_delta_0p_18_ext" || ch_name == "nc_delta_0p_19_ext" || ch_name == "nc_delta_0p_20_ext"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01_ext" || ch_name == "nc_delta_Np_02_ext" || ch_name == "nc_delta_Np_03_ext" || ch_name == "nc_delta_Np_04_ext"
            || ch_name == "nc_delta_Np_05_ext" || ch_name == "nc_delta_Np_06_ext" || ch_name == "nc_delta_Np_07_ext" || ch_name == "nc_delta_Np_08_ext"
            || ch_name == "nc_delta_Np_09_ext" || ch_name == "nc_delta_Np_10_ext" || ch_name == "nc_delta_Np_11_ext" || ch_name == "nc_delta_Np_12_ext"
            || ch_name == "nc_delta_Np_13_ext" || ch_name == "nc_delta_Np_14_ext" || ch_name == "nc_delta_Np_15_ext" || ch_name == "nc_delta_Np_16_ext"
            || ch_name == "nc_delta_Np_17_ext" || ch_name == "nc_delta_Np_18_ext" || ch_name == "nc_delta_Np_19_ext" || ch_name == "nc_delta_Np_20_ext"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_dirt" || ch_name == "nc_delta_0p_02_dirt" || ch_name == "nc_delta_0p_03_dirt" || ch_name == "nc_delta_0p_04_dirt"
            || ch_name == "nc_delta_0p_05_dirt" || ch_name == "nc_delta_0p_06_dirt" || ch_name == "nc_delta_0p_07_dirt" || ch_name == "nc_delta_0p_08_dirt"
            || ch_name == "nc_delta_0p_09_dirt" || ch_name == "nc_delta_0p_10_dirt" || ch_name == "nc_delta_0p_11_dirt" || ch_name == "nc_delta_0p_12_dirt"
            || ch_name == "nc_delta_0p_13_dirt" || ch_name == "nc_delta_0p_14_dirt" || ch_name == "nc_delta_0p_15_dirt" || ch_name == "nc_delta_0p_16_dirt"
            || ch_name == "nc_delta_0p_17_dirt" || ch_name == "nc_delta_0p_18_dirt" || ch_name == "nc_delta_0p_19_dirt" || ch_name == "nc_delta_0p_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01_dirt" || ch_name == "nc_delta_Np_02_dirt" || ch_name == "nc_delta_Np_03_dirt" || ch_name == "nc_delta_Np_04_dirt"
            || ch_name == "nc_delta_Np_05_dirt" || ch_name == "nc_delta_Np_06_dirt" || ch_name == "nc_delta_Np_07_dirt" || ch_name == "nc_delta_Np_08_dirt"
            || ch_name == "nc_delta_Np_09_dirt" || ch_name == "nc_delta_Np_10_dirt" || ch_name == "nc_delta_Np_11_dirt" || ch_name == "nc_delta_Np_12_dirt"
            || ch_name == "nc_delta_Np_13_dirt" || ch_name == "nc_delta_Np_14_dirt" || ch_name == "nc_delta_Np_15_dirt" || ch_name == "nc_delta_Np_16_dirt"
            || ch_name == "nc_delta_Np_17_dirt" || ch_name == "nc_delta_Np_18_dirt" || ch_name == "nc_delta_Np_19_dirt" || ch_name == "nc_delta_Np_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_nc_delta_overlay" || ch_name == "nc_delta_0p_01_nc_delta_overlay_add" || ch_name == "nc_delta_0p_02_nc_delta_overlay" || ch_name == "nc_delta_0p_03_nc_delta_overlay" || ch_name == "nc_delta_0p_04_nc_delta_overlay"
            || ch_name == "nc_delta_0p_05_nc_delta_overlay" || ch_name == "nc_delta_0p_06_nc_delta_overlay" || ch_name == "nc_delta_0p_07_nc_delta_overlay" || ch_name == "nc_delta_0p_08_nc_delta_overlay"
            || ch_name == "nc_delta_0p_09_nc_delta_overlay" || ch_name == "nc_delta_0p_10_nc_delta_overlay" || ch_name == "nc_delta_0p_11_nc_delta_overlay" || ch_name == "nc_delta_0p_12_nc_delta_overlay"
            || ch_name == "nc_delta_0p_13_nc_delta_overlay" || ch_name == "nc_delta_0p_14_nc_delta_overlay" || ch_name == "nc_delta_0p_15_nc_delta_overlay" || ch_name == "nc_delta_0p_16_nc_delta_overlay"
            || ch_name == "nc_delta_0p_17_nc_delta_overlay" || ch_name == "nc_delta_0p_18_nc_delta_overlay" || ch_name == "nc_delta_0p_19_nc_delta_overlay" || ch_name == "nc_delta_0p_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_delta_overlay" || ch_name == "nc_delta_Np_01_nc_delta_overlay_add" || ch_name == "nc_delta_Np_02_nc_delta_overlay" || ch_name == "nc_delta_Np_03_nc_delta_overlay" || ch_name == "nc_delta_Np_04_nc_delta_overlay"
            || ch_name == "nc_delta_Np_05_nc_delta_overlay" || ch_name == "nc_delta_Np_06_nc_delta_overlay" || ch_name == "nc_delta_Np_07_nc_delta_overlay" || ch_name == "nc_delta_Np_08_nc_delta_overlay"
            || ch_name == "nc_delta_Np_09_nc_delta_overlay" || ch_name == "nc_delta_Np_10_nc_delta_overlay" || ch_name == "nc_delta_Np_11_nc_delta_overlay" || ch_name == "nc_delta_Np_12_nc_delta_overlay"
            || ch_name == "nc_delta_Np_13_nc_delta_overlay" || ch_name == "nc_delta_Np_14_nc_delta_overlay" || ch_name == "nc_delta_Np_15_nc_delta_overlay" || ch_name == "nc_delta_Np_16_nc_delta_overlay"
            || ch_name == "nc_delta_Np_17_nc_delta_overlay" || ch_name == "nc_delta_Np_18_nc_delta_overlay" || ch_name == "nc_delta_Np_19_nc_delta_overlay" || ch_name == "nc_delta_Np_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_true_Np_nc_delta_overlay" || ch_name == "nc_delta_0p_01_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_delta_0p_01_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_delta_0p_01_true_Np_nc_delta_overlay"));
            }
            return false;

    }else if (ch_name == "nc_delta_Np_01_true_0p_nc_delta_overlay" || ch_name == "nc_delta_Np_01_true_Np_nc_delta_overlay"){
            if (flag_FC && flag_ncdelta_sel && !(flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_delta_Np_01_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_delta_Np_01_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_delta_0p_01_nc_pi0_overlay" || ch_name == "nc_delta_0p_02_nc_pi0_overlay" || ch_name == "nc_delta_0p_03_nc_pi0_overlay" || ch_name == "nc_delta_0p_04_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_05_nc_pi0_overlay" || ch_name == "nc_delta_0p_06_nc_pi0_overlay" || ch_name == "nc_delta_0p_07_nc_pi0_overlay" || ch_name == "nc_delta_0p_08_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_09_nc_pi0_overlay" || ch_name == "nc_delta_0p_10_nc_pi0_overlay" || ch_name == "nc_delta_0p_11_nc_pi0_overlay" || ch_name == "nc_delta_0p_12_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_13_nc_pi0_overlay" || ch_name == "nc_delta_0p_14_nc_pi0_overlay" || ch_name == "nc_delta_0p_15_nc_pi0_overlay" || ch_name == "nc_delta_0p_16_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_17_nc_pi0_overlay" || ch_name == "nc_delta_0p_18_nc_pi0_overlay" || ch_name == "nc_delta_0p_19_nc_pi0_overlay" || ch_name == "nc_delta_0p_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_pi0_overlay" || ch_name == "nc_delta_Np_02_nc_pi0_overlay" || ch_name == "nc_delta_Np_03_nc_pi0_overlay" || ch_name == "nc_delta_Np_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_05_nc_pi0_overlay" || ch_name == "nc_delta_Np_06_nc_pi0_overlay" || ch_name == "nc_delta_Np_07_nc_pi0_overlay" || ch_name == "nc_delta_Np_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_09_nc_pi0_overlay" || ch_name == "nc_delta_Np_10_nc_pi0_overlay" || ch_name == "nc_delta_Np_11_nc_pi0_overlay" || ch_name == "nc_delta_Np_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_13_nc_pi0_overlay" || ch_name == "nc_delta_Np_14_nc_pi0_overlay" || ch_name == "nc_delta_Np_15_nc_pi0_overlay" || ch_name == "nc_delta_Np_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_17_nc_pi0_overlay" || ch_name == "nc_delta_Np_18_nc_pi0_overlay" || ch_name == "nc_delta_Np_19_nc_pi0_overlay" || ch_name == "nc_delta_Np_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                        && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
   }else if (ch_name == "nc_delta_0p_01_overlay" || ch_name == "nc_delta_0p_02_overlay" || ch_name == "nc_delta_0p_03_overlay" || ch_name == "nc_delta_0p_04_overlay"
            || ch_name == "nc_delta_0p_05_overlay" || ch_name == "nc_delta_0p_06_overlay" || ch_name == "nc_delta_0p_07_overlay" || ch_name == "nc_delta_0p_08_overlay"
            || ch_name == "nc_delta_0p_09_overlay" || ch_name == "nc_delta_0p_10_overlay" || ch_name == "nc_delta_0p_11_overlay" || ch_name == "nc_delta_0p_12_overlay"
            || ch_name == "nc_delta_0p_13_overlay" || ch_name == "nc_delta_0p_14_overlay" || ch_name == "nc_delta_0p_15_overlay" || ch_name == "nc_delta_0p_16_overlay"
            || ch_name == "nc_delta_0p_17_overlay" || ch_name == "nc_delta_0p_18_overlay" || ch_name == "nc_delta_0p_19_overlay" || ch_name == "nc_delta_0p_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay" || ch_name == "nc_delta_Np_02_overlay" || ch_name == "nc_delta_Np_03_overlay" || ch_name == "nc_delta_Np_04_overlay"
            || ch_name == "nc_delta_Np_05_overlay" || ch_name == "nc_delta_Np_06_overlay" || ch_name == "nc_delta_Np_07_overlay" || ch_name == "nc_delta_Np_08_overlay"
            || ch_name == "nc_delta_Np_09_overlay" || ch_name == "nc_delta_Np_10_overlay" || ch_name == "nc_delta_Np_11_overlay" || ch_name == "nc_delta_Np_12_overlay"
            || ch_name == "nc_delta_Np_13_overlay" || ch_name == "nc_delta_Np_14_overlay" || ch_name == "nc_delta_Np_15_overlay" || ch_name == "nc_delta_Np_16_overlay"
            || ch_name == "nc_delta_Np_17_overlay" || ch_name == "nc_delta_Np_18_overlay" || ch_name == "nc_delta_Np_19_overlay" || ch_name == "nc_delta_Np_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_overlay_entire" || ch_name == "nc_delta_0p_02_overlay_entire" || ch_name == "nc_delta_0p_03_overlay_entire" || ch_name == "nc_delta_0p_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay_entire" || ch_name == "nc_delta_Np_02_overlay_entire" || ch_name == "nc_delta_Np_03_overlay_entire" || ch_name == "nc_delta_Np_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01" || ch_name == "nc_delta_Xp_02" || ch_name == "nc_delta_Xp_03" || ch_name == "nc_delta_Xp_04"
            || ch_name == "nc_delta_Xp_05" || ch_name == "nc_delta_Xp_06" || ch_name == "nc_delta_Xp_07" || ch_name == "nc_delta_Xp_08"
            || ch_name == "nc_delta_Xp_09" || ch_name == "nc_delta_Xp_10" || ch_name == "nc_delta_Xp_11" || ch_name == "nc_delta_Xp_12"
            || ch_name == "nc_delta_Xp_13" || ch_name == "nc_delta_Xp_14" || ch_name == "nc_delta_Xp_15" || ch_name == "nc_delta_Xp_16"
            || ch_name == "nc_delta_Xp_17" || ch_name == "nc_delta_Xp_18" || ch_name == "nc_delta_Xp_19" || ch_name == "nc_delta_Xp_20"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_ext" || ch_name == "nc_delta_Xp_02_ext" || ch_name == "nc_delta_Xp_03_ext" || ch_name == "nc_delta_Xp_04_ext"
            || ch_name == "nc_delta_Xp_05_ext" || ch_name == "nc_delta_Xp_06_ext" || ch_name == "nc_delta_Xp_07_ext" || ch_name == "nc_delta_Xp_08_ext"
            || ch_name == "nc_delta_Xp_09_ext" || ch_name == "nc_delta_Xp_10_ext" || ch_name == "nc_delta_Xp_11_ext" || ch_name == "nc_delta_Xp_12_ext"
            || ch_name == "nc_delta_Xp_13_ext" || ch_name == "nc_delta_Xp_14_ext" || ch_name == "nc_delta_Xp_15_ext" || ch_name == "nc_delta_Xp_16_ext"
            || ch_name == "nc_delta_Xp_17_ext" || ch_name == "nc_delta_Xp_18_ext" || ch_name == "nc_delta_Xp_19_ext" || ch_name == "nc_delta_Xp_20_ext"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_dirt" || ch_name == "nc_delta_Xp_02_dirt" || ch_name == "nc_delta_Xp_03_dirt" || ch_name == "nc_delta_Xp_04_dirt"
            || ch_name == "nc_delta_Xp_05_dirt" || ch_name == "nc_delta_Xp_06_dirt" || ch_name == "nc_delta_Xp_07_dirt" || ch_name == "nc_delta_Xp_08_dirt"
            || ch_name == "nc_delta_Xp_09_dirt" || ch_name == "nc_delta_Xp_10_dirt" || ch_name == "nc_delta_Xp_11_dirt" || ch_name == "nc_delta_Xp_12_dirt"
            || ch_name == "nc_delta_Xp_13_dirt" || ch_name == "nc_delta_Xp_14_dirt" || ch_name == "nc_delta_Xp_15_dirt" || ch_name == "nc_delta_Xp_16_dirt"
            || ch_name == "nc_delta_Xp_17_dirt" || ch_name == "nc_delta_Xp_18_dirt" || ch_name == "nc_delta_Xp_19_dirt" || ch_name == "nc_delta_Xp_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_delta_overlay" || ch_name == "nc_delta_Xp_02_nc_delta_overlay" || ch_name == "nc_delta_Xp_03_nc_delta_overlay" || ch_name == "nc_delta_Xp_04_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_05_nc_delta_overlay" || ch_name == "nc_delta_Xp_06_nc_delta_overlay" || ch_name == "nc_delta_Xp_07_nc_delta_overlay" || ch_name == "nc_delta_Xp_08_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_09_nc_delta_overlay" || ch_name == "nc_delta_Xp_10_nc_delta_overlay" || ch_name == "nc_delta_Xp_11_nc_delta_overlay" || ch_name == "nc_delta_Xp_12_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_13_nc_delta_overlay" || ch_name == "nc_delta_Xp_14_nc_delta_overlay" || ch_name == "nc_delta_Xp_15_nc_delta_overlay" || ch_name == "nc_delta_Xp_16_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_17_nc_delta_overlay" || ch_name == "nc_delta_Xp_18_nc_delta_overlay" || ch_name == "nc_delta_Xp_19_nc_delta_overlay" || ch_name == "nc_delta_Xp_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_pi0_overlay" || ch_name == "nc_delta_Xp_02_nc_pi0_overlay" || ch_name == "nc_delta_Xp_03_nc_pi0_overlay" || ch_name == "nc_delta_Xp_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_05_nc_pi0_overlay" || ch_name == "nc_delta_Xp_06_nc_pi0_overlay" || ch_name == "nc_delta_Xp_07_nc_pi0_overlay" || ch_name == "nc_delta_Xp_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_09_nc_pi0_overlay" || ch_name == "nc_delta_Xp_10_nc_pi0_overlay" || ch_name == "nc_delta_Xp_11_nc_pi0_overlay" || ch_name == "nc_delta_Xp_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_13_nc_pi0_overlay" || ch_name == "nc_delta_Xp_14_nc_pi0_overlay" || ch_name == "nc_delta_Xp_15_nc_pi0_overlay" || ch_name == "nc_delta_Xp_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_17_nc_pi0_overlay" || ch_name == "nc_delta_Xp_18_nc_pi0_overlay" || ch_name == "nc_delta_Xp_19_nc_pi0_overlay" || ch_name == "nc_delta_Xp_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_overlay" || ch_name == "nc_delta_Xp_02_overlay" || ch_name == "nc_delta_Xp_03_overlay" || ch_name == "nc_delta_Xp_04_overlay"
            || ch_name == "nc_delta_Xp_05_overlay" || ch_name == "nc_delta_Xp_06_overlay" || ch_name == "nc_delta_Xp_07_overlay" || ch_name == "nc_delta_Xp_08_overlay"
            || ch_name == "nc_delta_Xp_09_overlay" || ch_name == "nc_delta_Xp_10_overlay" || ch_name == "nc_delta_Xp_11_overlay" || ch_name == "nc_delta_Xp_12_overlay"
            || ch_name == "nc_delta_Xp_13_overlay" || ch_name == "nc_delta_Xp_14_overlay" || ch_name == "nc_delta_Xp_15_overlay" || ch_name == "nc_delta_Xp_16_overlay"
            || ch_name == "nc_delta_Xp_17_overlay" || ch_name == "nc_delta_Xp_18_overlay" || ch_name == "nc_delta_Xp_19_overlay" || ch_name == "nc_delta_Xp_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // NC Pi0 channels (with NC delta selected events removed):
    }else if (ch_name == "nc_pi0_0p" || ch_name == "nc_pi0_2_0p"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np" || ch_name == "nc_pi0_2_Np"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp" || ch_name == "nc_pi0_2_Xp"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_ext" || ch_name == "nc_pi0_2_0p_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_ext" || ch_name == "nc_pi0_2_Np_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_ext" || ch_name == "nc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_dirt" || ch_name == "nc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_dirt" || ch_name == "nc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_dirt" || ch_name == "nc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_nc_delta_overlay" || ch_name == "nc_pi0_0p_nc_delta_overlay_add" || ch_name == "nc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_delta_overlay" || ch_name == "nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_pi0_Np_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_pi0_Np_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_pi0_0p_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_0p_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_pi0_0p_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_pi0_0p_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_pi0_0p_nc_pi0_overlay" || ch_name == "nc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_pi0_overlay" || ch_name == "nc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_pi0_overlay" || ch_name == "nc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_overlay" || ch_name == "nc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_overlay" || ch_name == "nc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_overlay" || ch_name == "nc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // CC pi0 channels (with NC delta and NC pi0 selections removed) (not used for final NC delta fits):
     }else if (ch_name == "cc_pi0_0p" || ch_name == "cc_pi0_2_0p"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np" || ch_name == "cc_pi0_2_Np"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp" || ch_name == "cc_pi0_2_Xp"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_ext" || ch_name == "cc_pi0_2_0p_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_ext" || ch_name == "cc_pi0_2_Np_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_ext" || ch_name == "cc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_dirt" || ch_name == "cc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_dirt" || ch_name == "cc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_dirt" || ch_name == "cc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_delta_overlay" || ch_name == "cc_pi0_0p_nc_delta_overlay_add" || ch_name == "cc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_delta_overlay" || ch_name == "cc_pi0_Np_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_delta_overlay" || ch_name == "cc_pi0_Xp_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_pi0_overlay" || ch_name == "cc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_pi0_overlay" || ch_name == "cc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_pi0_overlay" || ch_name == "cc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_overlay" || ch_name == "cc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_overlay" || ch_name == "cc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_overlay" || ch_name == "cc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // numuCC channels (with NC delta, NC pi0, and CC pi0 events removed) (used for some constraint tests, not used for final fits):
    }else if (ch_name == "numuCC_noCCpi0_0p"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)  && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_0p_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_Np_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_Xp_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)  && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // numuCC channels (with NC delta and NC pi0 events removed):
    }else if (ch_name == "numuCC_0p"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_nc_delta_overlay" || ch_name == "numuCC_0p_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_delta_overlay" || ch_name == "numuCC_Np_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_delta_overlay" || ch_name == "numuCC_Xp_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_true_Np_nc_delta_overlay" || ch_name == "numuCC_Np_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="numuCC_Np_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="numuCC_Np_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "numuCC_0p_true_Np_nc_delta_overlay" || ch_name == "numuCC_0p_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="numuCC_0p_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="numuCC_0p_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "numuCC_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    // Pi0 1p selections (used for proton-pi0 invariant mass test)
    }else if (ch_name == "nc_pi0_1p" || ch_name == "nc_pi0_1p_ext" || ch_name == "nc_pi0_1p_dirt" || ch_name == "nc_pi0_1p_overlay_entire" || ch_name == "nc_pi0_2_1p" || ch_name == "nc_pi0_2_1p_ext" || ch_name == "nc_pi0_2_1p_dirt" || ch_name == "nc_pi0_2_1p_overlay_entire"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_1p) return true;
                  return false;

    }else if (ch_name == "cc_pi0_1p" || ch_name == "cc_pi0_1p_ext" || ch_name == "cc_pi0_1p_dirt" || ch_name == "cc_pi0_1p_overlay_entire" || ch_name == "cc_pi0_2_1p" || ch_name == "cc_pi0_2_1p_ext" || ch_name == "cc_pi0_2_1p_dirt" || ch_name == "cc_pi0_2_1p_overlay_entire"){
            if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_1p) return true;
                  return false;


  // end cuts from Lee's cuts.h
  }else if (ch_name == "nc_delta_energy_FC_0p" || ch_name == "nc_delta_score_FC_0p"
	    || ch_name == "nc_delta_energy_FC_0p_ncpio_overlay" || ch_name == "nc_delta_score_FC_0p_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay" || ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_overlay" || ch_name == "nc_delta_score_FC_0p_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_ext"
	    || ch_name == "nc_delta_energy_FC_0p_dirt" || ch_name == "nc_delta_score_FC_0p_dirt"){

    if (ch_name == "nc_delta_energy_FC_0p" ||  ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_energy_FC_0p_dirt"){
      if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
						     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_FC_0p" ||  ch_name == "nc_delta_score_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_dirt"){
      if (flag_FC  && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncpio_overlay"){
      if (flag_FC  && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
				  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"){
      if (flag_FC  && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_overlay"){
      if (flag_FC  && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC  && flag_0p) return true;
    }

    return false;
  }else if (ch_name == "nc_delta_energy_FC_Np" || ch_name == "nc_delta_score_FC_Np"
	    || ch_name == "nc_delta_energy_FC_Np_ncpio_overlay" || ch_name == "nc_delta_score_FC_Np_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay" || ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_overlay" || ch_name == "nc_delta_score_FC_Np_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_ext"
	    || ch_name == "nc_delta_energy_FC_Np_dirt" || ch_name == "nc_delta_score_FC_Np_dirt"){

    if (ch_name == "nc_delta_energy_FC_Np" ||  ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_energy_FC_Np_dirt"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np" ||  ch_name == "nc_delta_score_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_dirt"){
      if (flag_FC  && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncpio_overlay"){
      if (flag_FC  && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"){
      if (flag_FC  && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_overlay"){
      if (flag_FC  && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      // if (flag_FC  && (!flag_0p)) return true;
    }

    return false;

  //Erin
  }else if (ch_name == "all_bnb" || ch_name == "all_ext" || ch_name == "all_dirt" || ch_name == "all_bnb_LEE"){
    return true;
  }else if (ch_name == "all_bnb_nsbeam"){
    if (flag_nsbeam) return true;
    else return false;
  }else if (ch_name == "nodata_bnb"){
    return false;
  }else if (ch_name == "all_spoverlay" || ch_name == "all_spoverlay_2" || ch_name == "all_spoverlay_3"){
            if ((eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1)) return true;
            return false;
  }else if (ch_name == "all_ncpi0overlay" || ch_name == "all_ncpi0overlay_2" || ch_name == "all_ncpi0overlay_3"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "all_overlay_sp_ncpi0_BG" || ch_name == "all_overlay_sp_ncpi0_BG_2" || ch_name == "all_overlay_sp_ncpi0_BG_3"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "all_overlay_sp_BG"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "generic_nu_ext" || ch_name == "generic_nu_dirt" ||ch_name == "generic_nu_bnb_LEE"){
    if (flag_generic) return true;
    else return false;
  }else if (ch_name == "generic_nu_bnb_nsbeam"){
    if (flag_generic && flag_nsbeam) return true;
    else return false;
  }else if (ch_name == "nodata_bnb"){
    return false;
  }else if (ch_name == "generic_nu_spoverlay" || ch_name == "generic_nu_spoverlay_2" || ch_name == "generic_nu_spoverlay_3"){
            if (flag_generic &&
              (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_ncpi0overlay" || ch_name == "generic_nu_ncpi0overlay_2" || ch_name == "generic_nu_ncpi0overlay_3"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_overlay_sp_ncpi0_BG" || ch_name == "generic_nu_overlay_sp_ncpi0_BG_2" || ch_name == "generic_nu_overlay_sp_ncpi0_BG_3"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_overlay_sp_BG"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb" || ch_name == "single_photon_ext"
    || ch_name == "single_photon_overlay" || ch_name == "single_photon_dirt" 
    || ch_name == "single_photon_LEE"){
            if (flag_singlephoton_sel) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb" || ch_name == "single_photon_eff_ext"
    || ch_name == "single_photon_eff_overlay" || ch_name == "single_photon_eff_dirt" 
    || ch_name == "single_photon_eff_LEE"){
            if (flag_singlephoton_eff_sel) return true;
            return false;
  }else if (ch_name == "single_shower_bnb" || ch_name == "single_shower_ext"
    || ch_name == "single_shower_overlay" || ch_name == "single_shower_dirt" 
    || ch_name == "single_shower_LEE"){
            if (flag_singleshower_sel) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb" || ch_name == "single_shower_eff_ext"
    || ch_name == "single_shower_eff_overlay" || ch_name == "single_shower_eff_dirt"
    || ch_name == "single_shower_eff_bnb_2" || ch_name == "single_shower_eff_ext_2"
      || ch_name == "single_shower_eff_overlay_2" || ch_name == "single_shower_eff_dirt_2" 
    || ch_name == "single_shower_eff_LEE"){
            if (flag_singleshower_eff_sel) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p" || ch_name == "single_photon_ext_0p"
    || ch_name == "single_photon_overlay_0p" || ch_name == "single_photon_dirt_0p" 
    || ch_name == "single_photon_LEE_0p"){
            if (flag_singlephoton_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p" || ch_name == "single_photon_eff_ext_0p"
    || ch_name == "single_photon_eff_overlay_0p" || ch_name == "single_photon_eff_dirt_0p" 
    || ch_name == "single_photon_eff_LEE_0p"){
            if (flag_singlephoton_eff_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p" || ch_name == "single_shower_ext_0p"
    || ch_name == "single_shower_overlay_0p" || ch_name == "single_shower_dirt_0p" 
    || ch_name == "single_shower_LEE_0p"){
            if (flag_singleshower_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p" || ch_name == "single_shower_eff_ext_0p"
    || ch_name == "single_shower_eff_overlay_0p" || ch_name == "single_shower_eff_dirt_0p" 
    || ch_name == "single_shower_eff_LEE_0p"){
            if (flag_singleshower_eff_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np" || ch_name == "single_photon_ext_Np"
    || ch_name == "single_photon_overlay_Np" || ch_name == "single_photon_dirt_Np"  
    || ch_name == "single_photon_LEE_Np"){
            if (flag_singlephoton_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np" || ch_name == "single_photon_eff_ext_Np"
    || ch_name == "single_photon_eff_overlay_Np" || ch_name == "single_photon_eff_dirt_Np" 
    || ch_name == "single_photon_eff_LEE_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np" || ch_name == "single_shower_ext_Np"
    || ch_name == "single_shower_overlay_Np" || ch_name == "single_shower_dirt_Np" 
    || ch_name == "single_shower_LEE_Np"){
            if (flag_singleshower_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np" || ch_name == "single_shower_eff_ext_Np"
    || ch_name == "single_shower_eff_overlay_Np" || ch_name == "single_shower_eff_dirt_Np" 
    || ch_name == "single_shower_eff_LEE_Np"){
            if (flag_singleshower_eff_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_nodelta"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_delta"){
            if (flag_singlephoton_sel &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay" || ch_name == "single_photon_spoverlay_lee"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_otpc" || ch_name == "single_photon_spoverlay_otpc_lee"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay"){
            if (flag_singlephoton_eff_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay" || ch_name == "single_shower_spoverlay_lee"){
            if (flag_singleshower_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay"){
            if (flag_singleshower_eff_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p" || ch_name == "single_photon_spoverlay_lee_0p"){
            if (flag_singlephoton_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p"){
            if (flag_singleshower_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p"){
            if (flag_singleshower_eff_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np" || ch_name == "single_photon_spoverlay_lee_Np"){
            if (flag_singlephoton_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np"){
            if (flag_singleshower_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay
  }else if (ch_name == "single_photon_ncpi0overlay"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_otpc"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              //map_cuts_flag["SPNCPi0Sig"] || 
              map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              //&& !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              //&& pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              //&& pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam"){
            if (flag_singlephoton_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam"){
            if (flag_singlephoton_eff_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam"){
            if (flag_singleshower_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam"){
            if (flag_singleshower_eff_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam"){
            if (flag_singlephoton_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam"){
            if (flag_singlephoton_eff_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam"){
            if (flag_singleshower_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam"){
            if (flag_singleshower_eff_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam"){
            if (flag_singlephoton_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam"){
            if (flag_singlephoton_eff_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam"){
            if (flag_singleshower_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam"){
            if (flag_singleshower_eff_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  //FC
  }else if (ch_name == "single_photon_bnb_FC" || ch_name == "single_photon_ext_FC"
    || ch_name == "single_photon_overlay_FC" || ch_name == "single_photon_dirt_FC"){
            if (flag_singlephoton_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_FC" || ch_name == "single_photon_eff_ext_FC"
    || ch_name == "single_photon_eff_overlay_FC" || ch_name == "single_photon_eff_dirt_FC"){
            if (flag_singlephoton_eff_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_FC" || ch_name == "single_shower_ext_FC"
    || ch_name == "single_shower_overlay_FC" || ch_name == "single_shower_dirt_FC"){
            if (flag_singleshower_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_FC" || ch_name == "single_shower_eff_ext_FC"
    || ch_name == "single_shower_eff_overlay_FC" || ch_name == "single_shower_eff_dirt_FC"
    || ch_name == "single_shower_eff_bnb_2_FC" || ch_name == "single_shower_eff_ext_2_FC"
      || ch_name == "single_shower_eff_overlay_2_FC" || ch_name == "single_shower_eff_dirt_2_FC"){
            if (flag_singleshower_eff_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_FC" || ch_name == "single_photon_ext_0p_FC"
    || ch_name == "single_photon_overlay_0p_FC" || ch_name == "single_photon_dirt_0p_FC"){
            if (flag_singlephoton_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_FC" || ch_name == "single_photon_eff_ext_0p_FC"
    || ch_name == "single_photon_eff_overlay_0p_FC" || ch_name == "single_photon_eff_dirt_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_FC" || ch_name == "single_shower_ext_0p_FC"
    || ch_name == "single_shower_overlay_0p_FC" || ch_name == "single_shower_dirt_0p_FC"){
            if (flag_singleshower_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_FC" || ch_name == "single_shower_eff_ext_0p_FC"
    || ch_name == "single_shower_eff_overlay_0p_FC" || ch_name == "single_shower_eff_dirt_0p_FC"){
            if (flag_singleshower_eff_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_FC" || ch_name == "single_photon_ext_Np_FC"
    || ch_name == "single_photon_overlay_Np_FC" || ch_name == "single_photon_dirt_Np_FC"){
            if (flag_singlephoton_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_FC" || ch_name == "single_photon_eff_ext_Np_FC"
    || ch_name == "single_photon_eff_overlay_Np_FC" || ch_name == "single_photon_eff_dirt_Np_FC"){
            if (flag_singlephoton_eff_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_FC" || ch_name == "single_shower_ext_Np_FC"
    || ch_name == "single_shower_overlay_Np_FC" || ch_name == "single_shower_dirt_Np_FC"){
            if (flag_singleshower_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_FC" || ch_name == "single_shower_eff_ext_Np_FC"
    || ch_name == "single_shower_eff_overlay_Np_FC" || ch_name == "single_shower_eff_dirt_Np_FC"){
            if (flag_singleshower_eff_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_FC"){
            if (flag_singleshower_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay FC
  }else if (ch_name == "single_photon_ncpi0overlay_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  //PC
  }else if (ch_name == "single_photon_bnb_PC" || ch_name == "single_photon_ext_PC"
    || ch_name == "single_photon_overlay_PC" || ch_name == "single_photon_dirt_PC"){
            if (flag_singlephoton_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_PC" || ch_name == "single_photon_eff_ext_PC"
    || ch_name == "single_photon_eff_overlay_PC" || ch_name == "single_photon_eff_dirt_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_PC" || ch_name == "single_shower_ext_PC"
    || ch_name == "single_shower_overlay_PC" || ch_name == "single_shower_dirt_PC"){
            if (flag_singleshower_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_PC" || ch_name == "single_shower_eff_ext_PC"
    || ch_name == "single_shower_eff_overlay_PC" || ch_name == "single_shower_eff_dirt_PC"
    || ch_name == "single_shower_eff_bnb_2_PC" || ch_name == "single_shower_eff_ext_2_PC"
      || ch_name == "single_shower_eff_overlay_2_PC" || ch_name == "single_shower_eff_dirt_2_PC"){
            if (flag_singleshower_eff_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_PC" || ch_name == "single_photon_ext_0p_PC"
    || ch_name == "single_photon_overlay_0p_PC" || ch_name == "single_photon_dirt_0p_PC"){
            if (flag_singlephoton_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_PC" || ch_name == "single_photon_eff_ext_0p_PC"
    || ch_name == "single_photon_eff_overlay_0p_PC" || ch_name == "single_photon_eff_dirt_0p_PC"){
            if (flag_singlephoton_eff_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_PC" || ch_name == "single_shower_ext_0p_PC"
    || ch_name == "single_shower_overlay_0p_PC" || ch_name == "single_shower_dirt_0p_PC"){
            if (flag_singleshower_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_PC" || ch_name == "single_shower_eff_ext_0p_PC"
    || ch_name == "single_shower_eff_overlay_0p_PC" || ch_name == "single_shower_eff_dirt_0p_PC"){
            if (flag_singleshower_eff_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_PC" || ch_name == "single_photon_ext_Np_PC"
    || ch_name == "single_photon_overlay_Np_PC" || ch_name == "single_photon_dirt_Np_PC"){
            if (flag_singlephoton_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_PC" || ch_name == "single_photon_eff_ext_Np_PC"
    || ch_name == "single_photon_eff_overlay_Np_PC" || ch_name == "single_photon_eff_dirt_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_PC" || ch_name == "single_shower_ext_Np_PC"
    || ch_name == "single_shower_overlay_Np_PC" || ch_name == "single_shower_dirt_Np_PC"){
            if (flag_singleshower_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_PC" || ch_name == "single_shower_eff_ext_Np_PC"
    || ch_name == "single_shower_eff_overlay_Np_PC" || ch_name == "single_shower_eff_dirt_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay PC
  }else if (ch_name == "single_photon_ncpi0overlay_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  //nue bdt cut eff
  }else if (ch_name == "single_photon_eff_nue_overlay" ||
            ch_name == "single_photon_eff_nue_dirt" ||
            ch_name == "single_photon_eff_nue_ext" ||
            ch_name == "single_photon_eff_nue_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  //ncpi0 bdt cut
  }else if (ch_name == "single_photon_eff_ncpi0_overlay" ||
            ch_name == "single_photon_eff_ncpi0_dirt" ||
            ch_name == "single_photon_eff_ncpi0_ext" ||
            ch_name == "single_photon_eff_ncpi0_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay" ||
            ch_name == "single_photon_eff_other_dirt" ||
            ch_name == "single_photon_eff_other_ext" ||
            ch_name == "single_photon_eff_other_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay" ||
            ch_name == "single_photon_eff_numu_dirt" ||
            ch_name == "single_photon_eff_numu_ext" ||
            ch_name == "single_photon_eff_numu_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
    }else if (ch_name == "single_photon_eff_pre_overlay" ||
              ch_name == "single_photon_eff_pre_dirt" ||
              ch_name == "single_photon_eff_pre_ext" ||
              ch_name == "single_photon_eff_pre_bnb"){
      if (flag_singlephoton_pre) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_bnb_nsbeam"){
      if (flag_singlephoton_pre && flag_nsbeam_photon) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_spoverlay"){
      if (flag_singlephoton_pre &&
          (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
          map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
          map_cuts_flag["SPNumuCCSig"])) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_overlay_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
           map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
           map_cuts_flag["SPNumuCCSig"])
           && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
           && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
           && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_ncpi0overlay"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_overlay_ncpi0_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;

  //nue bdt cut pure
  }else if (ch_name == "single_photon_nue_overlay" ||
            ch_name == "single_photon_nue_dirt" ||
            ch_name == "single_photon_nue_ext" ||
            ch_name == "single_photon_nue_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue) return true;
      return false;
  }else if (ch_name == "single_photon_nue_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_nue_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_nue_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_nue_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_nue_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  //ncpi0 bdt cut
  }else if (ch_name == "single_photon_ncpi0_overlay" ||
            ch_name == "single_photon_ncpi0_dirt" ||
            ch_name == "single_photon_ncpi0_ext" ||
            ch_name == "single_photon_ncpi0_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay" ||
            ch_name == "single_photon_other_dirt" ||
            ch_name == "single_photon_other_ext" ||
            ch_name == "single_photon_other_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other) return true;
      return false;
  }else if (ch_name == "single_photon_other_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_other_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_other_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay" ||
            ch_name == "single_photon_numu_dirt" ||
            ch_name == "single_photon_numu_ext" ||
            ch_name == "single_photon_numu_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu) return true;
      return false;
  }else if (ch_name == "single_photon_numu_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_numu_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
    }else if (ch_name == "single_photon_pre_overlay" ||
              ch_name == "single_photon_pre_dirt" ||
              ch_name == "single_photon_pre_ext" ||
              ch_name == "single_photon_pre_bnb"){
      if (flag_singlephoton_pre) return true;
        return false;
    }else if (ch_name == "single_photon_pre_bnb_nsbeam"){
      if (flag_singlephoton_pre && flag_nsbeam_photon) return true;
        return false;
    }else if (ch_name == "single_photon_pre_spoverlay"){
      if (flag_singlephoton_pre &&
          (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
          map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
          map_cuts_flag["SPNumuCCSig"])) return true;
        return false;
    }else if (ch_name == "single_photon_pre_overlay_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
           map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
           map_cuts_flag["SPNumuCCSig"])
           && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
           && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
           && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
        return false;
    }else if (ch_name == "single_photon_pre_ncpi0overlay"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;
    }else if (ch_name == "single_photon_pre_overlay_ncpi0_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;

  }else if (ch_name == "single_photon_ncdel_sig_overlay" || ch_name == "single_photon_ncpi0_sig_overlay" ||
    ch_name == "single_photon_ncother_sig_overlay" || ch_name == "single_photon_numucc_sig_overlay" ||
    ch_name == "single_photon_outfv_sig_overlay" || ch_name == "single_photon_bkg_overlay"){
      if (ch_name == "single_photon_ncdel_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCDeltaSig"]) return true;
      }else if (ch_name == "single_photon_ncpi0_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCPi0Sig"]) return true;
      }else if (ch_name == "single_photon_ncother_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCOtherSig"]) return true;
      }else if (ch_name == "single_photon_numucc_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNumuCCSig"]) return true;
      }else if (ch_name == "single_photon_outfv_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPOutFVSig"]) return true;
      }else if (ch_name == "single_photon_bkg_overlay"){
        if (flag_singlephoton_sel &&
          !map_cuts_flag["SPNCDeltaSig"] && !map_cuts_flag["SPOutFVSig"] &&
          !map_cuts_flag["SPNCPi0Sig"] && !map_cuts_flag["SPNCOtherSig"] &&
          !map_cuts_flag["SPNumuCCSig"]) return true;
      }
      return false;
  // NC Pi0 channels (with single photon selected events removed):
}else if (ch_name == "sp_nc_pi0_0p" || ch_name == "sp_nc_pi0_2_0p" || ch_name == "sp_nc_pi0_3_0p" || ch_name == "sp_nc_pi0_4_0p"
           || ch_name == "sp_nc_pi0_5_0p" || ch_name == "sp_nc_pi0_6_0p" || ch_name == "sp_nc_pi0_7_0p" || ch_name == "sp_nc_pi0_8_0p"
           || ch_name == "sp_nc_pi0_9_0p" || ch_name == "sp_nc_pi0_10_0p" || ch_name == "sp_nc_pi0_11_0p" || ch_name == "sp_nc_pi0_12_0p"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np" || ch_name == "sp_nc_pi0_2_Np" || ch_name == "sp_nc_pi0_3_Np" || ch_name == "sp_nc_pi0_4_Np"
            || ch_name == "sp_nc_pi0_5_Np" || ch_name == "sp_nc_pi0_6_Np" || ch_name == "sp_nc_pi0_7_Np" || ch_name == "sp_nc_pi0_8_Np"
            || ch_name == "sp_nc_pi0_9_Np" || ch_name == "sp_nc_pi0_10_Np" || ch_name == "sp_nc_pi0_11_Np" || ch_name == "sp_nc_pi0_12_Np"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp" || ch_name == "sp_nc_pi0_2_Xp" || ch_name == "sp_nc_pi0_3_Xp" || ch_name == "sp_nc_pi0_4_Xp"
            || ch_name == "sp_nc_pi0_5_Xp" || ch_name == "sp_nc_pi0_6_Xp" || ch_name == "sp_nc_pi0_7_Xp" || ch_name == "sp_nc_pi0_8_Xp"
            || ch_name == "sp_nc_pi0_9_Xp" || ch_name == "sp_nc_pi0_10_Xp" || ch_name == "sp_nc_pi0_11_Xp" || ch_name == "sp_nc_pi0_12_Xp"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_nsbeam" || ch_name == "sp_nc_pi0_2_0p_nsbeam" || ch_name == "sp_nc_pi0_3_0p_nsbeam" || ch_name == "sp_nc_pi0_4_0p_nsbeam"
            || ch_name == "sp_nc_pi0_5_0p_nsbeam" || ch_name == "sp_nc_pi0_6_0p_nsbeam" || ch_name == "sp_nc_pi0_7_0p_nsbeam" || ch_name == "sp_nc_pi0_8_0p_nsbeam"
            || ch_name == "sp_nc_pi0_9_0p_nsbeam" || ch_name == "sp_nc_pi0_10_0p_nsbeam" || ch_name == "sp_nc_pi0_11_0p_nsbeam" || ch_name == "sp_nc_pi0_12_0p_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_nsbeam" || ch_name == "sp_nc_pi0_2_Np_nsbeam" || ch_name == "sp_nc_pi0_3_Np_nsbeam" || ch_name == "sp_nc_pi0_4_Np_nsbeam"
            || ch_name == "sp_nc_pi0_5_Np_nsbeam" || ch_name == "sp_nc_pi0_6_Np_nsbeam" || ch_name == "sp_nc_pi0_7_Np_nsbeam" || ch_name == "sp_nc_pi0_8_Np_nsbeam"
            || ch_name == "sp_nc_pi0_9_Np_nsbeam" || ch_name == "sp_nc_pi0_10_Np_nsbeam" || ch_name == "sp_nc_pi0_11_Np_nsbeam" || ch_name == "sp_nc_pi0_12_Np_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_nsbeam" || ch_name == "sp_nc_pi0_2_Xp_nsbeam" || ch_name == "sp_nc_pi0_3_Xp_nsbeam" || ch_name == "sp_nc_pi0_4_Xp_nsbeam"
            || ch_name == "sp_nc_pi0_5_Xp_nsbeam" || ch_name == "sp_nc_pi0_6_Xp_nsbeam" || ch_name == "sp_nc_pi0_7_Xp_nsbeam" || ch_name == "sp_nc_pi0_8_Xp_nsbeam"
            || ch_name == "sp_nc_pi0_9_Xp_nsbeam" || ch_name == "sp_nc_pi0_10_Xp_nsbeam" || ch_name == "sp_nc_pi0_11_Xp_nsbeam" || ch_name == "sp_nc_pi0_12_Xp_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_ext" || ch_name == "sp_nc_pi0_2_0p_ext" || ch_name == "sp_nc_pi0_3_0p_ext" || ch_name == "sp_nc_pi0_4_0p_ext"
            || ch_name == "sp_nc_pi0_5_0p_ext" || ch_name == "sp_nc_pi0_6_0p_ext" || ch_name == "sp_nc_pi0_7_0p_ext" || ch_name == "sp_nc_pi0_8_0p_ext"
            || ch_name == "sp_nc_pi0_9_0p_ext" || ch_name == "sp_nc_pi0_10_0p_ext" || ch_name == "sp_nc_pi0_11_0p_ext" || ch_name == "sp_nc_pi0_12_0p_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_ext" || ch_name == "sp_nc_pi0_2_Np_ext" || ch_name == "sp_nc_pi0_3_Np_ext" || ch_name == "sp_nc_pi0_4_Np_ext"
            || ch_name == "sp_nc_pi0_5_Np_ext" || ch_name == "sp_nc_pi0_6_Np_ext" || ch_name == "sp_nc_pi0_7_Np_ext" || ch_name == "sp_nc_pi0_8_Np_ext"
            || ch_name == "sp_nc_pi0_9_Np_ext" || ch_name == "sp_nc_pi0_10_Np_ext" || ch_name == "sp_nc_pi0_11_Np_ext" || ch_name == "sp_nc_pi0_12_Np_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_ext" || ch_name == "sp_nc_pi0_2_Xp_ext" || ch_name == "sp_nc_pi0_3_Xp_ext" || ch_name == "sp_nc_pi0_4_Xp_ext"
            || ch_name == "sp_nc_pi0_5_Xp_ext" || ch_name == "sp_nc_pi0_6_Xp_ext" || ch_name == "sp_nc_pi0_7_Xp_ext" || ch_name == "sp_nc_pi0_8_Xp_ext"
            || ch_name == "sp_nc_pi0_9_Xp_ext" || ch_name == "sp_nc_pi0_10_Xp_ext" || ch_name == "sp_nc_pi0_11_Xp_ext" || ch_name == "sp_nc_pi0_12_Xp_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_dirt" || ch_name == "sp_nc_pi0_2_0p_dirt" || ch_name == "sp_nc_pi0_3_0p_dirt" || ch_name == "sp_nc_pi0_4_0p_dirt"
            || ch_name == "sp_nc_pi0_5_0p_dirt" || ch_name == "sp_nc_pi0_6_0p_dirt" || ch_name == "sp_nc_pi0_7_0p_dirt" || ch_name == "sp_nc_pi0_8_0p_dirt"
            || ch_name == "sp_nc_pi0_9_0p_dirt" || ch_name == "sp_nc_pi0_10_0p_dirt" || ch_name == "sp_nc_pi0_11_0p_dirt" || ch_name == "sp_nc_pi0_12_0p_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_dirt" || ch_name == "sp_nc_pi0_2_Np_dirt" || ch_name == "sp_nc_pi0_3_Np_dirt" || ch_name == "sp_nc_pi0_4_Np_dirt"
            || ch_name == "sp_nc_pi0_5_Np_dirt" || ch_name == "sp_nc_pi0_6_Np_dirt" || ch_name == "sp_nc_pi0_7_Np_dirt" || ch_name == "sp_nc_pi0_8_Np_dirt"
            || ch_name == "sp_nc_pi0_9_Np_dirt" || ch_name == "sp_nc_pi0_10_Np_dirt" || ch_name == "sp_nc_pi0_11_Np_dirt" || ch_name == "sp_nc_pi0_12_Np_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_dirt" || ch_name == "sp_nc_pi0_2_Xp_dirt" || ch_name == "sp_nc_pi0_3_Xp_dirt" || ch_name == "sp_nc_pi0_4_Xp_dirt"
            || ch_name == "sp_nc_pi0_5_Xp_dirt" || ch_name == "sp_nc_pi0_6_Xp_dirt" || ch_name == "sp_nc_pi0_7_Xp_dirt" || ch_name == "sp_nc_pi0_8_Xp_dirt"
            || ch_name == "sp_nc_pi0_9_Xp_dirt" || ch_name == "sp_nc_pi0_10_Xp_dirt" || ch_name == "sp_nc_pi0_11_Xp_dirt" || ch_name == "sp_nc_pi0_12_Xp_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_nc_delta_overlay" || ch_name == "sp_nc_pi0_0p_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_0p_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_nc_delta_overlay" || ch_name == "sp_nc_pi0_Np_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_Np_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_nc_delta_overlay" || ch_name == "sp_nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_Xp_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_overlay" || ch_name == "sp_nc_pi0_2_Xp_overlay" || ch_name == "sp_nc_pi0_3_Xp_overlay" || ch_name == "sp_nc_pi0_4_Xp_overlay"
            || ch_name == "sp_nc_pi0_5_Xp_overlay" || ch_name == "sp_nc_pi0_6_Xp_overlay" || ch_name == "sp_nc_pi0_7_Xp_overlay" || ch_name == "sp_nc_pi0_8_Xp_overlay"
            || ch_name == "sp_nc_pi0_9_Xp_overlay" || ch_name == "sp_nc_pi0_10_Xp_overlay" || ch_name == "sp_nc_pi0_11_Xp_overlay" || ch_name == "sp_nc_pi0_12_Xp_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel))
                      return true;
                  return false;
  //nc pi0 sp sideband
}else if (ch_name == "sp_bdt_nc_pi0_0p" || ch_name == "sp_bdt_nc_pi0_2_0p" || ch_name == "sp_bdt_nc_pi0_3_0p" || ch_name == "sp_bdt_nc_pi0_4_0p"
          || ch_name == "sp_bdt_nc_pi0_5_0p" || ch_name == "sp_bdt_nc_pi0_6_0p" || ch_name == "sp_bdt_nc_pi0_7_0p" || ch_name == "sp_bdt_nc_pi0_8_0p"
          || ch_name == "sp_bdt_nc_pi0_9_0p" || ch_name == "sp_bdt_nc_pi0_10_0p" || ch_name == "sp_bdt_nc_pi0_11_0p" || ch_name == "sp_bdt_nc_pi0_12_0p"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np" || ch_name == "sp_bdt_nc_pi0_2_Np" || ch_name == "sp_bdt_nc_pi0_3_Np" || ch_name == "sp_bdt_nc_pi0_4_Np"
            || ch_name == "sp_bdt_nc_pi0_5_Np" || ch_name == "sp_bdt_nc_pi0_6_Np" || ch_name == "sp_bdt_nc_pi0_7_Np" || ch_name == "sp_bdt_nc_pi0_8_Np"
            || ch_name == "sp_bdt_nc_pi0_9_Np" || ch_name == "sp_bdt_nc_pi0_10_Np" || ch_name == "sp_bdt_nc_pi0_11_Np" || ch_name == "sp_bdt_nc_pi0_12_Np"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp" || ch_name == "sp_bdt_nc_pi0_2_Xp" || ch_name == "sp_bdt_nc_pi0_3_Xp" || ch_name == "sp_bdt_nc_pi0_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_5_Xp" || ch_name == "sp_bdt_nc_pi0_6_Xp" || ch_name == "sp_bdt_nc_pi0_7_Xp" || ch_name == "sp_bdt_nc_pi0_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_9_Xp" || ch_name == "sp_bdt_nc_pi0_10_Xp" || ch_name == "sp_bdt_nc_pi0_11_Xp" || ch_name == "sp_bdt_nc_pi0_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_ext" || ch_name == "sp_bdt_nc_pi0_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_ext" || ch_name == "sp_bdt_nc_pi0_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_ext" || ch_name == "sp_bdt_nc_pi0_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_dirt" || ch_name == "sp_bdt_nc_pi0_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_dirt" || ch_name == "sp_bdt_nc_pi0_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  //nc pi0 sp sideband with one shower
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p"
          || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p"
          || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  //nc pi0 sp sideband with not one shower
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p"
          || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p"
          || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  // numuCC channels (with single photon and NC pi0 events removed):
}else if (ch_name == "sp_numuCC_0p" || ch_name == "sp_numuCC_2_0p" || ch_name == "sp_numuCC_3_0p" || ch_name == "sp_numuCC_4_0p"
          || ch_name == "sp_numuCC_5_0p" || ch_name == "sp_numuCC_6_0p" || ch_name == "sp_numuCC_7_0p" || ch_name == "sp_numuCC_8_0p"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np" || ch_name == "sp_numuCC_2_Np"|| ch_name == "sp_numuCC_3_Np" || ch_name == "sp_numuCC_4_Np"
            || ch_name == "sp_numuCC_5_Np" || ch_name == "sp_numuCC_6_Np" || ch_name == "sp_numuCC_7_Np" || ch_name == "sp_numuCC_8_Np"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp" || ch_name == "sp_numuCC_2_Xp"|| ch_name == "sp_numuCC_3_Xp" || ch_name == "sp_numuCC_4_Xp"
            || ch_name == "sp_numuCC_5_Xp" || ch_name == "sp_numuCC_6_Xp" || ch_name == "sp_numuCC_7_Xp" || ch_name == "sp_numuCC_8_Xp"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_nsbeam" || ch_name == "sp_numuCC_2_0p_nsbeam" || ch_name == "sp_numuCC_3_0p_nsbeam" || ch_name == "sp_numuCC_4_0p_nsbeam"
             || ch_name == "sp_numuCC_5_0p_nsbeam" || ch_name == "sp_numuCC_6_0p_nsbeam" || ch_name == "sp_numuCC_7_0p_nsbeam" || ch_name == "sp_numuCC_8_0p_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_nsbeam" || ch_name == "sp_numuCC_2_Np_nsbeam" || ch_name == "sp_numuCC_3_Np_nsbeam" || ch_name == "sp_numuCC_4_Np_nsbeam"
             || ch_name == "sp_numuCC_5_Np_nsbeam" || ch_name == "sp_numuCC_6_Np_nsbeam" || ch_name == "sp_numuCC_7_Np_nsbeam" || ch_name == "sp_numuCC_8_Np_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_nsbeam" || ch_name == "sp_numuCC_2_Xp_nsbeam" || ch_name == "sp_numuCC_3_Xp_nsbeam" || ch_name == "sp_numuCC_4_Xp_nsbeam"
             || ch_name == "sp_numuCC_5_Xp_nsbeam" || ch_name == "sp_numuCC_6_Xp_nsbeam" || ch_name == "sp_numuCC_7_Xp_nsbeam" || ch_name == "sp_numuCC_8_Xp_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_ext" || ch_name == "sp_numuCC_2_0p_ext" || ch_name == "sp_numuCC_3_0p_ext" || ch_name == "sp_numuCC_4_0p_ext"
             || ch_name == "sp_numuCC_5_0p_ext" || ch_name == "sp_numuCC_6_0p_ext" || ch_name == "sp_numuCC_7_0p_ext" || ch_name == "sp_numuCC_8_0p_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_ext" || ch_name == "sp_numuCC_2_Np_ext" || ch_name == "sp_numuCC_3_Np_ext" || ch_name == "sp_numuCC_4_Np_ext"
             || ch_name == "sp_numuCC_5_Np_ext" || ch_name == "sp_numuCC_6_Np_ext" || ch_name == "sp_numuCC_7_Np_ext" || ch_name == "sp_numuCC_8_Np_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_ext" || ch_name == "sp_numuCC_2_Xp_ext" || ch_name == "sp_numuCC_3_Xp_ext" || ch_name == "sp_numuCC_4_Xp_ext"
             || ch_name == "sp_numuCC_5_Xp_ext" || ch_name == "sp_numuCC_6_Xp_ext" || ch_name == "sp_numuCC_7_Xp_ext" || ch_name == "sp_numuCC_8_Xp_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_dirt" || ch_name == "sp_numuCC_2_0p_dirt" || ch_name == "sp_numuCC_3_0p_dirt" || ch_name == "sp_numuCC_4_0p_dirt"
             || ch_name == "sp_numuCC_5_0p_dirt" || ch_name == "sp_numuCC_6_0p_dirt" || ch_name == "sp_numuCC_7_0p_dirt" || ch_name == "sp_numuCC_8_0p_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_dirt" || ch_name == "sp_numuCC_2_Np_dirt" || ch_name == "sp_numuCC_3_Np_dirt" || ch_name == "sp_numuCC_4_Np_dirt"
             || ch_name == "sp_numuCC_5_Np_dirt" || ch_name == "sp_numuCC_6_Np_dirt" || ch_name == "sp_numuCC_7_Np_dirt" || ch_name == "sp_numuCC_8_Np_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_dirt" || ch_name == "sp_numuCC_2_Xp_dirt" || ch_name == "sp_numuCC_3_Xp_dirt" || ch_name == "sp_numuCC_4_Xp_dirt"
             || ch_name == "sp_numuCC_5_Xp_dirt" || ch_name == "sp_numuCC_6_Xp_dirt" || ch_name == "sp_numuCC_7_Xp_dirt" || ch_name == "sp_numuCC_8_Xp_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_nc_delta_overlay" || ch_name == "sp_numuCC_0p_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_0p_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_nc_delta_overlay" || ch_name == "sp_numuCC_Np_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_Np_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_nc_delta_overlay" || ch_name == "sp_numuCC_Xp_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_Xp_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_overlay" || ch_name == "sp_numuCC_2_Xp_overlay"){
                if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p" || ch_name == "sp_bdt_numuCC_2_0p" || ch_name == "sp_bdt_numuCC_3_0p" || ch_name == "sp_bdt_numuCC_4_0p"
             || ch_name == "sp_bdt_numuCC_5_0p" || ch_name == "sp_bdt_numuCC_6_0p" || ch_name == "sp_bdt_numuCC_7_0p" || ch_name == "sp_bdt_numuCC_8_0p"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np" || ch_name == "sp_bdt_numuCC_2_Np" || ch_name == "sp_bdt_numuCC_3_Np" || ch_name == "sp_bdt_numuCC_4_Np"
             || ch_name == "sp_bdt_numuCC_5_Np" || ch_name == "sp_bdt_numuCC_6_Np" || ch_name == "sp_bdt_numuCC_7_Np" || ch_name == "sp_bdt_numuCC_8_Np"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp" || ch_name == "sp_bdt_numuCC_2_Xp" || ch_name == "sp_bdt_numuCC_3_Xp" || ch_name == "sp_bdt_numuCC_4_Xp"
             || ch_name == "sp_bdt_numuCC_5_Xp" || ch_name == "sp_bdt_numuCC_6_Xp" || ch_name == "sp_bdt_numuCC_7_Xp" || ch_name == "sp_bdt_numuCC_8_Xp"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_ext" || ch_name == "sp_bdt_numuCC_2_0p_ext" || ch_name == "sp_bdt_numuCC_3_0p_ext" || ch_name == "sp_bdt_numuCC_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_5_0p_ext" || ch_name == "sp_bdt_numuCC_6_0p_ext" || ch_name == "sp_bdt_numuCC_7_0p_ext" || ch_name == "sp_bdt_numuCC_8_0p_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_ext" || ch_name == "sp_bdt_numuCC_2_Np_ext" || ch_name == "sp_bdt_numuCC_3_Np_ext" || ch_name == "sp_bdt_numuCC_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_5_Np_ext" || ch_name == "sp_bdt_numuCC_6_Np_ext" || ch_name == "sp_bdt_numuCC_7_Np_ext" || ch_name == "sp_bdt_numuCC_8_Np_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_ext" || ch_name == "sp_bdt_numuCC_2_Xp_ext" || ch_name == "sp_bdt_numuCC_3_Xp_ext" || ch_name == "sp_bdt_numuCC_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_5_Xp_ext" || ch_name == "sp_bdt_numuCC_6_Xp_ext" || ch_name == "sp_bdt_numuCC_7_Xp_ext" || ch_name == "sp_bdt_numuCC_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_dirt" || ch_name == "sp_bdt_numuCC_2_0p_dirt" || ch_name == "sp_bdt_numuCC_3_0p_dirt" || ch_name == "sp_bdt_numuCC_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_5_0p_dirt" || ch_name == "sp_bdt_numuCC_6_0p_dirt" || ch_name == "sp_bdt_numuCC_7_0p_dirt" || ch_name == "sp_bdt_numuCC_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_dirt" || ch_name == "sp_bdt_numuCC_2_Np_dirt" || ch_name == "sp_bdt_numuCC_3_Np_dirt" || ch_name == "sp_bdt_numuCC_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_5_Np_dirt" || ch_name == "sp_bdt_numuCC_6_Np_dirt" || ch_name == "sp_bdt_numuCC_7_Np_dirt" || ch_name == "sp_bdt_numuCC_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_dirt" || ch_name == "sp_bdt_numuCC_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay" || ch_name == "sp_bdt_numuCC_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay" || ch_name == "sp_bdt_numuCC_2_0p_overlay" || ch_name == "sp_bdt_numuCC_3_0p_overlay" || ch_name == "sp_bdt_numuCC_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay" || ch_name == "sp_bdt_numuCC_6_0p_overlay" || ch_name == "sp_bdt_numuCC_7_0p_overlay" || ch_name == "sp_bdt_numuCC_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay" || ch_name == "sp_bdt_numuCC_2_Np_overlay" || ch_name == "sp_bdt_numuCC_3_Np_overlay" || ch_name == "sp_bdt_numuCC_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay" || ch_name == "sp_bdt_numuCC_6_Np_overlay" || ch_name == "sp_bdt_numuCC_7_Np_overlay" || ch_name == "sp_bdt_numuCC_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_spoverlay" || ch_name == "sp_bdt_numuCC_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_spoverlay" || ch_name == "sp_bdt_numuCC_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  //one shower numu sideband
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p" || ch_name == "sp_bdt_numuCC_oneshw_2_0p" || ch_name == "sp_bdt_numuCC_oneshw_3_0p" || ch_name == "sp_bdt_numuCC_oneshw_4_0p"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p" || ch_name == "sp_bdt_numuCC_oneshw_6_0p" || ch_name == "sp_bdt_numuCC_oneshw_7_0p" || ch_name == "sp_bdt_numuCC_oneshw_8_0p"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np" || ch_name == "sp_bdt_numuCC_oneshw_2_Np" || ch_name == "sp_bdt_numuCC_oneshw_3_Np" || ch_name == "sp_bdt_numuCC_oneshw_4_Np"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np" || ch_name == "sp_bdt_numuCC_oneshw_6_Np" || ch_name == "sp_bdt_numuCC_oneshw_7_Np" || ch_name == "sp_bdt_numuCC_oneshw_8_Np"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

//not one shower numu sideband
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

  //single photon other sideband
}else if (ch_name == "sp_bdt_other_0p" || ch_name == "sp_bdt_other_2_0p" || ch_name == "sp_bdt_other_3_0p" || ch_name == "sp_bdt_other_4_0p"
           || ch_name == "sp_bdt_other_5_0p" || ch_name == "sp_bdt_other_6_0p" || ch_name == "sp_bdt_other_7_0p" || ch_name == "sp_bdt_other_8_0p"
          || ch_name == "sp_bdt_other_9_0p" || ch_name == "sp_bdt_other_10_0p" || ch_name == "sp_bdt_other_11_0p" || ch_name == "sp_bdt_other_12_0p"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np" || ch_name == "sp_bdt_other_2_Np" || ch_name == "sp_bdt_other_3_Np" || ch_name == "sp_bdt_other_4_Np"
           || ch_name == "sp_bdt_other_5_Np" || ch_name == "sp_bdt_other_6_Np" || ch_name == "sp_bdt_other_7_Np" || ch_name == "sp_bdt_other_8_Np"
          || ch_name == "sp_bdt_other_9_Np" || ch_name == "sp_bdt_other_10_Np" || ch_name == "sp_bdt_other_11_Np" || ch_name == "sp_bdt_other_12_Np"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp" || ch_name == "sp_bdt_other_2_Xp" || ch_name == "sp_bdt_other_3_Xp" || ch_name == "sp_bdt_other_4_Xp"
             || ch_name == "sp_bdt_other_5_Xp" || ch_name == "sp_bdt_other_6_Xp" || ch_name == "sp_bdt_other_7_Xp" || ch_name == "sp_bdt_other_8_Xp"
            || ch_name == "sp_bdt_other_9_Xp" || ch_name == "sp_bdt_other_10_Xp" || ch_name == "sp_bdt_other_11_Xp" || ch_name == "sp_bdt_other_12_Xp"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_nsbeam" || ch_name == "sp_bdt_other_2_0p_nsbeam" || ch_name == "sp_bdt_other_3_0p_nsbeam" || ch_name == "sp_bdt_other_4_0p_nsbeam"
             || ch_name == "sp_bdt_other_5_0p_nsbeam" || ch_name == "sp_bdt_other_6_0p_nsbeam" || ch_name == "sp_bdt_other_7_0p_nsbeam" || ch_name == "sp_bdt_other_8_0p_nsbeam"
             || ch_name == "sp_bdt_other_9_0p_nsbeam" || ch_name == "sp_bdt_other_10_0p_nsbeam" || ch_name == "sp_bdt_other_11_0p_nsbeam" || ch_name == "sp_bdt_other_12_0p_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_nsbeam" || ch_name == "sp_bdt_other_2_Np_nsbeam" || ch_name == "sp_bdt_other_3_Np_nsbeam" || ch_name == "sp_bdt_other_4_Np_nsbeam"
             || ch_name == "sp_bdt_other_5_Np_nsbeam" || ch_name == "sp_bdt_other_6_Np_nsbeam" || ch_name == "sp_bdt_other_7_Np_nsbeam" || ch_name == "sp_bdt_other_8_Np_nsbeam"
             || ch_name == "sp_bdt_other_9_Np_nsbeam" || ch_name == "sp_bdt_other_10_Np_nsbeam" || ch_name == "sp_bdt_other_11_Np_nsbeam" || ch_name == "sp_bdt_other_12_Np_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_nsbeam" || ch_name == "sp_bdt_other_2_Xp_nsbeam" || ch_name == "sp_bdt_other_3_Xp_nsbeam" || ch_name == "sp_bdt_other_4_Xp_nsbeam"
           || ch_name == "sp_bdt_other_5_Xp_nsbeam" || ch_name == "sp_bdt_other_6_Xp_nsbeam" || ch_name == "sp_bdt_other_7_Xp_nsbeam" || ch_name == "sp_bdt_other_8_Xp_nsbeam"
             || ch_name == "sp_bdt_other_9_Xp_nsbeam" || ch_name == "sp_bdt_other_10_Xp_nsbeam" || ch_name == "sp_bdt_other_11_Xp_nsbeam" || ch_name == "sp_bdt_other_12_Xp_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_ext" || ch_name == "sp_bdt_other_2_0p_ext" || ch_name == "sp_bdt_other_3_0p_ext" || ch_name == "sp_bdt_other_4_0p_ext"
             || ch_name == "sp_bdt_other_5_0p_ext" || ch_name == "sp_bdt_other_6_0p_ext" || ch_name == "sp_bdt_other_7_0p_ext" || ch_name == "sp_bdt_other_8_0p_ext"
               || ch_name == "sp_bdt_other_9_0p_ext" || ch_name == "sp_bdt_other_10_0p_ext" || ch_name == "sp_bdt_other_11_0p_ext" || ch_name == "sp_bdt_other_12_0p_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_ext" || ch_name == "sp_bdt_other_2_Np_ext" || ch_name == "sp_bdt_other_3_Np_ext" || ch_name == "sp_bdt_other_4_Np_ext"
             || ch_name == "sp_bdt_other_5_Np_ext" || ch_name == "sp_bdt_other_6_Np_ext" || ch_name == "sp_bdt_other_7_Np_ext" || ch_name == "sp_bdt_other_8_Np_ext"
             || ch_name == "sp_bdt_other_9_Np_ext" || ch_name == "sp_bdt_other_10_Np_ext" || ch_name == "sp_bdt_other_11_Np_ext" || ch_name == "sp_bdt_other_12_Np_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_ext" || ch_name == "sp_bdt_other_2_Xp_ext" || ch_name == "sp_bdt_other_3_Xp_ext" || ch_name == "sp_bdt_other_4_Xp_ext"
             || ch_name == "sp_bdt_other_5_Xp_ext" || ch_name == "sp_bdt_other_6_Xp_ext" || ch_name == "sp_bdt_other_7_Xp_ext" || ch_name == "sp_bdt_other_8_Xp_ext"
            || ch_name == "sp_bdt_other_9_Xp_ext" || ch_name == "sp_bdt_other_10_Xp_ext" || ch_name == "sp_bdt_other_11_Xp_ext" || ch_name == "sp_bdt_other_12_Xp_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_dirt" || ch_name == "sp_bdt_other_2_0p_dirt" || ch_name == "sp_bdt_other_3_0p_dirt" || ch_name == "sp_bdt_other_4_0p_dirt"
             || ch_name == "sp_bdt_other_5_0p_dirt" || ch_name == "sp_bdt_other_6_0p_dirt" || ch_name == "sp_bdt_other_7_0p_dirt" || ch_name == "sp_bdt_other_8_0p_dirt"
             || ch_name == "sp_bdt_other_9_0p_dirt" || ch_name == "sp_bdt_other_10_0p_dirt" || ch_name == "sp_bdt_other_11_0p_dirt" || ch_name == "sp_bdt_other_12_0p_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_dirt" || ch_name == "sp_bdt_other_2_Np_dirt" || ch_name == "sp_bdt_other_3_Np_dirt" || ch_name == "sp_bdt_other_4_Np_dirt"
             || ch_name == "sp_bdt_other_5_Np_dirt" || ch_name == "sp_bdt_other_6_Np_dirt" || ch_name == "sp_bdt_other_7_Np_dirt" || ch_name == "sp_bdt_other_8_Np_dirt"
            || ch_name == "sp_bdt_other_9_Np_dirt" || ch_name == "sp_bdt_other_10_Np_dirt" || ch_name == "sp_bdt_other_11_Np_dirt" || ch_name == "sp_bdt_other_12_Np_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_dirt" || ch_name == "sp_bdt_other_2_Xp_dirt" || ch_name == "sp_bdt_other_3_Xp_dirt" || ch_name == "sp_bdt_other_4_Xp_dirt"
             || ch_name == "sp_bdt_other_5_Xp_dirt" || ch_name == "sp_bdt_other_6_Xp_dirt" || ch_name == "sp_bdt_other_7_Xp_dirt" || ch_name == "sp_bdt_other_8_Xp_dirt"
            || ch_name == "sp_bdt_other_9_Xp_dirt" || ch_name == "sp_bdt_other_10_Xp_dirt" || ch_name == "sp_bdt_other_11_Xp_dirt" || ch_name == "sp_bdt_other_12_Xp_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_nc_delta_overlay" || ch_name == "sp_bdt_other_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_nc_delta_overlay" || ch_name == "sp_bdt_other_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_nc_delta_overlay" || ch_name == "sp_bdt_other_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay" || ch_name == "sp_bdt_other_2_Xp_overlay" || ch_name == "sp_bdt_other_3_Xp_overlay" || ch_name == "sp_bdt_other_4_Xp_overlay"
             || ch_name == "sp_bdt_other_5_Xp_overlay" || ch_name == "sp_bdt_other_6_Xp_overlay" || ch_name == "sp_bdt_other_7_Xp_overlay" || ch_name == "sp_bdt_other_8_Xp_overlay"
             || ch_name == "sp_bdt_other_9_Xp_overlay" || ch_name == "sp_bdt_other_10_Xp_overlay" || ch_name == "sp_bdt_other_11_Xp_overlay" || ch_name == "sp_bdt_other_2_Xp_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay" || ch_name == "sp_bdt_other_2_0p_overlay" || ch_name == "sp_bdt_other_3_0p_overlay" || ch_name == "sp_bdt_other_4_0p_overlay"
             || ch_name == "sp_bdt_other_5_0p_overlay" || ch_name == "sp_bdt_other_6_0p_overlay" || ch_name == "sp_bdt_other_7_0p_overlay" || ch_name == "sp_bdt_other_8_0p_overlay"
            || ch_name == "sp_bdt_other_9_0p_overlay"|| ch_name == "sp_bdt_other_10_0p_overlay" || ch_name == "sp_bdt_other_11_0p_overlay" || ch_name == "sp_bdt_other_12_0p_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay" || ch_name == "sp_bdt_other_2_Np_overlay" || ch_name == "sp_bdt_other_3_Np_overlay" || ch_name == "sp_bdt_other_4_Np_overlay"
           || ch_name == "sp_bdt_other_5_Np_overlay" || ch_name == "sp_bdt_other_6_Np_overlay" || ch_name == "sp_bdt_other_7_Np_overlay" || ch_name == "sp_bdt_other_8_Np_overlay"
           || ch_name == "sp_bdt_other_9_Np_overlay" || ch_name == "sp_bdt_other_10_Np_overlay" || ch_name == "sp_bdt_other_11_Np_overlay" || ch_name == "sp_bdt_other_12_Np_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay_BG" || ch_name == "sp_bdt_other_2_Xp_overlay_BG" || ch_name == "sp_bdt_other_3_Xp_overlay_BG" || ch_name == "sp_bdt_other_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_other_5_Xp_overlay_BG" || ch_name == "sp_bdt_other_6_Xp_overlay_BG" || ch_name == "sp_bdt_other_7_Xp_overlay_BG" || ch_name == "sp_bdt_other_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_other_9_Xp_overlay_BG" || ch_name == "sp_bdt_other_10_Xp_overlay_BG" || ch_name == "sp_bdt_other_11_Xp_overlay_BG" || ch_name == "sp_bdt_other_12_Xp_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay_BG" || ch_name == "sp_bdt_other_2_0p_overlay_BG" || ch_name == "sp_bdt_other_3_0p_overlay_BG" || ch_name == "sp_bdt_other_4_0p_overlay_BG"
             || ch_name == "sp_bdt_other_5_0p_overlay_BG" || ch_name == "sp_bdt_other_6_0p_overlay_BG" || ch_name == "sp_bdt_other_7_0p_overlay_BG" || ch_name == "sp_bdt_other_8_0p_overlay_BG"
             || ch_name == "sp_bdt_other_9_0p_overlay_BG" || ch_name == "sp_bdt_other_10_0p_overlay_BG" || ch_name == "sp_bdt_other_11_0p_overlay_BG" || ch_name == "sp_bdt_other_12_0p_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay_BG" || ch_name == "sp_bdt_other_2_Np_overlay_BG" || ch_name == "sp_bdt_other_3_Np_overlay_BG" || ch_name == "sp_bdt_other_4_Np_overlay_BG"
             || ch_name == "sp_bdt_other_5_Np_overlay_BG" || ch_name == "sp_bdt_other_6_Np_overlay_BG" || ch_name == "sp_bdt_other_7_Np_overlay_BG" || ch_name == "sp_bdt_other_8_Np_overlay_BG"
             || ch_name == "sp_bdt_other_9_Np_overlay_BG" || ch_name == "sp_bdt_other_10_Np_overlay_BG" || ch_name == "sp_bdt_other_11_Np_overlay_BG" || ch_name == "sp_bdt_other_12_Np_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_spoverlay" || ch_name == "sp_bdt_other_2_Xp_spoverlay" || ch_name == "sp_bdt_other_3_Xp_spoverlay" || ch_name == "sp_bdt_other_4_Xp_spoverlay"
             || ch_name == "sp_bdt_other_5_Xp_spoverlay" || ch_name == "sp_bdt_other_6_Xp_spoverlay" || ch_name == "sp_bdt_other_7_Xp_spoverlay" || ch_name == "sp_bdt_other_8_Xp_spoverlay"
             || ch_name == "sp_bdt_other_9_Xp_spoverlay" || ch_name == "sp_bdt_other_10_Xp_spoverlay" || ch_name == "sp_bdt_other_11_Xp_spoverlay" || ch_name == "sp_bdt_other_12_Xp_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_spoverlay" || ch_name == "sp_bdt_other_2_0p_spoverlay" || ch_name == "sp_bdt_other_3_0p_spoverlay" || ch_name == "sp_bdt_other_4_0p_spoverlay"
            || ch_name == "sp_bdt_other_5_0p_spoverlay" || ch_name == "sp_bdt_other_6_0p_spoverlay" || ch_name == "sp_bdt_other_7_0p_spoverlay" || ch_name == "sp_bdt_other_8_0p_spoverlay"
            || ch_name == "sp_bdt_other_9_0p_spoverlay" || ch_name == "sp_bdt_other_10_0p_spoverlay" || ch_name == "sp_bdt_other_11_0p_spoverlay" || ch_name == "sp_bdt_other_12_0p_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_spoverlay" || ch_name == "sp_bdt_other_2_Np_spoverlay" || ch_name == "sp_bdt_other_3_Np_spoverlay" || ch_name == "sp_bdt_other_4_Np_spoverlay"
            || ch_name == "sp_bdt_other_5_Np_spoverlay" || ch_name == "sp_bdt_other_6_Np_spoverlay" || ch_name == "sp_bdt_other_7_Np_spoverlay" || ch_name == "sp_bdt_other_8_Np_spoverlay"
            || ch_name == "sp_bdt_other_9_Np_spoverlay" || ch_name == "sp_bdt_other_10_Np_spoverlay" || ch_name == "sp_bdt_other_11_Np_spoverlay" || ch_name == "sp_bdt_other_12_Np_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_12_Xp_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                    !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                    !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_ncpi0overlay" || ch_name == "sp_bdt_other_2_0p_ncpi0overlay" || ch_name == "sp_bdt_other_3_0p_ncpi0overlay" || ch_name == "sp_bdt_other_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_0p_ncpi0overlay" || ch_name == "sp_bdt_other_6_0p_ncpi0overlay" || ch_name == "sp_bdt_other_7_0p_ncpi0overlay" || ch_name == "sp_bdt_other_8_0p_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_0p_ncpi0overlay" || ch_name == "sp_bdt_other_10_0p_ncpi0overlay" || ch_name == "sp_bdt_other_11_0p_ncpi0overlay" || ch_name == "sp_bdt_other_12_0p_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_ncpi0overlay" || ch_name == "sp_bdt_other_2_Np_ncpi0overlay" || ch_name == "sp_bdt_other_3_Np_ncpi0overlay" || ch_name == "sp_bdt_other_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_Np_ncpi0overlay" || ch_name == "sp_bdt_other_6_Np_ncpi0overlay" || ch_name == "sp_bdt_other_7_Np_ncpi0overlay" || ch_name == "sp_bdt_other_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_Np_ncpi0overlay" || ch_name == "sp_bdt_other_10_Np_ncpi0overlay" || ch_name == "sp_bdt_other_11_Np_ncpi0overlay" || ch_name == "sp_bdt_other_12_Np_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  //single photon nue cc sideband
  }else if (ch_name == "sp_bdt_nue_0p" || ch_name == "sp_bdt_nue_2_0p"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np" || ch_name == "sp_bdt_nue_2_Np"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp" || ch_name == "sp_bdt_nue_2_Xp"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_nsbeam" || ch_name == "sp_bdt_nue_2_0p_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_nsbeam" || ch_name == "sp_bdt_nue_2_Np_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_nsbeam" || ch_name == "sp_bdt_nue_2_Xp_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_ext" || ch_name == "sp_bdt_nue_2_0p_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_ext" || ch_name == "sp_bdt_nue_2_Np_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_ext" || ch_name == "sp_bdt_nue_2_Xp_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_dirt" || ch_name == "sp_bdt_nue_2_0p_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_dirt" || ch_name == "sp_bdt_nue_2_Np_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_dirt" || ch_name == "sp_bdt_nue_2_Xp_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_nc_delta_overlay" || ch_name == "sp_bdt_nue_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_nc_delta_overlay" || ch_name == "sp_bdt_nue_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nue_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay" || ch_name == "sp_bdt_nue_2_Xp_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay" || ch_name == "sp_bdt_nue_2_0p_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay" || ch_name == "sp_bdt_nue_2_Np_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay_BG" || ch_name == "sp_bdt_nue_2_Xp_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay_BG" || ch_name == "sp_bdt_nue_2_0p_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay_BG" || ch_name == "sp_bdt_nue_2_Np_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_spoverlay" || ch_name == "sp_bdt_nue_2_Xp_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_spoverlay" || ch_name == "sp_bdt_nue_2_0p_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_spoverlay" || ch_name == "sp_bdt_nue_2_Np_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;

  }else if (ch_name == "sp_bdt_nue_Xp_ncpi0overlay" || ch_name == "sp_bdt_nue_2_Xp_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_ncpi0overlay" || ch_name == "sp_bdt_nue_2_0p_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_ncpi0overlay" || ch_name == "sp_bdt_nue_2_Np_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  
  //single photon nue cc sideband - all showers
  }else if (ch_name == "sp_bdt_nue_allshw_0p" || ch_name == "sp_bdt_nue_allshw_2_0p"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np" || ch_name == "sp_bdt_nue_allshw_2_Np"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp" || ch_name == "sp_bdt_nue_allshw_2_Xp"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_0p_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_Np_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_Xp_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_ext" || ch_name == "sp_bdt_nue_allshw_2_0p_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_ext" || ch_name == "sp_bdt_nue_allshw_2_Np_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_ext" || ch_name == "sp_bdt_nue_allshw_2_Xp_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_dirt" || ch_name == "sp_bdt_nue_allshw_2_0p_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_dirt" || ch_name == "sp_bdt_nue_allshw_2_Np_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_dirt" || ch_name == "sp_bdt_nue_allshw_2_Xp_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_0p_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_Np_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;

  }else if (ch_name == "sp_bdt_nue_allshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_0p_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_Np_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

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




bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;

  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;

  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;

  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;

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


bool LEEana::is_NCpio_sel(TaggerInfo& tagger_info, KineInfo& kine){ // includes all cuts except FC
  bool flag = false;
  if (tagger_info.nc_pio_score > 1.816 && tagger_info.numu_cc_flag >=0 && kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.) flag = true;
  return flag;
}
bool LEEana::is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){ // includes all cuts except FC
  bool flag = false;
  if (tagger_info.nc_delta_score > 2.61 && tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}

//Erin
bool LEEana::is_singlephoton_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    tagger_info.single_photon_ncpi0_score > -0.05 && tagger_info.single_photon_nue_score > -1.0 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score > -3.0 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singleshower_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    tagger_info.single_photon_ncpi0_score > -0.05 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singleshower_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_numu_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    //tagger_info.single_photon_numu_score < 0.4 && //pure
    tagger_info.single_photon_numu_score < 0.1 && //eff
    //tagger_info.single_photon_numu_score < -3.0 && //numu pure
    tagger_info.single_photon_numu_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_other_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score < -0.4 &&
    tagger_info.single_photon_other_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_ncpi0_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    //tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    //tagger_info.single_photon_ncpi0_score < -0.05 && //pure
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score < -0.4 && //eff
    //tagger_info.single_photon_numu_score > -1.0 && tagger_info.single_photon_other_score > -0.4 &&
    //tagger_info.single_photon_ncpi0_score < -0.4 && //nc pi0 eff
    tagger_info.single_photon_ncpi0_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_nue_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score < -3.0 &&
    tagger_info.shw_sp_n_20br1_showers==1 &&
    tagger_info.single_photon_nue_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_nue_sel_allshw(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score < -3.0 &&
    tagger_info.single_photon_nue_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_nsbeam(PFevalInfo& pfeval, EvalInfo& eval){
  bool flag = false;

  double delta_time_calc = -9999.;
  //Merge Peaks
  double gap=18.936;
  double Shift=0;
  if (pfeval.run >= 17380){ Shift=2916.0; }
  else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
  else if (pfeval.run >= 10812){ Shift = 3568.5; }
  else if (pfeval.run >= 8321){ Shift = 3610.7;}
  else if (pfeval.run >= 5800){ Shift = 3164.4;}
  else if (pfeval.run >= 0){ Shift = 3168.9;}
  //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
  //if(run>8000 && run<10812){Shift=3610.7; }
  //if(run>=10812 && run <12500){Shift=3568.5; }
  double TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
  double TT_merged = -9999.;

  //merge peaks
  if(TThelp>=0 && TThelp<gap*81.0){
    TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
  }

  delta_time_calc = TT_merged;

  if (abs(delta_time_calc) < 5.0) {flag = true;}
  return flag;
}

bool LEEana::is_nsbeam_photon(PFevalInfo& pfeval, EvalInfo& eval){
  bool flag = false;

  double delta_time_calc = -9999.;
  //Merge Peaks
  double gap=18.936;
  double Shift=0;
  if (pfeval.run >= 17380){ Shift=2916.0; }
  else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
  else if (pfeval.run >= 10812){ Shift = 3568.5; }
  else if (pfeval.run >= 8321){ Shift = 3610.7;}
  else if (pfeval.run >= 5800){ Shift = 3164.4;}
  else if (pfeval.run >= 0){ Shift = 3168.9;}
  //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
  //if(run>8000 && run<10812){Shift=3610.7; }
  //if(run>=10812 && run <12500){Shift=3568.5; }
  double TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
  double TT_merged = -9999.;

  //merge peaks
  if(TThelp>=0 && TThelp<gap*81.0){
    TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
  }

  delta_time_calc = TT_merged;

  if (delta_time_calc > -6.6 && delta_time_calc < 3.4) {flag = true;}
  return flag;
}

//break down selection
bool LEEana::is_singlephoton_pre(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_numu_score > 0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_other(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_other_score > 0.2) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_ncpi0_score > -0.05 ) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_nue_score > -1.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_numu_score > 0.1) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_other(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_other_score > -0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_ncpi0_score > -0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_nue_score > -3.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_oneshw(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
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

bool LEEana::is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
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
    if(Nproton==1) flag = true;
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


bool LEEana::is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data){
    bool flag = false;

    if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){

      double reco_Enu = get_reco_Enu_corr(kine, flag_data);

      Float_t Ehadron = reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
      if(Ehadron<200) // MeV
        {
	  flag = true;
        }
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

int LEEana::mcc8_pmuon_costheta_bin(float pmuon, float costh){

    if (costh>=-1 and costh<-0.5) {
      if (pmuon>=0 and pmuon<0.18) return 1;
      else if (pmuon>=0.18 and pmuon<0.30) return 2;
      else if (pmuon>=0.30 and pmuon<0.45) return 3;
      else if (pmuon>=0.45 and pmuon<0.77) return 4;
      else if (pmuon>=0.77 and pmuon<2.5) return 5;
      else return -10000;
    }
    else if (costh>=-0.5 and costh<0){
      if (pmuon>=0 and pmuon<0.18) return 6;
      else if (pmuon>=0.18 and pmuon<0.30) return 7;
      else if (pmuon>=0.30 and pmuon<0.45) return 8;
      else if (pmuon>=0.45 and pmuon<0.77) return 9;
      else if (pmuon>=0.77 and pmuon<2.5) return 10;
      else return -10000;
    }
    else if (costh>0 and costh<0.27){
      if (pmuon>=0 and pmuon<0.18) return 11;
      else if (pmuon>=0.18 and pmuon<0.30) return 12;
      else if (pmuon>=0.30 and pmuon<0.45) return 13;
      else if (pmuon>=0.45 and pmuon<0.77) return 14;
      else if (pmuon>=0.77 and pmuon<2.5) return 15;
      else return -10000;
    }
    else if (costh>=0.27 and costh<0.45){
      if (pmuon>=0 and pmuon<0.30) return 16;
      else if (pmuon>=0.30 and pmuon<0.45) return 17;
      else if (pmuon>=0.45 and pmuon<0.77) return 18;
      else if (pmuon>=0.77 and pmuon<2.5) return 19;
      else return -10000;
    }
    else if (costh>=0.45 and costh<0.62){
      if (pmuon>=0 and pmuon<0.30) return 20;
      else if (pmuon>=0.30 and pmuon<0.45) return 21;
      else if (pmuon>=0.45 and pmuon<0.77) return 22;
      else if (pmuon>=0.77 and pmuon<2.5) return 23;
      else return -10000;
    }
    else if (costh>=0.62 and costh<0.76){
      if (pmuon>=0 and pmuon<0.30) return 24;
      else if (pmuon>=0.30 and pmuon<0.45) return 25;
      else if (pmuon>=0.45 and pmuon<0.77) return 26;
      else if (pmuon>=0.77 and pmuon<2.5) return 27;
      else return -10000;
    }
    else if (costh>=0.76 and costh<0.86){
      if (pmuon>=0 and pmuon<0.30) return 28;
      else if (pmuon>=0.30 and pmuon<0.45) return 29;
      else if (pmuon>=0.45 and pmuon<0.77) return 30;
      else if (pmuon>=0.77 and pmuon<1.28) return 31;
      else if (pmuon>=1.28 and pmuon<2.5) return 32;
      else return -10000;
    }
    else if (costh>=0.86 and costh<0.94){
      if (pmuon>=0 and pmuon<0.30) return 33;
      else if (pmuon>=0.30 and pmuon<0.45) return 34;
      else if (pmuon>=0.45 and pmuon<0.77) return 35;
      else if (pmuon>=0.77 and pmuon<1.28) return 36;
      else if (pmuon>=1.28 and pmuon<2.5) return 37;
      else return -10000;
    }
    else if (costh>=0.94 and costh<1.00){
      if (pmuon>=0 and pmuon<0.30) return 38;
      else if (pmuon>=0.30 and pmuon<0.45) return 39;
      else if (pmuon>=0.45 and pmuon<0.77) return 40;
      else if (pmuon>=0.77 and pmuon<1.28) return 41;
      else if (pmuon>=1.28 and pmuon<2.5) return 42;
      else return -10000;
    }
    else return -10000;
}

// return bin number of a variable in alternative binning choice
// bin edges are defined in config (default: ./configurations/alt_var_xbins.txt)
int LEEana::alt_var_index(std::string var1, float val1, std::string var2, float val2, std::string config){

  if (map_var_hist.empty()) {
    // load all variables from configuration
    std::vector<float> xbins;
    float val;
    std::string varname, line;
    std::ifstream in(config);
    while ( std::getline(in, line) ) {
      std::istringstream ss(line);
      ss >> varname;
      xbins.clear();
      while (ss >> val) {
        xbins.push_back(val);
      }
      map_var_hist[varname] = TH1F(varname.c_str(), varname.c_str(), xbins.size()-1, xbins.data());
    }
    in.close();
  }

  auto varhist1 = map_var_hist[var1];
  auto varhist2 = map_var_hist[var2];
  int bin1 = varhist1.FindBin(val1);
  int bin2 = varhist2.FindBin(val2);
  int nBins1 = varhist1.GetNbinsX();
  int nBins2 = varhist2.GetNbinsX();
  if( bin1>0 and bin1<=nBins1 and bin2>0 and bin2<=nBins2){
    return bin1 + (bin2-1)*nBins1;
  }

  return -1;
}


//option:1 costheta, 2 theta in rad, 3 theta in deg, 4 cosphi, 5 phi in rad, 6 phi in deg
double LEEana::get_angle(double momentum_0, double momentum_1, double momentum_2, int option, bool to_numi){

  double theta=0;
  double phi=0;

  double momentum_0_new = momentum_0;
  double momentum_1_new = momentum_1;
  double momentum_2_new = momentum_2;

  if (to_numi){
    //define rotation
    double R11 = 9.21033958e-01;
    double R12 = 6.18338739e-05; 
    double R13 = -3.89544856e-01;
    double R21 = 2.26996774e-02;  
    double R22 = 9.98285581e-01;  
    double R23 = 5.38031249e-02;
    double R31 = 3.88820884e-01; 
    double R32 = -5.84161267e-02;  
    double R33 = 9.19386684e-01;
    
    momentum_0_new = momentum_0*R11 + momentum_1*R12 + momentum_2*R13;
    momentum_1_new = momentum_0*R21 + momentum_1*R22 + momentum_2*R23;
    momentum_2_new = momentum_0*R31 + momentum_1*R32 + momentum_2*R33;
  }

  double momentum_perp = sqrt(momentum_0_new * momentum_0_new + momentum_1_new * momentum_1_new);
  theta = atan2(momentum_perp, momentum_2_new);
  phi = atan2(momentum_1_new, momentum_0_new);
  if (option==1) return cos(theta);
  if (option==2) return theta;
  if (option==3) return theta*180/3.14159;
  if (option==4) return cos(phi);
  if (option==5) return phi;
  if (option==6) return phi*180/3.14159;
  std::cout<<"LEEana::get_angle: Warning, option not included, defaulting to returning cos(theta)"<<std::endl;
  return cos(theta);//just defualt to returning costheta
}

double LEEana::get_angle_to_absorber(double momentum_0, double momentum_1, double momentum_2){
    double absorber_dir_0 = 0.33;
    double absorber_dir_1 = 0.75;
    double absorber_dir_2 = -0.59;
    double cos = absorber_dir_0*momentum_0+absorber_dir_1*momentum_1+absorber_dir_2*momentum_2;
    cos = cos / sqrt(momentum_0*momentum_0 + momentum_1*momentum_1 + momentum_2*momentum_2) / sqrt(absorber_dir_0*absorber_dir_0 + absorber_dir_1*absorber_dir_1 + absorber_dir_2*absorber_dir_2); 
    double angle = TMath::ACos(cos);
    return angle;
}

//option:1 costheta, 2 theta in rad, 3 theta in deg, 4 cosphi, 5 phi in rad, 6 phi in deg
double LEEana::get_true_nu_angle(PFevalInfo& pfeval, int option, bool to_numi){

  double momentum_0 = pfeval.truth_nu_momentum[0];
  double momentum_1 = pfeval.truth_nu_momentum[1];
  double momentum_2 = pfeval.truth_nu_momentum[2];  

  return get_angle(momentum_0, momentum_1, momentum_2, option, to_numi);

} 

double LEEana::get_true_nu_angle_to_absorber(PFevalInfo& pfeval){

  double momentum_0 = pfeval.truth_nu_momentum[0];
  double momentum_1 = pfeval.truth_nu_momentum[1];
  double momentum_2 = pfeval.truth_nu_momentum[2];

  return get_angle_to_absorber(momentum_0, momentum_1, momentum_2);

}

bool LEEana::is_true_kdar_all(EvalInfo& eval, PFevalInfo& pfeval){
    bool flag=false;
    if(eval.truth_nuEnergy>237) return flag;
    if(eval.truth_nuEnergy<234) return flag;
    if(eval.truth_nuPdg!=14)  return flag;
    if(eval.truth_isCC!=1) return flag;
    //if(eval.truth_vtxInside!=1)  return flag;
    double truth_nu_angle_to_absorber = get_true_nu_angle_to_absorber(pfeval)*180/3.14159;
    if(truth_nu_angle_to_absorber>8) return flag;
    flag=true;
    return flag;
}

bool LEEana::is_true_kdar(EvalInfo& eval, PFevalInfo& pfeval){
    bool flag=false;
    if(eval.truth_vtxInside!=1)  return flag;
    if(!is_true_kdar_all(eval, pfeval)) return flag;
    flag=true;
    return flag;
}

bool LEEana::is_true_kdar_outFV(EvalInfo& eval, PFevalInfo& pfeval){
    bool flag=false;
    if(eval.truth_vtxInside!=0)  return flag;
    if(!is_true_kdar_all(eval, pfeval)) return flag;
    flag=true;
    return flag;
}

double LEEana::get_ssmE(TaggerInfo& tagger_info, bool flag_2track){
  if(tagger_info.ssm_kine_energy<0) return -999;
  double ssm_E=0;  
  ssm_E+=tagger_info.ssm_kine_energy+105.7;
  if(tagger_info.ssm_prim_track1_kine_energy_range>0) ssm_E+=tagger_info.ssm_prim_track1_kine_energy_range;
  if(tagger_info.ssm_prim_track2_kine_energy_range>0 && !flag_2track) ssm_E+=tagger_info.ssm_prim_track2_kine_energy_range;
  if(tagger_info.ssm_prim_shw1_kine_energy_range>0 && !flag_2track) ssm_E+=tagger_info.ssm_prim_shw1_kine_energy_best;
  if(tagger_info.ssm_prim_shw2_kine_energy_range>0 && !flag_2track) ssm_E+=tagger_info.ssm_prim_shw2_kine_energy_best;
  return ssm_E;
}


bool LEEana::is_kdar_presel(TaggerInfo& tagger_info, EvalInfo& eval){
  bool flag=false;
  if(eval.match_isFC==0) return flag;
  if(tagger_info.ssm_kine_energy<0) return flag;
  if(tagger_info.ssm_kine_reco_Enu<10) return flag;
  if(tagger_info.ssm_kine_reco_Enu>350) return flag;
  if(tagger_info.ssm_kine_pio_mass>50) return flag;
  if(tagger_info.ssm_cosmict_flag_9==1) return flag;
  if(get_ssmE(tagger_info)>270) return flag;
  flag=true;
  return flag;
}

bool LEEana::is_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut, double hiE_cut){
  bool flag=false;
  if(tagger_info.ssm_kdar_score_lowE<lowE_cut) return flag;
  if(tagger_info.ssm_kdar_score_hiE<hiE_cut) return flag;
  flag=true;
  return flag;
}

bool LEEana::is_lowE_kdar_bdtsel(TaggerInfo& tagger_info, double lowE_cut){
  bool flag=false;
  if(tagger_info.ssm_kdar_score_lowE<lowE_cut) return flag;
  flag=true;
  return flag;
}

bool LEEana::is_hiE_kdar_bdtsel(TaggerInfo& tagger_info, double hiE_cut){
  bool flag=false;
  if(tagger_info.ssm_kdar_score_hiE<hiE_cut) return flag;
  flag=true;
  return flag;
}

double LEEana::get_oldflux_weight(EvalInfo& eval, PFevalInfo& pfeval){
std::vector<std::vector<double>> ratio_numu = {
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.071, 1.044, 1.106, 1.126, 1.078, 1.108, 1.13, 1.098, 1.089, 1.026, 1.006, 0.99, 1.003, 1.04, 1.046, 1.009, 1.005, 0.976, 1.021, 1.014, 1.023, 1.036, 1.043, 1.046, 1.03, 1.04, 1.018, 1.013, 1.003, 1.15, 1.214, 1.22, 1.251, 1.241, 1.219, 1.222, 1.457, 2.138, 1.683, 1.508, 1.183, 1.47, 1.054, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.802, 0.799, 0.808, 0.808, 0.801, 0.8, 0.807, 0.819, 0.952, 0.998, 0.997, 0.997, 0.997, 0.998, 0.998, 0.995, 0.996, 0.996, 0.995, 0.996, 0.996, 0.995, 1.016, 1.027, 1.026, 1.026, 1.026, 1.027, 1.026, 1.115, 1.18, 1.179, 1.179, 1.178, 1.181, 1.177, 1.43, 1.646, 1.425, 1.546, 1.384, 1.419, 1.133, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.604, 0.599, 0.668, 0.677, 0.705, 0.69, 0.702, 0.752, 0.879, 0.863, 0.879, 0.899, 0.879, 0.893, 0.89, 0.922, 0.93, 0.937, 0.916, 0.917, 0.935, 0.932, 0.958, 0.959, 0.958, 0.957, 0.958, 0.96, 0.96, 1.039, 1.083, 1.083, 1.083, 1.083, 1.083, 1.083, 1.247, 1.471, 1.416, 1.819, 1.806, 1.531, 1.055, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.754, 0.749, 0.749, 0.771, 0.752, 0.733, 0.799, 0.717, 0.889, 0.913, 0.916, 0.913, 0.912, 0.912, 0.915, 0.921, 0.923, 0.925, 0.924, 0.922, 0.925, 0.923, 0.949, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 1.061, 1.102, 1.103, 1.103, 1.103, 1.105, 1.103, 1.196, 1.412, 1.383, 1.915, 1.842, 1.747, 1.094, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.983, 0.954, 0.946, 1.006, 0.973, 0.963, 0.988, 0.979, 0.964, 1.023, 1.009, 0.917, 0.944, 0.99, 0.949, 0.935, 0.928, 0.91, 0.902, 0.918, 0.921, 0.884, 0.92, 0.923, 0.937, 0.92, 0.924, 0.937, 0.925, 0.97, 1.0, 1.004, 1.001, 0.999, 1.002, 1.0, 1.102, 1.375, 1.386, 1.91, 1.798, 1.855, 1.086, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.593, 0.643, 0.633, 0.647, 0.637, 0.629, 0.623, 0.659, 0.681, 0.685, 0.699, 0.685, 0.675, 0.693, 0.686, 0.751, 0.763, 0.76, 0.755, 0.756, 0.758, 0.757, 0.875, 0.914, 0.921, 0.909, 0.916, 0.916, 0.916, 0.971, 1.004, 1.007, 1.004, 1.002, 1.004, 1.003, 1.173, 1.591, 1.244, 1.893, 1.675, 1.991, 1.116, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.804, 0.735, 0.786, 0.728, 0.81, 0.761, 0.775, 0.748, 0.804, 0.835, 0.777, 0.812, 0.811, 0.814, 0.82, 0.753, 0.732, 0.731, 0.736, 0.724, 0.73, 0.73, 0.797, 0.819, 0.819, 0.819, 0.819, 0.819, 0.819, 0.977, 1.063, 1.064, 1.062, 1.063, 1.062, 1.063, 1.244, 1.549, 1.222, 1.762, 1.706, 1.973, 1.118, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.662, 0.678, 0.679, 0.68, 0.68, 0.675, 0.673, 0.675, 0.681, 0.68, 0.681, 0.68, 0.682, 0.68, 0.681, 0.767, 0.789, 0.781, 0.791, 0.787, 0.78, 0.769, 0.934, 0.957, 0.954, 0.962, 0.955, 0.96, 0.959, 0.95, 0.939, 0.931, 0.933, 0.941, 0.932, 0.932, 1.159, 1.367, 1.638, 1.727, 1.709, 2.052, 1.142, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.735, 0.735, 0.692, 0.423, 0.635, 0.462, 0.44, 0.514, 0.979, 0.958, 0.995, 1.009, 0.96, 1.013, 0.989, 0.794, 0.751, 0.751, 0.75, 0.75, 0.751, 0.75, 0.873, 0.926, 0.927, 0.926, 0.929, 0.93, 0.927, 0.926, 0.938, 0.934, 0.933, 0.931, 0.931, 0.935, 1.025, 1.469, 1.158, 1.832, 1.513, 2.162, 1.093, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.77, 0.799, 0.814, 0.819, 0.823, 0.881, 0.76, 0.797, 0.792, 0.786, 0.774, 0.783, 0.78, 0.741, 0.75, 0.749, 0.752, 0.753, 0.746, 0.751, 0.723, 0.715, 0.715, 0.715, 0.715, 0.715, 0.715, 0.963, 1.064, 1.064, 1.064, 1.064, 1.062, 1.064, 1.329, 1.656, 1.295, 1.464, 1.352, 2.417, 1.208, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.456, 0.467, 0.522, 0.504, 0.47, 0.47, 0.484, 0.477, 0.498, 0.49, 0.486, 0.49, 0.489, 0.695, 0.757, 0.749, 0.755, 0.751, 0.753, 0.766, 0.813, 0.835, 0.84, 0.845, 0.843, 0.832, 0.862, 0.986, 1.059, 1.059, 1.062, 1.058, 1.061, 1.062, 1.131, 1.466, 1.356, 1.016, 1.418, 2.324, 1.162, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.439, 0.839, 1.439, 0.482, 0.556, 0.516, 0.671, 0.542, 0.347, 0.354, 0.362, 0.384, 0.402, 0.686, 0.743, 0.743, 0.742, 0.741, 0.742, 0.742, 0.785, 0.808, 0.808, 0.801, 0.807, 0.833, 0.79, 0.95, 1.064, 1.076, 1.059, 1.064, 1.068, 1.066, 0.937, 1.215, 1.219, 1.541, 1.447, 2.447, 1.223, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.66, 1.0, 1.0, 1.607, 1.041, 0.869, 0.612, 0.667, 0.621, 0.836, 0.7, 0.704, 0.719, 0.881, 0.909, 0.919, 0.897, 0.916, 0.909, 0.923, 0.78, 0.726, 0.732, 0.732, 0.741, 0.729, 0.727, 0.906, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.236, 1.843, 1.104, 1.164, 1.382, 2.662, 1.173, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.81, 1.0, 0.763, 0.777, 0.797, 0.773, 0.455, 0.368, 0.351, 0.376, 0.345, 0.325, 0.348, 0.806, 0.904, 0.912, 0.91, 0.904, 0.905, 0.922, 0.879, 0.863, 0.864, 0.856, 0.871, 0.868, 0.856, 0.906, 0.911, 0.922, 0.926, 0.921, 0.93, 0.927, 1.239, 1.187, 1.104, 1.215, 1.246, 2.486, 1.216, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.871, 1.0, 0.939, 1.7976931348623157e+308, 0.884, 0.728, 0.704, 0.625, 0.657, 0.639, 0.652, 0.757, 0.839, 0.849, 0.846, 0.85, 0.848, 0.843, 0.805, 0.794, 0.802, 0.801, 0.799, 0.798, 0.802, 0.8, 0.798, 0.818, 0.821, 0.807, 0.807, 0.799, 0.979, 1.244, 0.902, 1.086, 1.054, 2.212, 1.216, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.213, 0.425, 0.755, 0.893, 0.726, 0.772, 0.753, 0.748, 0.765, 0.766, 0.734, 0.733, 0.718, 0.738, 0.739, 0.739, 0.737, 0.768, 0.758, 0.776, 0.767, 0.764, 0.759, 0.77, 0.964, 1.145, 1.145, 1.145, 1.146, 1.145, 1.145, 1.215, 1.023, 1.134, 1.319, 1.185, 2.616, 1.223, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.476, 1.7976931348623157e+308, 1.0, 0.883, 0.853, 0.909, 0.858, 0.916, 0.87, 0.856, 0.756, 0.683, 0.697, 0.678, 0.681, 0.709, 0.711, 0.797, 0.848, 0.857, 0.839, 0.853, 0.832, 0.839, 0.863, 0.868, 0.869, 0.869, 0.868, 0.868, 0.868, 1.306, 1.154, 0.994, 1.27, 1.344, 2.178, 1.176, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 0.865, 0.805, 0.528, 0.658, 0.638, 0.729, 0.581, 0.847, 0.821, 0.841, 0.823, 0.816, 0.835, 0.832, 0.794, 0.783, 0.782, 0.783, 0.782, 0.782, 0.782, 0.797, 0.793, 0.8, 0.808, 0.805, 0.806, 0.813, 1.078, 1.297, 1.109, 1.451, 1.499, 2.389, 1.152, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 0.888, 1.194, 1.289, 1.125, 1.09, 0.766, 0.718, 0.717, 0.719, 0.714, 0.718, 0.721, 0.717, 0.841, 0.896, 0.904, 0.892, 0.899, 0.902, 0.9, 0.884, 0.873, 0.872, 0.873, 0.874, 0.875, 0.874, 1.018, 1.197, 1.159, 1.509, 1.268, 2.169, 1.112, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.337, 1.868, 1.308, 1.608, 1.724, 1.032, 1.197, 0.717, 0.717, 0.725, 0.734, 0.712, 0.72, 0.718, 0.713, 0.709, 0.714, 0.711, 0.708, 0.712, 0.838, 0.881, 0.882, 0.881, 0.88, 0.882, 0.88, 1.366, 1.307, 1.158, 1.507, 1.341, 2.035, 1.207, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.018, 3.37, 3.194, 2.43, 3.151, 2.604, 1.248, 1.082, 1.033, 1.068, 1.081, 1.094, 1.058, 0.827, 0.796, 0.811, 0.814, 0.798, 0.806, 0.801, 0.879, 0.909, 0.913, 0.917, 0.912, 0.904, 0.912, 0.888, 1.214, 1.31, 1.8, 1.503, 2.261, 1.056, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 0.903, 0.903, 0.903, 1.0, 1.411, 2.456, 0.814, 1.009, 1.019, 0.854, 0.825, 0.827, 0.952, 0.787, 0.739, 0.732, 0.732, 0.737, 0.737, 0.726, 0.834, 0.887, 0.887, 0.887, 0.887, 0.887, 0.887, 1.199, 0.9, 1.034, 1.674, 1.507, 1.992, 1.056, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.707, 1.0, 0.746, 0.844, 0.668, 0.847, 0.692, 0.554, 0.558, 0.534, 0.585, 0.523, 0.515, 0.523, 0.824, 0.858, 0.82, 0.847, 0.866, 0.849, 0.854, 0.8, 0.791, 0.781, 0.788, 0.777, 0.784, 0.783, 1.031, 1.462, 1.052, 1.713, 1.687, 2.153, 1.214, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.222, 1.111, 1.111, 1.7976931348623157e+308, 1.0, 1.739, 1.7976931348623157e+308, 0.714, 0.682, 0.687, 0.683, 0.711, 0.689, 0.684, 0.786, 0.813, 0.799, 0.776, 0.799, 0.805, 0.81, 0.931, 0.958, 0.965, 0.978, 0.967, 0.975, 0.961, 1.089, 1.747, 1.272, 1.562, 1.544, 2.231, 1.126, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.428, 0.768, 1.0, 0.768, 0.768, 0.333, 0.576, 0.74, 0.579, 0.562, 0.636, 0.609, 0.633, 0.919, 0.977, 0.995, 0.987, 0.966, 0.949, 1.012, 0.925, 0.925, 0.925, 0.906, 0.917, 0.913, 0.934, 1.241, 0.87, 1.342, 1.648, 1.295, 2.258, 1.051, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.363, 0.301, 1.0, 0.301, 0.51, 0.961, 0.812, 0.81, 1.373, 0.833, 0.9, 1.021, 0.79, 0.775, 0.813, 0.806, 0.815, 0.808, 0.806, 0.807, 0.902, 0.93, 0.929, 0.924, 0.927, 0.931, 0.929, 1.103, 1.132, 1.388, 1.587, 1.977, 1.976, 1.045, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.862, 1.448, 1.446, 1.36, 1.511, 1.419, 1.426, 0.856, 0.682, 0.727, 0.787, 0.744, 0.798, 0.811, 0.885, 0.928, 0.935, 0.923, 0.914, 0.924, 0.922, 0.915, 1.253, 0.913, 2.165, 1.414, 2.05, 1.154, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 9.72, 1.29, 1.293, 1.328, 1.021, 1.13, 1.08, 1.193, 0.672, 0.648, 0.648, 0.646, 0.648, 0.648, 0.647, 0.787, 0.844, 0.863, 0.852, 0.852, 0.838, 0.85, 1.33, 1.592, 1.294, 2.678, 1.751, 2.151, 1.023, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 0.37, 0.37, 0.37, 0.644, 0.462, 0.476, 0.51, 0.651, 0.974, 0.935, 0.949, 0.929, 1.013, 0.918, 0.95, 0.942, 0.938, 0.936, 0.945, 0.944, 0.934, 1.153, 1.018, 1.241, 1.911, 1.937, 1.711, 0.96, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.061, 1.0, 1.0, 1.187, 0.194, 0.328, 0.194, 0.277, 0.275, 0.218, 0.214, 0.448, 0.799, 0.83, 0.811, 0.833, 0.793, 0.81, 0.935, 1.084, 0.991, 1.082, 1.064, 1.011, 0.989, 1.297, 0.909, 1.269, 2.282, 1.517, 2.148, 1.151, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 2.349, 1.696, 2.984, 2.06, 2.241, 1.614, 2.178, 0.546, 0.554, 0.534, 0.561, 0.559, 0.57, 0.527, 0.778, 0.83, 0.896, 0.85, 0.805, 0.851, 0.851, 0.834, 0.986, 1.014, 1.583, 1.78, 2.488, 0.994, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 0.506, 0.413, 0.366, 0.299, 0.389, 0.449, 0.299, 0.77, 0.645, 0.683, 0.704, 0.685, 0.676, 0.697, 0.818, 0.923, 0.911, 0.91, 0.886, 0.903, 0.898, 0.613, 1.341, 1.068, 2.149, 2.021, 2.029, 1.194, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.121, 0.696, 2.241, 1.0, 0.726, 0.886, 0.887, 1.005, 1.004, 1.034, 1.094, 1.073, 1.12, 1.14, 0.865, 0.841, 0.837, 0.84, 0.842, 0.838, 0.839, 0.871, 0.959, 1.028, 1.366, 2.087, 2.136, 1.09, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.547, 3.524, 5.103, 3.619, 1.602, 1.206, 1.64, 3.277, 0.918, 0.827, 0.853, 0.958, 1.029, 0.839, 0.97, 1.005, 1.003, 1.005, 1.004, 1.006, 1.004, 1.091, 1.4, 1.162, 1.741, 2.15, 2.425, 1.185, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.037, 0.475, 0.552, 0.575, 0.512, 0.581, 0.586, 0.962, 1.144, 1.177, 1.147, 1.1, 1.143, 1.154, 1.294, 2.628, 1.115, 2.009, 1.475, 2.021, 1.089, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 14.844, 1.0, 1.0, 5.6, 4.338, 11.044, 1.004, 0.637, 0.562, 0.625, 0.595, 0.549, 0.696, 0.897, 0.94, 1.088, 0.939, 1.167, 0.998, 1.01, 0.786, 1.818, 1.514, 1.979, 1.996, 2.349, 1.136, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.701, 1.0, 0.804, 0.779, 0.799, 0.775, 0.791, 0.857, 0.876, 0.868, 0.865, 0.852, 0.892, 0.839, 0.875, 2.486, 0.829, 2.182, 2.018, 2.663, 1.158, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.572, 1.0, 1.0, 1.0, 0.409, 1.7976931348623157e+308, 1.0, 1.567, 1.082, 1.082, 1.2, 1.1, 1.065, 1.178, 1.116, 1.084, 1.099, 1.142, 1.104, 1.133, 3.066, 0.655, 1.058, 1.471, 2.482, 3.013, 1.108, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 0.773, 0.815, 0.799, 0.823, 0.791, 0.799, 1.007, 0.857, 1.066, 0.877, 0.966, 0.972, 0.989, 1.745, 1.408, 1.468, 1.062, 1.123, 3.532, 1.098, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.142, 1.446, 1.203, 1.287, 1.319, 1.172, 0.66, 0.702, 0.677, 0.729, 0.66, 0.684, 0.693, 1.57, 1.349, 2.201, 2.752, 1.812, 3.094, 1.244, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 0.596, 1.0, 4.621, 4.621, 1.422, 1.798, 4.621, 1.966, 2.192, 1.146, 1.058, 1.075, 1.024, 1.103, 1.132, 1.07, 1.933, 1.296, 1.8, 1.464, 2.666, 2.667, 1.296, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 0.349, 0.361, 0.358, 0.363, 0.361, 0.369, 0.796, 0.796, 0.594, 0.653, 0.704, 0.741, 0.686, 1.167, 1.193, 1.196, 1.352, 2.746, 3.223, 1.368, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.456, 0.747, 1.025, 0.747, 1.494, 0.94, 0.642, 0.645, 0.59, 0.657, 0.61, 0.655, 0.42, 1.003, 0.969, 1.759, 1.75, 2.99, 1.087, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 0.448, 0.618, 0.385, 0.981, 0.887, 0.845, 0.703, 0.918, 0.854, 0.723, 1.003, 3.827, 3.258, 1.468, 1.048, 3.234, 1.216, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 7.323, 7.258, 10.773, 4.513, 1.521, 0.67, 0.718, 0.78, 0.614, 0.614, 0.62, 0.577, 1.03, 1.06, 4.491, 2.007, 3.277, 1.316, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 0.152, 1.0, 1.0, 0.282, 0.192, 0.301, 0.295, 0.492, 1.289, 1.269, 1.202, 1.148, 1.22, 1.321, 0.803, 4.194, 1.868, 0.794, 2.874, 3.292, 1.333, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.6, 3.112, 1.914, 4.076, 1.827, 1.957, 1.581, 0.7, 1.162, 1.469, 1.689, 3.404, 3.572, 3.876, 1.954, 3.389, 1.214, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.514, 0.875, 1.142, 1.286, 1.725, 1.776, 1.936, 4.0, 2.098, 1.576, 2.396, 2.618, 1.7976931348623157e+308, 4.838, 6.468, 2.4, 3.373, 1.342, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.228, 0.988, 1.009, 1.223, 1.062, 1.014, 1.058, 2.614, 3.96, 3.339, 3.047, 3.37, 3.239, 1.41, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.018, 1.052, 1.026, 1.022, 1.008, 1.016, 1.744, 12.647, 1.017, 1.08, 2.375, 4.456, 1.29, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 0.715, 0.646, 0.588, 0.651, 0.596, 0.66, 0.711, 1.041, 1.7976931348623157e+308, 0.872, 0.942, 6.092, 3.7, 1.511, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 0.459, 0.276, 0.264, 0.279, 0.264, 0.297, 0.332, 0.524, 1.7976931348623157e+308, 1.896, 1.567, 0.866, 3.304, 1.424, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 2.297, 2.672, 3.005, 4.625, 5.202, 4.693, 5.701, 1.7976931348623157e+308, 2.761, 3.955, 11.248, 1.209, 3.686, 1.703, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.496, 1.0, 1.994, 2.992, 1.0, 1.833, 1.922, 1.892, 2.011, 1.904, 1.867, 4.927, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.38, 5.238, 1.43, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.1, 2.922, 3.827, 6.254, 1.276, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.008, 0.958, 1.008, 0.869, 0.903, 0.88, 1.197, 1.7976931348623157e+308, 1.7976931348623157e+308, 6.451, 7.43, 2.627, 1.587, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 2.202, 2.718, 1.7976931348623157e+308, 3.2, 14.265, 1.7976931348623157e+308, 1.7976931348623157e+308, 32.438, 1.943, 4.475, 1.717, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 12.824, 15.934, 4.275, 9.57, 11.376, 7.798, 38.251, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 19.564, 4.918, 1.691, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.898, 1.0, 1.0, 1.347, 1.0, 1.0, 1.356, 1.7976931348623157e+308, 1.247, 1.199, 1.7976931348623157e+308, 1.873, 1.417, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 12.967, 40.94, 4.118, 5.821, 1.478, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 7.165, 5.891, 1.162, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.437, 1.033, 0.958, 0.958, 0.992, 1.002, 0.998, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.64, 10.308, 39.996, 3.216, 1.487, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.216, 1.7976931348623157e+308, 0.635, 1.7976931348623157e+308, 1.7976931348623157e+308, 2.792, 1.801, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 5.472, 4.653, 2.022, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 3.291, 7.222, 1.678, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 2.244, 1.7976931348623157e+308, 1.7976931348623157e+308, 8.952, 4.285, 2.039, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 8.99, 4.522, 2.292, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 7.283, 4.817, 1.552, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 4.901, 15.028, 4.53, 1.539, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.0, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 7.3, 29.598, 1.7976931348623157e+308, 4.794, 1.546, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.0, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.0, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 1.7976931348623157e+308, 2.503, 1.665, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
  int nbins_Enu = ratio_numu.size();
  double bin_width_Enu = 7000/double(nbins_Enu);
  int nbins_dist = ratio_numu.at(0).size();
  double bin_width_dist = 850/double(nbins_dist);
  //std::cout<<"nbins_Enu="<<nbins_Enu<<" bin_width_Enu="<<bin_width_Enu<<" nbins_dist="<<nbins_dist<<" bin_width_dist="<<bin_width_dist<<std::endl; 
  int this_Enu_bin = std::floor(eval.truth_nuEnergy/bin_width_Enu);
  if (this_Enu_bin >= nbins_Enu) this_Enu_bin = nbins_Enu-1;

  double dist = pfeval.mcflux_dk2gen + pfeval.mcflux_gen2vtx;
  int this_dist_bin = std::floor(dist/bin_width_dist);
  if (this_dist_bin >= nbins_dist) this_dist_bin = nbins_dist-1;
  //std::cout<<"eval.truth_nuEnergy="<<eval.truth_nuEnergy<<" this_Enu_bin="<<this_Enu_bin<<" dist="<<dist<<" this_dist_bin="<<this_dist_bin<<std::endl;
  //std::cout<<ratio_numu.at(this_Enu_bin).at(this_dist_bin)<<" "<<1.0/ratio_numu.at(this_Enu_bin).at(this_dist_bin)<<std::endl;
  double new_weight = ratio_numu.at(this_Enu_bin).at(this_dist_bin);
  if (new_weight>10) new_weight=1;
  return new_weight;
 

}

double LEEana::recover_nsbeamtiming(PFevalInfo& pfeval, bool flag_data, bool flag_kdar){

  if(pfeval.PMT_Amp->size()>2){
    if(!flag_data && !flag_kdar) return pfeval.evtTimeNS_cor;
    return pfeval.evtTimeNS;
  }
  if(pfeval.PMT_Amp->size()==0){
    if(!flag_data && !flag_kdar) return pfeval.evtTimeNS_cor;
    return pfeval.evtTimeNS;
  }

  double f_evtTimeNS=-99999;

  float f_shiftoffset=0;
  float f_ccnd1_a=0.4343;
  float f_ccnd1_b=6.2884;
  float f_ccnd2_a=0.0637;
  float f_ccnd2_b=1.489;
  float f_ccnd3_a=0;
  float f_ccnd3_b=0;
  float f_ccnd3_c=0;
  float f_ccnd3_d=0;
  float f_ccnd4_a=0.0125;
  float f_ccnd4_b=2.3152;
  float f_ccnd4_2_a=0;
  float f_ccnd4_2_b=0;
  float f_dist_cut_x_cor=999999999;

  if(!flag_data){
    f_shiftoffset=0;
    f_ccnd1_a=0;
    f_ccnd1_b=0;
    f_ccnd2_a=0;
    f_ccnd2_b=0;
    f_ccnd3_a=0;
    f_ccnd3_b=0;
    f_ccnd3_c=0;
    f_ccnd3_d=0;
    f_ccnd4_a=0;
    f_ccnd4_b=0;
    f_ccnd4_2_a=0;
    f_ccnd4_2_b=0;
  }

  double RWM_T = pfeval.RWM_Time;
  
  Float_t x = pfeval.reco_nuvtxX;
  Float_t y = pfeval.reco_nuvtxY;
  Float_t z = pfeval.reco_nuvtxZ;
  
  std::vector<float> *sps_x = new std::vector<float>;
  std::vector<float> *sps_y = new std::vector<float>;
  std::vector<float> *sps_z = new std::vector<float>;
  std::vector<float> *sps_t = new std::vector<float>;
  //method using linear extrapolation and PID

  std::map<int,std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > >my_particle_times;
    //<id,<x,y,z,t>>

    //only do primary particles first
    for (int particle=0; particle<pfeval.reco_Ntrack; particle++){
      int mother = pfeval.reco_mother[particle];
      if(mother!=0) {continue;} 
      int id = pfeval.reco_id[particle];
      auto this_id = my_particle_times.find(id);      
      if (this_id==my_particle_times.end()){
        my_particle_times[id] = get_extrapolated_times(pfeval,particle, 0);
      }

      //now find the direct daughters of this primary particle
      std::vector< std::tuple<int,int,double> > daughters; //<id,particle,mother_time>
      for(int daughter_particle=0; daughter_particle<pfeval.reco_Ntrack; daughter_particle++){
        int mother_of_daughter = pfeval.reco_mother[daughter_particle];
        int daughter_id = pfeval.reco_id[daughter_particle];
        if (mother_of_daughter == id) {
          double mother_time = std::get<3>(my_particle_times[id])->back();
          daughters.push_back(std::make_tuple(daughter_id,daughter_particle,mother_time));
        }
      }

      //Now add this daughter and find the daughters of the daughter
      //keep going untill we exahust all daughters of daughters
      while (daughters.size() > 0){
        int daughter_id = std::get<0>(daughters.front());
        auto this_daughter_id = my_particle_times.find(daughter_id);
        //skip is we have already added it
	if (this_daughter_id!=my_particle_times.end()){
          daughters.erase(daughters.begin());
          continue;
        }
        my_particle_times[daughter_id] = get_extrapolated_times(pfeval, std::get<1>(daughters.front()), std::get<2>(daughters.front()));
        for(int daughter_daughter_particle=0; daughter_daughter_particle<pfeval.reco_Ntrack; daughter_daughter_particle++){
          int mother_of_daughter_daughter_id = pfeval.reco_mother[daughter_daughter_particle];
          int daughter_daughter_id = pfeval.reco_id[daughter_daughter_particle];
          auto this_daughter_daughter_id = my_particle_times.find(daughter_daughter_id);
	  if (daughter_id == mother_of_daughter_daughter_id && this_daughter_daughter_id==my_particle_times.end()) {
            //new daughter that we need to check on, its a daughter of this daughter.
	    //Set its mother time according to the time we just added for the daughter we were working on.
            daughters.push_back(std::make_tuple(daughter_daughter_id,daughter_daughter_particle,std::get<3>(my_particle_times[daughter_id])->back())); 
	  }
        } //no more daugters of the current daughter
        daughters.erase(daughters.begin());
      }//no more daughter or duaghters of daughters
    }//end loop over all primary particles

    //unpack this for the rest of the code
    for (auto my_particle = my_particle_times.begin(); my_particle != my_particle_times.end(); my_particle++) {
      std::vector<float>* _x = std::get<0>(my_particle->second);
      std::vector<float>* _y = std::get<1>(my_particle->second);
      std::vector<float>* _z = std::get<2>(my_particle->second);
      std::vector<float>* _t = std::get<3>(my_particle->second);
      for(uint point=0; point < _x->size(); point++){
        sps_x->push_back(_x->at(point));
        sps_y->push_back(_y->at(point));
        sps_z->push_back(_z->at(point));
        sps_t->push_back(_t->at(point));	 
	//if(pfeval.evtTimeNS>=0) std::cout<<point<<" "<<sps_x->at(point)<<" "<<sps_y->at(point)<<" "<<sps_z->at(point)<<" "<<sps_t->at(point)<<std::endl;
      }
    }
  //}
  
  double PMT0[3]={-11.4545, -28.625, 990.356};  double PMT1[3]={-11.4175, 27.607, 989.712};
  double PMT2[3]={-11.7755, -56.514, 951.865};  double PMT3[3]={-11.6415, 55.313, 951.861};
  double PMT4[3]={-12.0585, -56.309, 911.939};  double PMT5[3]={-11.8345, 55.822, 911.065};
  double PMT6[3]={-12.1765, -0.722, 865.599};   double PMT7[3]={-12.3045, -0.502, 796.208};
  double PMT8[3]={-12.6045, -56.284, 751.905};  double PMT9[3]={-12.5405, 55.625, 751.884};
  double PMT10[3]={-12.6125, -56.408, 711.274}; double PMT11[3]={-12.6615, 55.8, 711.073};
  double PMT12[3]={-12.6245, -0.051, 664.203};  double PMT13[3]={-12.6515, -0.549, 585.284};
  double PMT14[3]={-12.8735, 55.822, 540.929};  double PMT15[3]={-12.6205, -56.205, 540.616};
  double PMT16[3]={-12.5945, -56.323, 500.221}; double PMT17[3]={-12.9835, 55.771, 500.134};
  double PMT18[3]={-12.6185, -0.875, 453.096};  double PMT19[3]={-13.0855, -0.706, 373.839};
  double PMT20[3]={-12.6485, -57.022, 328.341}; double PMT21[3]={-13.1865, 54.693, 328.212};
  double PMT22[3]={-13.4175, 54.646, 287.976};  double PMT23[3]={-13.0075, -56.261, 287.639};
  double PMT24[3]={-13.1505, -0.829, 242.014};  double PMT25[3]={-13.4415, -0.303, 173.743};
  double PMT26[3]={-13.3965, 55.249, 128.354};  double PMT27[3]={-13.2784, -56.203, 128.18};
  double PMT28[3]={-13.2375, -56.615, 87.8695}; double PMT29[3]={-13.5415, 55.249, 87.7605};
  double PMT30[3]={-13.4345, 27.431, 51.1015};  double PMT31[3]={-13.1525, -28.576, 50.4745};
  double PMT[32][3];    for(int j=0; j<3; j++){ PMT[30][j]=PMT30[j]; PMT[31][j]=PMT31[j];
  PMT[0][j]=PMT0[j];   PMT[10][j]=PMT10[j]; PMT[20][j]=PMT20[j]; PMT[1][j]=PMT1[j];   PMT[11][j]=PMT11[j];
  PMT[21][j]=PMT21[j]; PMT[2][j]=PMT2[j];   PMT[12][j]=PMT12[j]; PMT[22][j]=PMT22[j]; PMT[3][j]=PMT3[j];
  PMT[13][j]=PMT13[j]; PMT[23][j]=PMT23[j]; PMT[4][j]=PMT4[j];   PMT[14][j]=PMT14[j]; PMT[24][j]=PMT24[j];
  PMT[5][j]=PMT5[j];   PMT[15][j]=PMT15[j]; PMT[25][j]=PMT25[j]; PMT[6][j]=PMT6[j];   PMT[16][j]=PMT16[j];
  PMT[18][j]=PMT18[j]; PMT[28][j]=PMT28[j]; PMT[9][j]=PMT9[j];   PMT[19][j]=PMT19[j]; PMT[29][j]=PMT29[j];
  PMT[26][j]=PMT26[j]; PMT[7][j]=PMT7[j];   PMT[17][j]=PMT17[j]; PMT[27][j]=PMT27[j]; PMT[8][j]=PMT8[j];}
  double offset[32]={1.03002, -5.18104, -2.11164, -5.99395, -1.25798, 0.633079, 2.87666, 2.21969, 0.885092, 2.35423,
    -1.63039, -1.83775, -0.859883, 3.4741, 1.84833, 1.58233, -2.71783, 0, 3.18776, 0.982666, 0.728438, 0.280592, -5.27068,
    -3.27857, -1.41196, 1.59643, 1.41425, -1.62682, -2.55772, 1.49136, -0.522791, 0.974533};
  if(!flag_data){
    for(int i=0; i<32; i++){offset[i]=0;}//no need to apply the additional pmt calibration to the MC
    for(int j=0; j<3; j++){//do the PMT remapping for the MC
      double temp = PMT[31][j];
      PMT[31][j] = PMT[30][j];
      PMT[30][j] = PMT[29][j];
      PMT[29][j] = PMT[28][j];
      PMT[28][j] = PMT[27][j];
      PMT[27][j] = PMT[26][j];
      PMT[26][j] = temp;
    }
  }
  //================================================================================================================
  double ccnd1, ccnd2,ccnd3, ccnd4;
  double Ph_Tot, nuToF, DPh,DLh, tPhelp,tp, tDPhelp,tDP, tDLhelp,tDL;
  double Med_TT3=-9999.;
  double TT_merged = -9999.;
  nuToF=0;
  //===================================================================================================================
  //===================================================================================================================
  //do not use a time cut in the NuMI case
  std::vector<double> timeProp = std::vector<double>(pfeval.PMT_Amp->size(),0);
  std::vector<double> timeDP = std::vector<double>(pfeval.PMT_Amp->size(),0);
  std::vector<double> timeDL = std::vector<double>(pfeval.PMT_Amp->size(),0);
  //--------------------------------------------------------------------------------------------------------------------
  
    double dist = z; //in BNB correct to front face of TPC, in NuMI correct to plane perpendicular to the beam
    //if(fIsNuMI) {
      TVector3 target_dir(-0.46, -0.05, -0.885);
      double min_a = -122.86902944472968;  
      double min_b = 80.60659897339974; 
      double min_c = 59.34119182916038;
      dist = ( (min_a-x)*target_dir[0] + (min_b-y)*target_dir[1] + (min_c-z)*target_dir[2] ) / sqrt(target_dir[0]*target_dir[0] + target_dir[1]*target_dir[1] + target_dir[2]*target_dir[2] );
    //}
    nuToF=dist*0.033356;
    for(uint i=0; i<pfeval.PMT_Amp->size(); i++){
        Ph_Tot += pfeval.PMT_Amp->at(i);
        tp=5000000000.0;
        tDL=5000000000.0;
        tDP=5000000000.0;
	for(uint j=0; j<sps_x->size(); j++){
          DPh=abs(sqrt(TMath::Power(x-sps_x->at(j),2)+TMath::Power(y-sps_y->at(j),2)+TMath::Power(z-sps_z->at(j),2)));
          DLh=abs(sqrt(TMath::Power(PMT[pfeval.PMT_ID->at(i)][0]-sps_x->at(j),2)+TMath::Power(PMT[pfeval.PMT_ID->at(i)][1]-sps_y->at(j),2)+TMath::Power(PMT[pfeval.PMT_ID->at(i)][2]-sps_z->at(j),2)));
          tPhelp=sps_t->at(j)+(DLh*0.0746);
          tDPhelp=sps_t->at(j);
          tDLhelp=(DLh*0.0746);
  	  if(tPhelp<tp){
            tp=tPhelp;
            tDP=tDPhelp;
            tDL=tDLhelp;
          }
	}
        timeProp[i]=tp;
        timeDP[i]=tDP;
        timeDL[i]=tDL;
    }
    double TT3_array[32];
    float RWM_offset = 5700.0 - f_shiftoffset;
    for(uint i=0; i<pfeval.PMT_Amp->size(); i++){
      ccnd1= timeProp[i]*(f_ccnd1_a)-(f_ccnd1_b);
      ccnd2= pfeval.PMT_Amp->at(i)*(f_ccnd2_a)-(f_ccnd2_b);
      ccnd4= x*(f_ccnd4_a)-(f_ccnd4_b);//for x dependent correction
      if(x>f_dist_cut_x_cor) ccnd4= x*(f_ccnd4_2_a)-(f_ccnd4_2_b);//for x dependent correction at high x in the MC
      if(Ph_Tot>150){ccnd3=f_ccnd3_a-f_ccnd3_b*Ph_Tot+f_ccnd3_c*Ph_Tot*Ph_Tot;}
      else{ccnd3=f_ccnd3_d;}

      
      //all the corrections
      TT3_array[i]=pfeval.PMT_Time->at(i)-RWM_T+RWM_offset-nuToF-timeProp[i]-offset[pfeval.PMT_ID->at(i)]+ccnd1+ccnd2+ccnd3+ccnd4;
    }
    Med_TT3=TMath::Median((Long64_t)pfeval.PMT_Amp->size(),TT3_array);
  
  f_evtTimeNS = Med_TT3;
  if(!flag_data && !flag_kdar){ 
    f_evtTimeNS = Med_TT3+pfeval.cor_nu_deltatime; 
  }
  //if(f_evtTimeNS>0) std::cout<<"f_evtTimeNS="<<f_evtTimeNS<<" pfeval.evtTimeNS_cor="<<pfeval.evtTimeNS_cor<<" pfeval.cor_nu_deltatime="<<pfeval.cor_nu_deltatime<<" pfeval.evtTimeNS="<<pfeval.evtTimeNS<<std::endl;
  return f_evtTimeNS;

}



std::tuple< std::vector<float>*,std::vector<float>*,std::vector<float>*,std::vector<float>* > LEEana::get_extrapolated_times(PFevalInfo& pfeval, int particle, double mother_time){

          double fsol = 0.033356;

  	  //<x,y,z,t>
          std::vector<float> *_x = new std::vector<float>;
          std::vector<float> *_y = new std::vector<float>;
          std::vector<float> *_z = new std::vector<float>;
          std::vector<float> *_t = new std::vector<float>;

          double dx = 0.5;
          int pdg = pfeval.reco_pdg[particle];
          double start_x_pos = pfeval.reco_startXYZT[particle][0];
          double start_y_pos = pfeval.reco_startXYZT[particle][1];
          double start_z_pos = pfeval.reco_startXYZT[particle][2];
          double end_x_pos = pfeval.reco_endXYZT[particle][0];
          double end_y_pos = pfeval.reco_endXYZT[particle][1];
          double end_z_pos = pfeval.reco_endXYZT[particle][2];
          double length = sqrt( pow(start_x_pos-end_x_pos,2) + pow(start_y_pos-end_y_pos,2) + pow(start_z_pos-end_z_pos,2) );
          // just assign the start and end points of neutrons
          if(pdg==2112){
             _x->push_back(start_x_pos);
             _y->push_back(start_y_pos);
             _z->push_back(start_z_pos);
             _t->push_back(mother_time);   
             _x->push_back(end_x_pos);
             _y->push_back(end_y_pos);
             _z->push_back(end_z_pos);
             _t->push_back(mother_time+length*fsol);        
             return std::make_tuple(_x,_y,_z,_t);
          }
          if(!isfinite(length)){length=0.3;}
          double residual_range = length;
          double x_pos = start_x_pos;
          double y_pos = start_y_pos;
          double z_pos = start_z_pos;
          double t_pos = mother_time;
          double mass = 0;
          if (pdg == 13){ mass = 0.1057;}
          if (pdg == 2212){ mass = 0.9397933;}//Don't do neutrons, KE is not assigned well so just assume c
          if (pdg == 211){ mass = 0.13982067;}
          double KE =  pfeval.reco_startMomentum[particle][3] - mass;
          double v = fsol*1/sqrt( 1-pow(mass/(mass+KE),2) );
          if(!isfinite(v)){v=fsol;}
	  if (pdg==22 || pdg==2112 || pdg==11){ v = fsol;}
          double gamma=0;//used for extrapolating the line from start to end position
          while (residual_range>=0){
            _x->push_back(x_pos);
            _y->push_back(y_pos);
            _z->push_back(z_pos);
            _t->push_back(t_pos);
            t_pos += v*dx;
            double dedx = get_dE_dx_range(residual_range,pdg)/1000;
            KE = KE-dedx*dx;
            v = fsol*1/sqrt( 1-pow(mass/(mass+KE),2) );
	    if (!isfinite(v)){v=fsol;}
            if (pdg==22 || pdg==2112 || pdg==11){ v = fsol;}
            gamma+=(dx/length);
            x_pos = start_x_pos + gamma*(end_x_pos-start_x_pos);
            y_pos = start_y_pos + gamma*(end_y_pos-start_y_pos);
            z_pos = start_z_pos + gamma*(end_z_pos-start_z_pos);
            residual_range = length - sqrt( pow(start_x_pos-x_pos,2) + pow(start_y_pos-y_pos,2) + pow(start_z_pos-z_pos,2));
            if (!isfinite(residual_range)){break;}
	  }

          return std::make_tuple(_x,_y,_z,_t);
}

double LEEana::get_dE_dx_range(double R, int pdg){
    if (pdg==22 || pdg==11 || pdg==2112){ return 0; }
    double A = 8; 
    double b = -0.37; 
    if (pdg==2212){    
        A = 17;
        b = -0.42;
    }
    double dedx = A*pow(R,b);
    return dedx;
}

#endif
