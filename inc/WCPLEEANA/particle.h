#ifndef UBOONE_LEE_PARTICLE
#define UBOONE_LEE_PARTICLE

#include <algorithm>
#include <numeric>

namespace LEEana{
struct ParticleInfo{
    std::vector<double> *spacepoints_x;
    std::vector<double> *spacepoints_y;
    std::vector<double> *spacepoints_z;
    std::vector<double> *spacepoints_q;
    
    double spacepoints_q_0;
    double spacepoints_q_1;
    double spacepoints_q_2;
    double spacepoints_q_3;
    double spacepoints_q_4;
    double spacepoints_q_5;
    double spacepoints_q_6;
    double spacepoints_q_7;
    double spacepoints_q_8;
    double spacepoints_q_9;
    double spacepoints_q_10;
    double spacepoints_q_11;
    double spacepoints_q_12;
    double spacepoints_q_13;
    double spacepoints_q_14;
    double spacepoints_q_15;
    double spacepoints_q_16;
    double spacepoints_q_17;
    double spacepoints_q_18;
    double spacepoints_q_19;
    double spacepoints_q_20;
    double spacepoints_q_21;
    double spacepoints_q_22;
    double spacepoints_q_23;
    double spacepoints_q_24;

    double spacepoints_q_bck_0;
    double spacepoints_q_bck_1;
    double spacepoints_q_bck_2;
    double spacepoints_q_bck_3;
    double spacepoints_q_bck_4;
    double spacepoints_q_bck_5;
    double spacepoints_q_bck_6;
    double spacepoints_q_bck_7;
    double spacepoints_q_bck_8;
    double spacepoints_q_bck_9;
    double spacepoints_q_bck_10;
    double spacepoints_q_bck_11;
    double spacepoints_q_bck_12;
    double spacepoints_q_bck_13;
    double spacepoints_q_bck_14;
    double spacepoints_q_bck_15;
    double spacepoints_q_bck_16;
    double spacepoints_q_bck_17;
    double spacepoints_q_bck_18;
    double spacepoints_q_bck_19;
    double spacepoints_q_bck_20;
    double spacepoints_q_bck_21;
    double spacepoints_q_bck_22;
    double spacepoints_q_bck_23;
    double spacepoints_q_bck_24;

    double spacepoints_q_med;

    int flag_prim_mu;

    int flag_is_contained;

    int flag_has_daught;
    int flag_has_daught_p;
    int flag_has_daught_el;
    int flag_has_daught_pi;
    
    int reco_truthMatch_pdg;
    int reco_truthMatch_id;
    int reco_truthMatch_mother;
    double reco_truthMatch_energy;
    
    double reco_momentum_0;
    double reco_momentum_1;
    double reco_momentum_2;
    double reco_momentum_3;
    int reco_pdg;

    int reco_larpid_pdg;
    double reco_larpid_pidScore_el;
    double reco_larpid_pidScore_ph;
    double reco_larpid_pidScore_mu;
    double reco_larpid_pidScore_pr;
    double reco_larpid_pidScore_pi;
    double reco_larpid_proccess;

    int true_is_n_induced; 
    int reco_is_n_induced; 
    int reco_is_g_induced; 

    double dist_to_vtx; 
    double cos_theta; 
    double proximity; 

    double track_len;
    double direct_track_len;
    double track_len_ratio;
};

void create_particle(SpaceInfo& space_info, PFeval&, pfeval, ParticleInfo& particle_info, int index, double tolerance_sp=0.5, double tolerance=0.0001);
//void create_particles(SpaceInfo& space_info, PFeval&, pfeval, std::vector<ParticleInfo& particle_info>,double tolerance_sp=0.5, double tolerance=0.0001);
void reset_particle();
}


void LEEana::create_particle(SpaceInfo& space_info, PFeval&, pfeval, ParticleInfo& particle_info, int index, double tolerance_sp, double tolerance){

  reset_particle();

  // Loop over all particles in this event 
  for(int reco_part=0; reco_part<pfeval.reco_Ntrack; reco_part++){

    // Find the one associated with this id, skip if its a pseudo particle
    if(pfeval.reco_id[reco_part]!=pfeval.reco_id[index]) continue;
    if(pfeval.reco_pdg[reco_part]==2112 || pfeval.reco_pdg[reco_part]==22 || pfeval.reco_pdg[reco_part]==111) continue;
    if(pfeval.reco_truthMatch_pdg[event][reco_part]<=0) continue;

    double part_x = pfeval.reco_startXYZT[reco_part][0];
    double part_y = pfeval.reco_startXYZT[reco_part][1];
    double part_z = pfeval.reco_startXYZT[reco_part][2];
    double part_end_x = pfeval.reco_endXYZT[reco_part][0];
    double part_end_y = pfeval.reco_endXYZT[reco_part][1];
    double part_end_z = pfeval.reco_endXYZT[reco_part][2];

    // Save the spacepoints with the same id as the current particle
    std::vector<double> *temp_spacepoints_x;
    std::vector<double> *temp_spacepoints_y;
    std::vector<double> *temp_spacepoints_z;
    std::vector<double> *temp_spacepoints_q;    
    for(int sp=0; sp<sapce_info.spacepoints_id.size(); sp++){ 
      if(Trecchargeblob_spacepoints_real_cluster_id.at(sp)==id){ 
        temp_spacepoints_x->push_back(spacepoints_x->at(sp));
        temp_spacepoints_y->push_back(spacepoints_y->at(sp));
        temp_spacepoints_z->push_back(spacepoints_z->at(sp));
        temp_spacepoints_q->push_back(spacepoints_q->at(sp));
      }
    }
    int n_spacepoints = temp_spacepoints_x->size();
    if(n_spacepoints==0) continue;

    //Get the length of the proton by adding up the distance between each pair of spacepoints
    for(int sp=0; sp<n_spacepoints-1; sp++){
      double dx = temp_spacepoints_x->at(sp) - temp_spacepoints_x->at(sp+1);
      double dy = temp_spacepoints_y->at(sp) - temp_spacepoints_y->at(sp+1);
      double dz = temp_spacepoints_z->at(sp) - temp_spacepoints_z->at(sp+1);
      double dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
      particle_info.track_len+=dist;
    }


    // Check the start and end of the proton track to see if one matches the first spacepoint. 
    // This tells us if the spacepoints for this track were saved track start-to-end or track end-to-start
    // If neither matches, throw an error, I think this can happen sometimes for mouns (vertex sometimes gets redefined in the reco)? But have not seen if for protons
    if( !(temp_spacepoints_x->at(0)>part_x-tolerance_sp && temp_spacepoints_x->at(0)<part_x+tolerance_sp && temp_spacepoints_y->at(0)>part_y-tolerance_sp && temp_spacepoints_y->at(0)<part_y+tolerance_sp && temp_spacepoints_z->at(0)>part_z-tolerance_sp && temp_spacepoints_z->at(0)<part_z+tolerance_sp) ){
      if(!(temp_spacepoints_x->back()>part_x-tolerance_sp && temp_spacepoints_x->back()<part_x+tolerance_sp && temp_spacepoints_y->back()>part_y-tolerance_sp && temp_spacepoints_y->back()<part_y+tolerance_sp && temp_spacepoints_z->back()>part_z-tolerance_sp && temp_spacepoints_z->back()<part_z+tolerance_sp) ){
        particle_info.spacepoints_x = temp_spacepoints_x;
        particle_info.spacepoints_y = temp_spacepoints_y;
        particle_info.spacepoints_z = temp_spacepoints_z;
        particle_info.spacepoints_q = temp_spacepoints_q; 
      }else{
        std::vector<double> *rev_temp_spacepoints_x;
        std::vector<double> *rev_temp_spacepoints_y;
        std::vector<double> *rev_temp_spacepoints_z;
        std::vector<double> *rev_temp_spacepoints_q;
        for(int sp=0; sp<n_spacepoints; sp++){
          rev_temp_spacepoints_x->push_back(temp_spacepoints_x->at(n_spacepoints-1-sp));
          rev_temp_spacepoints_y->push_back(temp_spacepoints_y->at(n_spacepoints-1-sp));
          rev_temp_spacepoints_z->push_back(temp_spacepoints_z->at(n_spacepoints-1-sp));
          rev_temp_spacepoints_q->push_back(temp_spacepoints_q->at(n_spacepoints-1-sp));
        }
        particle_info.spacepoints_x = rev_temp_spacepoints_x;
        particle_info.spacepoints_y = rev_temp_spacepoints_y;
        particle_info.spacepoints_z = rev_temp_spacepoints_z;
        particle_info.spacepoints_q = rev_temp_spacepoints_q;
      }else{ 
        particle_info.spacepoints_x = temp_spacepoints_x;
        particle_info.spacepoints_y = temp_spacepoints_y;
        particle_info.spacepoints_z = temp_spacepoints_z;
        particle_info.spacepoints_q = temp_spacepoints_q;
      }
    }

    // Save the median dqdx
    std::vector<double> *sorted_spacepoints_q; 
    *sorted_spacepoints_q = *(particle_info.spacepoints_q);
    std::sort(sorted_spacepoints_q->begin(), sorted_spacepoints_q->end());
    size_t size = sorted_spacepoints_q.size();
    if (size % 2 == 0) {
      particle_info.spacepoints_q_med = (sorted_spacepoints_q->at(size / 2 - 1) + sorted_spacepoints_q->at(size / 2)) / 2.0;
    }else{
      particle_info.spacepoints_q_med  = sorted_spacepoints_q->at(size / 2);
    }

    // Find the mother that was larpid matched and add some extra info on it
    int temp_truth_mother_id=-1;
    for(int truth_part=0; truth_part<pfeval.truth_Ntrack; truth_part++){         
      if(pfeval.truth_id[truth_part]==pfeval.reco_truthMatch_id[reco_part]){
        double mass = 0;
         if(pfeval.truth_pdg[truth_part]==13) mass = 105.7;
         if(pfeval.truth_pdg[truth_part]==211) mass = 138;
         if(pfeval.truth_pdg[truth_part]==2212 || pfeval.truth_pdg[truth_part]==2112) mass = 938;
           space_info.reco_truth_match_mother = pfeval.truth_mother[truth_part];
           temp_truth_mother_id=pfeval.truth_mother[truth_part];
           space_info.reco_truth_match_energy = pfeval.truth_momentum[truth_part][3]*1000-mass;
           break;
      }
    }
    space_info.reco_truth_match_pdg = pfeval.reco_truthMatch_pdg[reco_part];
    space_info.reco_truth_match_id = pfeval.reco_truthMatch_id[reco_part];
    int truth_is_n_induced=0;
    for(int truth_mother_part=0; truth_mother_part<pfeval.truth_Ntrack; truth_mother_part++){
      if(pfeval.truth_id[truth_mother_part]!=temp_truth_mother_id) continue;
      if(pfeval.truth_pdg[truth_mother_part]==2112) truth_is_n_induced = 1;
    }

    // Add some daughter information
    particle_info.flag_has_daught=0;
    particle_info.flag_has_daught_p=0;
    particle_info.flag_has_daught_el=0;
    particle_info.flag_has_daught_pi=0;
    for(int reco_daught_part=0; reco_daught_part<reco_Ntrack; reco_daught_part++){
      if(pfeval.reco_mother[reco_daught_part]!=pfeval.reco_id[reco_part]) continue;
      particle_info.flag_has_daught+=1;
      if(pfeval.reco_pdg[reco_daught_part]==2212) particle_info.flag_has_daught_p+=1;
      if(pfeval.reco_pdg[reco_daught_part]==11) particle_info.flag_has_daught_el+=1;
      if(pfeval.reco_pdg[reco_daught_part]==211) particle_info.flag_has_daught_pi+=1;
    }
	    
    particle_info.reco_is_n_induced = 0;
    particle_info.reco_is_g_induced = 0;
    for(int reco_mother_part=0; reco_mother_part<reco_Ntrack; reco_mother_part++){
      if(pfeval.reco_id[reco_mother_part]!=pfeval.reco_mother[reco_part]) continue;
      if(pfeval.reco_pdg[reco_mother_part]==2112) particle_info.reco_is_n_induced=1;
      else if(pfeval.reco_pdg[reco_mother_part]==22) particle_info.reco_is_g_induced=1;
    }


    // Save more general variables about the particel
    particle_info.dist_to_vtx = sqrt(pow(part_x-pfeval.reco_nuvtxX,2)+pow(part_y-pfeval.reco_nuvtxY,2)+pow(part_z-pfeval.reco_nuvtxZ,2));

    // Momentum and pid
    double mass = 0;
    if(pfeval.reco_pdg[reco_part]==13) mass = 0.1057;
    if(pfeval.reco_pdg[reco_part]==211) mass = 0.138;
    if(pfeval.reco_pdg[reco_part]==2212) mass = 0.938;
    particle_info.reco_momentum_0 = pfeval.reco_momentum[reco_part][0];
    particle_info.reco_momentum_1 = pfeval.reco_momentum[reco_part][1];
    particle_info.reco_momentum_2 = pfeval.reco_momentum[reco_part][2];
    particle_info.reco_momentum_3 = pfeval.reco_momentum[reco_part][3]-mass;
    particle_info.reco_pdg_list = pfeval.reco_pdg[reco_part];

    // Check if its the leading muon
    if(pfeval.reco_pdg[reco_part]==13 && pfeval.reco_startMomentum[reco_part][3]>pfeval.reco_muonMomentum[3]-tolerance && pfeval.reco_startMomentum[reco_part][3]<pfeval.reco_muonMomentum+tolerance) particle_info.flag_prim_mu = 1;
    else if(pfeval.reco_pdg[reco_part]==13) particle_info.flag_prim_mu = 0;
    else particle_info.flag_prim_mu = -1;

    // Scattering angle
    particle_info.cos_theta = pfeval.reco_startMomentum[reco_part][2] / sqrt( (pfeval.reco_startMomentum[reco_part][0],2) + pow(pfeval.reco_startMomentum[reco_part][1],2) + pow(pfeval.reco_startMomentum[reco_part][2],2) );

    // Check the containment
    particle_info.flag_contained_part = 1
    if(part_end_x < 3 || part_end_x > 250) particle_infoflag_contained_part = 0;
    else if(part_end_y < -113 || part_end_y > 113) particle_infoflag_contained_part = 0;
    else if(part_end_z < 3 || part_end_z > 1035) particle_infoflag_contained_part = 0;

    // Save the larpid vars
    particle_info.reco_larpid_pdg_part = pfeval.reco_larpid_pdg[reco_part];
    particle_info.reco_larpid_pidScore_el_part = pfeval.reco_larpid_pidScore_el[reco_part];
    particle_info.reco_larpid_pidScore_ph_part = pfeval.reco_larpid_pidScore_ph[reco_part];
    particle_info.reco_larpid_pidScore_mu_part = pfeval.reco_larpid_pidScore_mu[reco_part];
    particle_info.reco_larpid_pidScore_pr_part = pfeval.reco_larpid_pidScore_pr[reco_part];
    particle_info.reco_larpid_pidScore_pi_part = pfeval.reco_larpid_pidScore_pi[reco_part];
    particle_info.reco_larpid_proccess_part = pfeval.reco_larpid_proccess[reco_part];


    // Save the individual spacepoints for easier use in the BDT
    particle_info.spacepoints_q_0 = particle_info.temp_spacepoints_q->at(0);
    particle_info.spacepoints_q_bck_0 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-0);
    if(n_spacepoints>1){
      particle_info.spacepoints_q_1 = particle_info.temp_spacepoints_q->at(1);
      particle_info.spacepoints_q_bck_1 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-1);
    }else{
      particle_info.spacepoints_q_1 = -999;
      particle_info.spacepoints_q_bck_1 = -999;
    }
    if(n_spacepoints>2){
      particle_info.spacepoints_q_2 = particle_info.temp_spacepoints_q->at(2);
      particle_info.spacepoints_q_bck_2 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-2);
    }else{
      particle_info.spacepoints_q_2 = -999;
      particle_info.spacepoints_q_bck_2 = -999;
    }
    if(n_spacepoints>3){
      particle_info.spacepoints_q_3 = particle_info.temp_spacepoints_q->at(3);
      particle_info.spacepoints_q_bck_3 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-3);
    }else{
      particle_info.spacepoints_q_3 = -999;
      particle_info.spacepoints_q_bck_3 = -999;
    }
    if(n_spacepoints>4){
      particle_info.spacepoints_q_4 = particle_info.temp_spacepoints_q->at(4);
      particle_info.spacepoints_q_bck_4 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-4);
    }else{
      particle_info.spacepoints_q_4 = -999;
      particle_info.spacepoints_q_bck_4 = -999;
    }
    if(n_spacepoints>5){
      particle_info.spacepoints_q_5 = particle_info.temp_spacepoints_q->at(5);
      particle_info.spacepoints_q_bck_5 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-5);
    }else{
      particle_info.spacepoints_q_5 = -999;
      particle_info.spacepoints_q_bck_5 = -999;
    }
    if(n_spacepoints>6){
      particle_info.spacepoints_q_6 = particle_info.temp_spacepoints_q->at(6);
      particle_info.spacepoints_q_bck_6 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-6);
    }else{
      particle_info.spacepoints_q_6 = -999;
      particle_info.spacepoints_q_bck_6 = -999;
    }
    if(n_spacepoints>7){
      particle_info.spacepoints_q_7 = particle_info.temp_spacepoints_q->at(7);
      particle_info.spacepoints_q_bck_7 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-7);
    }else{
      particle_info.spacepoints_q_7 = -999;
      particle_info.spacepoints_q_bck_7 = -999;
    }
    if(n_spacepoints>8){
      particle_info.spacepoints_q_8 = particle_info.temp_spacepoints_q->at(8);
      particle_info.spacepoints_q_bck_8 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-8);
    }else{
      particle_info.spacepoints_q_8 = -999;
      particle_info.spacepoints_q_bck_8 = -999;
    }
    if(n_spacepoints>9){
      particle_info.spacepoints_q_9 = particle_info.temp_spacepoints_q->at(9);
      particle_info.spacepoints_q_bck_9 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-9);
    }else{
      particle_info.spacepoints_q_9 = -999;
      particle_info.spacepoints_q_bck_9 = -999;
    }
    if(n_spacepoints>10){
      particle_info.spacepoints_q_10 = particle_info.temp_spacepoints_q->at(10);
      particle_info.spacepoints_q_bck_10 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-10);
    }else{
      particle_info.spacepoints_q_10 = -999;
      particle_info.spacepoints_q_bck_10 = -999;
    }
    if(n_spacepoints>11){
      particle_info.spacepoints_q_11 = particle_info.temp_spacepoints_q->at(11);
      particle_info.spacepoints_q_bck_11 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-11);
    }else{
      particle_info.spacepoints_q_11 = -999;
      particle_info.spacepoints_q_bck_11 = -999;
    }
    if(n_spacepoints>12){
      particle_info.spacepoints_q_12 = particle_info.temp_spacepoints_q->at(12);
      particle_info.spacepoints_q_bck_12 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-12);
    }else{
      particle_info.spacepoints_q_12 = -999;
      particle_info.spacepoints_q_bck_12 = -999;
    }
    if(n_spacepoints>13){
      particle_info.spacepoints_q_13 = particle_info.temp_spacepoints_q->at(13);
      particle_info.spacepoints_q_bck_13 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-13);
    }else{
      particle_info.spacepoints_q_13 = -999;
      particle_info.spacepoints_q_bck_13 = -999;
    }
    if(n_spacepoints>14){
      particle_info.spacepoints_q_14 = particle_info.temp_spacepoints_q->at(14);
      particle_info.spacepoints_q_bck_14 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-14);
    }else{
      particle_info.spacepoints_q_14 = -999;
      particle_info.spacepoints_q_bck_14 = -999;
    }
    if(n_spacepoints>15){
      particle_info.spacepoints_q_15 = particle_info.temp_spacepoints_q->at(15);
      particle_info.spacepoints_q_bck_15 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-15);
    }else{
      particle_info.spacepoints_q_15 = -999;
      particle_info.spacepoints_q_bck_15 = -999;
    }
    if(n_spacepoints>16){
      particle_info.spacepoints_q_16 = particle_info.temp_spacepoints_q->at(16);
      particle_info.spacepoints_q_bck_16 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-16);
    }else{
      particle_info.spacepoints_q_16 = -999;
      particle_info.spacepoints_q_bck_16 = -999;
    }
    if(n_spacepoints>17){
      particle_info.spacepoints_q_17 = particle_info.temp_spacepoints_q->at(17);
      particle_info.spacepoints_q_bck_17 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-17);
    }else{
      particle_info.spacepoints_q_17 = -999;
      particle_info.spacepoints_q_bck_17 = -999;
    }
    if(n_spacepoints>18){
      particle_info.spacepoints_q_18 = particle_info.temp_spacepoints_q->at(18);
      particle_info.spacepoints_q_bck_18 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-18);
    }else{
      particle_info.spacepoints_q_18 = -999;
      particle_info.spacepoints_q_bck_18 = -999;
    }
    if(n_spacepoints>19){
      particle_info.spacepoints_q_19 = particle_info.temp_spacepoints_q->at(19);
      particle_info.spacepoints_q_bck_19 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-19);
    }else{
      particle_info.spacepoints_q_19 = -999;
      particle_info.spacepoints_q_bck_19 = -999;
    }
    if(n_spacepoints>20){
      particle_info.spacepoints_q_20 = particle_info.temp_spacepoints_q->at(20);
      particle_info.spacepoints_q_bck_20 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-20);
    }else{
      particle_info.spacepoints_q_20 = -999;
      particle_info.spacepoints_q_bck_20 = -999;
    }
    if(n_spacepoints>21){
      particle_info.spacepoints_q_21 = particle_info.temp_spacepoints_q->at(21);
      particle_info.spacepoints_q_bck_21 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-21);
    }else{
      particle_info.spacepoints_q_21 = -999;
      particle_info.spacepoints_q_bck_21 = -999;
    }
    if(n_spacepoints>22){
      particle_info.spacepoints_q_22 = particle_info.temp_spacepoints_q->at(22);
      particle_info.spacepoints_q_bck_22 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-22);
    }else{
      particle_info.spacepoints_q_22 = -999;
      particle_info.spacepoints_q_bck_22 = -999;
    }
    if(n_spacepoints>23){
      particle_info.spacepoints_q_23 = particle_info.temp_spacepoints_q->at(23);
      particle_info.spacepoints_q_bck_23 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-23);
    }else{
      particle_info.spacepoints_q_23 = -999;
      particle_info.spacepoints_q_bck_23 = -999;
    }
    if(n_spacepoints>24){
      particle_info.spacepoints_q_24 = particle_info.temp_spacepoints_q->at(24);
      particle_info.spacepoints_q_bck_24 = particle_info.temp_spacepoints_q->at(n_spacepoints-1-24);
    }else{
      particle_info.spacepoints_q_24 = -999;
      particle_info.spacepoints_q_bck_24 = -999;
    }

  }
}

void LEEana::reset_particle{
    std::vector<double> *spacepoints_x->clear();
    std::vector<double> *spacepoints_y->clear();
    std::vector<double> *spacepoints_z->clear();
    std::vector<double> *spacepoints_q->clear();

    double spacepoints_q_0=-999;
    double spacepoints_q_1=-999;
    double spacepoints_q_2=-999;
    double spacepoints_q_3=-999;
    double spacepoints_q_4=-999;
    double spacepoints_q_5=-999;
    double spacepoints_q_6=-999;
    double spacepoints_q_7=-999;
    double spacepoints_q_8=-999;
    double spacepoints_q_9=-999;
    double spacepoints_q_10=-999;
    double spacepoints_q_11=-999;
    double spacepoints_q_12=-999;
    double spacepoints_q_13=-999;
    double spacepoints_q_14=-999;
    double spacepoints_q_15=-999;
    double spacepoints_q_16=-999;
    double spacepoints_q_17=-999;
    double spacepoints_q_18=-999;
    double spacepoints_q_19=-999;
    double spacepoints_q_20=-999;
    double spacepoints_q_21=-999;
    double spacepoints_q_22=-999;
    double spacepoints_q_23=-999;
    double spacepoints_q_24=-999;

    double spacepoints_q_bck_0=-999;
    double spacepoints_q_bck_1=-999;
    double spacepoints_q_bck_2=-999;
    double spacepoints_q_bck_3=-999;
    double spacepoints_q_bck_4=-999;
    double spacepoints_q_bck_5=-999;
    double spacepoints_q_bck_6=-999;
    double spacepoints_q_bck_7=-999;
    double spacepoints_q_bck_8=-999;
    double spacepoints_q_bck_9=-999;
    double spacepoints_q_bck_10=-999;
    double spacepoints_q_bck_11=-999;
    double spacepoints_q_bck_12=-999;
    double spacepoints_q_bck_13=-999;
    double spacepoints_q_bck_14=-999;
    double spacepoints_q_bck_15=-999;
    double spacepoints_q_bck_16=-999;
    double spacepoints_q_bck_17=-999;
    double spacepoints_q_bck_18=-999;
    double spacepoints_q_bck_19=-999;
    double spacepoints_q_bck_20=-999;
    double spacepoints_q_bck_21=-999;
    double spacepoints_q_bck_22=-999;
    double spacepoints_q_bck_23=-999;
    double spacepoints_q_bck_24=-999;

    double spacepoints_q_med=-999;

    int flag_prim_mu=-999;

    int flag_is_contained=-999;

    int flag_has_daught=-999;
    int flag_has_daught_p=-999;
    int flag_has_daught_el=-999;
    int flag_has_daught_pi=-999;

    int reco_truthMatch_pdg=-999;
    int reco_truthMatch_id=-999;
    int reco_truthMatch_mother=-999;
    double reco_truthMatch_energy=-999;

    double reco_momentum_0=-999;
    double reco_momentum_1=-999;
    double reco_momentum_2=-999;
    double reco_momentum_3=-999;
    int reco_pdg=-999;

    int reco_larpid_pdg=-999;
    double reco_larpid_pidScore_el=-999;
    double reco_larpid_pidScore_ph=-999;
    double reco_larpid_pidScore_mu=-999;
    double reco_larpid_pidScore_pr=-999;
    double reco_larpid_pidScore_pi=-999;
    double reco_larpid_proccess=-999;

    int true_is_n_induced=-999;
    int reco_is_n_induced=-999;
    int reco_is_g_induced=-999;

    double dist_to_vtx=-999;
    double cos_theta=-999;
    double proximity=-999;

    double track_len=-999;
    double direct_track_len=-999;
    double track_len_ratio=-999;
}
