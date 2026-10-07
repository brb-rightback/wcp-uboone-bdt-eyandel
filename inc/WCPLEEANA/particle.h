#ifndef UBOONE_LEE_PARTICLE
#define UBOONE_LEE_PARTICLE

#include <algorithm>
#include <numeric>
#include <deque>
#include <map>
#include <set>
#include <tuple>
#include <array>
#include <cmath>
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/space.h"

namespace LEEana{
struct ParticleInfo{
    std::vector<float> *spacepoints_x;
    std::vector<float> *spacepoints_y;
    std::vector<float> *spacepoints_z;
    std::vector<float> *spacepoints_q;

    float spacepoints_q_0;
    float spacepoints_q_1;
    float spacepoints_q_2;
    float spacepoints_q_3;
    float spacepoints_q_4;
    float spacepoints_q_5;
    float spacepoints_q_6;
    float spacepoints_q_7;
    float spacepoints_q_8;
    float spacepoints_q_9;
    float spacepoints_q_10;
    float spacepoints_q_11;
    float spacepoints_q_12;
    float spacepoints_q_13;
    float spacepoints_q_14;
    float spacepoints_q_15;
    float spacepoints_q_16;
    float spacepoints_q_17;
    float spacepoints_q_18;
    float spacepoints_q_19;
    float spacepoints_q_20;
    float spacepoints_q_21;
    float spacepoints_q_22;
    float spacepoints_q_23;
    float spacepoints_q_24;

    float spacepoints_q_bck_0;
    float spacepoints_q_bck_1;
    float spacepoints_q_bck_2;
    float spacepoints_q_bck_3;
    float spacepoints_q_bck_4;
    float spacepoints_q_bck_5;
    float spacepoints_q_bck_6;
    float spacepoints_q_bck_7;
    float spacepoints_q_bck_8;
    float spacepoints_q_bck_9;
    float spacepoints_q_bck_10;
    float spacepoints_q_bck_11;
    float spacepoints_q_bck_12;
    float spacepoints_q_bck_13;
    float spacepoints_q_bck_14;
    float spacepoints_q_bck_15;
    float spacepoints_q_bck_16;
    float spacepoints_q_bck_17;
    float spacepoints_q_bck_18;
    float spacepoints_q_bck_19;
    float spacepoints_q_bck_20;
    float spacepoints_q_bck_21;
    float spacepoints_q_bck_22;
    float spacepoints_q_bck_23;
    float spacepoints_q_bck_24;

    float spacepoints_q_med;

    float flag_prim_mu;

    float flag_is_contained;

    float flag_has_daught;
    float flag_has_daught_p;
    float flag_has_daught_el;
    float flag_has_daught_pi;

    float reco_truthMatch_pdg;
    float reco_truthMatch_id;
    float reco_truthMatch_mother;
    float reco_truthMatch_energy;

    float reco_momentum_0;
    float reco_momentum_1;
    float reco_momentum_2;
    float reco_momentum_3;
    float reco_pdg;

    float reco_larpid_pdg;
    float reco_larpid_pidScore_el;
    float reco_larpid_pidScore_ph;
    float reco_larpid_pidScore_mu;
    float reco_larpid_pidScore_pr;
    float reco_larpid_pidScore_pi;
    float reco_larpid_proccess;

    float true_is_n_induced;
    float reco_is_n_induced;
    float reco_is_g_induced;

    float dist_to_vtx;
    float cos_theta;
    float proximity;

    float track_len;
    float direct_track_len;
    float track_len_ratio;

    // Secondary segment: the longest stored block of spacepoints that is not part of the main (start -> end) sequence,
    // ordered from its end nearest to the main sequence. (q+10000)*10, -999 if there is no secondary segment / fewer points.
    float spacepoints_q_sec_0;
    float spacepoints_q_sec_1;
    float spacepoints_q_sec_2;
    float spacepoints_q_sec_3;
    float spacepoints_q_sec_4;

    float spacepoints_q_sec_bck_0;
    float spacepoints_q_sec_bck_1;
    float spacepoints_q_sec_bck_2;
    float spacepoints_q_sec_bck_3;
    float spacepoints_q_sec_bck_4;

    // Summary of the pieces the particle's spacepoints are stored in (see build_particle_segments)
    float n_segments;             // 1 (main sequence) + number of blocks not in the main sequence
    float seg_total_len;          // path length of the main sequence + those blocks [cm]
    float seg_total_q;            // summed raw charge of the main sequence + those blocks
    float sec_seg_len;            // path length of the secondary segment [cm], -999 if none
    float sec_seg_direct_len;     // start-to-end distance of the secondary segment [cm], -999 if none
    float sec_seg_dist_to_main;   // closest distance between the secondary segment and the main sequence [cm], -999 if none
    float flag_start_end_swapped; // 1 if the reco start/end were swapped for the ordering (diagnostic, not meant as a BDT input)
};

// Result of build_particle_segments: indices into the particle's spacepoints (in stored order)
struct ParticleSegments{
    std::vector<int> main;        // main sequence, ordered from the particle start to its end
    std::vector<int> sec;         // secondary segment, ordered from the end nearest to the main sequence (empty if none)
    int n_segments;
    double seg_total_len;
    double seg_total_q;
    double sec_seg_len;
    double sec_seg_direct_len;
    double sec_seg_dist_to_main;
    bool flipped;
};

ParticleSegments build_particle_segments(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& q,
                                         const double start[3], const double end[3], const double* ref,
                                         double jump=2.0, double junction=2.5, double attach=2.0, double anchor_tol=0.05, double flip_tol=2.0);

void create_particle(SpaceInfo& space_info, PFevalInfo& pfeval, ParticleInfo& particle_info, int index, bool flag_data, double tolerance_sp=0.5, double tolerance=0.0001);
void reset_particle(ParticleInfo& particle_info);
void print_particle(ParticleInfo& particle_info);
}


void LEEana::create_particle(SpaceInfo& space_info, PFevalInfo& pfeval, ParticleInfo& particle_info, int index, bool flag_data, double tolerance_sp, double tolerance){

  //std::cout<<"Clearing particle"<<std::endl;
  reset_particle(particle_info);

  //std::cout<<"Searching for particle"<<std::endl;
  // Loop over all particles in this event 
  for(int reco_part=0; reco_part<pfeval.reco_Ntrack; reco_part++){

    // Find the one associated with this id, skip if its a pseudo particle
    if(pfeval.reco_id[reco_part]!=pfeval.reco_id[index]) continue;
    //std::cout<<"Found ID"<<std::endl;
    if(pfeval.reco_pdg[reco_part]==2112 || pfeval.reco_pdg[reco_part]==22 || pfeval.reco_pdg[reco_part]==111) continue;
    //if(pfeval.reco_truthMatch_pdg[reco_part]<=0 && !flag_data) continue;
    //std::cout<<"Real particle and found amtch"<<std::endl;

    double part_x = pfeval.reco_startXYZT[reco_part][0];
    double part_y = pfeval.reco_startXYZT[reco_part][1];
    double part_z = pfeval.reco_startXYZT[reco_part][2];
    double part_end_x = pfeval.reco_endXYZT[reco_part][0];
    double part_end_y = pfeval.reco_endXYZT[reco_part][1];
    double part_end_z = pfeval.reco_endXYZT[reco_part][2];

    //std::cout<<"Setting spacepoints"<<std::endl;
    // Save the spacepoints with the same id as the current particle
    std::vector<float> temp_spacepoints_x;
    std::vector<float> temp_spacepoints_y;
    std::vector<float> temp_spacepoints_z;
    std::vector<float> temp_spacepoints_q;
    //std::cout<<"Number of spacepoints: "<<space_info.Trecchargeblob_spacepoints_real_cluster_id->size()<<std::endl;    
    for(int sp=0; sp<space_info.Trecchargeblob_spacepoints_real_cluster_id->size(); sp++){
      //std::cout<<"Checking sp "<<sp<<std::endl; 
      if(space_info.Trecchargeblob_spacepoints_real_cluster_id->at(sp)==pfeval.reco_id[reco_part]){
        //std::cout<<"Adding sp "<<sp<<" x="<<space_info.Trecchargeblob_spacepoints_x->at(sp)<<" q="<<space_info.Trecchargeblob_spacepoints_q->at(sp)<<std::endl; 
        temp_spacepoints_x.push_back(space_info.Trecchargeblob_spacepoints_x->at(sp));
        temp_spacepoints_y.push_back(space_info.Trecchargeblob_spacepoints_y->at(sp));
        temp_spacepoints_z.push_back(space_info.Trecchargeblob_spacepoints_z->at(sp));
        temp_spacepoints_q.push_back(space_info.Trecchargeblob_spacepoints_q->at(sp));
      }
    }
    int n_spacepoints = temp_spacepoints_x.size();
    if(n_spacepoints==0) continue;

    //std::cout<<"Computing length"<<std::endl;
    //Get the length of the proton by adding up the distance between each pair of spacepoints, skipping jumps >= 2 cm.
    //Particles can be stored as separate pieces (e.g. the start/end points first, then the trajectory).
    particle_info.track_len=0;
    std::vector<bool> on_trajectory(n_spacepoints, false);
    for(int sp=0; sp<n_spacepoints-1; sp++){
      double dx = temp_spacepoints_x.at(sp) - temp_spacepoints_x.at(sp+1);
      double dy = temp_spacepoints_y.at(sp) - temp_spacepoints_y.at(sp+1);
      double dz = temp_spacepoints_z.at(sp) - temp_spacepoints_z.at(sp+1);
      double dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
      if(dist<2.0){
        particle_info.track_len+=dist;
        on_trajectory.at(sp) = true;
        on_trajectory.at(sp+1) = true;
      }
    }

    //Connect the reco start and end to the nearest trajectory point if it is within 2 cm.
    //If no steps were kept (e.g. only the start and end points are stored), fall back to the straight distance.
    bool has_trajectory = false;
    for(int sp=0; sp<n_spacepoints; sp++) if(on_trajectory.at(sp)) has_trajectory = true;
    if(has_trajectory){
      for(int ep=0; ep<2; ep++){
        double ex = (ep==0) ? part_x : part_end_x;
        double ey = (ep==0) ? part_y : part_end_y;
        double ez = (ep==0) ? part_z : part_end_z;
        double dmin = 1e9;
        for(int sp=0; sp<n_spacepoints; sp++){
          if(!on_trajectory.at(sp)) continue;
          double d = sqrt(pow(ex-temp_spacepoints_x.at(sp),2)+pow(ey-temp_spacepoints_y.at(sp),2)+pow(ez-temp_spacepoints_z.at(sp),2));
          if(d<dmin) dmin = d;
        }
        if(dmin<2.0) particle_info.track_len+=dmin;
      }
    }else{
      particle_info.track_len = sqrt(pow(part_x-part_end_x,2)+pow(part_y-part_end_y,2)+pow(part_z-part_end_z,2));
    }


    // Order the spacepoints from the particle start to its end (see build_particle_segments). WireCell stores a particle's
    // spacepoints as separate pieces (the start/end and other node points first, then blocks of trajectory points, each
    // block in its own direction), so the stored order is not the order along the particle.
    // Reference point for the start/end orientation: the neutrino vertex for primaries, the nearer of the mother's reco
    // start/end for secondaries. If the reco end is closer to it than the reco start, WireCell's direction is reversed.
    double part_start[3] = {part_x, part_y, part_z};
    double part_end[3] = {part_end_x, part_end_y, part_end_z};
    double ref_point[3] = {0, 0, 0};
    bool has_ref = false;
    if(pfeval.reco_mother[reco_part]==0){
      ref_point[0] = pfeval.reco_nuvtxX; ref_point[1] = pfeval.reco_nuvtxY; ref_point[2] = pfeval.reco_nuvtxZ;
      has_ref = true;
    }else{
      for(int mother_part=0; mother_part<pfeval.reco_Ntrack; mother_part++){
        if(pfeval.reco_id[mother_part]!=pfeval.reco_mother[reco_part]) continue;
        double d_mother_end = sqrt(pow(part_x-pfeval.reco_endXYZT[mother_part][0],2)+pow(part_y-pfeval.reco_endXYZT[mother_part][1],2)+pow(part_z-pfeval.reco_endXYZT[mother_part][2],2));
        double d_mother_start = sqrt(pow(part_x-pfeval.reco_startXYZT[mother_part][0],2)+pow(part_y-pfeval.reco_startXYZT[mother_part][1],2)+pow(part_z-pfeval.reco_startXYZT[mother_part][2],2));
        for(int c=0; c<3; c++) ref_point[c] = (d_mother_end<=d_mother_start) ? pfeval.reco_endXYZT[mother_part][c] : pfeval.reco_startXYZT[mother_part][c];
        has_ref = true;
        break;
      }
    }
    std::vector<double> dbl_x(temp_spacepoints_x.begin(), temp_spacepoints_x.end());
    std::vector<double> dbl_y(temp_spacepoints_y.begin(), temp_spacepoints_y.end());
    std::vector<double> dbl_z(temp_spacepoints_z.begin(), temp_spacepoints_z.end());
    std::vector<double> dbl_q(temp_spacepoints_q.begin(), temp_spacepoints_q.end());
    ParticleSegments segments = build_particle_segments(dbl_x, dbl_y, dbl_z, dbl_q, part_start, part_end, has_ref ? ref_point : nullptr);

    // The main sequence fills the spacepoint vectors, and so the spacepoints_q_* / _bck_* / _med variables below
    particle_info.spacepoints_x->clear();
    particle_info.spacepoints_y->clear();
    particle_info.spacepoints_z->clear();
    particle_info.spacepoints_q->clear();
    for(size_t k=0; k<segments.main.size(); k++){
      particle_info.spacepoints_x->push_back(temp_spacepoints_x.at(segments.main.at(k)));
      particle_info.spacepoints_y->push_back(temp_spacepoints_y.at(segments.main.at(k)));
      particle_info.spacepoints_z->push_back(temp_spacepoints_z.at(segments.main.at(k)));
      particle_info.spacepoints_q->push_back(temp_spacepoints_q.at(segments.main.at(k)));
    }

    // Secondary segment: first and last five charges, ordered from the end nearest to the main sequence
    float* q_sec[5] = {&particle_info.spacepoints_q_sec_0, &particle_info.spacepoints_q_sec_1, &particle_info.spacepoints_q_sec_2,
                       &particle_info.spacepoints_q_sec_3, &particle_info.spacepoints_q_sec_4};
    float* q_sec_bck[5] = {&particle_info.spacepoints_q_sec_bck_0, &particle_info.spacepoints_q_sec_bck_1, &particle_info.spacepoints_q_sec_bck_2,
                           &particle_info.spacepoints_q_sec_bck_3, &particle_info.spacepoints_q_sec_bck_4};
    int n_sec = segments.sec.size();
    for(int k=0; k<5; k++){
      *q_sec[k] = (k<n_sec) ? (temp_spacepoints_q.at(segments.sec.at(k))+10000)*10 : -999;
      *q_sec_bck[k] = (k<n_sec) ? (temp_spacepoints_q.at(segments.sec.at(n_sec-1-k))+10000)*10 : -999;
    }

    // Summary of the pieces
    particle_info.n_segments = segments.n_segments;
    particle_info.seg_total_len = segments.seg_total_len;
    particle_info.seg_total_q = segments.seg_total_q;
    particle_info.sec_seg_len = segments.sec_seg_len;
    particle_info.sec_seg_direct_len = segments.sec_seg_direct_len;
    particle_info.sec_seg_dist_to_main = segments.sec_seg_dist_to_main;
    particle_info.flag_start_end_swapped = segments.flipped ? 1 : 0;

    // Save the median dqdx
    std::vector<float> sorted_spacepoints_q = *(particle_info.spacepoints_q);
    std::sort(sorted_spacepoints_q.begin(), sorted_spacepoints_q.end());
    size_t size = sorted_spacepoints_q.size();
    if (size % 2 == 0) {
      particle_info.spacepoints_q_med = (((sorted_spacepoints_q.at(size / 2 - 1) + sorted_spacepoints_q.at(size / 2)) / 2.0)+10000)*10;
    }else{
      particle_info.spacepoints_q_med  = ((sorted_spacepoints_q.at(size / 2))+10000)*10;
    }

    // Find the mother that was larpid matched and add some extra info on it
    int temp_truth_mother_id=-1;
    for(int truth_part=0; truth_part<pfeval.truth_Ntrack; truth_part++){         
      if(pfeval.truth_id[truth_part]==pfeval.reco_truthMatch_id[reco_part]){
        double mass = 0;
         if(pfeval.truth_pdg[truth_part]==13) mass = 105.7;
         if(pfeval.truth_pdg[truth_part]==211) mass = 138;
         if(pfeval.truth_pdg[truth_part]==2212 || pfeval.truth_pdg[truth_part]==2112) mass = 938;
           particle_info.reco_truthMatch_mother = pfeval.truth_mother[truth_part];
           temp_truth_mother_id=pfeval.truth_mother[truth_part];
           particle_info.reco_truthMatch_energy = pfeval.truth_startMomentum[truth_part][3]*1000-mass;
           break;
      }
    }
    particle_info.reco_truthMatch_pdg = pfeval.reco_truthMatch_pdg[reco_part];
    particle_info.reco_truthMatch_id = pfeval.reco_truthMatch_id[reco_part];
    particle_info.true_is_n_induced=0;
    for(int truth_mother_part=0; truth_mother_part<pfeval.truth_Ntrack; truth_mother_part++){
      if(pfeval.truth_id[truth_mother_part]!=temp_truth_mother_id) continue;
      if(pfeval.truth_pdg[truth_mother_part]==2112) particle_info.true_is_n_induced = 1;
    }

    // Add some daughter information
    particle_info.flag_has_daught=0;
    particle_info.flag_has_daught_p=0;
    particle_info.flag_has_daught_el=0;
    particle_info.flag_has_daught_pi=0;
    for(int reco_daught_part=0; reco_daught_part<pfeval.reco_Ntrack; reco_daught_part++){
      if(pfeval.reco_mother[reco_daught_part]!=pfeval.reco_id[reco_part]) continue;
      particle_info.flag_has_daught+=1;
      if(pfeval.reco_pdg[reco_daught_part]==2212) particle_info.flag_has_daught_p+=1;
      if(pfeval.reco_pdg[reco_daught_part]==11) particle_info.flag_has_daught_el+=1;
      if(pfeval.reco_pdg[reco_daught_part]==211) particle_info.flag_has_daught_pi+=1;
    }
	    
    particle_info.reco_is_n_induced = 0;
    particle_info.reco_is_g_induced = 0;
    for(int reco_mother_part=0; reco_mother_part<pfeval.reco_Ntrack; reco_mother_part++){
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
    particle_info.reco_momentum_0 = pfeval.reco_startMomentum[reco_part][0];
    particle_info.reco_momentum_1 = pfeval.reco_startMomentum[reco_part][1];
    particle_info.reco_momentum_2 = pfeval.reco_startMomentum[reco_part][2];
    particle_info.reco_momentum_3 = pfeval.reco_startMomentum[reco_part][3]-mass;
    particle_info.reco_pdg = pfeval.reco_pdg[reco_part];

    // Check if its the leading muon
    if(pfeval.reco_pdg[reco_part]==13 && pfeval.reco_startMomentum[reco_part][3]>pfeval.reco_muonMomentum[3]-tolerance && pfeval.reco_startMomentum[reco_part][3]<pfeval.reco_muonMomentum[3]+tolerance) particle_info.flag_prim_mu = 1;
    else if(pfeval.reco_pdg[reco_part]==13) particle_info.flag_prim_mu = 0;
    else particle_info.flag_prim_mu = -1;

    // Scattering angle
    particle_info.cos_theta = pfeval.reco_startMomentum[reco_part][2] / sqrt( pow(pfeval.reco_startMomentum[reco_part][0],2) + pow(pfeval.reco_startMomentum[reco_part][1],2) + pow(pfeval.reco_startMomentum[reco_part][2],2) );

    // The oriximity to all other particles based off endpoints
    particle_info.proximity=99999999;
    for(int other_reco_part=0; other_reco_part<pfeval.reco_Ntrack; other_reco_part++){
      if (pfeval.reco_id[other_reco_part] == pfeval.reco_id[reco_part]) continue;
      if (pfeval.reco_pdg[other_reco_part] == 2112 || pfeval.reco_pdg[other_reco_part] == 22) continue;
      double temp_proximity_1 =  sqrt(pow(part_x-pfeval.reco_startXYZT[other_reco_part][0],2)+pow(part_y-pfeval.reco_startXYZT[other_reco_part][1],2)+pow(part_z-pfeval.reco_startXYZT[other_reco_part][2],2));
      double temp_proximity_2 = sqrt(pow(part_x-pfeval.reco_endXYZT[other_reco_part][0],2)+pow(part_y-pfeval.reco_endXYZT[other_reco_part][1],2)+pow(part_z-pfeval.reco_endXYZT[other_reco_part][2],2));
      double temp_proximity_3 = sqrt(pow(part_end_x-pfeval.reco_startXYZT[other_reco_part][0],2)+pow(part_end_y-pfeval.reco_startXYZT[other_reco_part][1],2)+pow(part_end_z-pfeval.reco_startXYZT[other_reco_part][2],2));
      double temp_proximity_4 = sqrt(pow(part_end_x-pfeval.reco_endXYZT[other_reco_part][0],2)+pow(part_end_y-pfeval.reco_endXYZT[other_reco_part][1],2)+pow(part_end_z-pfeval.reco_endXYZT[other_reco_part][2],2));
      double temp_proximity = std::min({temp_proximity_1,temp_proximity_2,temp_proximity_3,temp_proximity_4});
      if(temp_proximity<particle_info.proximity && temp_proximity<99999) particle_info.proximity=temp_proximity;
    }

    // Track lengths and ratio
    particle_info.direct_track_len = sqrt(pow(part_x-part_end_x,2)+pow(part_y-part_end_y,2)+pow(part_z-part_end_z,2));
    particle_info.track_len_ratio = particle_info.direct_track_len/particle_info.track_len;

    // Check the containment
    particle_info.flag_is_contained = 1;
    if(part_end_x < 3 || part_end_x > 250) particle_info.flag_is_contained = 0;
    else if(part_end_y < -113 || part_end_y > 113) particle_info.flag_is_contained = 0;
    else if(part_end_z < 3 || part_end_z > 1035) particle_info.flag_is_contained = 0;

    // Save the larpid vars
    particle_info.reco_larpid_pdg = pfeval.reco_larpid_pdg[reco_part];
    particle_info.reco_larpid_pidScore_el = pfeval.reco_larpid_pidScore_el[reco_part];
    particle_info.reco_larpid_pidScore_ph = pfeval.reco_larpid_pidScore_ph[reco_part];
    particle_info.reco_larpid_pidScore_mu = pfeval.reco_larpid_pidScore_mu[reco_part];
    particle_info.reco_larpid_pidScore_pr = pfeval.reco_larpid_pidScore_pr[reco_part];
    particle_info.reco_larpid_pidScore_pi = pfeval.reco_larpid_pidScore_pi[reco_part];
    particle_info.reco_larpid_proccess = pfeval.reco_larpid_proccess[reco_part];


    // Save the individual spacepoints for easier use in the BDT (from the ordered main sequence)
    int n_main = particle_info.spacepoints_q->size();
    particle_info.spacepoints_q_0 = (particle_info.spacepoints_q->at(0)+10000)*10;
    particle_info.spacepoints_q_bck_0 = (particle_info.spacepoints_q->at(n_main-1-0)+10000)*10;
    if(n_main>1){
      particle_info.spacepoints_q_1 = (particle_info.spacepoints_q->at(1)+10000)*10;
      particle_info.spacepoints_q_bck_1 = (particle_info.spacepoints_q->at(n_main-1-1)+10000)*10;
    }else{
      particle_info.spacepoints_q_1 = -999;
      particle_info.spacepoints_q_bck_1 = -999;
    }
    if(n_main>2){
      particle_info.spacepoints_q_2 = (particle_info.spacepoints_q->at(2)+10000)*10;
      particle_info.spacepoints_q_bck_2 = (particle_info.spacepoints_q->at(n_main-1-2)+10000)*10;
    }else{
      particle_info.spacepoints_q_2 = -999;
      particle_info.spacepoints_q_bck_2 = -999;
    }
    if(n_main>3){
      particle_info.spacepoints_q_3 = (particle_info.spacepoints_q->at(3)+10000)*10;
      particle_info.spacepoints_q_bck_3 = (particle_info.spacepoints_q->at(n_main-1-3)+10000)*10;
    }else{
      particle_info.spacepoints_q_3 = -999;
      particle_info.spacepoints_q_bck_3 = -999;
    }
    if(n_main>4){
      particle_info.spacepoints_q_4 = (particle_info.spacepoints_q->at(4)+10000)*10;
      particle_info.spacepoints_q_bck_4 = (particle_info.spacepoints_q->at(n_main-1-4)+10000)*10;
    }else{
      particle_info.spacepoints_q_4 = -999;
      particle_info.spacepoints_q_bck_4 = -999;
    }
    if(n_main>5){
      particle_info.spacepoints_q_5 = (particle_info.spacepoints_q->at(5)+10000)*10;
      particle_info.spacepoints_q_bck_5 = (particle_info.spacepoints_q->at(n_main-1-5)+10000)*10;
    }else{
      particle_info.spacepoints_q_5 = -999;
      particle_info.spacepoints_q_bck_5 = -999;
    }
    if(n_main>6){
      particle_info.spacepoints_q_6 = (particle_info.spacepoints_q->at(6)+10000)*10;
      particle_info.spacepoints_q_bck_6 = (particle_info.spacepoints_q->at(n_main-1-6)+10000)*10;
    }else{
      particle_info.spacepoints_q_6 = -999;
      particle_info.spacepoints_q_bck_6 = -999;
    }
    if(n_main>7){
      particle_info.spacepoints_q_7 = (particle_info.spacepoints_q->at(7)+10000)*10;
      particle_info.spacepoints_q_bck_7 = (particle_info.spacepoints_q->at(n_main-1-7)+10000)*10;
    }else{
      particle_info.spacepoints_q_7 = -999;
      particle_info.spacepoints_q_bck_7 = -999;
    }
    if(n_main>8){
      particle_info.spacepoints_q_8 = (particle_info.spacepoints_q->at(8)+10000)*10;
      particle_info.spacepoints_q_bck_8 = (particle_info.spacepoints_q->at(n_main-1-8)+10000)*10;
    }else{
      particle_info.spacepoints_q_8 = -999;
      particle_info.spacepoints_q_bck_8 = -999;
    }
    if(n_main>9){
      particle_info.spacepoints_q_9 = (particle_info.spacepoints_q->at(9)+10000)*10;
      particle_info.spacepoints_q_bck_9 = (particle_info.spacepoints_q->at(n_main-1-9)+10000)*10;
    }else{
      particle_info.spacepoints_q_9 = -999;
      particle_info.spacepoints_q_bck_9 = -999;
    }
    if(n_main>10){
      particle_info.spacepoints_q_10 = (particle_info.spacepoints_q->at(10)+10000)*10;
      particle_info.spacepoints_q_bck_10 = (particle_info.spacepoints_q->at(n_main-1-10)+10000)*10;
    }else{
      particle_info.spacepoints_q_10 = -999;
      particle_info.spacepoints_q_bck_10 = -999;
    }
    if(n_main>11){
      particle_info.spacepoints_q_11 = (particle_info.spacepoints_q->at(11)+10000)*10;
      particle_info.spacepoints_q_bck_11 = (particle_info.spacepoints_q->at(n_main-1-11)+10000)*10;
    }else{
      particle_info.spacepoints_q_11 = -999;
      particle_info.spacepoints_q_bck_11 = -999;
    }
    if(n_main>12){
      particle_info.spacepoints_q_12 = (particle_info.spacepoints_q->at(12)+10000)*10;
      particle_info.spacepoints_q_bck_12 = (particle_info.spacepoints_q->at(n_main-1-12)+10000)*10;
    }else{
      particle_info.spacepoints_q_12 = -999;
      particle_info.spacepoints_q_bck_12 = -999;
    }
    if(n_main>13){
      particle_info.spacepoints_q_13 = (particle_info.spacepoints_q->at(13)+10000)*10;
      particle_info.spacepoints_q_bck_13 = (particle_info.spacepoints_q->at(n_main-1-13)+10000)*10;
    }else{
      particle_info.spacepoints_q_13 = -999;
      particle_info.spacepoints_q_bck_13 = -999;
    }
    if(n_main>14){
      particle_info.spacepoints_q_14 = (particle_info.spacepoints_q->at(14)+10000)*10;
      particle_info.spacepoints_q_bck_14 = (particle_info.spacepoints_q->at(n_main-1-14)+10000)*10;
    }else{
      particle_info.spacepoints_q_14 = -999;
      particle_info.spacepoints_q_bck_14 = -999;
    }
    if(n_main>15){
      particle_info.spacepoints_q_15 = (particle_info.spacepoints_q->at(15)+10000)*10;
      particle_info.spacepoints_q_bck_15 = (particle_info.spacepoints_q->at(n_main-1-15)+10000)*10;
    }else{
      particle_info.spacepoints_q_15 = -999;
      particle_info.spacepoints_q_bck_15 = -999;
    }
    if(n_main>16){
      particle_info.spacepoints_q_16 = (particle_info.spacepoints_q->at(16)+10000)*10;
      particle_info.spacepoints_q_bck_16 = (particle_info.spacepoints_q->at(n_main-1-16)+10000)*10;
    }else{
      particle_info.spacepoints_q_16 = -999;
      particle_info.spacepoints_q_bck_16 = -999;
    }
    if(n_main>17){
      particle_info.spacepoints_q_17 = (particle_info.spacepoints_q->at(17)+10000)*10;
      particle_info.spacepoints_q_bck_17 = (particle_info.spacepoints_q->at(n_main-1-17)+10000)*10;
    }else{
      particle_info.spacepoints_q_17 = -999;
      particle_info.spacepoints_q_bck_17 = -999;
    }
    if(n_main>18){
      particle_info.spacepoints_q_18 = (particle_info.spacepoints_q->at(18)+10000)*10;
      particle_info.spacepoints_q_bck_18 = (particle_info.spacepoints_q->at(n_main-1-18)+10000)*10;
    }else{
      particle_info.spacepoints_q_18 = -999;
      particle_info.spacepoints_q_bck_18 = -999;
    }
    if(n_main>19){
      particle_info.spacepoints_q_19 = (particle_info.spacepoints_q->at(19)+10000)*10;
      particle_info.spacepoints_q_bck_19 = (particle_info.spacepoints_q->at(n_main-1-19)+10000)*10;
    }else{
      particle_info.spacepoints_q_19 = -999;
      particle_info.spacepoints_q_bck_19 = -999;
    }
    if(n_main>20){
      particle_info.spacepoints_q_20 = (particle_info.spacepoints_q->at(20)+10000)*10;
      particle_info.spacepoints_q_bck_20 = (particle_info.spacepoints_q->at(n_main-1-20)+10000)*10;
    }else{
      particle_info.spacepoints_q_20 = -999;
      particle_info.spacepoints_q_bck_20 = -999;
    }
    if(n_main>21){
      particle_info.spacepoints_q_21 = (particle_info.spacepoints_q->at(21)+10000)*10;
      particle_info.spacepoints_q_bck_21 = (particle_info.spacepoints_q->at(n_main-1-21)+10000)*10;
    }else{
      particle_info.spacepoints_q_21 = -999;
      particle_info.spacepoints_q_bck_21 = -999;
    }
    if(n_main>22){
      particle_info.spacepoints_q_22 = (particle_info.spacepoints_q->at(22)+10000)*10;
      particle_info.spacepoints_q_bck_22 = (particle_info.spacepoints_q->at(n_main-1-22)+10000)*10;
    }else{
      particle_info.spacepoints_q_22 = -999;
      particle_info.spacepoints_q_bck_22 = -999;
    }
    if(n_main>23){
      particle_info.spacepoints_q_23 = (particle_info.spacepoints_q->at(23)+10000)*10;
      particle_info.spacepoints_q_bck_23 = (particle_info.spacepoints_q->at(n_main-1-23)+10000)*10;
    }else{
      particle_info.spacepoints_q_23 = -999;
      particle_info.spacepoints_q_bck_23 = -999;
    }
    if(n_main>24){
      particle_info.spacepoints_q_24 = (particle_info.spacepoints_q->at(24)+10000)*10;
      particle_info.spacepoints_q_bck_24 = (particle_info.spacepoints_q->at(n_main-1-24)+10000)*10;
    }else{
      particle_info.spacepoints_q_24 = -999;
      particle_info.spacepoints_q_bck_24 = -999;
    }


  }
  //if(particle_info.reco_pdg>0) print_particle(particle_info);
}

// Order a particle's spacepoints along its trajectory and describe any extra pieces (same logic as the python
// build_particle_segments used for the BDT training).
//
// WireCell stores a particle's spacepoints as separate pieces: the reco start/end (and other node) points first, then one
// or more blocks of trajectory points (0.6 cm spacing), each block stored in its own direction.
//
// 0. ref (optional, nullptr if none): the point the particle should start from (the neutrino vertex for primaries, the
//    nearer of the mother's reco start/end for secondaries). If the reco end is closer to ref than the reco start by more
//    than min(flip_tol, 0.5 x start-end length), the reco start and end are swapped before ordering.
// 1. Anchors: the stored points at the reco start and reco end (within anchor_tol) are taken out first.
// 2. The remaining points are split, in stored order, into blocks (runs of steps < jump, >= 2 points) and isolated points.
// 3. <= 1 block (almost all protons/pions, and all short particles): main = [start anchor] + block + [end anchor], where
//    isolated points within attach of a block end are attached to that end, the block is oriented from the reco start,
//    and each anchor is only added if it is within attach of that end of the sequence. With no block, the isolated
//    points are ordered by distance from the reco start.
// 4. >= 2 blocks: endpoints (anchors, isolated points, block ends) closer than junction are merged into junctions (never
//    merging the two ends of the same block, or the start and end anchors). The blocks are edges between junctions.
//    main = path from the start junction to the end junction (isolated points at junctions on the path are kept in
//    order between the blocks). If there is no such path, the longest path reachable from the start junction is used;
//    if the start junction is not found, the longest block.
// 5. Secondary = the longest block not in main, oriented so the end closest to main comes first.
LEEana::ParticleSegments LEEana::build_particle_segments(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& q,
                                                 const double start[3], const double end[3], const double* ref,
                                                 double jump, double junction, double attach, double anchor_tol, double flip_tol){

  ParticleSegments res;
  res.n_segments = 0;
  res.seg_total_len = 0;
  res.seg_total_q = 0;
  res.sec_seg_len = -999;
  res.sec_seg_direct_len = -999;
  res.sec_seg_dist_to_main = -999;
  res.flipped = false;

  auto dist3 = [](const double* a, const double* b){ return sqrt(pow(a[0]-b[0],2)+pow(a[1]-b[1],2)+pow(a[2]-b[2],2)); };

  double S[3] = {start[0], start[1], start[2]};
  double E[3] = {end[0], end[1], end[2]};
  if(ref!=nullptr){
    if(dist3(E,ref) < dist3(S,ref) - std::min(flip_tol, 0.5*dist3(E,S))){
      for(int c=0; c<3; c++) std::swap(S[c], E[c]);
      res.flipped = true;
    }
  }

  int n = x.size();
  if(n==0) return res;
  std::vector<std::array<double,3>> P(n);
  for(int k=0; k<n; k++) P[k] = {x.at(k), y.at(k), z.at(k)};
  auto dpp = [&](int a, int b){ return dist3(P[a].data(), P[b].data()); };
  auto plen = [&](const std::vector<int>& idx){
    double len = 0;
    for(size_t k=1; k<idx.size(); k++) len += dpp(idx[k-1], idx[k]);
    return len;
  };

  // 1. anchors
  int s_anchor = -1, e_anchor = -1;
  double ds_min = 1e9, de_min = 1e9;
  for(int k=0; k<n; k++){
    double ds = dist3(P[k].data(), S), de = dist3(P[k].data(), E);
    if(ds<ds_min){ ds_min = ds; s_anchor = k; }
    if(de<de_min){ de_min = de; e_anchor = k; }
  }
  if(ds_min>=anchor_tol) s_anchor = -1;
  if(de_min>=anchor_tol) e_anchor = -1;
  if(e_anchor==s_anchor) e_anchor = -1;
  std::vector<int> rest;
  for(int k=0; k<n; k++) if(k!=s_anchor && k!=e_anchor) rest.push_back(k);

  // 2. blocks and isolated points among the remaining points (stored order)
  std::vector<std::vector<int>> blocks;
  std::vector<int> isolated;
  std::vector<int> run;
  for(size_t k=0; k<rest.size(); k++){
    if(k>0 && dpp(rest[k-1], rest[k])>=jump){
      if(run.size()>=2) blocks.push_back(run); else isolated.push_back(run[0]);
      run.clear();
    }
    run.push_back(rest[k]);
  }
  if(run.size()>=2) blocks.push_back(run); else if(run.size()==1) isolated.push_back(run[0]);

  // seq is ordered from the start side; with no other points, keep whichever anchors exist
  auto add_anchors = [&](const std::vector<int>& seq){
    std::vector<int> out;
    if(seq.empty()){
      if(s_anchor>=0) out.push_back(s_anchor);
      if(e_anchor>=0) out.push_back(e_anchor);
      return out;
    }
    if(s_anchor>=0 && dpp(seq.front(), s_anchor)<attach) out.push_back(s_anchor);
    out.insert(out.end(), seq.begin(), seq.end());
    if(e_anchor>=0 && dpp(out.back(), e_anchor)<attach) out.push_back(e_anchor);
    return out;
  };
  auto orient_from_start = [&](std::vector<int>& b){
    if(dist3(P[b.back()].data(), S) < dist3(P[b.front()].data(), S)) std::reverse(b.begin(), b.end());
  };

  std::vector<int> main;
  std::vector<int> main_blocks; // indices into blocks that are part of main
  if(blocks.size()<=1){
    // 3. simple path
    if(!blocks.empty()){
      std::vector<int> b = blocks[0];
      // attach isolated points to the nearest block end
      std::vector<std::pair<double,int>> front, back;
      for(int k : isolated){
        double d0 = dpp(k, b.front()), d1 = dpp(k, b.back());
        if(std::min(d0,d1)<attach){
          if(d0<=d1) front.push_back({d0,k}); else back.push_back({d1,k});
        }
      }
      std::sort(front.begin(), front.end(), std::greater<std::pair<double,int>>());
      std::sort(back.begin(), back.end());
      std::vector<int> bb;
      for(auto& f : front) bb.push_back(f.second);
      bb.insert(bb.end(), b.begin(), b.end());
      for(auto& f : back) bb.push_back(f.second);
      orient_from_start(bb);
      main = add_anchors(bb);
      main_blocks.push_back(0);
    }else{
      std::vector<int> iso = isolated;
      std::stable_sort(iso.begin(), iso.end(), [&](int a, int b){ return dist3(P[a].data(), S) < dist3(P[b].data(), S); });
      main = add_anchors(iso);
    }
  }else{
    // 4. junction path
    struct EndPoint{ int pt; char kind; int block; int side; };
    std::vector<EndPoint> ends;
    if(s_anchor>=0) ends.push_back({s_anchor, 'S', -1, -1});
    if(e_anchor>=0) ends.push_back({e_anchor, 'E', -1, -1});
    for(int k : isolated) ends.push_back({k, 'I', -1, -1});
    for(size_t bi=0; bi<blocks.size(); bi++){
      ends.push_back({blocks[bi].front(), 'B', (int)bi, 0});
      ends.push_back({blocks[bi].back(), 'B', (int)bi, 1});
    }
    int m = ends.size();
    std::vector<int> parent(m);
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<std::set<int>> mem_blocks(m);
    std::vector<bool> mem_S(m), mem_E(m);
    for(int a=0; a<m; a++){
      if(ends[a].kind=='B') mem_blocks[a].insert(ends[a].block);
      mem_S[a] = ends[a].kind=='S';
      mem_E[a] = ends[a].kind=='E';
    }
    auto find = [&](int a){
      while(parent[a]!=a){ parent[a] = parent[parent[a]]; a = parent[a]; }
      return a;
    };
    std::vector<std::tuple<double,int,int>> pairs;
    for(int a=0; a<m; a++) for(int b=a+1; b<m; b++){
      double d = dpp(ends[a].pt, ends[b].pt);
      if(d<junction) pairs.push_back(std::make_tuple(d,a,b));
    }
    std::sort(pairs.begin(), pairs.end());
    for(auto& pr : pairs){
      int ra = find(std::get<1>(pr)), rb = find(std::get<2>(pr));
      if(ra==rb) continue;
      bool shared_block = false;
      for(int bi : mem_blocks[rb]) if(mem_blocks[ra].count(bi)) shared_block = true;
      if(shared_block) continue;                                                 // would join the two ends of one block
      if((mem_S[ra] && mem_E[rb]) || (mem_E[ra] && mem_S[rb])) continue;
      parent[rb] = ra;
      mem_blocks[ra].insert(mem_blocks[rb].begin(), mem_blocks[rb].end());
      mem_S[ra] = mem_S[ra] || mem_S[rb];
      mem_E[ra] = mem_E[ra] || mem_E[rb];
    }
    std::vector<int> junc_of(m);
    for(int a=0; a<m; a++) junc_of[a] = find(a);
    std::map<std::pair<int,int>,int> bend_j;
    std::map<int,std::vector<int>> iso_at;                                       // junction -> isolated points there
    for(int a=0; a<m; a++){
      if(ends[a].kind=='B') bend_j[std::make_pair(ends[a].block, ends[a].side)] = junc_of[a];
      if(ends[a].kind=='I') iso_at[junc_of[a]].push_back(ends[a].pt);
    }
    std::map<int,std::vector<std::array<int,3>>> adj;                            // junction -> (other junction, block, side)
    for(size_t bi=0; bi<blocks.size(); bi++){
      int a = bend_j[std::make_pair((int)bi,0)], b = bend_j[std::make_pair((int)bi,1)];
      adj[a].push_back({b, (int)bi, 0});
      adj[b].push_back({a, (int)bi, 1});
    }

    auto junction_near = [&](const double* pt){
      int best = -1; double dbest = 1e9;
      for(int a=0; a<m; a++){
        double d = dist3(P[ends[a].pt].data(), pt);
        if(d<dbest){ dbest = d; best = a; }
      }
      return (best>=0 && dbest<attach) ? junc_of[best] : -1;
    };
    int js = -1, je = -1;
    for(int a=0; a<m; a++){
      if(ends[a].kind=='S' && js<0) js = junc_of[a];
      if(ends[a].kind=='E' && je<0) je = junc_of[a];
    }
    if(s_anchor<0) js = junction_near(S);
    if(e_anchor<0) je = junction_near(E);

    // BFS over junctions; prev[j] = (previous junction, block, side), with -1 for the starting junction
    std::map<int,std::array<int,3>> prev;
    std::vector<int> visit_order;
    auto walk = [&](int jt){
      std::vector<std::array<int,4>> steps;                                      // (from junction, block, side, to junction)
      int u = jt;
      while(prev[u][0]>=0){
        std::array<int,3> p = prev[u];
        steps.push_back({p[0], p[1], p[2], u});
        u = p[0];
      }
      std::reverse(steps.begin(), steps.end());                                  // from js to jt
      return steps;
    };
    std::vector<std::array<int,4>> steps;
    if(js>=0){
      prev[js] = {-1,-1,-1};
      visit_order.push_back(js);
      std::deque<int> dq;
      dq.push_back(js);
      while(!dq.empty()){
        int u = dq.front(); dq.pop_front();
        for(auto& w : adj[u]){
          if(prev.count(w[0])) continue;
          prev[w[0]] = {u, w[1], w[2]};
          visit_order.push_back(w[0]);
          dq.push_back(w[0]);
        }
      }
      if(je>=0 && prev.count(je) && je!=js){
        steps = walk(je);
      }else{
        // longest reachable path (by number of points)
        int best_npts = -1;
        for(int jt : visit_order){
          if(jt==js) continue;
          std::vector<std::array<int,4>> w = walk(jt);
          int npts = 0;
          for(auto& st : w) npts += blocks[st[1]].size();
          if(npts>best_npts){ best_npts = npts; steps = w; }
        }
      }
    }
    if(!steps.empty()){
      std::vector<int> seq;
      auto add_iso = [&](int j, const double* r){
        if(!iso_at.count(j)) return;
        std::vector<int> iso = iso_at[j];
        std::stable_sort(iso.begin(), iso.end(), [&](int a, int b){ return dist3(P[a].data(), r) < dist3(P[b].data(), r); });
        seq.insert(seq.end(), iso.begin(), iso.end());
      };
      add_iso(steps.front()[0], S);
      int j_last = steps.back()[3];
      for(auto& st : steps){
        std::vector<int> b = blocks[st[1]];
        if(st[2]==1) std::reverse(b.begin(), b.end());                           // side 0: traversed from its end 0
        seq.insert(seq.end(), b.begin(), b.end());
        main_blocks.push_back(st[1]);
        if(st[3]!=j_last){ double r[3] = {P[seq.back()][0], P[seq.back()][1], P[seq.back()][2]}; add_iso(st[3], r); }
      }
      double r[3] = {P[seq.back()][0], P[seq.back()][1], P[seq.back()][2]};
      add_iso(j_last, r);
      main = add_anchors(seq);
    }else{
      int bi_best = 0;
      for(size_t bi=1; bi<blocks.size(); bi++) if(plen(blocks[bi])>plen(blocks[bi_best])) bi_best = bi;
      std::vector<int> b = blocks[bi_best];
      orient_from_start(b);
      main = add_anchors(b);
      main_blocks.push_back(bi_best);
    }
  }
  if(main.empty()){ main.resize(n); std::iota(main.begin(), main.end(), 0); } // safety, should not happen
  res.main = main;

  std::vector<int> others;
  for(size_t bi=0; bi<blocks.size(); bi++) if(std::find(main_blocks.begin(), main_blocks.end(), (int)bi)==main_blocks.end()) others.push_back(bi);
  res.n_segments = 1 + others.size();
  res.seg_total_len = plen(main);
  for(int k : main) res.seg_total_q += q.at(k);
  for(int bi : others){
    res.seg_total_len += plen(blocks[bi]);
    for(int k : blocks[bi]) res.seg_total_q += q.at(k);
  }
  if(!others.empty()){
    int bi_sec = others[0];
    for(int bi : others) if(plen(blocks[bi])>plen(blocks[bi_sec])) bi_sec = bi;
    std::vector<int> sec = blocks[bi_sec];
    auto dist_to_main = [&](int k){
      double d = 1e9;
      for(int k2 : main) d = std::min(d, dpp(k, k2));
      return d;
    };
    if(dist_to_main(sec.back())<dist_to_main(sec.front())) std::reverse(sec.begin(), sec.end());
    res.sec = sec;
    res.sec_seg_len = plen(sec);
    res.sec_seg_direct_len = dpp(sec.front(), sec.back());
    res.sec_seg_dist_to_main = 1e9;
    for(int k : sec) res.sec_seg_dist_to_main = std::min(res.sec_seg_dist_to_main, dist_to_main(k));
  }
  return res;
}

void LEEana::reset_particle(ParticleInfo& particle_info){
    particle_info.spacepoints_x->clear();
    particle_info.spacepoints_y->clear();
    particle_info.spacepoints_z->clear();
    particle_info.spacepoints_q->clear();

    particle_info.spacepoints_q_0=-999;
    particle_info.spacepoints_q_1=-999;
    particle_info.spacepoints_q_2=-999;
    particle_info.spacepoints_q_3=-999;
    particle_info.spacepoints_q_4=-999;
    particle_info.spacepoints_q_5=-999;
    particle_info.spacepoints_q_6=-999;
    particle_info.spacepoints_q_7=-999;
    particle_info.spacepoints_q_8=-999;
    particle_info.spacepoints_q_9=-999;
    particle_info.spacepoints_q_10=-999;
    particle_info.spacepoints_q_11=-999;
    particle_info.spacepoints_q_12=-999;
    particle_info.spacepoints_q_13=-999;
    particle_info.spacepoints_q_14=-999;
    particle_info.spacepoints_q_15=-999;
    particle_info.spacepoints_q_16=-999;
    particle_info.spacepoints_q_17=-999;
    particle_info.spacepoints_q_18=-999;
    particle_info.spacepoints_q_19=-999;
    particle_info.spacepoints_q_20=-999;
    particle_info.spacepoints_q_21=-999;
    particle_info.spacepoints_q_22=-999;
    particle_info.spacepoints_q_23=-999;
    particle_info.spacepoints_q_24=-999;

    particle_info.spacepoints_q_bck_0=-999;
    particle_info.spacepoints_q_bck_1=-999;
    particle_info.spacepoints_q_bck_2=-999;
    particle_info.spacepoints_q_bck_3=-999;
    particle_info.spacepoints_q_bck_4=-999;
    particle_info.spacepoints_q_bck_5=-999;
    particle_info.spacepoints_q_bck_6=-999;
    particle_info.spacepoints_q_bck_7=-999;
    particle_info.spacepoints_q_bck_8=-999;
    particle_info.spacepoints_q_bck_9=-999;
    particle_info.spacepoints_q_bck_10=-999;
    particle_info.spacepoints_q_bck_11=-999;
    particle_info.spacepoints_q_bck_12=-999;
    particle_info.spacepoints_q_bck_13=-999;
    particle_info.spacepoints_q_bck_14=-999;
    particle_info.spacepoints_q_bck_15=-999;
    particle_info.spacepoints_q_bck_16=-999;
    particle_info.spacepoints_q_bck_17=-999;
    particle_info.spacepoints_q_bck_18=-999;
    particle_info.spacepoints_q_bck_19=-999;
    particle_info.spacepoints_q_bck_20=-999;
    particle_info.spacepoints_q_bck_21=-999;
    particle_info.spacepoints_q_bck_22=-999;
    particle_info.spacepoints_q_bck_23=-999;
    particle_info.spacepoints_q_bck_24=-999;

    particle_info.spacepoints_q_med=-999;

    particle_info.flag_prim_mu=-999;

    particle_info.flag_is_contained=-999;

    particle_info.flag_has_daught=-999;
    particle_info.flag_has_daught_p=-999;
    particle_info.flag_has_daught_el=-999;
    particle_info.flag_has_daught_pi=-999;

    particle_info.reco_truthMatch_pdg=-999;
    particle_info.reco_truthMatch_id=-999;
    particle_info.reco_truthMatch_mother=-999;
    particle_info.reco_truthMatch_energy=-999;

    particle_info.reco_momentum_0=-999;
    particle_info.reco_momentum_1=-999;
    particle_info.reco_momentum_2=-999;
    particle_info.reco_momentum_3=-999;
    particle_info.reco_pdg=-999;

    particle_info.reco_larpid_pdg=-999;
    particle_info.reco_larpid_pidScore_el=-999;
    particle_info.reco_larpid_pidScore_ph=-999;
    particle_info.reco_larpid_pidScore_mu=-999;
    particle_info.reco_larpid_pidScore_pr=-999;
    particle_info.reco_larpid_pidScore_pi=-999;
    particle_info.reco_larpid_proccess=-999;

    particle_info.true_is_n_induced=-999;
    particle_info.reco_is_n_induced=-999;
    particle_info.reco_is_g_induced=-999;

    particle_info.dist_to_vtx=-999;
    particle_info.cos_theta=-999;
    particle_info.proximity=-999;

    particle_info.track_len=-999;
    particle_info.direct_track_len=-999;
    particle_info.track_len_ratio=-999;

    particle_info.spacepoints_q_sec_0=-999;
    particle_info.spacepoints_q_sec_1=-999;
    particle_info.spacepoints_q_sec_2=-999;
    particle_info.spacepoints_q_sec_3=-999;
    particle_info.spacepoints_q_sec_4=-999;
    particle_info.spacepoints_q_sec_bck_0=-999;
    particle_info.spacepoints_q_sec_bck_1=-999;
    particle_info.spacepoints_q_sec_bck_2=-999;
    particle_info.spacepoints_q_sec_bck_3=-999;
    particle_info.spacepoints_q_sec_bck_4=-999;

    particle_info.n_segments=0;
    particle_info.seg_total_len=0;
    particle_info.seg_total_q=0;
    particle_info.sec_seg_len=-999;
    particle_info.sec_seg_direct_len=-999;
    particle_info.sec_seg_dist_to_main=-999;
    particle_info.flag_start_end_swapped=0;
}

void LEEana::print_particle(ParticleInfo& particle_info){


    std::cout<<"particle_info.spacepoints_q_0 "<<particle_info.spacepoints_q_0<<std::endl;
    std::cout<<"particle_info.spacepoints_q_1 "<<particle_info.spacepoints_q_1<<std::endl;
    std::cout<<"particle_info.spacepoints_q_2 "<<particle_info.spacepoints_q_2<<std::endl;
    std::cout<<"particle_info.spacepoints_q_3 "<<particle_info.spacepoints_q_3<<std::endl;
    std::cout<<"particle_info.spacepoints_q_4 "<<particle_info.spacepoints_q_4<<std::endl;
    std::cout<<"particle_info.spacepoints_q_5 "<<particle_info.spacepoints_q_5<<std::endl;
    std::cout<<"particle_info.spacepoints_q_6 "<<particle_info.spacepoints_q_6<<std::endl;
    std::cout<<"particle_info.spacepoints_q_7 "<<particle_info.spacepoints_q_7<<std::endl;
    std::cout<<"particle_info.spacepoints_q_8 "<<particle_info.spacepoints_q_8<<std::endl;
    std::cout<<"particle_info.spacepoints_q_9 "<<particle_info.spacepoints_q_9<<std::endl;
    std::cout<<"particle_info.spacepoints_q_10 "<<particle_info.spacepoints_q_10<<std::endl;
    std::cout<<"particle_info.spacepoints_q_11 "<<particle_info.spacepoints_q_11<<std::endl;
    std::cout<<"particle_info.spacepoints_q_12 "<<particle_info.spacepoints_q_12<<std::endl;
    std::cout<<"particle_info.spacepoints_q_13 "<<particle_info.spacepoints_q_13<<std::endl;
    std::cout<<"particle_info.spacepoints_q_14 "<<particle_info.spacepoints_q_14<<std::endl;
    std::cout<<"particle_info.spacepoints_q_15 "<<particle_info.spacepoints_q_15<<std::endl;
    std::cout<<"particle_info.spacepoints_q_16 "<<particle_info.spacepoints_q_16<<std::endl;
    std::cout<<"particle_info.spacepoints_q_17 "<<particle_info.spacepoints_q_17<<std::endl;
    std::cout<<"particle_info.spacepoints_q_18 "<<particle_info.spacepoints_q_18<<std::endl;
    std::cout<<"particle_info.spacepoints_q_19 "<<particle_info.spacepoints_q_19<<std::endl;
    std::cout<<"particle_info.spacepoints_q_20 "<<particle_info.spacepoints_q_20<<std::endl;
    std::cout<<"particle_info.spacepoints_q_21 "<<particle_info.spacepoints_q_21<<std::endl;
    std::cout<<"particle_info.spacepoints_q_22 "<<particle_info.spacepoints_q_22<<std::endl;
    std::cout<<"particle_info.spacepoints_q_23 "<<particle_info.spacepoints_q_23<<std::endl;
    std::cout<<"particle_info.spacepoints_q_24 "<<particle_info.spacepoints_q_24<<std::endl;

    std::cout<<"particle_info.spacepoints_q_bck_0 "<<particle_info.spacepoints_q_bck_0<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_1 "<<particle_info.spacepoints_q_bck_1<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_2 "<<particle_info.spacepoints_q_bck_2<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_3 "<<particle_info.spacepoints_q_bck_3<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_4 "<<particle_info.spacepoints_q_bck_4<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_5 "<<particle_info.spacepoints_q_bck_5<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_6 "<<particle_info.spacepoints_q_bck_6<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_7 "<<particle_info.spacepoints_q_bck_7<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_8 "<<particle_info.spacepoints_q_bck_8<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_9 "<<particle_info.spacepoints_q_bck_9<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_10 "<<particle_info.spacepoints_q_bck_10<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_11 "<<particle_info.spacepoints_q_bck_11<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_12 "<<particle_info.spacepoints_q_bck_12<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_13 "<<particle_info.spacepoints_q_bck_13<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_14 "<<particle_info.spacepoints_q_bck_14<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_15 "<<particle_info.spacepoints_q_bck_15<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_16 "<<particle_info.spacepoints_q_bck_16<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_17 "<<particle_info.spacepoints_q_bck_17<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_18 "<<particle_info.spacepoints_q_bck_18<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_19 "<<particle_info.spacepoints_q_bck_19<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_20 "<<particle_info.spacepoints_q_bck_20<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_21 "<<particle_info.spacepoints_q_bck_21<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_22 "<<particle_info.spacepoints_q_bck_22<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_23 "<<particle_info.spacepoints_q_bck_23<<std::endl;
    std::cout<<"particle_info.spacepoints_q_bck_24 "<<particle_info.spacepoints_q_bck_24<<std::endl;

    std::cout<<"particle_info.spacepoints_q_med "<<particle_info.spacepoints_q_med<<std::endl;

    std::cout<<"particle_info.flag_prim_mu "<<particle_info.flag_prim_mu<<std::endl;

    std::cout<<"particle_info.flag_is_contained "<<particle_info.flag_is_contained<<std::endl;

    std::cout<<"particle_info.flag_has_daught "<<particle_info.flag_has_daught<<std::endl;
    std::cout<<"particle_info.flag_has_daught_p "<<particle_info.flag_has_daught_p<<std::endl;
    std::cout<<"particle_info.flag_has_daught_el "<<particle_info.flag_has_daught_el<<std::endl;
    std::cout<<"particle_info.flag_has_daught_pi "<<particle_info.flag_has_daught_pi<<std::endl;

    std::cout<<"particle_info.reco_truthMatch_pdg "<<particle_info.reco_truthMatch_pdg<<std::endl;
    std::cout<<"particle_info.reco_truthMatch_id "<<particle_info.reco_truthMatch_id<<std::endl;
    std::cout<<"particle_info.reco_truthMatch_mother "<<particle_info.reco_truthMatch_mother<<std::endl;
    std::cout<<"particle_info.reco_truthMatch_energy "<<particle_info.reco_truthMatch_energy<<std::endl;

    std::cout<<"particle_info.reco_momentum_0 "<<particle_info.reco_momentum_0<<std::endl;
    std::cout<<"particle_info.reco_momentum_1 "<<particle_info.reco_momentum_1<<std::endl;
    std::cout<<"particle_info.reco_momentum_2 "<<particle_info.reco_momentum_2<<std::endl;
    std::cout<<"particle_info.reco_momentum_3 "<<particle_info.reco_momentum_3<<std::endl;
    std::cout<<"particle_info.reco_pdg "<<particle_info.reco_pdg<<std::endl;

    std::cout<<"particle_info.reco_larpid_pdg "<<particle_info.reco_larpid_pdg<<std::endl;
    std::cout<<"particle_info.reco_larpid_pidScore_el "<<particle_info.reco_larpid_pidScore_el<<std::endl;
    std::cout<<"particle_info.reco_larpid_pidScore_ph "<<particle_info.reco_larpid_pidScore_ph<<std::endl;
    std::cout<<"particle_info.reco_larpid_pidScore_mu "<<particle_info.reco_larpid_pidScore_mu<<std::endl;
    std::cout<<"particle_info.reco_larpid_pidScore_pr "<<particle_info.reco_larpid_pidScore_pr<<std::endl;
    std::cout<<"particle_info.reco_larpid_pidScore_pi "<<particle_info.reco_larpid_pidScore_pi<<std::endl;
    std::cout<<"particle_info.reco_larpid_proccess "<<particle_info.reco_larpid_proccess<<std::endl;

    std::cout<<"particle_info.true_is_n_induced "<<particle_info.true_is_n_induced<<std::endl;
    std::cout<<"particle_info.reco_is_n_induced "<<particle_info.reco_is_n_induced<<std::endl;
    std::cout<<"particle_info.reco_is_g_induced "<<particle_info.reco_is_g_induced<<std::endl;

    std::cout<<"particle_info.dist_to_vtx "<<particle_info.dist_to_vtx<<std::endl;
    std::cout<<"particle_info.cos_theta "<<particle_info.cos_theta<<std::endl;
    std::cout<<"particle_info.proximity "<<particle_info.proximity<<std::endl;

    std::cout<<"particle_info.track_len "<<particle_info.track_len<<std::endl;
    std::cout<<"particle_info.direct_track_len "<<particle_info.direct_track_len<<std::endl;
    std::cout<<"particle_info.track_len_ratio "<<particle_info.track_len_ratio<<std::endl;

    std::cout<<"particle_info.spacepoints_q_sec_0 "<<particle_info.spacepoints_q_sec_0<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_1 "<<particle_info.spacepoints_q_sec_1<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_2 "<<particle_info.spacepoints_q_sec_2<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_3 "<<particle_info.spacepoints_q_sec_3<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_4 "<<particle_info.spacepoints_q_sec_4<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_bck_0 "<<particle_info.spacepoints_q_sec_bck_0<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_bck_1 "<<particle_info.spacepoints_q_sec_bck_1<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_bck_2 "<<particle_info.spacepoints_q_sec_bck_2<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_bck_3 "<<particle_info.spacepoints_q_sec_bck_3<<std::endl;
    std::cout<<"particle_info.spacepoints_q_sec_bck_4 "<<particle_info.spacepoints_q_sec_bck_4<<std::endl;

    std::cout<<"particle_info.n_segments "<<particle_info.n_segments<<std::endl;
    std::cout<<"particle_info.seg_total_len "<<particle_info.seg_total_len<<std::endl;
    std::cout<<"particle_info.seg_total_q "<<particle_info.seg_total_q<<std::endl;
    std::cout<<"particle_info.sec_seg_len "<<particle_info.sec_seg_len<<std::endl;
    std::cout<<"particle_info.sec_seg_direct_len "<<particle_info.sec_seg_direct_len<<std::endl;
    std::cout<<"particle_info.sec_seg_dist_to_main "<<particle_info.sec_seg_dist_to_main<<std::endl;
    std::cout<<"particle_info.flag_start_end_swapped "<<particle_info.flag_start_end_swapped<<std::endl;
}

#endif
