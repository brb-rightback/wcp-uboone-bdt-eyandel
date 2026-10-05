#ifndef UBOONE_LEE_GLEE
#define UBOONE_LEE_GLEE

// gLEE (singlephotonana) GENIE truth information, read from singlephotonana/eventweight_tree.
// The tree is entry-aligned with wcpselection/T_eval.

namespace LEEana{
struct GleeInfo{
  Int_t run;
  Int_t subrun;
  Int_t event;

  Int_t GTruth_ResNum; // GENIE resonance number (-1 if not a resonant interaction)
};

 void set_tree_address(TTree *tree0, GleeInfo& glee_info);
 void put_tree_address(TTree *tree0, GleeInfo& glee_info);

}


void LEEana::set_tree_address(TTree *tree0, GleeInfo& glee_info) {
  glee_info.GTruth_ResNum = -1;

  tree0->SetBranchAddress("run", &glee_info.run);
  tree0->SetBranchAddress("subrun", &glee_info.subrun);
  tree0->SetBranchAddress("event", &glee_info.event);

  tree0->SetBranchAddress("GTruth_ResNum", &glee_info.GTruth_ResNum);
}

void LEEana::put_tree_address(TTree *tree0, GleeInfo& glee_info) {
  tree0->Branch("run", &glee_info.run, "run/I");
  tree0->Branch("subrun", &glee_info.subrun, "subrun/I");
  tree0->Branch("event", &glee_info.event, "event/I");

  tree0->Branch("GTruth_ResNum", &glee_info.GTruth_ResNum, "GTruth_ResNum/I");
}

#endif
