#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <set>
#include <limits>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "WCPLEEANA/Util.h"

#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

#include "WCPLEEANA/tree_wrangler.h"


using TreeKey = std::pair<std::string,std::string>; // (dir, tree)

std::map<TreeKey, TTree*> build_tree_map(tree_wrangler& w)
{
  std::map<TreeKey, TTree*> out;

  for (const auto& kv : w.names_wi_directories_and_trees) {
    const std::string& dir = kv.first;
    auto trees = std::get<1>(kv.second);

    for (auto* t : *trees) {
      out[std::make_pair(dir, t->GetName())] = t;
    }
  }

  return out;
}


void build_output_trees(tree_wrangler& w1, tree_wrangler& w2, std::vector<TTree*>* new_trees, TFile* out_file){

  auto map1 = build_tree_map(w1);
  auto map2 = build_tree_map(w2);

  std::set<TreeKey> all_keys;
  for (auto& kv : map1) all_keys.insert(kv.first);
  for (auto& kv : map2) all_keys.insert(kv.first);

  for (const auto& key : all_keys) {

    const std::string& dir  = key.first;
    const std::string& name = key.second;

    TTree* t1 = map1.count(key) ? map1[key] : nullptr;
    TTree* t2 = map2.count(key) ? map2[key] : nullptr;

    // Ensure directory exists and go there
    if (!out_file->GetDirectory(dir.c_str())){
      out_file->mkdir(dir.c_str());
    }
    out_file->cd(dir.c_str());

    TTree* out_tree = nullptr;

    // Get the trees only in one file or another
    if (t1 && !t2) {
      out_tree = t1->CloneTree(0);
    } else if (t2) {
      out_tree = t2->CloneTree(0);
    }

    // Merge branches from t2 into t1 if both exist
    if (t1 && t2) {

      // Enable all branches
      t1->SetBranchStatus("*",1);
      t2->SetBranchStatus("*",1);

      // Clone t1
      out_tree = t1->CloneTree(0);

      // Now add t2 branches safely
      TObjArray* b2 = t2->GetListOfBranches();
      for (int i=0;i<b2->GetEntries();i++) {
        TBranch* br = (TBranch*)b2->At(i);
        std::string name = br->GetName();
        if (out_tree->GetBranch(name.c_str())) continue;
        t2->SetBranchAddress(name.c_str(), br->GetAddress());
        out_tree->Branch(name.c_str(), br->GetAddress(), br->GetTitle());
      }

    }

    out_tree->Write("",TTree::kOverwrite);
    new_trees->push_back(out_tree);

  }

}


int main( int argc, char** argv )
{
  if(argc==2 && argv[1][1]=='h'){
    std::cout<<"TODO"<<std::endl;
    return 0;
  }
  else if(argc==2 && argv[1][1]=='H'){
    std::cout<<"TODO"<<std::endl;
    return 0;
  }
  else if (argc < 6) {
    std::cout << "tree_merger #input_file1 #input_file2 #output_file #config1.txt #config2.txt" << std::endl;
    std::cout << "tree_merger -h for further help and instructions." << std::endl;
    return -1;
  }


  TString input_file1 = argv[1];
  TString input_file2 = argv[2];
  TString out_file = argv[3];

  std::string config_file1_name=argv[4];
  std::string config_file2_name=argv[5];

  char delimiter = ',';

  int flag_overwrite=0;

  int flag_mismatch_events=0;

  int set_verbose=10000;

  for (Int_t i=4;i!=argc;i++){
    switch(argv[i][1]){
    case 'd':
      delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
      break;
    case 'o':
      flag_overwrite = atoi(&argv[i][2]);
      break;
    case 'm':
      flag_mismatch_events = atoi(&argv[i][2]);
      if(flag_mismatch_events==0){ std::cout<<"Fail if file1 and file2 have a different number of events."<<'\n'<<std::endl; }
      else if(flag_mismatch_events==1){ std::cout<<"Try to drop events from file2 if it has more events than file1"<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -m option, setting to default flag_mismatch_events=0"<<'\n'<<std::endl;
        flag_mismatch_events=0;
      }
      break;
    case 'v':
      if(atoi(&argv[i][2])>0) set_verbose = atoi(&argv[i][2]);
      break;
    }
  }


  // Check if the output file exists if overwrite is not set.
  if(flag_overwrite!=1){
    TFile *f = TFile::Open(out_file, "READ");
    if (f && !f->IsZombie()) {
        std::cout<<'\n'<< "File exists. Exiting." << std::endl;
        f->Close();
        std::cout<<"Outputfile file "<<out_file<<" already exists and overwrite not set."<<std::endl;
        std::cout<<"Pick a new file name or use -o to force and overwrite of existing file."<<'\n'<<std::endl;
        return 1;
    }
  }


  // Initiate the tree wranglers
  tree_wrangler wrangler1(true, config_file1_name, delimiter);
  tree_wrangler wrangler_ex1(true, config_file1_name, delimiter,2);
  tree_wrangler wrangler_pot1(true, config_file1_name, delimiter,1);

  tree_wrangler wrangler2(true, config_file2_name, delimiter);
  tree_wrangler wrangler_ex2(true, config_file2_name, delimiter,2);
  tree_wrangler wrangler_pot2(true, config_file2_name, delimiter,1);


  TFile *file1 = new TFile(input_file1);
  TFile *file2 = new TFile(input_file2);


  // Load other trees from directories as specified by the config file
  wrangler1.get_old_trees(file1);
  wrangler_ex1.get_old_trees(file1);
  wrangler_pot1.get_old_trees(file1);
  wrangler2.get_old_trees(file2);
  wrangler_ex2.get_old_trees(file2);
  wrangler_pot2.get_old_trees(file2);


  // Figure out which tree we can load RSE from for file 1
  int run1;
  int subrun1;
  int event1;
  bool found_rse_tree1;
  TTree *T_rse1;
  while(1){
  if(file1->GetDirectory("wcpselection")){
    T_rse1=(TTree*)file1->Get("wcpselection/T_eval");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("subrun") && T_rse1->GetBranch("event")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("subrun",&subrun1);
        T_rse1->SetBranchAddress("event",&event1);
        std::cout<<'\n'<<"Using T_eval for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse1=(TTree*)file1->Get("wcpselection/T_PFeval");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("subrun") && T_rse1->GetBranch("event")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("subrun",&subrun1);
        T_rse1->SetBranchAddress("event",&event1);
        std::cout<<'\n'<<"Using T_PFeval for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse1=(TTree*)file1->Get("wcpselection/T_KINEvars");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("subrun") && T_rse1->GetBranch("event")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("subrun",&subrun1);
        T_rse1->SetBranchAddress("event",&event1);
        std::cout<<'\n'<<"Using T_KINEvars for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse1=(TTree*)file1->Get("wcpselection/T_BDTvars");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("subrun") && T_rse1->GetBranch("event")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("subrun",&subrun1);
        T_rse1->SetBranchAddress("event",&event1);
        std::cout<<'\n'<<"Using T_BDTvars for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("nuselection")){
    T_rse1=(TTree*)file1->Get("nuselection/NeutrinoSelectionFilter");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("sub") && T_rse1->GetBranch("evt")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("sub",&subrun1);
        T_rse1->SetBranchAddress("evt",&event1);
        std::cout<<'\n'<<"Using NeutrinoSelectionFilter for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("singlephotonana")){
    T_rse1=(TTree*)file1->Get("singlephotonana/vertex_tree");
    if(T_rse1){
      if(T_rse1->GetBranch("run_number") && T_rse1->GetBranch("subrun_number") && T_rse1->GetBranch("event_number")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run_number",&run1);
        T_rse1->SetBranchAddress("subrun_number",&subrun1);
        T_rse1->SetBranchAddress("event_number",&event1);
        std::cout<<'\n'<<"Using vertex_tree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("lantern")){
    T_rse1=(TTree*)file1->Get("lantern/EventTree");
    if(T_rse1){
      if(T_rse1->GetBranch("run") && T_rse1->GetBranch("subrun") && T_rse1->GetBranch("event")){
        found_rse_tree1=true;
        T_rse1->SetBranchAddress("run",&run1);
        T_rse1->SetBranchAddress("subrun",&subrun1);
        T_rse1->SetBranchAddress("event",&event1);
        std::cout<<'\n'<<"Using EventTree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  std::cout<<'\n'<<"Could not find RS tree in file1. Exiting."<<std::endl;
  return 1;
  }

  // Figure out which tree we can load RSE from for file2
  int run2;
  int subrun2;
  int event2;
  bool found_rse_tree2;
  TTree *T_rse2;
  while(1){
  if(file2->GetDirectory("wcpselection")){
    T_rse2=(TTree*)file2->Get("wcpselection/T_eval");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("subrun") && T_rse2->GetBranch("event")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("subrun",&subrun2);
        T_rse2->SetBranchAddress("event",&event2);
        std::cout<<'\n'<<"Using T_eval for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse2=(TTree*)file2->Get("wcpselection/T_PFeval");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("subrun") && T_rse2->GetBranch("event")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("subrun",&subrun2);
        T_rse2->SetBranchAddress("event",&event2);
        std::cout<<'\n'<<"Using T_PFeval for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse2=(TTree*)file2->Get("wcpselection/T_KINEvars");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("subrun") && T_rse2->GetBranch("event")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("subrun",&subrun2);
        T_rse2->SetBranchAddress("event",&event2);
        std::cout<<'\n'<<"Using T_KINEvars for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
    T_rse2=(TTree*)file2->Get("wcpselection/T_BDTvars");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("subrun") && T_rse2->GetBranch("event")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("subrun",&subrun2);
        T_rse2->SetBranchAddress("event",&event2);
        std::cout<<'\n'<<"Using T_BDTvars for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file2->GetDirectory("nuselection")){
    T_rse2=(TTree*)file2->Get("nuselection/NeutrinoSelectionFilter");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("sub") && T_rse2->GetBranch("evt")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("sub",&subrun2);
        T_rse2->SetBranchAddress("evt",&event2);
        std::cout<<'\n'<<"Using NeutrinoSelectionFilter for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file2->GetDirectory("singlephotonana")){
    T_rse2=(TTree*)file2->Get("singlephotonana/vertex_tree");
    if(T_rse2){
      if(T_rse2->GetBranch("run_number") && T_rse2->GetBranch("subrun_number") && T_rse2->GetBranch("event_number")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run_number",&run2);
        T_rse2->SetBranchAddress("subrun_number",&subrun2);
        T_rse2->SetBranchAddress("event_number",&event2);
        std::cout<<'\n'<<"Using vertex_tree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("lantern")){
    T_rse2=(TTree*)file1->Get("lantern/EventTree");
    if(T_rse2){
      if(T_rse2->GetBranch("run") && T_rse2->GetBranch("subrun") && T_rse2->GetBranch("event")){
        found_rse_tree2=true;
        T_rse2->SetBranchAddress("run",&run2);
        T_rse2->SetBranchAddress("subrun",&subrun2);
        T_rse2->SetBranchAddress("event",&event2);
        std::cout<<'\n'<<"Using EventTree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  std::cout<<'\n'<<"Could not find RS tree in file1. Exiting."<<std::endl;
  return 1;
  }

  // Check entry numbers make sense.
  int nentry1 = T_rse1->GetEntries();
  int nentry2 = T_rse2->GetEntries();
  if(nentry1!=nentry2){
    std::cout<<"Number of events in each file mismatch. Exiting."<<std::endl;
    if(flag_mismatch_events==0){
      std::cout<<"Exiting. You can overide this behavior with -m1"<<std::endl;
      return 1;
    }
    if(flag_mismatch_events==1 && nentry1<nentry2){
      std::cout<<"OK: file2 has more events than file1, will drop file2 events not in file1."<<std::endl;
      std::cout<<"The program will still fail if file1 is not a subset of file2."<<std::endl;
    }
    else{
      std::cout<<"BAD: file1 has more events than file2, but file1 must be a subset of file2."<<std::endl;
      std::cout<<"Try switching the file order if you expect this is the case."<<std::endl;
      std::cout<<"Exiting"<<std::endl;
      return 1;
    }
  }


  TFile *fileout = new TFile(out_file,"RECREATE");

  std::vector<TTree*>* new_trees = new std::vector<TTree*>;

  build_output_trees(wrangler1, wrangler2, new_trees, fileout);
  build_output_trees(wrangler_ex1, wrangler_ex2, new_trees, fileout);


  // Build RSE map for file 1
  // Map out the relation between run -> subrun -> (event,index) for file 2
  int verbose_counter=0;
  int index_counter=0;
  int event_counter=0;
  std::cout<<"Starting first pass loop over file1. Will pass over "<<nentry1<<" entries."<<std::endl;
  std::unordered_map<int, std::unordered_map<int, std::unordered_map<int,int > > > run_subrun_event_index_file1;
  for (Int_t i=0;i!=nentry1;i++){
    if ((i-verbose_counter)%set_verbose == 0) {
      std::cout << "    seen: "<<i<<"    passed: "  << event_counter << std::endl;
      verbose_counter=(int(int(i)/int(set_verbose)))*set_verbose;
    }
    T_rse1->GetEntry(i);
    run_subrun_event_index_file1[run1][subrun1][event1] = i;
    index_counter=index_counter+1;
    event_counter=event_counter+1;
  }

  // Build RSE map for file 2
  // Map out the relation between run -> subrun -> (event,index) for file 2
  verbose_counter=0;
  index_counter=0;
  event_counter=0;
  std::cout<<"Starting first pass loop over file1. Will pass over "<<nentry2<<" entries."<<std::endl;
  std::unordered_map<int, std::unordered_map<int, std::unordered_map<int,int > > > run_subrun_event_index_file2;
  for (Int_t i=0;i!=nentry2;i++){
    if ((i-verbose_counter)%set_verbose == 0) {
      std::cout << "    seen: "<<i<<"    passed: "  << event_counter << std::endl;
      verbose_counter=(int(int(i)/int(set_verbose)))*set_verbose;
    }
    T_rse2->GetEntry(i);
    run_subrun_event_index_file2[run2][subrun2][event2] = i;
    index_counter=index_counter+1;
    event_counter=event_counter+1;
  }


  fileout->cd();


  // Loop over file 1 and do the filling

  std::cout<<'\n'<<'\n'<<"Begin looping over "<<nentry1<<" events to fill trees."<<std::endl;
  verbose_counter=0;
  index_counter=0;
  event_counter=0;

  for (auto r_it = run_subrun_event_index_file1.begin(); r_it != run_subrun_event_index_file1.end(); r_it++){

    int this_run = (*r_it).first;
    if (run_subrun_event_index_file2.find(this_run) == run_subrun_event_index_file2.end()) {
      std::cout<<"file 2 is missing run = "<<this_run<<". Exiting."<<std::endl;
      return 1;
    }

    for (auto s_it = (*r_it).second.begin(); s_it != (*r_it).second.end(); s_it++){

      int this_subrun = (*s_it).first;
      if (run_subrun_event_index_file2[this_run].find(this_subrun) == run_subrun_event_index_file2[this_run].end()) {
        std::cout<<"file 2 is missing run,subrun = "<<this_run<<", "<<this_subrun<<". Exiting."<<std::endl;
        return 1;
      }

      for (auto e_it = (*s_it).second.begin(); e_it != (*s_it).second.end(); e_it++){

      	int this_event = (*e_it).first;
        if (run_subrun_event_index_file2[this_run][this_subrun].find(this_event) == run_subrun_event_index_file2[this_run][this_subrun].end()) {
          std::cout<<"file 2 is missing run,subrun,event = "<<this_run<<", "<<this_subrun<<", "<<this_event<<". Exiting."<<std::endl;
          return 1;
        }

        int this_index1 = (*e_it).second;
        int this_index2 = run_subrun_event_index_file2[this_run][this_subrun][this_event];

        // Now fill all the trees.
        for(auto tree_it=wrangler1.old_trees->begin(); tree_it!=wrangler1.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index1);
        }
        for(auto tree_it=wrangler2.old_trees->begin(); tree_it!=wrangler2.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index2);
        }
        for(auto tree_it=wrangler_ex1.old_trees->begin(); tree_it!=wrangler_ex1.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index1);
        }
        for(auto tree_it=wrangler_ex2.old_trees->begin(); tree_it!=wrangler_ex2.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index2);
        }
        for(auto tree_it=new_trees->begin(); tree_it!=new_trees->end(); tree_it++){
          (*tree_it)->Fill();
        }
      
        if ((index_counter-verbose_counter)%set_verbose == 0) {
          std::cout << "    seen: "<<index_counter<<"    saved: "  << event_counter<< std::endl;
          verbose_counter=(int(int(index_counter)/int(set_verbose)))*set_verbose;
        }
	index_counter=index_counter+1;
	event_counter=event_counter+1;

      } // e_it

    } // s_it

  } // r_it


  // Directly copy POT trees from file1, then file2.
  // Prioritize file1, if both trees contain a tree, it is taken from file1

  auto map1 = build_tree_map(wrangler_pot1);
  auto map2 = build_tree_map(wrangler_pot2);

  std::set<TreeKey> all_keys;
  for (auto& kv : map1) all_keys.insert(kv.first);
  for (auto& kv : map2) all_keys.insert(kv.first);

  for (const auto& key : all_keys) {

    const std::string& dir  = key.first;
    const std::string& name = key.second;

    TTree* t1 = map1.count(key) ? map1[key] : nullptr;
    TTree* t2 = map2.count(key) ? map2[key] : nullptr;

    // Ensure directory exists and go there
    if (!fileout->GetDirectory(dir.c_str())){
      fileout->mkdir(dir.c_str());
    }
    fileout->cd(dir.c_str());

    TTree* out_tree = nullptr;

    // Get the trees only in one file or another
    if (t1) {
      out_tree = t1->CloneTree(-1);
    } else if (t2) {
      out_tree = t2->CloneTree(-1);
    }

    out_tree->Write("",TTree::kOverwrite);

  }


  fileout->Write("",TFile::kOverwrite);
  fileout->Close();


  std::cout<<'\n'<<"Saving output file: "<<out_file<<'\n'<<std::endl;


  return 0;

}
