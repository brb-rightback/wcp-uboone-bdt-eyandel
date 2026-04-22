// cz: code modified from tutorials/tmva/TMVAClassification.C

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

int main( int argc, char** argv )
{
  if (argc < 3) {
    std::cout << "tree_trimmer #input_file #output_file #config.txt" << std::endl;
    return -1;
  }


  TString input_file = argv[1];
  TString out_file = argv[2];

  std::string config_file_name=argv[3];

  char delimiter = ',';

  int flag_overwrite=0;

  int flag_kill_duplicates=1;

  int max_events=std::numeric_limits<int>::max();
  int start_events=0;
  int start_run=0;
  int stop_run=std::numeric_limits<int>::max();

  int skip_cut = 1;
  int flag_numi = 0;
  int flag_data = 0;

  int remove_lantern_fails = 1;

  TString training_list = "";
  string global_file_type = "";
  int flag_keep_only_bdt_train = -1;

  int flag_set_samdef = 0;
  TString samdef="";

  for (Int_t i=4;i!=argc;i++){
    switch(argv[i][1]){
    case 'd':
      delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
      break;
    case 'o':
      flag_overwrite = atoi(&argv[i][2]);
      break;
    case 'k':
      flag_kill_duplicates = atoi(&argv[i][2]);
      if(flag_kill_duplicates==0) { std::cout<<"Allowing duplicates if found."<<'\n'<<std::endl; }
      else if(flag_kill_duplicates==1) { std::cout<<"Will remove duplicates if found."<<'\n'<<std::endl; }
      else if(flag_kill_duplicates==2) { std::cout<<"Will exit if duplicates are found."<<'\n'<<std::endl; }
      else {
        flag_kill_duplicates=1;
        std::cout<<"Unknown -k option, setting to default flag_kill_duplicates=1"<<'\n'<<std::endl; 
      }
      break;
    case 'm':
      max_events = atoi(&argv[i][2]);
      if(max_events<0){
        max_events=0;
      }
      std::cout<<"Will be saving at most "<<max_events<<" events to the output file."<<std::endl;
      std::cout<<"Note that this will ''round up'' to the nearest subrun in order to ensure the POT is correct."<<'\n'<<std::endl;
      break;
    case 'e':
      start_events = atoi(&argv[i][2]);
      std::cout<<"Will start saving to the output file at event "<<start_events<<std::endl;
      std::cout<<"Note that this will ''round down'' to the nearest subrun in order to ensure the POT is correct."<<'\n'<<std::endl;
      break;
    case 'y':
      start_run = atoi(&argv[i][2]);
      std::cout<<"Will start saving to the output file at run "<<start_run<<'\n'<<std::endl;
      break;
    case 'z':
      stop_run = atoi(&argv[i][2]);
      std::cout<<"Will stop saving to the output file at run "<<stop_run<<'\n'<<std::endl;
      break;
    case 's':
      skip_cut = atoi(&argv[i][2]);
      if(skip_cut==0) { std::cout<<"Will keep runs not in the good runs list/."<<'\n'<<std::endl; }
      else if(skip_cut==1) { std::cout<<"Will remove runs not in the good runs list."<<'\n'<<std::endl; }
      else{
        skip_cut=1;
        std::cout<<"Unknown -s option, setting to default skip_cut=1"<<'\n'<<std::endl; 
      }
      break;
    case 'n':
      flag_numi = atoi(&argv[i][2]);
      if(flag_numi==0){ std::cout<<"Good runs list will be the one for BNB."<<'\n'<<std::endl; }
      else if(flag_numi==1){ std::cout<<"Good runs list will be the one for NuMI."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -n option, setting to default flag_numi=0"<<'\n'<<std::endl; 
      }
      break;
    case 'c':
      flag_data = atoi(&argv[i][2]);
      if(flag_data==0){ std::cout<<"Good runs for list will be the one for MC."<<'\n'<<std::endl; }
      else if(flag_data==1){ std::cout<<"Good runs for list will be the one for Data."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -c option, setting to default flag_data=0"<<'\n'<<std::endl; 
      }
      break;
    case 'r':
      remove_lantern_fails = atoi(&argv[i][2]);
      if(remove_lantern_fails==0){ std::cout<<"Will Keep subruns where Lantern container failed."<<'\n'<<std::endl; }
      else if(remove_lantern_fails==1){ std::cout<<"Removing subruns where Lantern container failed."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -r option, setting to default remove_lantern_fails=1"<<'\n'<<std::endl; 
      }
      break;
    case 'l':
      training_list = &argv[i][2];
      std::cout<<"Loading Wire-Cell training list from "<<training_list<<'\n'<<std::endl;
      break;
    case 'g':
      global_file_type = &argv[i][2];
      std::cout<<"Setting Wire-Cell BDT training file type to "<<global_file_type<<'\n'<<std::endl;
      break;
    case 'b':
      flag_keep_only_bdt_train = atoi(&argv[i][2]);
      if(flag_keep_only_bdt_train==-1){ std::cout<<"Saving all runs regardless of Wire-Cell BDT training status."<<'\n'<<std::endl; }
      else if(flag_keep_only_bdt_train==0){ std::cout<<"Only saving subruns that were not used for Wire-Cell BDT training."<<'\n'<<std::endl; }
      else if(flag_keep_only_bdt_train==1){ std::cout<<"Only saving subruns that WERE used for Wire-Cell BDT training."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -b option, setting to default lag_keep_only_bdt_train=-1"<<'\n'<<std::endl; 
      }
      break;
    case 'a':
      flag_set_samdef = 1;
      samdef = &argv[i][2];
      std::cout<<"Saving the following samdef to trees: "<<samdef<<'\n'<<std::endl;
      break;
    }
  }


  if( (training_list=="" || global_file_type=="") && flag_keep_only_bdt_train>=0){
    std::cout<<"Have flag_keep_only_bdt_train>=0, but no file list of file type set."<<std::endl; 
    std::cout<<"Please set -b-1 or add a file list and file type with -l and -g"<<std::endl;
    std::cout<<"Exiting."<<std::endl;
    return 1;
  }


  // Initiate the tree wranglers
  tree_wrangler wrangler(true, config_file_name, delimiter);
  if(flag_set_samdef) wrangler.set_samdef(flag_set_samdef, samdef);
  tree_wrangler wrangler_ex(true, config_file_name, delimiter,2);
  if(flag_set_samdef) wrangler_ex.set_samdef(flag_set_samdef, samdef);
  tree_wrangler wrangler_pot(true, config_file_name, delimiter,1);


  TFile *file1 = new TFile(input_file);


  // Load other trees from directories as specified by the config file
  wrangler.get_old_trees(file1);
  wrangler_ex.get_old_trees(file1);
  wrangler_pot.get_old_trees(file1);


  // Load good runs lists

  std::vector<int>good_run_list_vec = wrangler.get_good_run_list();
  std::set<int> good_runlist_set(good_run_list_vec.begin(), good_run_list_vec.end());

  std::vector<int> low_lifetime_runs = wrangler.get_low_lifetime_runs();
  std::set<int> low_lifetime_set(low_lifetime_runs.begin(), low_lifetime_runs.end());

  std::vector<int> low_neutrino_count_numi_run2RHC = wrangler.get_low_neutrino_count_numi_run2RHC();
  std::set<int> low_neutrino_count_numi_run2RHC_set(low_neutrino_count_numi_run2RHC.begin(), low_neutrino_count_numi_run2RHC.end());


  // Load the BDT subruns if that flag is set
  std::map<string, std::set<std::pair<int, int> > > map_type_run_subrun;
  if(flag_keep_only_bdt_train>=0){
    if (training_list != ""){
      ifstream infile(training_list);
      if (!infile.good()) {
        std::cout<<"Unable to open Wire-Cell BDT training list. Exiting"<<std::endl;
        return 1;
      }
      string tmp_type;
      int run, subrun;
      while(!infile.eof()){
        infile >> tmp_type >> run >> subrun;
        map_type_run_subrun[tmp_type].insert(std::make_pair(run, subrun));
      }
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


  TFile *file2 = new TFile(out_file,"RECREATE");


  //Setup the directories specified in the config file
  wrangler.set_new_trees(file2);
  wrangler_ex.set_new_trees(file2);
  wrangler_pot.set_new_trees(file2);


  // Build the pairs of pot trees
  wrangler_pot.grow_pot_arboretum();


  // Load Lantern, used when dropping subruns where container failed.
  int haveReco;
  TTree *T_lantern = (TTree*)file1->Get("lantern/EventTree");
  if(T_lantern && remove_lantern_fails==1) T_lantern->SetBranchAddress("haveReco",&haveReco);
  if(!T_lantern && remove_lantern_fails==1){
    std::cout<<"WARNING: remove_lantern_fails==1, but Lantern tree not found."<<'\n'<<std::endl;
  }

  // Figure out which tree you can load for the rse map
  int run;
  int subrun;
  int event;
  bool found_rse_tree;
  TTree *T_rse; 
  while(1){
  if(file1->GetDirectory("wcpselection")){ 
    T_rse=(TTree*)file1->Get("wcpselection/T_eval");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun") && T_rse->GetBranch("event")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
	T_rse->SetBranchAddress("subrun",&subrun);
        T_rse->SetBranchAddress("event",&event);
        std::cout<<'\n'<<"Using T_eval for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    } 
    T_rse=(TTree*)file1->Get("wcpselection/T_PFeval");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun") && T_rse->GetBranch("event")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        T_rse->SetBranchAddress("event",&event);
        std::cout<<'\n'<<"Using T_PFeval for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
    T_rse=(TTree*)file1->Get("wcpselection/T_KINEvars");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun") && T_rse->GetBranch("event")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        T_rse->SetBranchAddress("event",&event);
        std::cout<<'\n'<<"Using T_KINEvars for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
    T_rse=(TTree*)file1->Get("wcpselection/T_BDTvars");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun") && T_rse->GetBranch("event")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        T_rse->SetBranchAddress("event",&event);
        std::cout<<'\n'<<"Using T_BDTvars for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
  }
  if(file1->GetDirectory("nuselection")){
    T_rse=(TTree*)file1->Get("nuselection/NeutrinoSelectionFilter");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("sub") && T_rse->GetBranch("evt")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("sub",&subrun);
        T_rse->SetBranchAddress("evt",&event);
        std::cout<<'\n'<<"Using NeutrinoSelectionFilter for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("singlephotonana")){
    T_rse=(TTree*)file1->Get("singlephotonana/vertex_tree");
    if(T_rse){
      if(T_rse->GetBranch("run_number") && T_rse->GetBranch("subrun_number") && T_rse->GetBranch("event_number")){
        found_rse_tree=true;
        T_rse->SetBranchAddress("run_number",&run);
        T_rse->SetBranchAddress("subrun_number",&subrun);
        T_rse->SetBranchAddress("event_number",&event);
        std::cout<<'\n'<<"Using vertex_tree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("lantern")){
    T_rse=(TTree*)file1->Get("lantern/EventTree");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun") && T_rse->GetBranch("event")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        T_rse->SetBranchAddress("event",&event);
        std::cout<<'\n'<<"Using EventTree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  std::cout<<'\n'<<"Could not find RS tree. Exiting."<<std::endl;
  return 1;
  }


  // Map out the relation between index and run-subrun for POT trees
  wrangler_pot.map_rs_to_entry();


  // Create run-subrun-index map and exlude and bad subruns.

  std::set<std::tuple<int, int, int, int>> run_sub_event_entry;

  int nentry = T_rse->GetEntries();

  std::cout<<"Starting first pass loop to order run-subruns. Will pass over "<<nentry<<" entries."<<std::endl;

  for (Int_t i=0;i!=nentry;i++){

    if (i%10000 == 0) std::cout << i/1000 << " k " << double(i)/nentry*100. << " %"<< std::endl;

    T_rse->GetEntry(i);

    if(run>=stop_run) continue;
    if(run<start_run) continue;

    // Check if the Lantern container failed on this event, if so throw out the subrun if the flag is set.
    if(T_lantern) T_lantern->GetEntry(i);
    if(remove_lantern_fails==1 && haveReco==0){
      continue;
    }

    // Cutting or including WireCell BDT training subruns, if the flag is set.
    if (flag_keep_only_bdt_train>=0){
      auto it1 = map_type_run_subrun.find(global_file_type);
      if (it1 != map_type_run_subrun.end()){
        // removing run-subruns used to train the BDTs
        if ( it1->second.find(std::make_pair(run, subrun)) != it1->second.end() && flag_keep_only_bdt_train==0 ) {
          continue;
	}
        // removing run-subruns NOT used to train the BDTs
        if ( it1->second.find(std::make_pair(run, subrun)) == it1->second.end() && flag_keep_only_bdt_train==1 ) {
          continue;
	}
      }
      // removing run-subruns NOT used to train the BDTs
      else if(it1 == map_type_run_subrun.end() && flag_keep_only_bdt_train==1){
        continue;
      }
    }

    run_sub_event_entry.insert(std::make_tuple(run, subrun, event, i));

  }


  // Copy all the event level trees to the new file. 

  std::set<std::tuple<int, int, int, int>> run_sub_event_entry_drain = run_sub_event_entry;

  if(start_events>>T_rse->GetEntries()){
    std::cout<<'\n'<<"start_events>T_rse->GetEntries(). Exiting."<<std::endl;
    return 1;
  }

  std::cout<<'\n'<<'\n'<<"Begin looping over "<<T_rse->GetEntries()<<" events"<<std::endl;

  int index_counter=0;
  int event_counter=0;
  int first_run=std::numeric_limits<int>::max();
  int first_subrun=std::numeric_limits<int>::max();
  int last_run=0;
  int last_subrun=0;

  int final_run=-1;
  int final_subrun=-1;

  std::map<std::tuple<int,int,int>, int> counts;

  for (const auto& t : run_sub_event_entry) {
    counts[std::make_tuple(std::get<0>(t), std::get<1>(t), std::get<2>(t))]++;
  }

  for (auto it = run_sub_event_entry.begin(); it != run_sub_event_entry.end(); it++){

    // Get the right tree entry from the map
    int index = std::get<3>(*it);

    T_rse->GetEntry(index);

    // Checking for duplicates.
    if( counts[std::make_tuple(run, subrun, event)] !=1 ){
      if(counts[std::make_tuple(run, subrun, event)]<1){
        std::cout<<"Found duplicate event: run,subrun,event = "<<std::get<0>(*it)<<", "<<std::get<1>(*it)<<", "<<std::get<2>(*it)<<std::endl;
        if(flag_kill_duplicates==1){
          std::cout<<"Removing Duplicate."<<std::endl;
          continue;
        }
      }
      if(flag_kill_duplicates==2){
        std::cout<<"Exiting. Can overide this by re-running with -k0"<<std::endl;
        return 1;
      }
      // First time seeing the event, set the counts to -1 so we know about it when we come back around.
      counts[std::make_tuple(run, subrun, event)]=-1;
    }


    if (index_counter%10000 == 0) std::cout << " seen: "<<index_counter<<"    saved: "  << event_counter<< std::endl;


    // Exit when we have finished the last subrun
    if(final_run>=0 && final_subrun>=0 && (final_run!=run || final_subrun!=subrun)) break;

    if(run>=stop_run){ 
      std::cout<<"Reached the specified maximum run of "<<stop_run<<". Exiting event loop having saved "<<event_counter<<" events."<<'\n'<<std::endl;
      break;
    }

    if(run<start_run) continue;

    // Remove bad run-subruns if the flag is set.
    if (flag_data && skip_cut == 0){
      if (good_runlist_set.find(run) == good_runlist_set.end()) continue;
      if (low_lifetime_set.find(run) != low_lifetime_set.end()) continue;
      if (flag_numi && low_neutrino_count_numi_run2RHC_set.find(run) != low_neutrino_count_numi_run2RHC_set.end()) continue;
      // ext bnb in run 1, high rate
      if (run>=7004 && run <=7070) continue;
      // bnb run 2 high rate
      if (run >= 8321 && run <=8404) continue;
      // bnb run 3 high rate
      if (run >=15369 && run <= 15402) continue;
    }
    if (skip_cut == 0){
      // low lifetime, docdb 39787
      if (run >= 19753 && run <= 19850) continue;
      // low lifetime, docdb 40093
      if (run >= 25447 && run <= 25512) continue;
    }

    index_counter+=1;
    if(index_counter<=start_events) continue;

    // Now fill all the trees.
    for(auto tree_it=wrangler.old_trees->begin(); tree_it!=wrangler.old_trees->end(); tree_it++){
      (*tree_it)->GetEntry(index);
    }
    for(auto tree_it=wrangler.new_trees->begin(); tree_it!=wrangler.new_trees->end(); tree_it++){
      (*tree_it)->Fill();
    }
    for(auto tree_it=wrangler_ex.old_trees->begin(); tree_it!=wrangler_ex.old_trees->end(); tree_it++){
      (*tree_it)->GetEntry(index);
    }
    for(auto tree_it=wrangler_ex.new_trees->begin(); tree_it!=wrangler_ex.new_trees->end(); tree_it++){
      (*tree_it)->Fill();
    }

    event_counter+=1;

    run_sub_event_entry_drain.erase(std::make_tuple(run, subrun, event, index));

    if( (first_run>run) || (first_run==run && first_subrun>subrun) ){
      first_run=run;
      first_subrun=subrun;
    }
    if( (last_run<run) || (last_run==run && last_subrun<subrun) ){
      last_run=run;
      last_subrun=subrun;
    }

    // Set loop to exit when you finish the subrun
    if(event_counter>=max_events) {
      final_run=run;
      final_subrun=subrun;
    }

  }

  // If saving the whole file overwrite limits, otherwise recover the events from the first subrun we started at if that is not complete.
  if(max_events>T_rse->GetEntries() && stop_run>last_run && start_events==0 && start_run==0){
    first_run=-1;
    first_subrun=-1;
    last_run=std::numeric_limits<int>::max();
    last_subrun=std::numeric_limits<int>::max();
  } 
  else{
    auto lower = run_sub_event_entry_drain.lower_bound({first_run, first_subrun, std::numeric_limits<int>::min(),std::numeric_limits<int>::min()});
    auto upper = run_sub_event_entry_drain.upper_bound({first_run, first_subrun, std::numeric_limits<int>::max(),std::numeric_limits<int>::max()});
    for (auto it = lower; it != upper; ++it) {
      std::cout<<"Adding run,subrun,event = "<<std::get<0>(*it)<<", "<<std::get<1>(*it)<<", "<<std::get<2>(*it)<<" from the starting subrun."<<std::endl;
      if( count_matches(run_sub_event_entry_drain,std::get<0>(*it),std::get<1>(*it),std::get<2>(*it)) > 1 ){
        std::cout<<"Found duplicate event: run,subrun,event = "<<std::get<0>(*it)<<", "<<std::get<1>(*it)<<", "<<std::get<2>(*it)<<std::endl;
        if(flag_kill_duplicates==1){
          std::cout<<"Removing Duplicate"<<std::endl;
          continue;
        }
        if(flag_kill_duplicates==2){
          std::cout<<"Exiting. Can overide this by re-running with -k0"<<std::endl;
          return 1;
        }
      }
      // Get the right tree entry from the map
      int index = std::get<3>(*it);
      // Now fill all the trees.
      for(auto tree_it=wrangler.old_trees->begin(); tree_it!=wrangler.old_trees->end(); tree_it++){
        (*tree_it)->GetEntry(index);
      }
      for(auto tree_it=wrangler.new_trees->begin(); tree_it!=wrangler.new_trees->end(); tree_it++){
        (*tree_it)->Fill();
      }
      for(auto tree_it=wrangler_ex.old_trees->begin(); tree_it!=wrangler_ex.old_trees->end(); tree_it++){
        (*tree_it)->GetEntry(index);
      }
      for(auto tree_it=wrangler_ex.new_trees->begin(); tree_it!=wrangler_ex.new_trees->end(); tree_it++){
        (*tree_it)->Fill();
      }
      event_counter+=1;
    }
  }

  // Loop over each POT tree seperatly

  for(auto pot_tree_it=wrangler_pot.pot_arboretum->begin(); pot_tree_it!=wrangler_pot.pot_arboretum->end(); pot_tree_it++){

    int nentries = (*pot_tree_it)->old_pot_tree->GetEntries();
    std::cout<<'\n'<<"Begin looping over "<<(*pot_tree_it)->old_pot_tree->GetName()<<" tree with "<<nentries<<" entries"<<std::endl;
    for (Int_t i=0;i!=nentries;i++){

      if (i%10000 == 0) std::cout << " seen: "<<i<<"    saved: "  << (*pot_tree_it)->new_pot_tree->GetEntries()<< std::endl; 

      (*pot_tree_it)->old_pot_tree->GetEntry(i);

      if( ( (*pot_tree_it)->runNo<first_run ) || ( (*pot_tree_it)->runNo==first_run && (*pot_tree_it)->subRunNo<first_subrun) ) continue;
      if( ( (*pot_tree_it)->runNo>last_run  ) || ( (*pot_tree_it)->runNo==last_run  && (*pot_tree_it)->subRunNo>last_subrun ) ) continue;

      // Remove bad run-subruns if the flag is set.
      if (flag_data && skip_cut == 0){
        if (good_runlist_set.find((*pot_tree_it)->runNo) == good_runlist_set.end()) continue;
        if (low_lifetime_set.find((*pot_tree_it)->runNo) != low_lifetime_set.end()) continue;
        if (flag_numi && low_neutrino_count_numi_run2RHC_set.find((*pot_tree_it)->runNo) != low_neutrino_count_numi_run2RHC_set.end()) continue;
        // ext bnb in run 1, high rate
        if ((*pot_tree_it)->runNo >=7004 && (*pot_tree_it)->runNo <=7070) continue;
        // bnb run 2 high rate
        if ((*pot_tree_it)->runNo >= 8321 && (*pot_tree_it)->runNo <=8404) continue;
        // bnb run 3 high rate
        if ((*pot_tree_it)->runNo >=15369 && (*pot_tree_it)->runNo <= 15402) continue;
      }
      if (skip_cut == 0){
        // low lifetime, docdb 39787
        if ((*pot_tree_it)->runNo >= 19753 && (*pot_tree_it)->runNo <= 19850) continue;
        // low lifetime, docdb 40093
        if ((*pot_tree_it)->runNo >= 25447 && (*pot_tree_it)->runNo <= 25512) continue;
      }

      (*pot_tree_it)->new_pot_tree->Fill();

    }//i, loop over events in a given pot tree set

  }//pot_trees_it, loop over sets of pot trees


  file2->Write("",TFile::kOverwrite);
  file2->Close();


  return 0;


}
