// cz: code modified from tutorials/tmva/TMVAClassification.C

//#include <cstdlib>
//#include <iomanip>
#include <iostream>
//#include <fstream>
#include <map>
#include <string>
#include <set>
#include <limits>

//#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
//#include "TMath.h"
//#include "TKey.h"

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
      std::cout<<"Will remove bad runs."<<'\n'<<std::endl;
      break;
    case 'n':
      flag_numi = atoi(&argv[i][2]);
      std::cout<<"Will remove bad runs for NuMI."<<'\n'<<std::endl;
      break;
    case 'c':
      flag_numi = atoi(&argv[i][2]);
      std::cout<<"Will remove bad runs for Data."<<'\n'<<std::endl;
      break;
    case 'r':
      remove_lantern_fails = atoi(&argv[i][2]);
      break;
    case 'l':
      training_list = &argv[i][2];
      break;
    case 'g':
      global_file_type = &argv[i][2];
      break;
    case 'b':
      flag_keep_only_bdt_train = atoi(&argv[i][2]);
      break;
    case 'a':
      flag_set_samdef = 1;
      samdef = &argv[i][2];
      break;
    }
  }


  // Initiate the tree wranglers
  tree_wrangler wrangler(true, config_file_name, delimiter);
  tree_wrangler wrangler_ex(true, config_file_name, delimiter,2);
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
        std::cout<<"Pick a new file name or use -o to force and overwrite of existing file."<<std::endl;
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


  // Figure out which tree you can load for the rse map
  int run;
  int subrun;
  bool found_rse_tree;
  TTree *T_rse; 
  while(1){
  if(file1->GetDirectory("wcpselection")){ 
    T_rse=(TTree*)file1->Get("wcpselection/T_eval");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
	T_rse->SetBranchAddress("subrun",&subrun);
        std::cout<<'\n'<<"Using T_eval for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    } 
    T_rse=(TTree*)file1->Get("wcpselection/T_PFeval");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        std::cout<<'\n'<<"Using T_PFeval for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
    T_rse=(TTree*)file1->Get("wcpselection/T_KINEvars");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        std::cout<<'\n'<<"Using T_KINEvars for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
    T_rse=(TTree*)file1->Get("wcpselection/T_BDTvars");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        std::cout<<'\n'<<"Using T_BDTvars for the run-subrun tree"<<'\n'<<std::endl;
	break;
      }
    }
  }
  if(file1->GetDirectory("nuselection")){
    T_rse=(TTree*)file1->Get("nuselection/NeutrinoSelectionFilter");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("sub")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("sub",&subrun);
        std::cout<<'\n'<<"Using NeutrinoSelectionFilter for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("lantern")){
    T_rse=(TTree*)file1->Get("lantern/EventTree");
    if(T_rse){
      if(T_rse->GetBranch("run") && T_rse->GetBranch("subrun")){
        found_rse_tree=true;
	T_rse->SetBranchAddress("run",&run);
        T_rse->SetBranchAddress("subrun",&subrun);
        std::cout<<'\n'<<"Using EventTree for the run-subrun tree"<<'\n'<<std::endl;
        break;
      }
    }
  }
  if(file1->GetDirectory("singlephotonana")){
    T_rse=(TTree*)file1->Get("singlephotonana/vertex_tree");
    if(T_rse){
      if(T_rse->GetBranch("run_number") && T_rse->GetBranch("subrun_number")){
        found_rse_tree=true;
        T_rse->SetBranchAddress("run_number",&run);
        T_rse->SetBranchAddress("subrun_number",&subrun);
        std::cout<<'\n'<<"Using vertex_tree for the run-subrun tree"<<'\n'<<std::endl;
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

  std::set<std::tuple<int, int, int>> run_sub_entry;

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

    run_sub_entry.insert(std::make_tuple(run, subrun, i));

  }


  // Copy all the event level trees to the new file. 

  if(start_events>>T_rse->GetEntries()){
    std::cout<<'\n'<<"start_events>T_rse->GetEntries(). Exiting."<<std::endl;
    return 1;
  }

  std::cout<<'\n'<<'\n'<<"Begin looping over "<<T_rse->GetEntries()<<" events"<<std::endl;

  int index_counter=0;
  int event_counter=0;
  int first_run=-1;
  int first_subrun=-1;
  int last_run=0;
  int last_subrun=0;

  for (auto it = run_sub_entry.begin(); it != run_sub_entry.end(); it++){

    if (event_counter%10000 == 0) std::cout << event_counter/1000 << " k " << std::endl;

    // Get the right tree entry from the map
    int index = std::get<2>(*it); 

    T_rse->GetEntry(index);

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

    if( (first_run>run) || (first_run==run && first_subrun>subrun) ){
      first_run=run;
      first_subrun=subrun;
    }
    if( (last_run<run) || (last_run==run && last_subrun<subrun) ){
      last_run=run;
      last_subrun=subrun;
    }

    if(event_counter>=max_events) break;

  }


  // Loop over each POT tree seperatly

  for(auto pot_tree_it=wrangler_pot.pot_arboretum->begin(); pot_tree_it!=wrangler_pot.pot_arboretum->end(); pot_tree_it++){

    int nentries = (*pot_tree_it)->old_pot_tree->GetEntries();
    std::cout<<'\n'<<"Begin looping over "<<(*pot_tree_it)->old_pot_tree->GetName()<<" tree with "<<nentries<<" entries"<<std::endl;
    for (Int_t i=0;i!=nentries;i++){

      if (i%10000 == 0) std::cout << i/1000 << " k " << std::endl;

      (*pot_tree_it)->old_pot_tree->GetEntry(i);

      if( ( (*pot_tree_it)->runNo<first_run ) || ( (*pot_tree_it)->runNo==first_run && (*pot_tree_it)->subRunNo<first_subrun) ) continue;
      if( ( (*pot_tree_it)->runNo>last_run  ) || ( (*pot_tree_it)->runNo==last_run  && (*pot_tree_it)->subRunNo>last_subrun ) ) continue;
//	std::cout<<"Passed the last run-subrun seen by the event tree. Done saving from this POT tree."<<std::endl;
      //  break;
      //}

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
