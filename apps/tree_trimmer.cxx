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


struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};
using Key = std::pair<int,int>;
using Value = std::pair< bool, std::unordered_map<int,std::vector<int>> >;


void print_help() {
  std::cout << R"(

========================================
 tree_trimmer : Help
========================================

Overview:
---------
tree_trimmer reads a ROOT file containing multiple directories and TTrees and produces a
filtered output file. It organizes events by (run, subrun, event), removes
duplicates (optionally), applies data-quality and selection cuts, and ensures
POT consistency by operating at the subrun level.

The program:
  - Loads trees defined in a config file (run tree_trimmer -H for more info)
  - Optionally filters based on Wire-Cell BDT training lists
  - Removes bad runs and detector-quality failures
  - Writes a trimmed ROOT file with consistent POT accounting

Usage:
------
  tree_trimmer <input_file> <output_file> <config.txt> [options]

Required arguments:
-------------------
  input_file     Input ROOT file
  output_file    Output ROOT file
  config.txt     Tree configuration file

Options:
--------

  -h
      Show this help message and exit

  -H
      Show help message for configuration file and exit

  -d<char>
      Set delimiter used in config file (default: ',')

  -o<int>
      Overwrite output file if it exists
        0 = do not overwrite (default)
        1 = overwrite

  -k<int>
      Duplicate handling:
        0 = allow duplicates
        1 = remove duplicates (default)
        2 = exit if duplicates are found

  -m<int>
      Maximum number of events to save
      (rounded up to nearest subrun for correct POT)

  -e<int>
      Start saving from this event index
      (rounded down to nearest subrun)

  -y<int>
      Start saving events at this run number

  -z<int>
      Stop saving events at this run number

  -w<int>
      Start saving events at this subrun number (within the starting run)

  -x<int>
      Stop saving events at this subrun number (within the stopping run)

  -s<int>
      Skip data-quality cuts:
        0 = keep all runs (ignore good run list)
        1 = apply good run list cuts (default)

  -n<int>
      Beam selection:
        0 = BNB (default)
        1 = NuMI

  -c<int>
      Data type:
        0 = MC (default)
        1 = Data

  -r<int>
      Lantern failure handling:
        0 = keep subruns where Lantern failed
        1 = remove subruns where Lantern failed (default)

  -l<string>
      Path to Wire-Cell BDT training list file

  -g<string>
      Global file type label used with training list

  -b<int>
      Filter based on BDT training usage:
        -1 = keep all subruns (default)
         0 = keep only subruns NOT used in training
         1 = keep only subruns used in training

      NOTE: Requires both -l and -g

  -a<string>
      Set SAM definition string to be stored in output trees

  -v<int>
      Verbosity level for event loop (default: 10000)

  -u<int>
      Verbosity level for POT loop (default: 1000)

Configuration File:
-------------------
run tree_trimmer -H for more info

Examples:
---------

  Basic usage:
    tree_trimmer input.root output.root config.txt

  Overwrite output and limit to 10000 events:
    tree_trimmer input.root output.root config.txt -o1 -m10000

  Keep only BDT training subruns:
    tree_trimmer input.root output.root config.txt -ltrain.txt -gmytype -b1

  Run on data with quality cuts:
    tree_trimmer input.root output.root config.txt -c1 -s1

Notes:
------
- Event selection is applied at the subrun level to preserve POT consistency.
- Duplicate events are identified via (run, subrun, event).
- Some run ranges are hard-coded as bad and will always be removed unless -s0.

)";
}


void print_help_config() {
  std::cout << R"(

Configuration File:
-------------------

The configuration file controls:
  - Which directories and TTrees are processed
  - Which trees are skipped
  - Which branches (variables) are copied

The file is read line-by-line and divided into sections. Each line
generally has the format:

  <directory> <tree_list> [additional_fields...]

Where:
  directory   = ROOT directory name
  tree_list   = comma-separated list (delimiter configurable with -d)

------------------------------------------------------------
Sections:
------------------------------------------------------------

The configuration file can contain up to three sections:

  (1) Default section (top of file)
  (2) exclusive / pot section
  (3) pick section

The parser stops or switches behavior when encountering:

  exclusive / Exclusive / pot / POT
  pick / Pick
  end / End

----------------------------------------
1. Default section
----------------------------------------

Copies all trees in the specified directory. 
The only ones that are skipped are the ones specified in the file.
These trees should all have the same number of events across reconstructions.

Format:
  <directory> <trees_to_skip>

Example:
  nuselection SubRun

Meaning:
  - Load all trees in "nuselection"
  - Skip the tree named "SubRun"

----------------------------------------
2. "exclusive" (or "POT") section
----------------------------------------

This is used for POT trees. These will be looped over seprate by the code.
This accounts for the fact that POT trees can have different number of entries across reconstructions.
Able to specify which barnches will be used to load run/subrun/event information.

Format:
  <directory> <trees> <pot_variables>

Where:
  trees          = list of tree names
  pot_variables  = grouped in triples:
                   (run, subrun, POT) per tree

Example:
  nuselection SubRun run,subRun,POT

Meaning:
  - Use only the specified trees
  - Associate POT information via the listed branches
  - Each tree can have its own (run, subrun, POT) triplet

Special case:
  If "None,None,None" is provided, POT variables are ignored for that tree.

----------------------------------------
3. "pick" section
----------------------------------------

Copies only the specified trees and branches in the specified directory. 
These are otherwise reated the same as the nominal section.


This section is used when:
  tree_wrangler(..., set_flag_exclusive=2)

Format:
  <directory> <trees> <branch_list_per_tree...>

Example:
  wcpselection T_BDTvars,T_eval kine_reco_Enu,event run,subrun,event

Meaning:
  - Only the listed trees are processed
  - For each tree, only selected branches are copied

Important:
  - Each tree must have a corresponding branch list
  - If no branch list is provided, ALL branches are copied (with warning)

----------------------------------------
General Notes:
----------------------------------------

- The delimiter for lists (trees, branches, variables) defaults to ','
  but can be changed with the -d option.

- Empty lines are ignored.

- Duplicate tree names in a line will trigger a warning.

- Parsing stops when "end" or "End" is encountered.

------------------------------------------------------------
Example Configuration:
------------------------------------------------------------

nuselection SubRun

exclusive
nuselection SubRun run,subRun,pot

pick
wcpselection T_BDTvars,T_eval kine_reco_Enu,event run,subrun,event
end

Explanation:

1. Default section:
   - In directory "nuselection", skip the "SubRun" tree, but copy all other trees and branches

2. Exclusive section:
   - For "nuselection", explicitly process "SubRun" as a POT tree
   - Use branches (run, subRun, pot) for POT accounting

3. Pick section:
   - In "wcpselection", process only:
       T_BDTvars and T_eval
   - From T_BDTvars keep on the following branches: kine_reco_Enu, event
   - From T_eval keep on the following branche: run, subrun, event

------------------------------------------------------------


========================================

)";
}


bool keep_subrun(int run, int subrun,
                 int start_run, int start_subrun, int stop_run, int stop_subrun,
                 int remove_lantern_fails, bool haveReco, 
                 int flag_keep_only_bdt_train, const std::string& global_file_type, const std::map<string, std::set<std::pair<int, int>>>& map_type_run_subrun,
                 int skip_cut, int flag_data, int flag_numi,
                 const std::set<int>& good_runlist_set, const std::set<int>& low_lifetime_set, const std::set<int>& low_neutrino_count_numi_run2RHC_set) {

  bool flag_keep_subrun = false;

  if ( ( run<start_run ) || ( run==start_run && subrun<start_subrun) ) {
      return flag_keep_subrun;
  } 
  if ( ( run>stop_run  ) || ( run==stop_run  && subrun>stop_subrun ) ) {
      return flag_keep_subrun;
  } 

  // check if the lantern container failed on this subrun
  if (remove_lantern_fails == 1 && haveReco == 0) {
      return flag_keep_subrun;
  }

  // Cutting or including WireCell BDT training subruns
  if (flag_keep_only_bdt_train >= 0) {
      auto filetype_it = map_type_run_subrun.find(global_file_type);
      if (filetype_it != map_type_run_subrun.end()) {
          auto rs_it = filetype_it->second.find(std::make_pair(run, subrun));
          if (rs_it != filetype_it->second.end() && flag_keep_only_bdt_train == 0) {
              return flag_keep_subrun;
          }
          if (rs_it == filetype_it->second.end() && flag_keep_only_bdt_train == 1) {
              return flag_keep_subrun;
          }
      } else if (flag_keep_only_bdt_train == 1) {
          return flag_keep_subrun;
      }
  }

  // Remove bad run-subruns if the flag is set.
  if (flag_data == 1 && skip_cut == 0) {
      if (good_runlist_set.find(run) == good_runlist_set.end()) {
          return flag_keep_subrun;
      }
      if (low_lifetime_set.find(run) != low_lifetime_set.end()) {
          return flag_keep_subrun;
      }
      if (flag_numi == 1 && low_neutrino_count_numi_run2RHC_set.find(run) != low_neutrino_count_numi_run2RHC_set.end()) {
          return flag_keep_subrun;
      }
      if (run >= 7004 && run <= 7070) {
          return flag_keep_subrun;
      }
      if (run >= 8321 && run <= 8404) {
          return flag_keep_subrun;
      }
      if (run >= 15369 && run <= 15402) {
          return flag_keep_subrun;
      }
  }

  if (skip_cut == 0) {
      if (run >= 19753 && run <= 19850) {
          return flag_keep_subrun;
      }
      if (run >= 25447 && run <= 25512) {
          return flag_keep_subrun;
      }
  }


  flag_keep_subrun = true;
  return flag_keep_subrun;

}


int main( int argc, char** argv )
{
  if(argc==2 && argv[1][1]=='h'){
    print_help();
    return 0;
  }
  else if(argc==2 && argv[1][1]=='H'){
    print_help_config();
    return 0;
  }
  else if (argc < 4) {
    std::cout << "tree_trimmer #input_file #output_file #config.txt" << std::endl;
    std::cout << "tree_trimmer -h for further help and instructions." << std::endl;
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
  int start_subrun=0;
  int stop_run=std::numeric_limits<int>::max();
  int stop_subrun=std::numeric_limits<int>::max();

  int skip_cut = 1;
  int flag_numi = 0;
  int flag_data = 0;

  int remove_lantern_fails = 1;

  TString training_list = "";
  string global_file_type = "";
  int flag_keep_only_bdt_train = -1;

  int flag_set_samdef = 0;
  TString samdef="";

  int set_verbose=10000;
  int set_verbose_pot=1000;

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
        std::cout<<"Unknown -k option, setting to default flag_kill_duplicates=1"<<'\n'<<std::endl; 
        flag_kill_duplicates=1;
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
    case 'w':
      start_subrun = atoi(&argv[i][2]);
      std::cout<<"Will start saving to the output file at subrun "<<start_subrun<<'\n'<<std::endl;
      break;
    case 'x':
      stop_subrun = atoi(&argv[i][2]);
      std::cout<<"Will stop saving to the output file at subrun "<<stop_subrun<<'\n'<<std::endl;
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
        std::cout<<"Unknown -s option, setting to default skip_cut=1"<<'\n'<<std::endl; 
        skip_cut=1;
      }
      break;
    case 'n':
      flag_numi = atoi(&argv[i][2]);
      if(flag_numi==0){ std::cout<<"Good runs list will be the one for BNB."<<'\n'<<std::endl; }
      else if(flag_numi==1){ std::cout<<"Good runs list will be the one for NuMI."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -n option, setting to default flag_numi=0"<<'\n'<<std::endl; 
        flag_numi=0;
      }
      break;
    case 'c':
      flag_data = atoi(&argv[i][2]);
      if(flag_data==0){ std::cout<<"Good runs for list will be the one for MC."<<'\n'<<std::endl; }
      else if(flag_data==1){ std::cout<<"Good runs for list will be the one for Data."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -c option, setting to default flag_data=0"<<'\n'<<std::endl; 
        flag_data=0;
      }
      break;
    case 'r':
      remove_lantern_fails = atoi(&argv[i][2]);
      if(remove_lantern_fails==0){ std::cout<<"Will Keep subruns where Lantern container failed."<<'\n'<<std::endl; }
      else if(remove_lantern_fails==1){ std::cout<<"Removing subruns where Lantern container failed."<<'\n'<<std::endl; }
      else{
        std::cout<<"Unknown -r option, setting to default remove_lantern_fails=1"<<'\n'<<std::endl; 
        remove_lantern_fails=1;
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
        flag_keep_only_bdt_train=-1;
      }
      break;
    case 'a':
      flag_set_samdef = 1;
      samdef = &argv[i][2];
      std::cout<<"Saving the following samdef to trees: "<<samdef<<'\n'<<std::endl;
      break;
    case 'v':
      if(atoi(&argv[i][2])>0) set_verbose = atoi(&argv[i][2]);
      else{
        std::cout<<"Bad -v optioion, must be greater than 0. Leaving at 10000."<<'\n'<<std::endl;
      }
      break;
    case 'u':
      if(atoi(&argv[i][2])>0) set_verbose_pot = atoi(&argv[i][2]);      
      else{
        std::cout<<"Bad -u optioion, must be greater than 0. Leaving at 1000."<<'\n'<<std::endl;
      } 
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
      while(infile >> tmp_type >> run >> subrun){
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

  int nentry = T_rse->GetEntries();
  if(start_events>nentry){
    std::cout<<'\n'<<"start_events>T_rse->GetEntries(). Exiting."<<std::endl;
    return 1;
  }


  // Map out the relation between index and run-subrun for POT trees
  wrangler_pot.map_rs_to_entry();


  // List of runs where Lantern container failed.
  std::unordered_set<std::pair<int,int>, PairHash> lantern_fail;


  // Create run-subrun-index map and exlude and bad subruns.
  // Key is run, subrun pair, each holds a pair, first being if its flagged as a good run, second being a vecotr of pairs of event,tree index.
  std::unordered_map<Key, Value, PairHash> run_sub_event_entry;

  std::cout<<"Starting first pass loop to order run-subruns. Will pass over "<<nentry<<" entries."<<std::endl;

  int verbose_counter=0;
  int event_counter=0;

  for (Int_t i=0;i!=nentry;i++){

    if ((i-verbose_counter)%set_verbose == 0) {
      std::cout << "    seen: "<<i<<"    passed: "  << event_counter << std::endl;
      verbose_counter=(int(int(i)/int(set_verbose)))*set_verbose;
    }

    T_rse->GetEntry(i);

    std::pair<int, int> run_subrun_pair = std::make_pair(run,subrun);

    // Already seen this run-subrun pair, no need to re-check it.
    if(run_sub_event_entry.find(run_subrun_pair)!=run_sub_event_entry.end()){
      auto& event_index_map = run_sub_event_entry[run_subrun_pair].second;
      // Event already exists, duplicate so append the index.
      if(event_index_map.find(event)!=event_index_map.end()){   
        std::cout<<"Found duplicate event: run,subrun,event = "<<run<<", "<<subrun<<", "<<event<<std::endl;
        if(flag_kill_duplicates==2){
          std::cout<<"Exiting. Can overide this by re-running with -k0"<<std::endl;
          return 1;
        }
        run_sub_event_entry[run_subrun_pair].second[event].push_back(i);
        if(run_sub_event_entry[run_subrun_pair].first) event_counter++;
        continue;
      }        
      std::vector<int> index_vector = {i};
      run_sub_event_entry[run_subrun_pair].second[event] = index_vector;
      if(run_sub_event_entry[run_subrun_pair].first) event_counter++;
      continue;
    }

    //                                      std::vector<int>                                                                                  indcies = {i};
    //               std::unordered_map<int,std::vector<int>>                                                       event_index_map = {{event,indcies}};
    //std::pair<bool,std::unordered_map<int,std::vector<int> >> pair_goodrun_event_index_map = std::make_pair(false,event_index_map);

    if(T_lantern) T_lantern->GetEntry(i);
    else haveReco=1;
    if (haveReco==0){
      lantern_fail.insert(std::make_pair(run,subrun));
    }
    bool temp_keep_subrun = keep_subrun(run, subrun, 
                                        start_run, start_subrun, stop_run, stop_subrun,
                                        remove_lantern_fails, haveReco, 
                                        flag_keep_only_bdt_train, global_file_type, map_type_run_subrun,
                                        skip_cut, flag_data, flag_numi,
                                        good_runlist_set, low_lifetime_set, low_neutrino_count_numi_run2RHC_set);

    auto& entry = run_sub_event_entry[run_subrun_pair];
    entry.first = temp_keep_subrun;
    entry.second[event].push_back(i);

    if(temp_keep_subrun) event_counter++;

  }
  
  std::cout << "    seen: "<<nentry<<"    passed: "  << event_counter << std::endl;


  // Turn the map into an ordered map
  std::map<Key, Value> run_sub_event_entry_first_subrun(run_sub_event_entry.begin(), run_sub_event_entry.end());


  // Copy all the event level trees to the new file. 

  std::cout<<'\n'<<'\n'<<"Begin looping over "<<nentry<<" events to fill trees."<<std::endl;

  int first_run=std::numeric_limits<int>::max();
  int first_subrun=std::numeric_limits<int>::max();
  int last_run=0;
  int last_subrun=0;

  int index_counter=0;
  event_counter=0;
  verbose_counter=0;

  for (auto rs_it = run_sub_event_entry_first_subrun.begin(); rs_it != run_sub_event_entry_first_subrun.end(); rs_it++){

    if ((index_counter-verbose_counter)%set_verbose == 0) {
      std::cout << "    seen: "<<index_counter<<"    saved: "  << event_counter<< std::endl;
      verbose_counter=(int(int(index_counter)/int(set_verbose)))*set_verbose;
    }

    int this_run = (*rs_it).first.first;
    int this_subrun = (*rs_it).first.second;
    bool good_subrun = (*rs_it).second.first;
    auto& event_index_map = (*rs_it).second.second;

    index_counter+=event_index_map.size();
    if(index_counter<=start_events) continue;

    // Skip event if the run-subrun got flagged as bad.
    if(!good_subrun){
      continue;
    }

    // Begin loop over events in this subrun
    for (auto e_it = event_index_map.begin(); e_it != event_index_map.end(); e_it++){

      int this_event = (*e_it).first;
      auto& index_vector = (*e_it).second; 
     
      // If allowing duplicates, loop over all duplicates of this event.
      for (int i_it=0; i_it<index_vector.size(); i_it++){
        if(flag_kill_duplicates==1 && i_it>0){
          std::cout<<"Removing Duplicat: run,subrun,event = "<<this_run<<", "<<this_subrun<<", "<<this_event<<"."<<std::endl;
          continue;
        }
        event_counter++;
        int this_index = index_vector.at(i_it);
        // Now fill all the trees.
        for(auto tree_it=wrangler.old_trees->begin(); tree_it!=wrangler.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index);
        }
        for(auto tree_it=wrangler.new_trees->begin(); tree_it!=wrangler.new_trees->end(); tree_it++){
          (*tree_it)->Fill();
        }
        for(auto tree_it=wrangler_ex.old_trees->begin(); tree_it!=wrangler_ex.old_trees->end(); tree_it++){
          (*tree_it)->GetEntry(this_index);
        }
        for(auto tree_it=wrangler_ex.new_trees->begin(); tree_it!=wrangler_ex.new_trees->end(); tree_it++){
          (*tree_it)->Fill();
        }
      } // i_it, loop over duplicates of the event.

    } // e_it, loop over all events in the run-subrun.

    // Make sure you have the first and last run-subrun set as such
    if( (first_run>this_run) || (first_run==this_run && first_subrun>this_subrun) ){
      first_run=this_run;
      first_subrun=this_subrun;
    }
    if( (last_run<this_run) || (last_run==this_run && last_subrun<this_subrun) ){
      last_run=this_run;
      last_subrun=this_subrun;
    }

    // End the loop when you have enough events
    if(event_counter>=max_events){
      break;
    }

  } // rs_it, loop over all run-subruns

  std::cout << "    seen: "<<index_counter<<"    saved: "  << event_counter<< std::endl;


  // If saving the whole file overwrite limits, otherwise recover the events from the first subrun we started at if that is not complete.
  if(max_events>T_rse->GetEntries() && stop_run>last_run && start_events==0 && start_run==0){
    first_run=-1;
    first_subrun=-1;
    last_run=std::numeric_limits<int>::max();
    last_subrun=std::numeric_limits<int>::max();
  } 


  // Loop over each POT tree seperatly

  int rs_counter=0;

  for(auto pot_tree_it=wrangler_pot.pot_arboretum->begin(); pot_tree_it!=wrangler_pot.pot_arboretum->end(); pot_tree_it++){

    int nentry_pot = (*pot_tree_it)->old_pot_tree->GetEntries();
    std::cout<<'\n'<<"Begin looping over "<<(*pot_tree_it)->old_pot_tree->GetName()<<" tree with "<<nentry_pot<<" entries"<<std::endl;

    verbose_counter=0;
    rs_counter=0;

    for (Int_t i=0;i!=nentry_pot;i++){

      if ((i-verbose_counter)%set_verbose_pot == 0) {
        std::cout << "    seen: "<<i<<"    saved: "  << rs_counter<< std::endl;
        verbose_counter=(int(int(i)/int(set_verbose)))*set_verbose_pot;
      }

      (*pot_tree_it)->old_pot_tree->GetEntry(i);
      
      bool temp_haveReco=1;
      if (lantern_fail.find(std::make_pair((*pot_tree_it)->runNo,(*pot_tree_it)->subRunNo)) != lantern_fail.end()) {
        temp_haveReco=0;
      } 

      bool flag_keep_run =   keep_subrun((*pot_tree_it)->runNo, (*pot_tree_it)->subRunNo,
                                                         first_run, first_subrun, last_run, last_subrun, 
                                                         remove_lantern_fails, temp_haveReco, 
                                                         flag_keep_only_bdt_train, global_file_type, map_type_run_subrun,
                                                         skip_cut, flag_data, flag_numi,
                                                         good_runlist_set, low_lifetime_set, low_neutrino_count_numi_run2RHC_set);

      if(!flag_keep_run) continue;

      (*pot_tree_it)->new_pot_tree->Fill();

      rs_counter+=1;

    }//i, loop over events in a given pot tree set

    std::cout << "    seen: "<<nentry_pot<<"    saved: "  << rs_counter << std::endl;

  }//pot_trees_it, loop over sets of pot trees


  file2->Write("",TFile::kOverwrite);
  file2->Close();


  std::cout<<'\n'<<"Saving output file: "<<out_file<<'\n'<<std::endl;


  return 0;


}
