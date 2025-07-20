#include "WCPLEEANA/tree_wrangler.h"
#include "WCPLEEANA/Util.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TROOT.h"

#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

using namespace LEEana;

LEEana::tree_wrangler::tree_wrangler(bool configure, std::string config_file_name, char delimiter, bool set_verbose){
  verbose = set_verbose;
  if(configure){
    std::ifstream config_file(config_file_name);
    if(config_file.is_open()){
      std::cout<<"Loading directories and trees based off the configuration in "<<config_file_name<<std::endl;
      while(!config_file.eof()){
        std::string line;
        std::getline(config_file, line);
        while(line.empty()) {std::getline(config_file, line);}//Skip empty lines
        if(line == "end" || line == "End") break;
        std::istringstream iss(line);
        std::string temp_dir;
        std::string temp_all_trees_to_skip;
        std::vector<std::string> temp_trees_to_skip;
        iss >> temp_dir >> temp_all_trees_to_skip;//Read the line
        temp_trees_to_skip = splitString(temp_all_trees_to_skip, delimiter);//Break up the trees
        if (directories_wi_trees_to_skip_names.find(temp_dir) != directories_wi_trees_to_skip_names.end()) {std::cout<<"WARNING: Failed to open the config file, just loading WC"<<std::endl;}
        sort( temp_trees_to_skip.begin(), temp_trees_to_skip.end() );
        const bool hasDuplicates = std::adjacent_find(temp_trees_to_skip.begin(), temp_trees_to_skip.end()) != temp_trees_to_skip.end();
        if(hasDuplicates){
          std::cout<<"Duplicates in trees to skip for "<<temp_dir<<". You should check the config file."<<std::endl;
          temp_trees_to_skip.erase( unique( temp_trees_to_skip.begin(), temp_trees_to_skip.end() ), temp_trees_to_skip.end() );
        }
        directories_wi_trees_to_skip_names[temp_dir] = temp_trees_to_skip;
      }
    }else {std::cout<<"WARNING: Failed to open the config file: " <<config_file_name<<". Just loading WC"<<std::endl;}
  }//configure
}

LEEana::tree_wrangler::~tree_wrangler(){}


void LEEana::tree_wrangler::CopyDir(TDirectory *source, bool blank_tree, std::vector<std::string> to_skip) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  if(verbose) source->ls();
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
     const char *classname = key->GetClassName();
     TClass *cl = gROOT->GetClass(classname);
     if (!cl) continue;
     if (cl->InheritsFrom(TDirectory::Class())) {
        source->cd(key->GetName());
        TDirectory *subdir = gDirectory;
        adir->cd();
        CopyDir(subdir);
        adir->cd();
     } else if (cl->InheritsFrom(TTree::Class())) {
        TTree *T = (TTree*)source->Get(key->GetName());
        std::string temp_name = T->GetName();
        if (std::find(to_skip.begin(), to_skip.end(), temp_name) != to_skip.end()) continue;
        adir->cd();
        int nentry = -1;
        if (blank_tree){nentry=0;}
        TTree *newT = T->CloneTree(nentry,"fast");
        newT->Write();
     } else {
        source->cd();
        TObject *obj = key->ReadObj();
        adir->cd();
        obj->Write();
        delete obj;
    }
 }
 adir->SaveSelf(kTRUE);
 savdir->cd();
}

void LEEana::tree_wrangler::CopyDir(TDirectory *source, TString TDirectory_name, bool blank_tree, std::vector<std::string> to_skip) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  if(verbose) source->ls();
  TDirectory *savdir = gDirectory;
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
     const char *classname = key->GetClassName();
     TClass *cl = gROOT->GetClass(classname);
     if (!cl) continue;
     if (cl->InheritsFrom(TTree::Class())) {
        TTree *T = (TTree*)source->Get(key->GetName());
        std::string temp_name = T->GetName();
        if (std::find(to_skip.begin(), to_skip.end(), temp_name) != to_skip.end()) continue;
        savdir->cd();
        int nentry = -1;
        if (blank_tree){nentry=0;}
        TTree *newT = T->CloneTree(nentry,"fast");
        newT->SetObject(key->GetName()+TDirectory_name,key->GetName()+TDirectory_name);
        newT->Write();
    }
 }
 savdir->cd();
}


std::vector<TTree*>* LEEana::tree_wrangler::CopyTrees(TDirectory *source, bool blank_tree,  bool rename, TString TDirectory_name, std::vector<std::string> to_skip) {
  if(verbose) source->ls();
  TDirectory *savdir = gDirectory;
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  std::vector<TTree*>* ttree_vec = new std::vector<TTree*>();
  while ((key = (TKey*)nextkey())) {
     const char *classname = key->GetClassName();
     TClass *cl = gROOT->GetClass(classname);
     if (!cl) continue;
     if (cl->InheritsFrom(TTree::Class())) {
        TTree *T = (TTree*)source->Get(key->GetName());
        std::string temp_name = T->GetName();
        if (std::find(to_skip.begin(), to_skip.end(), temp_name) != to_skip.end()) continue;
       savdir->cd();
        int nentry = -1;
        if (blank_tree) {nentry=0;}
        TTree *newT = T->CloneTree(nentry,"fast");
        if (rename) {newT->SetObject(key->GetName()+TDirectory_name,key->GetName()+TDirectory_name);}
        newT->Write();
        ttree_vec->push_back(newT);
    }
  }
  savdir->cd();
  return ttree_vec;
}

std::vector<TTree*>* LEEana::tree_wrangler::GetTrees(TDirectory *source, std::vector<std::string> to_skip) {
  if(verbose) source->ls();
  TDirectory *savdir = gDirectory;
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  std::vector<TTree*>* ttree_vec = new std::vector<TTree*>();
  while ((key = (TKey*)nextkey())) {
     const char *classname = key->GetClassName();
     TClass *cl = gROOT->GetClass(classname);
     if (!cl) continue;
     if (cl->InheritsFrom(TTree::Class())) {
        TTree *T = (TTree*)source->Get(key->GetName());
        std::string temp_name = T->GetName();
        if (std::find(to_skip.begin(), to_skip.end(), temp_name) != to_skip.end()) continue;
        ttree_vec->push_back(T);
    }
  }
  savdir->cd();
  return ttree_vec;
}


std::vector<TTree*>* LEEana::tree_wrangler::get_old_trees(TFile* file){
   std::vector<TTree*>* old_trees = new std::vector<TTree*>;
   for (const auto& it_directories_trees_names : directories_wi_trees_to_skip_names) {
    std::string directory_name = it_directories_trees_names.first;
    std::vector<std::string> trees_to_skip_names = it_directories_trees_names.second;
    TDirectory *topdir = gDirectory;
    if (file->GetDirectory(directory_name.c_str())){
      file->cd(directory_name.c_str());
      TDirectory *temp_directory = new TDirectory;
      temp_directory = gDirectory;
      std::vector<TTree*>* temp_trees = new std::vector<TTree*>;
      temp_trees = GetTrees(temp_directory,trees_to_skip_names);
      old_trees->insert(old_trees->end(), temp_trees->begin(), temp_trees->end());
      std::tuple<TDirectory*,std::vector<TTree*>*> temp_directory_trees(temp_directory,temp_trees);
      names_wi_directories_and_trees[directory_name] = temp_directory_trees;
      topdir->cd();
    }else{std::cout<<"WARNING: can't find the directory "<<directory_name<<" in the input file. Check the config."<<std::endl;}
  }
  return old_trees;
}

std::vector<TTree*>* LEEana::tree_wrangler::set_new_trees(TFile* file){
  TDirectory *topdirout = gDirectory;
  std::vector<TTree*>* new_trees = new std::vector<TTree*>;
  for (const auto& it_directories_trees_names : directories_wi_trees_to_skip_names) {
    std::string directory_name = it_directories_trees_names.first;
    std::vector<std::string> trees_to_skip_names = it_directories_trees_names.second;
    TDirectory * temp_input_directory = std::get<0>(names_wi_directories_and_trees[directory_name]);
    file->mkdir(directory_name.c_str());
    file->cd(directory_name.c_str());
    std::vector<TTree*>* temp_new_trees = new std::vector<TTree*>;
    temp_new_trees = CopyTrees(temp_input_directory,true,false,"",trees_to_skip_names);
    new_trees->insert(new_trees->end(), temp_new_trees->begin(), temp_new_trees->end());
    topdirout->cd();
  }
  return new_trees;
}
                                                                                                                                                                                                                                                                                                                             
