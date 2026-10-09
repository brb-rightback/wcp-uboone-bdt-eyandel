#ifndef LEEANA_CMD_OPTIONS_H
#define LEEANA_CMD_OPTIONS_H

#include <cctype>
#include <cstring>

namespace LEEana{
  // Lets the apps take an option value either attached (-xVALUE) or as the next argument (-x VALUE): from argv[first] on,
  // an option with nothing attached ("-x") is joined with the next argument into "-xVALUE", unless the next argument is
  // itself an option (a negative number such as -1 is a value). argv is rewritten in place; returns the new argc.
  // Call it just before the option loop, with first = the index that loop starts at:
  //   argc = LEEana::join_option_values(argc, argv, 3);
  //   for (Int_t i=3;i!=argc;i++){ switch(argv[i][1]){ ... } }
  inline int join_option_values(int argc, char** argv, int first){
    int n = first;
    for (int i=first; i<argc; i++){
      char* opt = argv[i];
      if (opt[0] == '-' && opt[1] != '\0' && opt[2] == '\0' && i+1 < argc){
        const char* val = argv[i+1];
        bool is_value = val[0] != '-' || std::isdigit((unsigned char)val[1]) || (val[1] == '.' && std::isdigit((unsigned char)val[2]));
        if (is_value){
          char* joined = new char[std::strlen(opt) + std::strlen(val) + 1];
          std::strcpy(joined, opt);
          std::strcat(joined, val);
          argv[n++] = joined;
          i++;
          continue;
        }
      }
      argv[n++] = opt;
    }
    if (n < argc) argv[n] = nullptr;
    return n;
  }
}

#endif
