//
// Created by Misha on 3/7/2023.
//

#define QNTOOLS_INCLUDE_PATH "/lustre/home/user/p/parfenov/Soft/QnTools/install-bmn/lib/cmake/QnTools/../../../include/QnTools"

#include <TROOT.h>
#include <stdexcept>
#include <TSystem.h>
#include <TSystem.h>
#include "QnDataFrame.hpp"

int main(int n_args, char** args){
  if( n_args < 2 )
    throw std::runtime_error( "No argumets provided" );
  gInterpreter->AddIncludePath(QNTOOLS_INCLUDE_PATH);
  std::string macro{args[1]};
  std::string macro_full = macro+"(";
  for( int i=2; i<n_args; ++i ){
    macro_full+="\""+ std::string{args[i]}+"\",";
  }
  macro_full.back()=')';
  gROOT->Macro(macro_full.c_str());
  return 0;
}
