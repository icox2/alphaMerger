#include "alphaMerger.hpp"
#include "TApplication.h"
#include <chrono>
//using namespace std::chrono;

void runAlphaMerger(){
  TString rootfile="/SCRATCH/DScratch7/e21069rootfile/utkscan/devTesting/lyso210PoWithDiode21us_DD.root";
  TString configfile="tstconfig.txt";
  TString outroot = "/SCRATCH/DScratch7/e21069rootfile/utkscan/devTesting/lyso210PoNoDiodeMerge.root";
  alphaMerger am(rootfile,outroot);
  cout<<"Merging file: "<<rootfile<<endl;
  am.SetParameters(configfile);
  am.InitAlphaData();
  am.InitIonData();
  auto start = std::chrono::high_resolution_clock::now();
  am.Merge();
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  cout<<"Merged in "<<duration.count()/1000.<<" seconds."<<endl;
}
