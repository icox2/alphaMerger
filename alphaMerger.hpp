#ifndef ALPHAMERGER
#define ALPHAMERGER

#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "PaassRootStruct.hpp"

class alphaMerger{
  public:
    alphaMerger(TString in, TString out)
    {
      infilename=in;
      outfilename=out;
    };
    virtual ~alphaMerger(){};

    TFile *infile;
    TString infilename;

    int t_off=0; // in units of ns
    int t_window=500;// in units of ns
    double corrRad;
    double dT_;
    double dr_;
    double alphaEnHG_;
    double ionEnLG_;

    // map of {time,[entrynum,posx,posy,dynenergy]}  positions are from HG (or LG) for alphas and LG for ions    
    map<double,tuple<Long64_t,double,double,double> > mpIon;
    map<double,tuple<Long64_t,double,double,double> > mpAlpha;

    TString outfilename;
    TFile *outfile = nullptr;
    TTree *outtree = nullptr;

    TTreeReader fReader; //!reads the ttree
    TTreeReader fReaderIon; //!reads the ttree
    TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain
    TTreeReaderValue<std::vector<processor_struct::PSPMTSUMMARY> > pspmtsum_vec_ = {fReader, "pspmtsum_vec_"};
    TTreeReaderValue<std::vector<processor_struct::GAMMASCINT> > gammascint_vec_ = {fReader, "gammascint_vec_"};
    TTreeReaderValue<std::vector<processor_struct::PID> > pid_vec_ = {fReader,"pid_vec_"};
    TTreeReaderValue<std::vector<processor_struct::PSPMTSUMMARY> > ionpspmtsum_vec_ = {fReaderIon, "pspmtsum_vec_"};
    TTreeReaderValue<std::vector<processor_struct::PID> > ionpid_vec_ = {fReaderIon,"pid_vec_"};
    TTreeReaderArray<Double_t> pspmtsum_vec__dynEnergylow = {fReaderIon, "pspmtsum_vec_.dynEnergylow"};

    std::vector<processor_struct::PSPMTSUMMARY> pspmtsum_data_;
    std::vector<processor_struct::GAMMASCINT> gammascint_data_;
    std::vector<processor_struct::PID> pid_data_;
    std::vector<processor_struct::PSPMTSUMMARY> ionpspmtsum_data_;

    void SetParameters(TString fname);
    void InitAlphaData();
    void InitIonData();
    void Merge();

  private:
    // Parameters for alphas
    double alphaMinHGDynQDC_ = -999;
    double alphaMinHGDynEn_ = -999;
    double alphaMinLGDynQDC_ = -999;
    double alphaMinLGDynEn_ = -999;
    double alphaMaxHGDynQDC_ = 9999999;
    double alphaMaxHGDynEn_ = 9999999;
    double alphaMaxLGDynQDC_ = 9999999;
    double alphaMaxLGDynEn_ = 9999999;
    double alphaFIThresh_ = 99999;
    double alphaRIThres_ = 0;
    double ionMinHGDynQDC_ = -999;
    double ionMinHGDynEn_ = -999;
    double ionMinLGDynQDC_ = -999;
    double ionMinLGDynEn_ = -999;
    double ionMaxHGDynQDC_ = 9999999;
    double ionMaxHGDynEn_ = 9999999;
    double ionMaxLGDynQDC_ = 9999999;
    double ionMaxLGDynEn_ = 9999999;
};


#endif //ALPHAMERGER
