#include "alphaMerger.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

void alphaMerger::SetParameters(TString fname){
  string line;
  ifstream fin;
  fin.open(fname.Data(),std::ifstream::in);
  while(getline(fin,line)){
    stringstream sline(line);
    string index, value;
    sline>>index>>value;
    if(index.empty() || index[0]=='#') continue;
    else if(index=="AlphaMinHGDynQDC") alphaMinHGDynQDC_ = stod(value);
    else if(index=="AlphaMinHGDynEn") alphaMinHGDynEn_ = stod(value);
    else if(index=="AlphaMinLGDynQDC") alphaMinLGDynQDC_ = stod(value);
    else if(index=="AlphaMinLGDynEn") alphaMinLGDynEn_ = stod(value);
    else if(index=="AlphaMaxHGDynQDC") alphaMaxHGDynQDC_ = stod(value);
    else if(index=="AlphaMaxHGDynEn") alphaMaxHGDynEn_ = stod(value);
    else if(index=="AlphaMaxLGDynQDC") alphaMaxLGDynQDC_ = stod(value);
    else if(index=="AlphaMaxLGDynEn") alphaMaxLGDynEn_ = stod(value);
    else if(index=="AlphaFIThresh") alphaFIThresh_ = stod(value);
    else if(index=="AlphaRIThres") alphaRIThres_ = stod(value);
    else if(index=="IonMinHGDynQDC") ionMinHGDynQDC_ = stod(value);
    else if(index=="IonMinHGDynEn") ionMinHGDynEn_ = stod(value);
    else if(index=="IonMinLGDynQDC") ionMinLGDynQDC_ = stod(value);
    else if(index=="IonMinLGDynEn") ionMinLGDynEn_ = stod(value);
    else if(index=="IonMaxHGDynQDC") ionMaxHGDynQDC_ = stod(value);
    else if(index=="IonMaxHGDynEn") ionMaxHGDynEn_ = stod(value);
    else if(index=="IonMaxLGDynQDC") ionMaxLGDynQDC_ = stod(value);
    else if(index=="IonMaxLGDynEn") ionMaxLGDynEn_ = stod(value);
    else if(index=="CorrelationRadius") corrRad = stod(value);
    else if(index=="TimeOffset") t_off = stoi(value); // in units of ns
    else if(index=="TimeWindow") t_window = stoi(value); // in units of ns
  }
  fin.close();
  std::cout<<corrRad<<std::endl;
  return;
}

void alphaMerger::InitAlphaData(){
  infile = new TFile(infilename,"READ");
  TTree *intree = (TTree*)infile->Get("PixTree");
  fReader.SetTree(intree);
  
  Long64_t nentries = fReader.GetEntries();
  cout<<nentries<<endl;
  for(Long64_t i=0;i<nentries;i++){
    if(i%100000==0) cout<<"Entry: "<<i<<", "<<(double(i))/(double(nentries))*100<<"\% finished"<<endl;
    fReader.SetLocalEntry(i);
    if(pspmtsum_vec_->size()==0){
      i++;
      continue;
    }
    // Get variables
    auto pspmtsum_vec = pspmtsum_vec_.Get();
    double dynEnLow_ = pspmtsum_vec->at(0).dynEnergylow;
    double dynQDCLow_ = pspmtsum_vec->at(0).dynQdclow;
    double timeLow_ = pspmtsum_vec->at(0).timelow;
    bool validPosLow_ = pspmtsum_vec->at(0).validPoslow;
    double posXLow_ = pspmtsum_vec->at(0).posXlow;
    double posYLow_ = pspmtsum_vec->at(0).posYlow;
    double dynEnHigh_ = pspmtsum_vec->at(0).dynEnergyhigh;
    double dynQDCHigh_ = pspmtsum_vec->at(0).dynQdchigh;
    double timeHigh_ = pspmtsum_vec->at(0).timehigh;
    bool validPosHigh_ = pspmtsum_vec->at(0).validPoshigh;
    double posXHigh_ = pspmtsum_vec->at(0).posXhigh;
    double posYHigh_ = pspmtsum_vec->at(0).posYhigh;
    // First check if event qualifies as an alpha
    if(dynQDCHigh_>alphaMinHGDynQDC_ && dynQDCHigh_<alphaMaxHGDynQDC_ && dynQDCLow_>alphaMinLGDynQDC_ && dynQDCLow_<alphaMaxLGDynQDC_ && dynEnHigh_>alphaMinHGDynEn_ && dynEnHigh_<alphaMaxHGDynEn_ && dynEnLow_>alphaMinLGDynEn_ && dynEnLow_<alphaMaxLGDynEn_){
      if(validPosHigh_ && posXHigh_!=0.25 && posYHigh_!=0.25){
        mpAlpha[timeHigh_] = std::make_tuple(i,posXHigh_,posYHigh_,dynEnHigh_);
      }
      else if(validPosLow_ && posXLow_!=0.25 && posYLow_!=0.25){
        mpAlpha[timeHigh_] = std::make_tuple(i,posXLow_,posYLow_,dynEnHigh_);
      }
    }
  }
  cout<<"Number of Alphas: "<<mpAlpha.size()<<endl;
  return;
}

void alphaMerger::InitIonData(){
  infile = new TFile(infilename,"READ");
  TTree *intree = (TTree*)infile->Get("PixTree");
  fReaderIon.SetTree(intree);
  
  Long64_t nentries = fReaderIon.GetEntries();
  cout<<nentries<<endl;
  for(Long64_t i=0;i<nentries;i++){
    if(i%100000==0) cout<<"Entry: "<<i<<", "<<(double(i))/(double(nentries))*100<<"\% finished"<<endl;
    fReaderIon.SetLocalEntry(i);
    if(ionpspmtsum_vec_->size()==0){
      i++;
      continue;
    }
    // Get variables
    auto pspmtsum_vec = ionpspmtsum_vec_.Get();
    double dynEnLow_ = pspmtsum_vec->at(0).dynEnergylow;
    double dynQDCLow_ = pspmtsum_vec->at(0).dynQdclow;
    double timeLow_ = pspmtsum_vec->at(0).timelow;
    bool validPosLow_ = pspmtsum_vec->at(0).validPoslow;
    double posXLow_ = pspmtsum_vec->at(0).posXlow;
    double posYLow_ = pspmtsum_vec->at(0).posYlow;
    double dynEnHigh_ = pspmtsum_vec->at(0).dynEnergyhigh;
    double dynQDCHigh_ = pspmtsum_vec->at(0).dynQdchigh;
    double timeHigh_ = pspmtsum_vec->at(0).timehigh;
    bool validPosHigh_ = pspmtsum_vec->at(0).validPoshigh;
    double posXHigh_ = pspmtsum_vec->at(0).posXhigh;
    double posYHigh_ = pspmtsum_vec->at(0).posYhigh;
    // First check if event qualifies as an ion
    if(dynQDCHigh_>ionMinHGDynQDC_ && dynQDCHigh_<ionMaxHGDynQDC_ && dynQDCLow_>ionMinLGDynQDC_ && dynQDCLow_<ionMaxLGDynQDC_ && dynEnHigh_>ionMinHGDynEn_ && dynEnHigh_<ionMaxHGDynEn_ && dynEnLow_>ionMinLGDynEn_ && dynEnLow_<ionMaxLGDynEn_){
      if(validPosLow_ && posXLow_!=0.25 && posYLow_!=0.25){
        mpIon[timeLow_] = std::make_tuple(i,posXLow_,posYLow_,dynEnLow_);
      }
    }
  }
  cout<<"Number of Ions:   "<<mpIon.size()<<endl;
  return;
}


void alphaMerger::Merge(){
  double tmin = 1e-7;
  double tmax = 1e0;
  double lmin = log10(tmin);
  double lmax = log10(tmax);
  double range = lmax-lmin;
  const int numBins = 1000;
  double maxEn = 10000.;
  double ybins[numBins+1];
  double val = range/((double)numBins);

  for(int i=0;i<numBins+1;i++){
    ybins[i] = pow(10.,val*((double)i)+lmin);
  }

  outfile = new TFile(outfilename,"RECREATE");
  outtree = new TTree("mergedAlpha","mergedAlpha");
  //outtree->Branch("pspmtsum_vec_","std::vector<processor_struct::PSPMTSUMMARY",&pspmtsum_data_);
  //outtree->Branch("gammascint_vec_","std::vector<processor_struct::GAMMASCINT",&gammascint_data_);
  //outtree->Branch("pid_vec_","std::vector<processor_struct::PID",&pid_data_);
  outtree->Branch("dT",&dT_,"dT/D");
  outtree->Branch("dr",&dr_,"dr/D");
  outtree->Branch("alphaEnHG",&alphaEnHG_,"alphaEnHG/D");
  outtree->Branch("ionEnLG",&ionEnLG_,"ionEnLG/D");
  //outtree->Branch("alpha_vec_","std::vector<processor_struct::PSPMTSUMMARY",&pspmtsum_data_);
  //outtree->Branch("ion_vec_","std::vector<processor_struct::PSPMTSUMMARY",&ionpspmtsum_data_);

  TH2D *dten = new TH2D("dten","dten",numBins,0,maxEn,numBins,ybins);

  cout<<"Merging with time window: "<<t_window<<", and correlation radius: "<<corrRad<<endl;

  map<double,tuple<Long64_t,double,double,double> >::iterator impIon, impAlpha;

  Long64_t count=0;
  // Want to match all the corresponding ions to a single alpha
  for(impAlpha=mpAlpha.begin();impAlpha!=mpAlpha.end();impAlpha++){
    dT_ = 0.;
    dr_ = 999.;
    alphaEnHG_ = -999.;
    pspmtsum_data_.clear();
    gammascint_data_.clear();
    pid_data_.clear();

    Long64_t entrya = std::get<0>(impAlpha->second);
    double alphaX = std::get<1>(impAlpha->second);
    double alphaY = std::get<2>(impAlpha->second);

    fReader.SetLocalEntry(entrya);
    //auto pspmtsum_vec = pspmtsum_vec_.Get();
    //alphaEnHG_ = pspmtsum_vec->at(0).dynEnergyhigh;
    alphaEnHG_ = std::get<3>(impAlpha->second);

    double ts1 = impAlpha->first+t_off-t_window;
    double ts2 = ts1 + 2*t_window;

    //pspmtsum_data_ = *pspmtsum_vec_.Get();

    //impIon = mpIon.lower_bound(ts1);
    //if(impIon!=mpIon.end() && impIon->first<ts2){
    for(impIon=mpIon.lower_bound(ts1);(impIon!=mpIon.end() && impIon->first<ts2);impIon++){
      //ionpspmtsum_data_.clear();
      ionEnLG_ = -999.;
      Long64_t entryion = std::get<0>(impIon->second);
      double ionX = std::get<1>(impIon->second);
      double ionY = std::get<2>(impIon->second);
      dT_ = (impAlpha->first-impIon->first)*1e-9;
      dr_ = sqrt((alphaX-ionX)*(alphaX-ionX)+(alphaY-ionY)*(alphaY-ionY));

      ////// Need to be sure to add alpha energy to the output + maybe Ion Chamber results
      ////// Could also think about merging with BigRIPS data for PID
      if(dr_<corrRad && dT_!=0){
        if(dT_>0)
          dten->Fill(alphaEnHG_,dT_);
        fReaderIon.SetLocalEntry(entryion);
        ionEnLG_ = std::get<3>(impIon->second);
        //auto ionpspmt_vec = ionpspmtsum_vec_.Get();
        //auto ionpspmt_vec = pspmtsum_vec_.Get();
        //ionpspmtsum_data_ = *ionpspmtsum_vec_.Get();
        //ionEnLG_ = ionpspmt_vec->at(0).dynEnergylow;
        outtree->Fill();
      }
      ionEnLG_ = -999.;
    }
    alphaEnHG_ = -999.;
    if(count%10000==0) cout<<"Entry: "<<count<<", "<<(double(count))/(double(mpAlpha.size()))*100<<"\% finished"<<endl;
    count++;
  }
  dten->Write();
  outtree->Write();
  outfile->Close();
  return;
}
