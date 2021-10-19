////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        MyAnalysis_module.cc
//
// Generated at Sun Oct  4 21:09:21 2020 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TTimeStamp.h"

//Sim photons?
#include <cmath>
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"

//OpHits?
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"

//OpDetWaveForms
#include "lardataobj/RawData/OpDetWaveform.h"


class MyAnalysis;


class MyAnalysis : public art::EDAnalyzer {
public:
  explicit MyAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyAnalysis(MyAnalysis const&) = delete;
  MyAnalysis(MyAnalysis&&) = delete;
  MyAnalysis& operator=(MyAnalysis const&) = delete;
  MyAnalysis& operator=(MyAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  bool fdump_OpHits;
  bool fdump_G4_PE;         
  bool fdump_DE;            
  bool fdump_Op_Waveforms;  
  
  TH1F* fHist;  //!< Output histogram
  TNtuple* fNtuple_XArapucas;
  TNtuple* fNtuple_PMTs;
  TNtuple* fNtuple_IDE;
  TNtuple* fNtuple_OpHits_XArapuca;
  TNtuple* fNtuple_OpHits_PMT;
  TTree*   fTree_OpDetWaveforms;
  opdet::sbndPDMapAlg pdsMap;  //map for photon detector types
};


MyAnalysis::MyAnalysis(fhicl::ParameterSet const& p) : EDAnalyzer{p}
{
  float maxEnergy         = p.get<float>("MaxNuEnergy", 3.0);
  bool dump_OpHits        = p.get<bool >("Dump_OpHits", false);
  bool dump_G4_PE         = p.get<bool >("Dump_G4_PE", false);
  bool dump_DE            = p.get<bool >("Dump_DE", false);
  bool dump_Op_Waveforms  = p.get<bool >("Dump_Op_Waveforms", false);

  fdump_OpHits            = dump_OpHits;
  fdump_G4_PE             = dump_G4_PE;         
  fdump_DE                = dump_DE;            
  fdump_Op_Waveforms      = dump_Op_Waveforms;  

  art::ServiceHandle<art::TFileService> tfs;
  fHist                                          = tfs->make<TH1F>   ("enu", ";E_{#nu} (GeV);Events", 100, 0, maxEnergy);
  if (dump_DE)           fNtuple_IDE             = tfs->make<TNtuple>("IDEs","IDEs","ev:totalDE:totalPEarapucas:totalPEarapucasVUV:totalPEarapucasVis");
  if (dump_G4_PE)        fNtuple_XArapucas       = tfs->make<TNtuple>("xarapuca","xarapuca","ev:ch:meanphotons:t:arapucatype");
  if (dump_G4_PE)        fNtuple_PMTs            = tfs->make<TNtuple>("pmt","pmt","ev:ch:meanphotons:t");
  if (dump_OpHits)       fNtuple_OpHits_XArapuca = tfs->make<TNtuple>("ophits_xarapuca","ophits_xarapuca","ev:ch:peak_time_abs:peak_time:width:area:amplitude:PE");
  if (dump_OpHits)       fNtuple_OpHits_PMT      = tfs->make<TNtuple>("ophits_pmt","ophits_pmt"          ,"ev:ch:peak_time_abs:peak_time:width:area:amplitude:PE");
  if (dump_Op_Waveforms) fTree_OpDetWaveforms    = tfs->make<TTree>("OpDetWaveforms","OpDetWaveforms");
}

void MyAnalysis::analyze(art::Event const& e)
{
  //--------------------------------Variables--------------------------------

  float totalEdep = 0;
  float totalPEarapucas=0;
  float totalPEarapucasVUV=0;
  float totalPEarapucasVis=0;
  //~float totalne = 0;
  //~float x,y,z;
  int fEvNumber = e.id().event();
  std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> fPhotonLiteHandles;
  art::Handle< std::vector<sim::SimChannel> > simchannels;
  art::Handle<std::vector<simb::MCTruth> > mctruths;
  art::Handle<std::vector<recob::OpHit>> ophits_ara_h;
  art::Handle<std::vector<recob::OpHit>> ophits_pmt_h;
  art::Handle< std::vector< raw::OpDetWaveform > > waveHandle;


  //--------------------------------SimPhotons (G4)--------------------------------
  
  //Get *ALL* SimPhotonsCollectionLite from Event
  fPhotonLiteHandles.clear();
  fPhotonLiteHandles = e.getMany<std::vector<sim::SimPhotonsLite>>();
  const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = fPhotonLiteHandles;
  if (fdump_G4_PE){
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
      // this now tells you if light collection is reflected
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

      for (auto const& litesimphotons : (*opdetHandle)) {
        const unsigned ch = litesimphotons.OpChannel;
        const std::string pdtype = pdsMap.pdType(ch);
        std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;

        for (auto const& photonMember : photonMap) {
          auto meanPhotons = photonMember.second;
          auto tphoton = photonMember.first;

          if((pdtype == "xarapuca_vuv" && !Reflected) || (pdtype == "xarapuca_vis" && Reflected) ){
            //~if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
            int arapucatype=2;
            if(pdtype == "xarapuca_vuv"){
              totalPEarapucasVUV+=meanPhotons;
              arapucatype=1;
            }else{
              totalPEarapucasVis+=meanPhotons;
              arapucatype=0;
            } //1 VUV, 0 Visible
            fNtuple_XArapucas->Fill(fEvNumber,ch,meanPhotons,tphoton,arapucatype);
            totalPEarapucas+=meanPhotons;
          }//xarapucas

          if((pdtype == "pmt_coated") || (pdtype == "pmt_uncoated") ){
            //~if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
            fNtuple_PMTs->Fill(fEvNumber,ch,meanPhotons,tphoton);
          }
        }
      }
    }
  }
  //--------------------------------Deposited Energy--------------------------------
  if(fdump_DE){
    e.getByLabel("simdrift", simchannels); 
    for (auto const& channel : *simchannels) {
      //~const unsigned int ch = channel.Channel();
      //~std::cout<<ch<<std::endl;//verify handle works
          for(auto const &tdcide : channel.TDCIDEMap() ) {
            for(const auto& ide : tdcide.second) {//tdcide=pair   <  tick(time)    ,    vector of sim::IDEs    >
              totalEdep += ide.energy;
              //~totalne += ide.numElectrons;
              //~x = ide.x;
              //~y = ide.y;
              //~z = ide.z;

            }
      }
    }
    float_t kNplanes=3;
    totalEdep/=kNplanes  ;
    fNtuple_IDE->Fill(fEvNumber,totalEdep,totalPEarapucas,totalPEarapucasVUV,totalPEarapucasVis);
  }
  
  //--------------------------------OpHits---------------------------------
  if(fdump_OpHits){

    //-----XArapucas
    e.getByLabel("ophitarapuca", ophits_ara_h);

    for(auto const& oph : *ophits_ara_h) {
      auto ch        = oph.OpChannel();
      auto peak_abs  = oph.PeakTimeAbs();
      auto peak      = oph.PeakTime();
      auto width     = oph.Width();
      auto area      = oph.Area();
      auto amplitude = oph.Amplitude();
      auto pe        = oph.PE();

      fNtuple_OpHits_XArapuca->Fill(fEvNumber,
                               ch,
                               peak_abs,
                               peak,
                               width,
                               area,
                               amplitude,
                               pe);
    }

    //-----PMTs
    e.getByLabel("ophitpmt", ophits_pmt_h);

    for(auto const& oph : *ophits_pmt_h) {
      auto ch        = oph.OpChannel();
      auto peak_abs  = oph.PeakTimeAbs();
      auto peak      = oph.PeakTime();
      auto width     = oph.Width();
      auto area      = oph.Area();
      auto amplitude = oph.Amplitude();
      auto pe        = oph.PE();

      fNtuple_OpHits_PMT->Fill(fEvNumber,
                               ch,
                               peak_abs,
                               peak,
                               width,
                               area,
                               amplitude,
                               pe);
    }
  }


  //--------------------------------OpDetWaveforms-----------------------------
  if(fdump_Op_Waveforms){

    e.getByLabel("opdaq", waveHandle);
    unsigned int fChNumber;
    double fStartTime;
    std::vector<double> fwave={};
    // std::string fpdtype
    fTree_OpDetWaveforms->Branch("ev",&fEvNumber);
    fTree_OpDetWaveforms->Branch("ch",&fChNumber);
    fTree_OpDetWaveforms->Branch("timestamp",&fStartTime);
    fTree_OpDetWaveforms->Branch("waveform",&fwave);
    // fTree_OpDetWaveforms->Branch("pd_type",&fpdtype);

    for(auto const& wvf : (*waveHandle)) {
      fChNumber   = wvf.ChannelNumber();
      fStartTime  = wvf.TimeStamp(); //in us
      // fpdtype     = pdsMap.pdType(fChNumber);
      fwave       = {};
      for(unsigned int i = 0; i < wvf.size(); i++) {
        fwave.push_back((double)wvf[i]);
      }
      fTree_OpDetWaveforms->Fill();
    }
  }

  // //--------------------------------OpDetWaveforms----------------------------- as vectors
  // e.getByLabel("opdaq", waveHandle);
  
  // std::vector< unsigned int >        fChNumber;
  // std::vector< double >              fStartTime ;
  // std::vector< std::vector<double> > fwave;
  // std::vector< std::string >               fpdtype;

  // fTree_OpDetWaveforms->Branch("ch",&fChNumber);
  // fTree_OpDetWaveforms->Branch("timestamp",&fStartTime);
  // fTree_OpDetWaveforms->Branch("waveform",&fwave);
  // fTree_OpDetWaveforms->Branch("pd_type","",&fpdtype);

  // for(auto const& wvf : (*waveHandle)) {
  //   fChNumber.push_back( wvf.ChannelNumber());
  //   fStartTime.push_back(wvf.TimeStamp()); //in us
  //   fpdtype.push_back(   pdsMap.pdType( wvf.ChannelNumber() ) );
  //   std::vector<double> wave       = {};
  //   for(unsigned int i = 0; i < wvf.size(); i++) {
  //     wave.push_back((double)wvf[i]);
  //   }
  //   fwave.push_back(wave);
  // }
  // fTree_OpDetWaveforms->Fill();


  //--------------------------------MCTruth--------------------------------

  e.getByLabel("generator", mctruths); 

  //Andy Stuff
    //~for (auto const& truth : *mctruths) {
      //~const simb::MCNeutrino& mcnu = truth.GetNeutrino();
      //~const simb::MCParticle& nu = mcnu.Nu();
      //~float enu = nu.E();
      //~fHist->Fill(enu);
    //~}

}

DEFINE_ART_MODULE(MyAnalysis)

