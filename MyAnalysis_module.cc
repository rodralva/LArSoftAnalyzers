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

//Sim photons?
#include <cmath>
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"


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
  TH1F* fHist;  //!< Output histogram
  TNtuple* fNtuple_XArapucas;
  TNtuple* fNtuple_PMTs;
  TNtuple* fNtuple_IDE;
  opdet::sbndPDMapAlg pdsMap;  //map for photon detector types
};


MyAnalysis::MyAnalysis(fhicl::ParameterSet const& p) : EDAnalyzer{p}
{
  float maxEnergy = p.get<float>("MaxNuEnergy", 3.0);
  art::ServiceHandle<art::TFileService> tfs;
  fHist = tfs->make<TH1F>("enu", ";E_{#nu} (GeV);Events", 100, 0, maxEnergy);
  fNtuple_IDE =tfs->make<TNtuple>("IDEs","IDEs","ev:totalDE:totalPEarapucas:totalPEarapucasVUV:totalPEarapucasVis");
  fNtuple_XArapucas =tfs->make<TNtuple>("xarapuca","xarapuca","ev:ch:meanphotons:t:arapucatype");
  fNtuple_PMTs =tfs->make<TNtuple>("pmt","pmt","ev:ch:meanphotons:t");
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


  //--------------------------------SimPhotons--------------------------------

    e.getManyByType(fPhotonLiteHandles);

 const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = fPhotonLiteHandles;
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
		//~if (meanPhotons>2) std::cout<<fEvNumber<<"	"<<ch<<"	"<<meanPhotons<<"	"<<tphoton<<std::endl;
		int arapucatype=2;
		
		if(pdtype == "xarapuca_vuv"){
			totalPEarapucasVUV+=meanPhotons;
			arapucatype=1;}else{
			totalPEarapucasVis+=meanPhotons;
			arapucatype=0;} //1 VUV, 0 Visible
		fNtuple_XArapucas->Fill(fEvNumber,ch,meanPhotons,tphoton,arapucatype);
		totalPEarapucas+=meanPhotons;
		
		}//xarapucas
		
	  if((pdtype == "pmt_coated") || (pdtype == "pmt_uncoated") ){
		//~if (meanPhotons>2) std::cout<<fEvNumber<<"	"<<ch<<"	"<<meanPhotons<<"	"<<tphoton<<std::endl;
		fNtuple_PMTs->Fill(fEvNumber,ch,meanPhotons,tphoton);}
	  
	  }
	}
	}
	
	
  //--------------------------------Deposited Energy--------------------------------

  e.getByLabel("simdrift", simchannels); 

  
  for (auto const& channel : *simchannels) {
	  //~const unsigned int ch = channel.Channel();
	  //~std::cout<<ch<<std::endl;//verify handle works
        for(auto const &tdcide : channel.TDCIDEMap() ) {
          for(const auto& ide : tdcide.second) {//tdcide=pair	 <	tick(time)		,		vector of sim::IDEs		>
            totalEdep += ide.energy;
            //~totalne += ide.numElectrons;
            //~x = ide.x;
            //~y = ide.y;
            //~z = ide.z;
            
          }
		}
	}
	float_t kNplanes=3;
	totalEdep/=kNplanes	;
	fNtuple_IDE->Fill(fEvNumber,totalEdep,totalPEarapucas,totalPEarapucasVUV,totalPEarapucasVis);
	  
  
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

