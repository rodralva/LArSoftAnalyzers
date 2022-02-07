/////////////////////////////////////////////////////////////////////////
// Class:       MyAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        MyAnalysis_module.cc
//
// Generated at Sun Oct  4 21:09:21 2020 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_10_00.
/////////////////////////////////////////////////////////////////////////

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
#include "TString.h"

// Sim photons?
#include <cmath>
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"

// OpHits?
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"

// OpDetWaveForms
#include "lardataobj/RawData/OpDetWaveform.h"

class MyAnalysis;

class MyAnalysis : public art::EDAnalyzer
{
public:
  explicit MyAnalysis(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyAnalysis(MyAnalysis const &) = delete;
  MyAnalysis(MyAnalysis &&) = delete;
  MyAnalysis &operator=(MyAnalysis const &) = delete;
  MyAnalysis &operator=(MyAnalysis &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  bool fdump_G4_PE;
  bool fdump_DE;

  TTree *fMyAnaTree;
  TNtuple *fNtuple_XArapucas;
  TNtuple *fNtuple_PMTs;
  TNtuple *fNtuple_IDE;
  TTree *fTree_MC_Truth;
  opdet::sbndPDMapAlg pdsMap; // map for photon detector types
  std::vector<std::vector<std::string>> fProductsToDump;
  TTree *tree_array[100];
  TNtuple *ntuple_array[100];
  bool debug = true;

  std::vector<double> my_vector[10][10][10];
  // MCTruth
  int f_pdg;
  float f_E;

  double fx_i;
  double fy_i;
  double fz_i;

  double fx_f;
  double fy_f;
  double fz_f;
  int fEvNumber;

  // Raw Digitized Waveforms
  unsigned int fChNumber;
  double fStartTime;
  std::vector<double> fwave = {};
  std::vector<std::vector<double>> fwaves[10];
};

MyAnalysis::MyAnalysis(fhicl::ParameterSet const &p) : EDAnalyzer{p}
{
  debug = true;

  bool dump_G4_PE = p.get<bool>("Dump_G4_PE", false);
  bool dump_DE = p.get<bool>("Dump_DE", false);

  fdump_G4_PE = dump_G4_PE;
  fdump_DE = dump_DE;

  art::ServiceHandle<art::TFileService> tfs;
  if (dump_DE)
    fNtuple_IDE = tfs->make<TNtuple>("IDEs", "IDEs", "ev:totalDE:totalPEarapucas:totalPEarapucasVUV:totalPEarapucasVis");
  if (dump_G4_PE)
    fNtuple_XArapucas = tfs->make<TNtuple>("xarapuca", "xarapuca", "ev:ch:meanphotons:t:arapucatype");
  if (dump_G4_PE)
    fNtuple_PMTs = tfs->make<TNtuple>("pmt", "pmt", "ev:ch:meanphotons:t");

  // Vector/Mechanised approach
  fMyAnaTree = tfs->make<TTree>("AnaTree", "AnaTree");

  std::vector<std::vector<std::string>> ProductsToDump = p.get<std::vector<std::vector<std::string>>>("ProductsToDump");
  fProductsToDump = ProductsToDump;
  int counter = 0;
  for (auto i : fProductsToDump)
  {
    auto i_product = i[0];
    auto i_label = i[1];
    auto o_label = i[2];
    TString out_label = o_label;

    if (debug)
      std::cout << i_product << "\t" << i_label << "\t" << o_label << std::endl;

    // Set branches: Ifs over the different types of data
    if (i_product == "raw::OpDetWaveform")
    {
      //--------------------------------Raw Waveforms-----------------------------
      tree_array[counter] = tfs->make<TTree>(o_label.c_str(), o_label.c_str());
      tree_array[counter]->Branch("ev", &fEvNumber);
      tree_array[counter]->Branch("ch", &fChNumber);
      tree_array[counter]->Branch("timestamp", &fStartTime);
      tree_array[counter]->Branch("waveform", &fwave);

      fMyAnaTree->Branch(out_label + "_ch",        &my_vector[1][1][counter]);
      fMyAnaTree->Branch(out_label + "_timestamp", &my_vector[1][2][counter]);
      fMyAnaTree->Branch(out_label + "_waveforms", &fwaves[counter]);
    }
    else if (i_product == "recob::OpHit")
    {
      //--------------------------------Optical Hits-----------------------------
      ntuple_array[counter] = tfs->make<TNtuple>(o_label.c_str(), o_label.c_str(), "ev:ch:peak_time_abs:peak_time:width:area:amplitude:PE");
      fMyAnaTree->Branch(out_label + "_ch",            &my_vector[2][1][counter]);
      fMyAnaTree->Branch(out_label + "_peak_time_abs", &my_vector[2][2][counter]);
      fMyAnaTree->Branch(out_label + "_peak_time",     &my_vector[2][3][counter]);
      fMyAnaTree->Branch(out_label + "_width",         &my_vector[2][4][counter]);
      fMyAnaTree->Branch(out_label + "_area",          &my_vector[2][5][counter]);
      fMyAnaTree->Branch(out_label + "_amplitude",     &my_vector[2][6][counter]);
      fMyAnaTree->Branch(out_label + "_PE",            &my_vector[2][7][counter]);
    }
    else if (i_product == "sim::SimPhotonsLite")
    {
      //--------------------------------G4 Photons-----------------------------
      ntuple_array[counter] = tfs->make<TNtuple>(o_label.c_str(), o_label.c_str(), "ev:ch:meanphotons:t");
      fMyAnaTree->Branch(out_label + "_ch", &my_vector[3][1][counter]);
      fMyAnaTree->Branch(out_label + "_Photons", &my_vector[3][2][counter]);
      fMyAnaTree->Branch(out_label + "_t", &my_vector[3][3][counter]);
    }
    else if (i_product == "simb::MCTruth")
    {
      //--------------------------------Gen MCTruth-----------------------------
      tree_array[counter] = tfs->make<TTree>(o_label.c_str(), o_label.c_str());
      tree_array[counter]->Branch("ev", &fEvNumber);
      tree_array[counter]->Branch("pdg", &f_pdg);
      tree_array[counter]->Branch("E", &f_E);
      tree_array[counter]->Branch("x_i", &fx_i);
      tree_array[counter]->Branch("y_i", &fy_i);
      tree_array[counter]->Branch("z_i", &fz_i);
      tree_array[counter]->Branch("x_f", &fx_f);
      tree_array[counter]->Branch("y_f", &fy_f);
      tree_array[counter]->Branch("z_f", &fz_f);

      fMyAnaTree->Branch(out_label + "pdg", &my_vector[4][1][counter]);
      fMyAnaTree->Branch(out_label + "E",   &my_vector[4][2][counter]);
      fMyAnaTree->Branch(out_label + "x_i", &my_vector[4][3][counter]);
      fMyAnaTree->Branch(out_label + "y_i", &my_vector[4][4][counter]);
      fMyAnaTree->Branch(out_label + "z_i", &my_vector[4][5][counter]);
      fMyAnaTree->Branch(out_label + "x_f", &my_vector[4][6][counter]);
      fMyAnaTree->Branch(out_label + "y_f", &my_vector[4][7][counter]);
      fMyAnaTree->Branch(out_label + "z_f", &my_vector[4][8][counter]);
    }

    ++counter;
    // delete aux;
  };
}

void MyAnalysis::analyze(art::Event const &e)
{
  //--------------------------------Variables--------------------------------

  float totalEdep = 0;
  float totalPEarapucas = 0;
  float totalPEarapucasVUV = 0;
  float totalPEarapucasVis = 0;
  // float totalne = 0;
  // float x,y,z;
  fEvNumber = e.id().event();
  // std::cout<<"rodrigoa debug: "<<fEvNumber<<std::endl;

  std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> fPhotonLiteHandles;
  art::Handle<std::vector<sim::SimChannel>> simchannels;
  art::Handle<std::vector<simb::MCTruth>> mctruths;
  art::Handle<std::vector<recob::OpHit>> ophits_h;
  art::Handle<std::vector<raw::OpDetWaveform>> waveHandle;

  int counter = 0;
  for (auto i : fProductsToDump)
  {
    auto i_product = i[0];
    auto i_label = i[1];
    auto o_label = i[2];

    if (debug)
      std::cout << i_product << "\t" << i_label << "\t" << o_label << std::endl;

    // Fill branches: Ifs over the different types of data
    if (i_product == "raw::OpDetWaveform")
    {
      //--------------------------------Raw Waveforms-----------------------------
      e.getByLabel(i_label.c_str(), waveHandle);

      for (auto const &wvf : (*waveHandle))
      {
        fChNumber = wvf.ChannelNumber();
        fStartTime = wvf.TimeStamp(); // in us
        my_vector[1][1][counter].push_back(fChNumber);
        my_vector[1][2][counter].push_back(fStartTime);
        // fpdtype     = pdsMap.pdType(fChNumber);
        fwave = {};
        for (unsigned int i = 0; i < wvf.size(); i++)
        {
          fwave.push_back((double)wvf[i]);
        }
        // if (debug) std::cout<<"inside raw::OpDetWaveform loop"<<std::endl;
        tree_array[counter]->Fill();
        fwaves[counter].push_back(fwave);
      }
    }
    else if (i_product == "recob::OpHit")
    {
      //--------------------------------Optical Hits-----------------------------
      e.getByLabel(i_label.c_str(), ophits_h);

      for (auto const &oph : *ophits_h)
      {
        auto ch = oph.OpChannel();
        auto peak_abs = oph.PeakTimeAbs();
        auto peak = oph.PeakTime();
        auto width = oph.Width();
        auto area = oph.Area();
        auto amplitude = oph.Amplitude();
        auto pe = oph.PE();
        my_vector[2][1][counter].push_back(ch);
        my_vector[2][2][counter].push_back(peak_abs);
        my_vector[2][3][counter].push_back(peak);
        my_vector[2][4][counter].push_back(width);
        my_vector[2][5][counter].push_back(area);
        my_vector[2][6][counter].push_back(amplitude);
        my_vector[2][7][counter].push_back(pe);

        ntuple_array[counter]->Fill(fEvNumber,
                                    ch,
                                    peak_abs,
                                    peak,
                                    width,
                                    area,
                                    amplitude,
                                    pe);
      }
    }
    else if (i_product == "sim::SimPhotonsLite")
    {
      //--------------------------------SimPhotonsLite (G4)--------------------------------
      fPhotonLiteHandles.clear();
      fPhotonLiteHandles = e.getMany<std::vector<sim::SimPhotonsLite>>();
      const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = fPhotonLiteHandles;
      for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles)
      {
        // this now tells you if light collection is reflected
        const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

        for (auto const &litesimphotons : (*opdetHandle))
        {
          const unsigned ch = litesimphotons.OpChannel;
          const std::string pdtype = pdsMap.pdType(ch);
          std::map<int, int> const &photonMap = litesimphotons.DetectedPhotons;

          for (auto const &photonMember : photonMap)
          {
            auto meanPhotons = photonMember.second;
            auto tphoton = photonMember.first;

            if ((pdtype == "xarapuca_vuv" && !Reflected) || (pdtype == "xarapuca_vis" && Reflected))
            {
              // if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
              ntuple_array[counter]->Fill(fEvNumber, ch, meanPhotons, tphoton);
              totalPEarapucas += meanPhotons;

              my_vector[3][1][counter].push_back(ch);
              my_vector[3][2][counter].push_back(meanPhotons);
              my_vector[3][3][counter].push_back(tphoton);
            } // xarapucas

            if ((pdtype == "pmt_coated") || (pdtype == "pmt_uncoated"))
            {
              // if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
              ntuple_array[counter]->Fill(fEvNumber, ch, meanPhotons, tphoton);
            }
          }
        }
      }
    }
    else if (i_product == "simb::MCTruth")
    {
      //--------------------------------MCTruth--------------------------------
      e.getByLabel(i_label.c_str(), mctruths);
      for (auto const &truth : *mctruths)
      {
        int N = truth.NParticles();
        // std::cout<<"rodrigoa debug: "<<N<<std::endl;
        for (int i = 0; i < N; ++i)
        {
          const simb::MCParticle &nu = truth.GetParticle(i);
          float E = nu.E();
          const int pdg = nu.PdgCode();

          const TLorentzVector &v4_i = nu.Position();
          const TLorentzVector &v4_f = nu.EndPosition();

          auto x_i = v4_i.X();
          auto y_i = v4_i.Y();
          auto z_i = v4_i.Z();

          auto x_f = v4_f.X();
          auto y_f = v4_f.Y();
          auto z_f = v4_f.Z();

          f_E = E;
          f_pdg = pdg;
          fx_i = x_i;
          fy_i = y_i;
          fz_i = z_i;
          fx_f = x_f;
          fy_f = y_f;
          fz_f = z_f;

          my_vector[4][1][counter].push_back(pdg);
          my_vector[4][2][counter].push_back(E);
          my_vector[4][3][counter].push_back(x_i);
          my_vector[4][4][counter].push_back(y_i);
          my_vector[4][5][counter].push_back(z_i);
          my_vector[4][6][counter].push_back(x_f);
          my_vector[4][7][counter].push_back(y_f);
          my_vector[4][8][counter].push_back(z_f);
          tree_array[counter]->Fill();
        }
      }
    }

    ++counter;
  } // fProductsToDump loop
  //--------------------------------SimPhotons (G4)--------------------------------

  // Get *ALL* SimPhotonsCollectionLite from Event
  if (fdump_G4_PE)
  {
    fPhotonLiteHandles.clear();
    fPhotonLiteHandles = e.getMany<std::vector<sim::SimPhotonsLite>>();
    const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = fPhotonLiteHandles;
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles)
    {
      // this now tells you if light collection is reflected
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

      for (auto const &litesimphotons : (*opdetHandle))
      {
        const unsigned ch = litesimphotons.OpChannel;
        const std::string pdtype = pdsMap.pdType(ch);
        std::map<int, int> const &photonMap = litesimphotons.DetectedPhotons;

        for (auto const &photonMember : photonMap)
        {
          auto meanPhotons = photonMember.second;
          auto tphoton = photonMember.first;

          if ((pdtype == "xarapuca_vuv" && !Reflected) || (pdtype == "xarapuca_vis" && Reflected))
          {
            // if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
            int arapucatype = 2;
            if (pdtype == "xarapuca_vuv")
            {
              totalPEarapucasVUV += meanPhotons;
              arapucatype = 1;
            }
            else
            {
              totalPEarapucasVis += meanPhotons;
              arapucatype = 0;
            } // 1 VUV, 0 Visible
            fNtuple_XArapucas->Fill(fEvNumber, ch, meanPhotons, tphoton, arapucatype);
            totalPEarapucas += meanPhotons;
          } // xarapucas

          if ((pdtype == "pmt_coated") || (pdtype == "pmt_uncoated"))
          {
            // if (meanPhotons>2) std::cout<<fEvNumber<<"  "<<ch<<"  "<<meanPhotons<<"  "<<tphoton<<std::endl;
            fNtuple_PMTs->Fill(fEvNumber, ch, meanPhotons, tphoton);
          }
        }
      }
    }
  }
  //--------------------------------Deposited Energy--------------------------------
  if (fdump_DE)
  {
    e.getByLabel("simdrift", simchannels);
    for (auto const &channel : *simchannels)
    {
      // const unsigned int ch = channel.Channel();
      // std::cout<<ch<<std::endl;//verify handle works
      for (auto const &tdcide : channel.TDCIDEMap())
      {
        for (const auto &ide : tdcide.second)
        { // tdcide=pair   <  tick(time)    ,    vector of sim::IDEs    >
          totalEdep += ide.energy;
          // totalne += ide.numElectrons;
          // x = ide.x;
          // y = ide.y;
          // z = ide.z;
        }
      }
    }
    float_t kNplanes = 3;
    totalEdep /= kNplanes;
    fNtuple_IDE->Fill(fEvNumber, totalEdep, totalPEarapucas, totalPEarapucasVUV, totalPEarapucasVis);
  }

  fMyAnaTree->Fill();
  for (auto &v1 : my_vector)
  {
    for (auto &v2 : v1)
    {
      for (auto &v3 : v2)
      {
        v3.clear();
      }
    }
  }
}

DEFINE_ART_MODULE(MyAnalysis)