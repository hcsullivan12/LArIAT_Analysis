////////////////////////////////////////////////////////////////////////
// Class:       CalculatePionXS
// Module Type: analyzer
// File:        CalculatePionXS_module.cc
//
// Generated at Thu May 30 13:07:47 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h" 

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TVector3.h"
#include "TH1D.h"
#include "TTree.h"

// LArIAT includes
#include "LArIATFilterModule/InelasticSubClassifier.h"

namespace lariat {

class CalculatePionXS;

class CalculatePionXS : public art::EDAnalyzer {
public:
  explicit CalculatePionXS(fhicl::ParameterSet const & p);

  CalculatePionXS(CalculatePionXS const &) = delete;
  CalculatePionXS(CalculatePionXS &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  
private:

  void ResetVars();
  bool Filter(size_t& xsRecoTrkId, const std::vector<art::Ptr<recob::Track> >& tracklist);
  void FillTruthInfo(const sim::ParticleList& plist);
  bool InActiveRegion(const TVector3& thePos);
  void ComputeAngles(float& phi, float& theta, const TVector3& p0Hat);
  void TrackMatchedTrack(const size_t& xsRecoTrkId, 
                         const std::vector<art::Ptr<recob::Track> >& tracklist,
                         art::FindManyP<anab::Calorimetry>& fmcal);
  void DetermineInelasticity(const size_t& xsRecoTrkId,
                             const TVector3& furthestInZCaloPoint,
                             const std::vector<art::Ptr<recob::Track> >& tracklist);

  // labels
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;

  // =========================================================================
  // ==============================  VARIABLES  ==============================
  // =========================================================================  
  int    fTotalEvents = 0;
  int    fPrimariesEntered = 0;
  int    fGoodMCEvents = 0;
  int    fEventsFrontTpcTrk = 0;
  int    fEventsUpperTpcTrkCount = 0;
  int    fEventsWcTpcUniqueMatchAlpha = 0;
  int    fEventsDeltaMatch = 0;
  int    fEventsWcTpcUniqueMatch = 0;

  int    fRun;
  int    fSubrun;
  int    fEvent;
  int    fIsTrackInteracting;
  int    fIsTrackSignal;
  int    fDidDetermineInteracting;
  int    fDidDetermineSignal;
  double fTrueKEFF;
  double fTrueInteractingKE;
  double fRecoInteractingKE;

  TTree* fTree;

  // ====================================================================
  // ==============================  CUTS  ==============================
  // ====================================================================

  // maximum distance to allow track to be seperated from vertex
  float fVertexCut;

  // minimum length of secondary tracks
  float fSecondaryLengthCut;

  // minimum angle between primary and secondary (degrees)
  float fSecondaryAngleCut;

  // =====================================================================================
  // ==============================  VARIABLES FOR G4 INFO  ==============================
  // =====================================================================================
  std::vector<int>                      fG4PrimaryTrkId;     // track id 
  std::vector<std::string>              fG4PrimaryProcesses; // processes
  std::vector<std::string>              fG4PrimarySubProcesses; // subprocesses for inelastic
  std::vector<TVector3>                 fG4PrimaryPos0;      // start pos
  std::vector<TVector3>                 fG4PrimaryPosf;      // final pos 
  std::vector<TVector3>                 fG4PrimaryMom0;      // momentum
  std::vector<TVector3>                 fG4PrimaryProjPos0;  // projected position
  std::vector<std::vector<TVector3>>    fG4PrimaryTrTrjPos;  // trajectory sp
  std::vector<std::vector<TVector3>>    fG4PrimaryTrTrjMom;  // trajectory momentum 
  std::vector<int>                      fG4PrimaryVtx;       // primary interaction points

  
  
  // ==========================================================================================
  // ==================================== CONSTANTS  ==========================================
  // ==========================================================================================
  // tpc boundaries
  const float TPC_X_BOUND[2] = {   0.0, 47.0 };
  const float TPC_Y_BOUND[2] = { -20.0, 20.0 };
  const float TPC_Z_BOUND[2] = {   0.0, 90.0 };
  
  // fiducial volume definition
  const float FV_X_BOUND[2] = {   2.0, 45.0 };
  const float FV_Y_BOUND[2] = { -15.0, 15.0 };
  const float FV_Z_BOUND[2] = {   0.0, 86.0 };
   
  // mass of pion in MeV
  const float PARTICLE_MASS = 139.57; 
  
  // number of centimeters in Z we require a track to have a spacepoint
  const float FIRST_SP_Z_POS = 5.0;
   
  // portion of upstream TPC which we will restrict the number of tracks
  const float UPPER_PART_OF_TPC = 14.0;
   
  // number of upper tpc tracks allowed
  const size_t N_UPPER_TPC_TRACKS = 4;
   
  // Delta X Between Wire Chamber Track and TPC Track 
  const float DELTA_X_BOUND[2] = { -2.0, 6.0 };
   
  // Delta Y Between Wire Chamber Track and TPC Track 
  const float DELTA_Y_BOUND[2] = { -3.0, 6.0 };
   
  // the assumed energy loss between the cryostat and the TPC 
  const float ENTRY_TPC_ENERGY_LOSS = 36; //MeV
  
  // setting the global event weight based on
  //   open box WCTrack momentum spectrum     
  // 100 MeV < P < 200 MeV = 0.02
  // 200 MeV < P < 300 MeV = 0.10
  // 300 MeV < P < 400 MeV = 0.535
  // 400 MeV < P < 500 MeV = 0.84
  // 500 MeV < P < 600 MeV = 0.965
  // 600 MeV < P < 700 MeV = 1.0
  // 700 MeV < P < 800 MeV = 0.62
  // 800 MeV < P < 900 MeV = 0.225
  // 900 MeV < P < 1000MeV = 0.094
  // 1000MeV < P < 1100MeV = 0.0275
  // 1100MeV < P < 1500MeV = 0.01
  float EVENT_WEIGHT = 1.0;
  
  // True  = Use the momentum based weighting  
  // False = Don't weight events
  const bool USE_EVENT_WEIGHT = false;
   
  // constants for cross section calculation
  const double RHO            = 1395;          //kg/m^3
  const double MOLAR_MASS     = 39.948;        //g/mol
  const double G_PER_KG       = 1000;
  const double AVOGADRO       = 6.022140857e+23;        //number/mol
  const double NUMBER_DENSITY =  (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO;
  const double SLAB_WIDTH     = 0.0045;        //in m
  const double PI             = 3.141592654;
  const double BARN_PER_M2    = 1e28;
   
  // cut on angle between wc and tpc track (degrees)
  const float ALPHA_CUT = 10.; 

  // threshold for hadron dEdx
  double HIT_DEDX_THRESHOLD = 40.;



  // ============================================================================================
  // ====================================  HISTOGRAMS  ==========================================
  // ============================================================================================

  TH1D* hMCPrimaryMissedTpcX;
  TH1D* hMCPrimaryMissedTpcY;
  TH1D* hMCPrimaryMissedTpcZ;
  TH1D* hMCELossUpstream;
  TH1D* hDataUpstreamZPos;
  TH1D* hDeltaX;
  TH1D* hDeltaY;
  TH1D* hDeltaZ;
  TH1D* hAlpha;
  TH1D* hMCPrimaryPx;
  TH1D* hMCPrimaryPy;
  TH1D* hMCPrimaryPz;
  TH1D* hMCPrimaryP ;
  TH1D* hMCPrimaryProjX0;
  TH1D* hMCPrimaryProjY0;
  TH1D* hMCPrimaryProjZ0;
  TH1D* hTrueLength;
  TH1D* hFurthestInZCaloX;
  TH1D* hFurthestInZCaloY;
  TH1D* hFurthestInZCaloZ;
  TH1D* hRecoMCIncidentKE;
  TH1D* hRecoMCInteractingKE;
  TH1D* hVertexDiff;
  TH1D* hSecondaryLength;
};




//----------------------------------------------------------------------------------------------------
// Contructor
CalculatePionXS::CalculatePionXS(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);
}




//----------------------------------------------------------------------------------------------------
// Reconfigure
void CalculatePionXS::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel        = p.get< std::string >("TrackModuleLabel"      , "pmtrack");
  fCalorimetryModuleLabel  = p.get< std::string >("CalorimetryModuleLabel", "calo"     );

  fVertexCut          = p.get< float >("VertexCut",          2.0);
  fSecondaryLengthCut = p.get< float >("SecondaryLengthCut", 3.0); 
  fSecondaryAngleCut  = p.get< float >("SecondaryAngleCut",  10.0); 
}





//----------------------------------------------------------------------------------------------------
// Analyze
void CalculatePionXS::analyze(art::Event const & e)
{
  // ### Reset vars
  ResetVars();

  fRun    = e.run();
  fSubrun = e.subRun();
  fEvent  = e.event();

  // ### Handles
  art::Handle< std::vector<recob::Track> > trackListHandle; // Container of reco tracks
  std::vector<art::Ptr<recob::Track> >     tracklist;       // Vector of wc tracks

  // ### Get the track information
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << std::endl;
  if(!e.getByLabel(fTrackModuleLabel,trackListHandle)) return;
  fTotalEvents++; 
  if (fTotalEvents%1000 == 0) std::cout << "EVENT = " << fTotalEvents << std::endl; 
  art::fill_ptr_vector(tracklist, trackListHandle);
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, e, fCalorimetryModuleLabel);

  // ### Backtracker to recover truth information
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // =====================================================
  // =================  FILL TRUTH INFO  =================
  // =====================================================
  FillTruthInfo(plist);

  // ============================================
  // =================  FILTER  =================
  // ============================================
  size_t xsRecoTrkId;
  if (!Filter(xsRecoTrkId, tracklist)) return;

  // ==========================================================================
  // =================  BEGIN TRACKING OUR MATCHED TPC TRACK  =================
  // ==========================================================================
  // ### We are know tracking a uniquely WC/TPC matched track
  // ### The ID is in xsRecoTrkId variable

  std::cout << "///////////////////////////////////////////\n";
  // ### Get the sub process 
  piinelastic::InelasticSubClassifier subclassifier;
  std::cout << fG4PrimaryTrkId.size() << std::endl;
  if (!fG4PrimaryTrkId.size()) return;

  auto subprocess = subclassifier.Classify(plist, fG4PrimaryTrkId[0]);
  fG4PrimarySubProcesses.push_back(subprocess);

  // ### Track the matched track
  TrackMatchedTrack(xsRecoTrkId, tracklist, fmcal);

  // ### Fill our tree and get out of here
  fTree->Fill();
}





//----------------------------------------------------------------------------------------------------
// Reset vars
void CalculatePionXS::ResetVars()
{
  fRun    = -999;
  fSubrun = -999;
  fEvent  = -999;
  fIsTrackInteracting = 0;
  fIsTrackSignal      = 0;
  fDidDetermineInteracting = 0;
  fDidDetermineSignal = 0;
  fTrueKEFF = -999;
  fTrueInteractingKE = -999;
  fRecoInteractingKE = -999;

  fG4PrimaryTrkId.clear();
  fG4PrimaryProcesses.clear();
  fG4PrimarySubProcesses.clear();
  fG4PrimaryPos0.clear();
  fG4PrimaryPosf.clear();
  fG4PrimaryMom0.clear();
  fG4PrimaryProjPos0.clear();
  fG4PrimaryTrTrjPos.clear();
  fG4PrimaryTrTrjMom.clear();
  fG4PrimaryVtx.clear();
}





//----------------------------------------------------------------------------------------------------
// Begin job
void CalculatePionXS::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  
  hMCPrimaryMissedTpcX = tfs->make<TH1D>("hMCPrimaryMissedTpcX", "X Position of Missed Primary", 1000, -50, 100);
  hMCPrimaryMissedTpcY = tfs->make<TH1D>("hMCPrimaryMissedTpcY", "Y Position of Missed Primary", 1000, -40, 40);
  hMCPrimaryMissedTpcZ = tfs->make<TH1D>("hMCPrimaryMissedTpcZ", "Z Position of Missed Primary", 1000, -50, 150);
  hMCELossUpstream     = tfs->make<TH1D>("hMCELossUpstream",     "Primary Energy Loss Upstream", 1000, 0, 1100);
  hDataUpstreamZPos    = tfs->make<TH1D>("hDataUpstreamZPos",    "Minimum Z Position Upstream", 1000, 0, 100);
  hDeltaX              = tfs->make<TH1D>("hDeltaX",              "Primary Delta X", 1000, -10, 10);
  hDeltaY              = tfs->make<TH1D>("hDeltaY",              "Primary Delta Y", 1000, -10, 10);
  hDeltaZ              = tfs->make<TH1D>("hDeltaZ",              "Primary Delta Z", 1000, -10, 10);
  hAlpha               = tfs->make<TH1D>("hAlpha",               "Alpha", 1000, 0, 90);
  hMCPrimaryPx         = tfs->make<TH1D>("hMCPrimaryPx",         "Primary P_{x}", 1000, 0, 1100);
  hMCPrimaryPy         = tfs->make<TH1D>("hMCPrimaryPy",         "Primary P_{y}", 1000, 0, 1100);
  hMCPrimaryPz         = tfs->make<TH1D>("hMCPrimaryPz",         "Primary P_{z}", 1000, 0, 1100);
  hMCPrimaryP          = tfs->make<TH1D>("hMCPrimaryP",          "Primary P", 1000, 0, 1100);
  hMCPrimaryProjX0     = tfs->make<TH1D>("hMCPrimaryProjX0",     "Projected X Position on Front Face", 1000, 0, 47);
  hMCPrimaryProjY0     = tfs->make<TH1D>("hMCPrimaryProjY0",     "Projected Y Position on Front Face", 1000, -20, 40);
  hMCPrimaryProjZ0     = tfs->make<TH1D>("hMCPrimaryProjZ0",     "Projected Z Position on Front Face", 1000, 0, 100);
  hTrueLength          = tfs->make<TH1D>("hTrueLength",          "Primary True Length", 1000, 0, 100);
  hFurthestInZCaloX   = tfs->make<TH1D>("hFurthestInZCaloX",   "Furthest Calorimetry X Position", 1000, 0, 47);
  hFurthestInZCaloY   = tfs->make<TH1D>("hFurthestInZCaloY",   "Furthest Calorimetry Y Position", 1000, -20, 20);
  hFurthestInZCaloZ   = tfs->make<TH1D>("hFurthestInZCaloZ",   "Furthest Calorimetry Z Position", 1000, 0, 100);
  hRecoMCIncidentKE    = tfs->make<TH1D>("hRecoMCIncidentKE",    "Incident",    23, -50, 1100);
  hRecoMCInteractingKE = tfs->make<TH1D>("hRecoMCInteractingKE", "Interacting", 23, -50, 1100);
  hVertexDiff          = tfs->make<TH1D>("hVertexDiff",          "Distance of Candidate Endpoints", 100, 0, 50); 
  hSecondaryLength     = tfs->make<TH1D>("hSecondaryLength",     "Length of Secondary Tracks", 100, 0, 100);

  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&fRun      ,"run/I");
  fTree->Branch("subrun"   ,&fSubrun   ,"subrun/I");
  fTree->Branch("event"    ,&fEvent    ,"event/I");
  fTree->Branch("isTrackInteracting"  ,&fIsTrackInteracting   ,"isTrackInteracting/I");
  fTree->Branch("isTrackSignal"       ,&fIsTrackSignal        ,"isTrackSignal/I");
  fTree->Branch("trueProcess"         ,&fG4PrimaryProcesses);
  fTree->Branch("trueSubProcess"      ,&fG4PrimarySubProcesses);
  fTree->Branch("didDetermineInteracting" , &fDidDetermineInteracting   ,"didDetermineInteracting/I");
  fTree->Branch("didDetermineSignal" , &fDidDetermineSignal   ,"didDetermineSignal/I");
  fTree->Branch("trueKEFF" ,&fTrueKEFF ,"trueKEFF/D");
  fTree->Branch("trueInteractingKE" ,&fTrueInteractingKE ,"trueInteractingKE/D");
  fTree->Branch("recoInteractingKE" ,&fRecoInteractingKE ,"recoInteractingKE/D");
  double mass = PARTICLE_MASS;
  fTree->Branch("mass" , &mass ,"mass/D");
}





//----------------------------------------------------------------------------------------------------
// Fill truth info
void CalculatePionXS::FillTruthInfo(const sim::ParticleList& plist)
{
  // ### Loop over g4 particles 
  // ### Counter for the primary
  size_t iPrim(0);
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    auto mcParticle = plist.Particle(iG4);

    // ### If this is not a primary, skip it
    if (mcParticle->Process().find("primary") == std::string::npos) continue;
     
    // ### Store the track id 
    fG4PrimaryTrkId.push_back(mcParticle->TrackId());

    // ### Store the processes for this primary, only those that occur in TPC
    simb::MCTrajectory truetraj = mcParticle->Trajectory();
    auto thisTrjProcessMap = truetraj.TrajectoryProcesses();
  
    // ### If there's nothing, check the end of the track
    if (thisTrjProcessMap.size() == 0)
    {
      // ### Get the vertex
      if (!InActiveRegion(mcParticle->EndPosition().Vect())) continue;
      fG4PrimaryVtx.push_back( mcParticle->NumberTrajectoryPoints()-1 );

      auto tempMom = mcParticle->Momentum(mcParticle->NumberTrajectoryPoints()-2).Vect();
      fTrueInteractingKE = std::sqrt( tempMom*tempMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
 
      // ### If it doesn't have daughters, it's through going
      if (mcParticle->NumberDaughters())
      {
        auto theDauId = mcParticle->Daughter(0);
        for (size_t iD = 0; iD < plist.size(); iD++) 
        {
          if (plist.Particle(iD)->TrackId() == theDauId) fG4PrimaryProcesses.push_back(plist.Particle(iD)->Process());
        }
      }//<-- End if has daughters
      else fG4PrimaryProcesses.push_back("throughgoing");
    }//<-- End if map is zero
    else 
    {
      for (const auto& couple : thisTrjProcessMap)
      {
        // ### Get the vertex, only if it's in the TPC
        int interestingPoint = (int) couple.first;
        if (!InActiveRegion( truetraj.Position(interestingPoint).Vect()) ) continue;

        fG4PrimaryVtx.push_back( interestingPoint );
        fG4PrimaryProcesses.push_back( truetraj.KeyToProcess(couple.second) );

        auto tempMom = mcParticle->Momentum(interestingPoint-1).Vect();
        fTrueInteractingKE = std::sqrt( tempMom*tempMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
      }
    }//<-- End if map is not zero

    // ### Store the positions and momentum 
    fG4PrimaryPos0.push_back( mcParticle->Position().Vect() );
    fG4PrimaryPosf.push_back( mcParticle->EndPosition().Vect() );
    fG4PrimaryMom0.push_back( mcParticle->Momentum().Vect()*1000 ); // convert to MeV

    // ### Setting a global event weight 
    // ### Setting Event weight 
    if(USE_EVENT_WEIGHT)
    {
      if(fG4PrimaryMom0[iPrim].Z() > 0   && fG4PrimaryMom0[iPrim].Z() < 100) EVENT_WEIGHT = 0.010;
      if(fG4PrimaryMom0[iPrim].Z() > 100 && fG4PrimaryMom0[iPrim].Z() < 200) EVENT_WEIGHT = 0.020;
      if(fG4PrimaryMom0[iPrim].Z() > 200 && fG4PrimaryMom0[iPrim].Z() < 300) EVENT_WEIGHT = 0.100;
      if(fG4PrimaryMom0[iPrim].Z() > 300 && fG4PrimaryMom0[iPrim].Z() < 400) EVENT_WEIGHT = 0.535;
      if(fG4PrimaryMom0[iPrim].Z() > 400 && fG4PrimaryMom0[iPrim].Z() < 500) EVENT_WEIGHT = 0.840;
      if(fG4PrimaryMom0[iPrim].Z() > 500 && fG4PrimaryMom0[iPrim].Z() < 600) EVENT_WEIGHT = 0.965;
      if(fG4PrimaryMom0[iPrim].Z() > 600 && fG4PrimaryMom0[iPrim].Z() < 700) EVENT_WEIGHT = 1.000;
      if(fG4PrimaryMom0[iPrim].Z() > 700 && fG4PrimaryMom0[iPrim].Z() < 800) EVENT_WEIGHT = 0.620;
      if(fG4PrimaryMom0[iPrim].Z() > 800 && fG4PrimaryMom0[iPrim].Z() < 900) EVENT_WEIGHT = 0.225;
      if(fG4PrimaryMom0[iPrim].Z() > 900 && fG4PrimaryMom0[iPrim].Z() <1000) EVENT_WEIGHT = 0.094;
      if(fG4PrimaryMom0[iPrim].Z() >1000 && fG4PrimaryMom0[iPrim].Z() <1100) EVENT_WEIGHT = 0.0275;
      if(fG4PrimaryMom0[iPrim].Z() >1100)                             EVENT_WEIGHT = 0.010;
    }

    // ### Fill momentum histos
    hMCPrimaryPx->Fill(fG4PrimaryMom0[iPrim].X(),   EVENT_WEIGHT);
    hMCPrimaryPy->Fill(fG4PrimaryMom0[iPrim].Y(),   EVENT_WEIGHT);
    hMCPrimaryPz->Fill(fG4PrimaryMom0[iPrim].Z(),   EVENT_WEIGHT);
    hMCPrimaryP ->Fill(fG4PrimaryMom0[iPrim].Mag(), EVENT_WEIGHT);
    
    // ### Fill true length histo;
    hTrueLength->Fill( (fG4PrimaryPosf[iPrim] - fG4PrimaryPos0[iPrim]).Mag() );
    
    // ### Project onto tpc 
    TVector3 thisPosProjVec = fG4PrimaryPos0[iPrim] - ( fG4PrimaryPos0[iPrim].Z()/fG4PrimaryMom0[iPrim].Z() )*fG4PrimaryMom0[iPrim];
    fG4PrimaryProjPos0.push_back(thisPosProjVec);

    // ### Fill the proj histos
    hMCPrimaryProjX0->Fill( fG4PrimaryProjPos0[iPrim].X() );
    hMCPrimaryProjY0->Fill( fG4PrimaryProjPos0[iPrim].Y() );
    hMCPrimaryProjZ0->Fill( fG4PrimaryProjPos0[iPrim].Z() );

    // ### Store the trajectory points and momentum
    std::vector<TVector3> temp;
    fG4PrimaryTrTrjPos.push_back(temp);
    fG4PrimaryTrTrjMom.push_back(temp);
    std::map<float, float> keMap;
    for (size_t iPoint = 0; iPoint < mcParticle->NumberTrajectoryPoints(); iPoint++)
    {
      fG4PrimaryTrTrjPos[iPrim].push_back( mcParticle->Position(iPoint).Vect() );
      fG4PrimaryTrTrjMom[iPrim].push_back( mcParticle->Momentum(iPoint).Vect()*1000 ); // convert to MeV
      if (mcParticle->Position(iPoint).Vect().Z() < 0)
      {
        float ke = std::sqrt( fG4PrimaryTrTrjMom[iPrim].back()*fG4PrimaryTrTrjMom[iPrim].back() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; 
        keMap.emplace(mcParticle->Position(iPoint).Vect().Z(), ke);
      }
    }

    // ### Set the KE at z = 0
    if(!keMap.empty()) fTrueKEFF = (--keMap.end())->second;

    iPrim++;
  }//<--- End loop over G4 particles 

  // ### Set the interacting variables
  if (fG4PrimaryProcesses.size() > 0) {fIsTrackInteracting = 1;}
  for (const auto& p : fG4PrimaryProcesses)
  {
    if (p.find("pi") != std::string::npos && p.find("Inelastic") != std::string::npos) {fIsTrackSignal = 1;}
  }
}





//----------------------------------------------------------------------------------------------------
// Fill truth info
bool CalculatePionXS::Filter(size_t& xsRecoTrkId, const std::vector<art::Ptr<recob::Track> >& tracklist)
{
  // ==========================================================================
  // =================  LOOKING AT EVENTS THAT ENTER THE TPC  =================
  // ==========================================================================
  // ### Loop over primaries 
  bool isGoodEvent(false);
  for (size_t iPrim = 0; iPrim < fG4PrimaryTrkId.size(); iPrim++)
  {
    // only look if a primary entered the tpc
    if ( fG4PrimaryPosf[iPrim].Z() > 0 ) isGoodEvent = true;
    if ( !isGoodEvent ) 
    {
      hMCPrimaryMissedTpcX->Fill(fG4PrimaryPosf[iPrim].X());
      hMCPrimaryMissedTpcY->Fill(fG4PrimaryPosf[iPrim].Y());
      hMCPrimaryMissedTpcZ->Fill(fG4PrimaryPosf[iPrim].Z());
      continue;
    }
    // ####################################
    // ### This primary entered the TPC ###
    // ####################################
    fPrimariesEntered++;
    // calculate energy loss
    float energyLoss(0);
    // loop over trj points for this primary
    for (size_t iPoint = 0; iPoint < fG4PrimaryTrTrjPos[iPrim].size(); iPoint++)
    {
      // only look at the upstream portion
      if ( fG4PrimaryTrTrjPos[iPrim][iPoint].Z() > 0 ) continue;
      // ignore last point
      if ( (iPoint+1) >= fG4PrimaryTrTrjPos[iPrim].size() ) break;
      auto mom1Vec = fG4PrimaryTrTrjMom[iPrim][iPoint];
      auto mom2Vec = fG4PrimaryTrTrjMom[iPrim][iPoint+1];
      float energy1 = std::sqrt( mom1Vec.Mag()*mom1Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS );
      float energy2 = std::sqrt( mom2Vec.Mag()*mom2Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS ); 
      energyLoss += (energy1 - energy2);
    }//<--- End loop over true traj points
    hMCELossUpstream->Fill(energyLoss);
  }//<--- End loop over primaries
  if (!isGoodEvent) return false;
  std::cout << "[ x ] Good Event" << std::endl;
  fGoodMCEvents++;
  



  // ============================================================================
  // =================      APPLY FRONT FACE TPC TRACK CUT      =================
  // =================  APPLY CUT ON NUMBER OF UPSTREAM TRACKS  =================
  // =================          (Cuts on EM showers)            =================
  // ============================================================================
  std::vector<TVector3> frontFaceTrksPos;
  std::vector<TVector3> frontFaceTrksP0Hat;
  std::vector<double> frontFaceTrksPhi, frontFaceTrksTheta;
  std::vector<size_t> frontFaceTrksId;
  size_t nUpstreamTpcTrks(0);
  // ##########################################################
  // ### Only keeping events if there is a track within the ###
  // ###            the first x cm of the TPC               ###
  // ##########################################################
  for (size_t iTrk = 0; iTrk < tracklist.size(); iTrk++)
  {
    // looking for most upstream point still in FV of TPC
    // most upstream point in FV
    float tempMinZ(FV_Z_BOUND[1]);
    bool isFrontTpcTrk(false); 
    bool isUpstreamTpcTrk(false); 
    // dummy variables for observables attached to tempMinZ
    TVector3 dummyTrkPos(0,0,0);
    TVector3 dummyTrkP0Hat(0,0,0);
    size_t dummyTrkId;
 
    // ### Loop over trajectory points 
    auto thisTrk = tracklist[iTrk];
    for (size_t iTrjPoint = 0; iTrjPoint < thisTrk->NumberTrajectoryPoints(); iTrjPoint++)
    {
      // ### Fiducial check
      TVector3 theTrjPointVec( thisTrk->LocationAtPoint(iTrjPoint).X(), 
                               thisTrk->LocationAtPoint(iTrjPoint).Y(), 
                               thisTrk->LocationAtPoint(iTrjPoint).Z() );
      TVector3 theTrjP0HatVec( thisTrk->DirectionAtPoint(iTrjPoint).X(),
                               thisTrk->DirectionAtPoint(iTrjPoint).Y(),
                               thisTrk->DirectionAtPoint(iTrjPoint).Z() );
      if ( FV_X_BOUND[0] < theTrjPointVec.X() && theTrjPointVec.X() < FV_X_BOUND[1] &&
           FV_Y_BOUND[0] < theTrjPointVec.Y() && theTrjPointVec.Y() < FV_Y_BOUND[1] &&
           FV_Z_BOUND[0] < theTrjPointVec.Z() && theTrjPointVec.Z() < tempMinZ        )
      {
        // this point is in the fiducial volume and our new minimum
        tempMinZ = theTrjPointVec.Z();
        // we're storing all tracks that are in the front of the TPC
        if (tempMinZ < FIRST_SP_Z_POS)
        {
          isFrontTpcTrk = true;
          
          // Store these points for later
          dummyTrkPos   = theTrjPointVec;
          dummyTrkP0Hat = theTrjP0HatVec; 
          dummyTrkId    = iTrk;
        }//<--- End if within first few cm of TPC
        // check if this point is in the upstream portion
        if (tempMinZ < UPPER_PART_OF_TPC) isUpstreamTpcTrk = true;
      }//<--- End if in FV
    }//<--- End loop over trj points 
    
    // we should have the observables attached to the most upstream trj point for this track
    // if it is at the front, append to our vectors
    if (isFrontTpcTrk)
    {
      frontFaceTrksPos.  push_back( dummyTrkPos );
      frontFaceTrksP0Hat.push_back( dummyTrkP0Hat ); 
      frontFaceTrksId.   push_back( dummyTrkId );
    }
    // increment if in upstream portion
    if (isUpstreamTpcTrk) nUpstreamTpcTrks++;
    
    // fill histos
    hDataUpstreamZPos->Fill(tempMinZ);
  }//<-- End loop over reconstructed tracks
  
  // ### Skip events that do not have a trk at the front
  if (frontFaceTrksId.size() == 0) return false;
  std::cout << "[ x ] Front face tracks" << std::endl;
  fEventsFrontTpcTrk++;

  // ### Skip events that have too many tracks upstream        
  if(nUpstreamTpcTrks > N_UPPER_TPC_TRACKS || nUpstreamTpcTrks == 0) return false;
  std::cout << "[ x ] Upper TPC tracks" << std::endl;
  fEventsUpperTpcTrkCount++;





  // ========================================================================
  // =================  MATCHING MC TO RECONSTRUCTED TRACK  =================
  // ========================================================================
  size_t nMcTpcMatch(0);
  // ###########################################
  // ### Loop over all the front face Tracks ###
  // ###########################################
  for(size_t iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
  {
    auto deltaVec = frontFaceTrksPos[iFrFaTrk] - fG4PrimaryProjPos0[0];
    hDeltaX->Fill(deltaVec.X());
    hDeltaY->Fill(deltaVec.Y());
    hDeltaZ->Fill(deltaVec.Z());
    // matching in delta X and delta Y 
    if( DELTA_X_BOUND[0] < deltaVec.X() && deltaVec.X() < DELTA_X_BOUND[1] && 
        DELTA_Y_BOUND[0] < deltaVec.Y() && deltaVec.Y() < DELTA_Y_BOUND[1] )
    {
      nMcTpcMatch++;
    }
  }//<---End bb loop




  // ========================================================
  // =================  CALCULATING ANGLES  =================
  // ========================================================
  // #################################################
  // ### Calculating the angles for the G4 primary ###
  // #################################################
  TVector3 mcP0Hat = fG4PrimaryMom0[0]; 
  mcP0Hat = mcP0Hat.Unit();
  float    mcPhi(0), mcTheta(0);

  // compute the angles for primary
  ComputeAngles( mcPhi, mcTheta, mcP0Hat );

  // ##############################################################
  // ### Calculating the angles for the front face tracks (TPC) ###
  // ##############################################################
  for(size_t iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
  {
    // setting the TVector 
    auto tpcP0Hat = frontFaceTrksP0Hat[iFrFaTrk];
    tpcP0Hat = tpcP0Hat.Unit();
    float thisPhi(0), thisTheta(0);
    ComputeAngles( thisPhi, thisTheta, tpcP0Hat );
    frontFaceTrksPhi.push_back(thisPhi);
    frontFaceTrksTheta.push_back(thisTheta);
  }//<--- End iFrFaTrk loop
  // sanity check
  assert( frontFaceTrkPos.size()  == frontFaceTrkP0Hat.size()   == 
          frontFaceTrksPhi.size() == frontFaceTrksTheta.size()  == 
          frontFaceTrksId.size() );


  


  // ======================================================
  // =================  APPLY ANGLE CUTS  =================
  // ======================================================
  bool isAlphaMatch(false);
  for(size_t iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
  {
    // compute the g4 primary momentum unit vector
    TVector3 theMCUnit( std::sin(mcTheta)*std::cos(mcPhi),
                        std::sin(mcTheta)*std::sin(mcPhi),
                        std::cos(mcTheta) );
    // compute the front face tpc trk momontum unit vector
    TVector3 theTpcUnit( std::sin(frontFaceTrksTheta[iFrFaTrk])*std::cos(frontFaceTrksPhi[iFrFaTrk]),
                         std::sin(frontFaceTrksTheta[iFrFaTrk])*std::sin(frontFaceTrksPhi[iFrFaTrk]),
                         std::cos(frontFaceTrksTheta[iFrFaTrk]) );
    
    // ###########################################################
    // ### Calculating the angle between WCTrack and TPC Track ###
    // ###########################################################
    float alpha = std::acos( theMCUnit.Dot(theTpcUnit) ) * (180.0/PI); // convert to degrees
    // fill histo 
    hAlpha->Fill(alpha);
    // checking if this made our cut
    if (alpha < ALPHA_CUT)
    {
      isAlphaMatch = true;
      xsRecoTrkId = frontFaceTrksId[iFrFaTrk];
    }
  }
  // did we get a match?
  if (nMcTpcMatch > 0) fEventsDeltaMatch++;
  // force there to be one match
  if (nMcTpcMatch != 1) return false;
  std::cout << "[ x ] WC to TPC unique match" << std::endl;
  fEventsWcTpcUniqueMatch++;
  // force this unique match to pass our alpha cut
  if (!isAlphaMatch) return false;
  std::cout << "[ x ] Alpha cut" << std::endl;
  fEventsWcTpcUniqueMatchAlpha++;

  // ### We passed all cuts :)
  return true;
}





//----------------------------------------------------------------------------------------------------
// Compute angles
void CalculatePionXS::ComputeAngles(float& phi, float& theta, const TVector3& p0Hat)
{
  TVector3 zHat(0,0,1);

  // calculate theta 
  theta = std::acos( zHat.Dot(p0Hat) );
  // calculate phi
       if ( p0Hat.Y() > 0  && p0Hat.X() > 0 )  phi = atan(p0Hat.Y()/p0Hat.X());             
  else if ( p0Hat.Y() > 0  && p0Hat.X() < 0 )  phi = atan(p0Hat.Y()/p0Hat.X())+PI; 
  else if ( p0Hat.Y() < 0  && p0Hat.X() < 0 )  phi = atan(p0Hat.Y()/p0Hat.X())+PI; 
  else if ( p0Hat.Y() < 0  && p0Hat.X() > 0 )  phi = atan(p0Hat.Y()/p0Hat.X())+2*PI;     
  else if ( p0Hat.Y() == 0 && p0Hat.X() == 0 ) phi = 0; //defined by convention
  else if ( p0Hat.Y() == 0 )
  {
    if ( p0Hat.X() > 0 ) phi = 0; 
    else                 phi = PI; 
  }
  else if ( p0Hat.X() == 0 )
  {
    if ( p0Hat.Y() > 0 ) phi = PI/2; 
    else                 phi = 3*PI/2; 
  }
}






//----------------------------------------------------------------------------------------------------
// In active region
bool CalculatePionXS::InActiveRegion( const TVector3& thePos  )
{
  if ( FV_X_BOUND[0] < thePos.X() && thePos.X() < FV_X_BOUND[1] &&
       FV_Y_BOUND[0] < thePos.Y() && thePos.Y() < FV_Y_BOUND[1] &&
       FV_Z_BOUND[0] < thePos.Z() && thePos.Z() < FV_Z_BOUND[1] ) return true;

  return false;
}







//----------------------------------------------------------------------------------------------------
// Track matched track
void CalculatePionXS::TrackMatchedTrack(const size_t& xsRecoTrkId, 
                                        const std::vector<art::Ptr<recob::Track> >& tracklist,
                                        art::FindManyP<anab::Calorimetry>& fmcal)
{
  // ### Define some containers important for XS
  std::vector<double> pitchVec;
  std::vector<double> dEdXVec;
  std::vector<double> eDepVec; 
  std::vector<double> resRanVec;
  std::vector<double> zPosVec;

  int furthestInZCaloPointIndex  = -1;
  TVector3 furthestInZCaloPoint(0,0,0);

  // ### We need the calorimetry information for the track
  art::Ptr<recob::Track> theRecoTrk = tracklist[xsRecoTrkId];
  if (!fmcal.isValid()) return;
  std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(theRecoTrk.key());
  
  // ### Just in case tracking is backwards
  bool isInvertedTracking = false;
  auto realFirstValidPt = (theRecoTrk->TrajectoryPoint(theRecoTrk->FirstValidPoint())).position;
  auto realLastValidPt  = (theRecoTrk->TrajectoryPoint(theRecoTrk->LastValidPoint( ))).position;
  if ( realFirstValidPt.Z() - realLastValidPt.Z() > 0) isInvertedTracking = true;

  // ### Loop over calos
  for (size_t j = 0; j < calos.size(); j++)
  {
    // ### Skip Induction Plane
    if (!calos[j]->PlaneID().isValid)    continue;
    if (calos[j]->PlaneID().Plane == 0)  continue;

    // ### Loop over entries
    size_t nTrkHits = calos[j]->dEdx().size();
    for (size_t k = 0; k < nTrkHits; k++)
    {
      // ### Make sure we can actually see this
      TVector3 theXYZ( calos[j]->XYZ()[k].X(),
                       calos[j]->XYZ()[k].Y(), 
                       calos[j]->XYZ()[k].Z() );
      if ( !InActiveRegion(theXYZ) ) continue;

      // ### We're determining the most downstream point
      if (theXYZ.Z() > furthestInZCaloPoint.Z()) 
      {
        furthestInZCaloPointIndex = k;
        furthestInZCaloPoint = theXYZ;
      }

      // ### Check the dEdX value
      if (calos[j]->dEdx()[k] < 0 || calos[j]->dEdx()[k] > HIT_DEDX_THRESHOLD) continue;

      pitchVec. push_back( calos[j]->TrkPitchVec()[k] );
      dEdXVec.  push_back( calos[j]->dEdx()[k] );
      eDepVec.  push_back( calos[j]->dEdx()[k] * calos[j]->TrkPitchVec()[k] );
      resRanVec.push_back( calos[j]->ResidualRange()[k] );
      zPosVec.  push_back( calos[j]->XYZ()[k].Z() );
    }//<-- End loop over entries 
  }//<-- End loop over calos

  // ### Sanity check
  auto sc = pitchVec.size();
  if (pitchVec.size() != sc && dEdXVec.size()   != sc && 
      eDepVec.size()  != sc && resRanVec.size() != sc && 
      zPosVec.size()  != sc) throw cet::exception("CalculatePionXS") << "Calorimetry sizes are not the same!\n";

  std::cout << "Most downstream point from calorimetry is (" 
            << furthestInZCaloPoint.X() << ", "
            << furthestInZCaloPoint.Y() << ", "
            << furthestInZCaloPoint.Z() << ")\n\n";

  // ### Determine if the particle interacted 
  // ### This is simply if the furthest Z calo point is within active volume
  if (furthestInZCaloPointIndex >= 0)
  {
    // ### Fill histos
    hFurthestInZCaloX->Fill(furthestInZCaloPoint.X());
    hFurthestInZCaloY->Fill(furthestInZCaloPoint.Y());
    hFurthestInZCaloZ->Fill(furthestInZCaloPoint.Z());
    if (InActiveRegion(furthestInZCaloPoint)) fDidDetermineInteracting = 1;
  }

  // ### Make sure we have stuff here to work with
  if (pitchVec.size() == 0) return;

  // ### Check track inversion
  if (isInvertedTracking)
  {
    std::reverse(pitchVec.begin(),  pitchVec.end());
    std::reverse(dEdXVec.begin(),   dEdXVec.end());
    std::reverse(eDepVec.begin(),   eDepVec.end());
    std::reverse(resRanVec.begin(), resRanVec.end());
    std::reverse(zPosVec.begin(),   zPosVec.end());
  }

  // ### Determine if the interaction was inelastic
  if (fDidDetermineInteracting) DetermineInelasticity(xsRecoTrkId, furthestInZCaloPoint, tracklist);

  // ### Get the kinetic energy at "WC 4" 
  // (for single particle gun just use KE at z=0)
  TVector3 theWCMom(0,0,0);
  for (size_t iPt = 0; iPt < fG4PrimaryTrTrjMom[0].size(); iPt++)
  {
    if (fG4PrimaryTrTrjPos[0][iPt].Z() > 0) break;
    theWCMom = fG4PrimaryTrTrjMom[0][iPt];
  }

  // ### This is the first incident KE
  std::vector<double> incidentKEVec;
  double theWCKE = std::sqrt( theWCMom.Mag()*theWCMom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
  incidentKEVec.push_back(theWCKE);

  // ### Fill the energy depositions
  double totalEnDep(0);
  for (size_t iDep = 0; iDep < eDepVec.size(); iDep++)
  {
    // ### Exit if we've passed zero
    totalEnDep = totalEnDep + eDepVec[iDep];
    if ((theWCKE - totalEnDep) < 0) break;
    incidentKEVec.push_back(theWCKE - totalEnDep);
  }

  // ### Fill the Incident and interacting histograms
  for (auto iKE : incidentKEVec) hRecoMCIncidentKE->Fill(iKE);

  fRecoInteractingKE = incidentKEVec.back();
  if (fDidDetermineSignal && incidentKEVec.size()) hRecoMCInteractingKE->Fill(incidentKEVec.back());
}





//----------------------------------------------------------------------------------------------------
// Determine inelasticity
void CalculatePionXS::DetermineInelasticity(const size_t& xsRecoTrkId,
                                            const TVector3& furthestInZCaloPoint,
                                            const std::vector<art::Ptr<recob::Track> >& tracklist)
{
 /*
  * Cuts...
  *
  *   1) There must be another track eminating from vertex
  *        * Angle between secondary and primary > primSecAngle
  *
  *
  */ 


  // ### Look for any daughters around this furthest in Z calo point
  // ### Loop over tracks to find a starting point close to this
  std::cout << "There are " << tracklist.size() << " reconstructed tracks." << std::endl;
  std::map<size_t, bool> theTracksLeaving; // trk id, is tracking inverted

  for(size_t iTrk = 0; iTrk < tracklist.size(); iTrk++)
  {
    // ### Skip our primary track
    if (iTrk == xsRecoTrkId) continue;

    // ### What determines the start and ending point here??
    auto theTrkExtent = tracklist[iTrk]->Extent();
    TVector3 dXYZstart( theTrkExtent.first.X(),  theTrkExtent.first.Y(), theTrkExtent.first.Z() );
    TVector3 dXYZend  ( theTrkExtent.second.X(), theTrkExtent.second.Y(), theTrkExtent.second.Z() );
    double length = (dXYZstart-dXYZend).Mag();
    hSecondaryLength->Fill(length);

    // ### Length cut
    if (length < fSecondaryLengthCut) continue;
    
    // ### Which endpoint is connected to vertex
    double diff1 = (furthestInZCaloPoint - dXYZstart).Mag();
    double diff2 = (furthestInZCaloPoint - dXYZend  ).Mag();

    std::cout << "Distances from vertex " << diff1 << "  " << diff2 << std::endl;
    hVertexDiff->Fill(diff1); hVertexDiff->Fill(diff2);

    size_t minId    = diff1 < diff2 ? 0 : 1;
    std::cout << "MinId is " << minId << std::endl;
    if (diff1 < fVertexCut || diff2 < fVertexCut)
    {
      if (minId == 0) theTracksLeaving.emplace(iTrk, false);
      if (minId == 1) theTracksLeaving.emplace(iTrk, true);
    }

    std::cout << "Track " << iTrk << " reco info\n";
    std::cout << "Start Point "; dXYZstart.Print();
    std::cout << "End Point   "; dXYZend.Print();
    std::cout << "True info\n";
    for (const auto& v : fG4PrimaryVtx) fG4PrimaryTrTrjPos[0][v].Print();
    std::cout << std::endl;
  }//<-- End loop over trks

  std::cout << "Tracks leaving = " << theTracksLeaving.size() << std::endl;

  // ### If there are zero visible tracks leaving this vertex, 
  if (theTracksLeaving.size() == 0) 
  {
    // ### WHAT TO DO HERE
  }

  // ### If there is just 1
  if (theTracksLeaving.size() == 1)
  {
    // ### I'm still interested in supplementing this with 
    // ### small energy deposits around vertex
    auto primTrk = tracklist[xsRecoTrkId];
    auto secTrk  = tracklist[theTracksLeaving.begin()->first];
    bool isInverted = theTracksLeaving.begin()->second;
    
    // ### To compute the angle between tracks, I'll only use the first n
    // ### points, in case there's something funky going on in tracking

    // ### Grab the end points and the nth point from them
    std::vector<TVector3> primPoints, secPoints;
    int primNpointsBack = 10; int primNpointsTot = primTrk->NumberTrajectoryPoints();
    int secNpointsBack  = 10; int secNpointsTot  = secTrk->NumberTrajectoryPoints();

    primNpointsBack = primNpointsBack < primNpointsTot ? primNpointsBack : primNpointsTot;
    secNpointsBack  = secNpointsBack  < secNpointsTot  ? secNpointsBack  : secNpointsTot;

    // ### For primary, just in case, grab the most downstream
    if (primTrk->Start().Z() > primTrk->End().Z())
    {
      // ### It's inverted, grab the first few points
      primPoints.push_back( TVector3(primTrk->LocationAtPoint(primNpointsBack-1).X(), 
                                     primTrk->LocationAtPoint(primNpointsBack-1).Y(), 
                                     primTrk->LocationAtPoint(primNpointsBack-1).Z()) );      
      primPoints.push_back( TVector3(primTrk->LocationAtPoint(0).X(), 
                                     primTrk->LocationAtPoint(0).Y(), 
                                     primTrk->LocationAtPoint(0).Z()) );
    }
    else
    {
      // ### It's not inverted, grab the last few points
      primPoints.push_back( TVector3(primTrk->LocationAtPoint(primNpointsTot-primNpointsBack).X(), 
                                     primTrk->LocationAtPoint(primNpointsTot-primNpointsBack).Y(), 
                                     primTrk->LocationAtPoint(primNpointsTot-primNpointsBack).Z()) );
      primPoints.push_back( TVector3(primTrk->LocationAtPoint(primNpointsTot-1).X(), 
                                     primTrk->LocationAtPoint(primNpointsTot-1).Y(), 
                                     primTrk->LocationAtPoint(primNpointsTot-1).Z()) );                                     
    }

    // ### Now handle secondary
    if (!isInverted) 
    {
      // ### It's the start point
      secPoints.push_back( TVector3(secTrk->LocationAtPoint(0).X(), 
                                    secTrk->LocationAtPoint(0).Y(), 
                                    secTrk->LocationAtPoint(0).Z()) );
      secPoints.push_back( TVector3(secTrk->LocationAtPoint(secNpointsBack-1).X(), 
                                    secTrk->LocationAtPoint(secNpointsBack-1).Y(), 
                                    secTrk->LocationAtPoint(secNpointsBack-1).Z()) );      
    }
    else 
    {
      // ### It's the end point
      secPoints.push_back( TVector3(secTrk->LocationAtPoint(secNpointsTot-1).X(), 
                                    secTrk->LocationAtPoint(secNpointsTot-1).Y(), 
                                    secTrk->LocationAtPoint(secNpointsTot-1).Z()) );
      secPoints.push_back( TVector3(secTrk->LocationAtPoint(secNpointsTot-secNpointsBack).X(), 
                                    secTrk->LocationAtPoint(secNpointsTot-secNpointsBack).Y(), 
                                    secTrk->LocationAtPoint(secNpointsTot-secNpointsBack).Z()) );
    }

    // ### Now form our direction vectors
    TVector3 primDir(primPoints[1] - primPoints[0]);
    TVector3 secDir (secPoints[1]  - secPoints[0]);

    double theta = (180/TMath::Pi())*std::acos( primDir.Unit().Dot(secDir.Unit()) );
    if (theta > fSecondaryAngleCut) fDidDetermineSignal = 1;

    std::cout << "Primary Start Point "; primPoints[0].Print();
    std::cout << "Primary End Point   "; primPoints[1].Print();
    std::cout << "Second. Start Point "; secPoints[0].Print();
    std::cout << "Second. End Point   "; secPoints[1].Print();
    std::cout << "THETA " << theta << std::endl;
  }//<-- End if 1 visible secondary

  // ### If there are at least 2 visible tracks leaving this vertex, yes
  if (theTracksLeaving.size() >= 2) fDidDetermineSignal = 1;

  std::cout << fDidDetermineSignal << std::endl;
  if (fDidDetermineSignal) std::cout << "Determined as inelastic!" << std::endl;
  else std::cout << "Determined NOT inelastic!\n";
}


DEFINE_ART_MODULE(CalculatePionXS)
}
