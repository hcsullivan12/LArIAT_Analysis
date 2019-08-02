////////////////////////////////////////////////////////////////////////
// Class:       PionInelasticXS
// Module Type: analyzer
// File:        PionInelasticXS_module.cc
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
#include "LArIATAnaModule/InelasticSubClassifier.h"

namespace lariat {

class PionInelasticXS;

class PionInelasticXS : public art::EDAnalyzer {
public:
  explicit PionInelasticXS(fhicl::ParameterSet const & p);

  PionInelasticXS(PionInelasticXS const &) = delete;
  PionInelasticXS(PionInelasticXS &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  
private:

  void resetVars();
  inline void PrintVec(const TVector3& v){std::cout<<"("<<v.X()<<", "<<v.Y()<<", "<<v.Z()<<")\n";}
  bool filter(size_t& xsRecoTrkId, const std::vector<art::Ptr<recob::Track> >& tracklist);
  void fillTruthInfo(const sim::ParticleList& plist);
  bool inActiveRegion(const TVector3& thePos);
  void ComputeAngles(float& phi, float& theta, const TVector3& p0Hat);
  void trackMatchedTrack(const size_t& xsRecoTrkId, 
                         const std::vector<art::Ptr<recob::Track> >& tracklist,
                         art::FindManyP<anab::Calorimetry>& fmcal);
  void determineInelasticity(const size_t& xsRecoTrkId,
                             const TVector3& furthestInZCaloPoint,
                             const std::vector<art::Ptr<recob::Track> >& tracklist);

  // labels
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;

  // Useful variables   
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

  // Cuts

  // maximum distance to allow track to be seperated from vertex
  float fVertexCut;
  // minimum length of secondary tracks
  float fSecondaryLengthCut;
  // minimum angle between primary and secondary (degrees)
  float fSecondaryAngleCut;
  // Is data
  bool fIsData;

  // VARIABLES FOR G4 INFO  
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

  //  CONSTANTS  
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
   
   
  // the assumed energy loss between the cryostat and the TPC 
  const float ENTRY_TPC_ENERGY_LOSS = 36; //MeV
   
  // constants for cross section calculation
  const double RHO            = 1395;          //kg/m^3
  const double MOLAR_MASS     = 39.948;        //g/mol
  const double G_PER_KG       = 1000;
  const double AVOGADRO       = 6.022140857e+23;        //number/mol
  const double NUMBER_DENSITY =  (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO;
  const double SLAB_WIDTH     = 0.0045;        //in m
  const double PI             = 3.141592654;
  const double BARN_PER_M2    = 1e28;

  // threshold for hadron dEdx
  double HIT_DEDX_THRESHOLD = 40.;

  // Histograms
  TH1D* hRecoUpstreamZPos;
  TH1D* hRecoFurthestInZCaloX;
  TH1D* hRecoFurthestInZCaloY;
  TH1D* hRecoFurthestInZCaloZ;
  TH1D* hRecoIncidentKE;
  TH1D* hRecoInteractingKE;
  TH1D* hRecoVertexDiff;
  TH1D* hRecoSecLength;

  TH1D* hMCLength;
};

//----------------------------------------------------------------------------------------------------
// Contructor
PionInelasticXS::PionInelasticXS(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);
}

//----------------------------------------------------------------------------------------------------
// Reconfigure
void PionInelasticXS::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel        = p.get< std::string >("TrackModuleLabel"      , "pmtrack");
  fCalorimetryModuleLabel  = p.get< std::string >("CalorimetryModuleLabel", "calo"     );

  fVertexCut          = p.get< float >("VertexCut",          2.0);
  fSecondaryLengthCut = p.get< float >("SecondaryLengthCut", 3.0); 
  fSecondaryAngleCut  = p.get< float >("SecondaryAngleCut",  10.0); 

  fIsData = p.get<bool>("IsData", true);
}

//----------------------------------------------------------------------------------------------------
// Analyze
void PionInelasticXS::analyze(art::Event const & e)
{
  //Reset vars
  resetVars();

  fRun    = e.run();
  fSubrun = e.subRun();
  fEvent  = e.event();

  // FILL TRUTH INFO  
  if (!fIsData) 
  {
    // Backtracker to recover truth information
    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    fillTruthInfo(plist);
  }

  // Get the track information
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << std::endl;

  // Let's get the reco and wc stuff
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;   // Container of wc tracks  
  std::vector<art::Ptr<ldp::WCTrack> >     wctrack;         // Vector of wc tracks  
  art::Handle< std::vector<recob::Track> > trackListHandle; // Container of reco tracks
  std::vector<art::Ptr<recob::Track> >     tracklist;       // Vector of wc tracks  

  // Find which recontructed track is associated with the WC
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;
  if(!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return;

  // Fill the container of reco tracks
  art::fill_ptr_vector(wctrack, wctrackHandle);
  art::fill_ptr_vector(tracklist, trackListHandle);

  art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
  int matchedRecoTrkKey = -99999;
  if (fWC2TPC.isValid())
  {
    for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn)
    {
      cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));
      if (!trackWC2TPC) continue;
      recob::Track const &aTrack(trackWC2TPC.ref());
      matchedRecoTrkKey = aTrack.ID(); // This checks out OK
    }
  } // if there's an associated track
  // Now we know what ID identifies the reco track we're interested in.

  // If we didn't get a track :,(
  if (matchedRecoTrkKey == -99999) return;

  // Track the matched track
  trackMatchedTrack(matchedRecoTrkKey, tracklist, fmcal);

  // Fill our tree and get out of here
  fTree->Fill();
}

//----------------------------------------------------------------------------------------------------
// Reset vars
void PionInelasticXS::resetVars()
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
void PionInelasticXS::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  hRecoUpstreamZPos       = tfs->make<TH1D>("hRecoUpstreamZPos",    "Minimum Z Position Upstream", 1000, 0, 100);
  hMCLength               = tfs->make<TH1D>("hMCLength",          "Primary True Length", 1000, 0, 100);
  hRecoFurthestInZCaloX   = tfs->make<TH1D>("hRecoFurthestInZCaloX",   "Furthest Calorimetry X Position", 1000, 0, 47);
  hRecoFurthestInZCaloY   = tfs->make<TH1D>("hRecoFurthestInZCaloY",   "Furthest Calorimetry Y Position", 1000, -20, 20);
  hRecoFurthestInZCaloZ   = tfs->make<TH1D>("hRecoFurthestInZCaloZ",   "Furthest Calorimetry Z Position", 1000, 0, 100);
  hRecoIncidentKE    = tfs->make<TH1D>("hRecoIncidentKE",    "Incident",    23, -50, 1100);
  hRecoInteractingKE = tfs->make<TH1D>("hRecoInteractingKE", "Interacting", 23, -50, 1100);
  hRecoVertexDiff          = tfs->make<TH1D>("hRecoVertexDiff",          "Distance of Candidate Endpoints", 100, 0, 50); 
  hRecoSecLength     = tfs->make<TH1D>("hRecoSecLength",     "Length of Secondary Tracks", 100, 0, 100);

  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&fRun      ,"run/I");
  fTree->Branch("subrun"   ,&fSubrun   ,"subrun/I");
  fTree->Branch("event"    ,&fEvent    ,"event/I");
  fTree->Branch("isTrackInteracting"  ,&fIsTrackInteracting   ,"isTrackInteracting/I");
  fTree->Branch("isTrackSignal"       ,&fIsTrackSignal        ,"isTrackSignal/I");
  fTree->Branch("trueProcess"         ,&fG4PrimaryProcesses);
  fTree->Branch("didDetermineInteracting" , &fDidDetermineInteracting   ,"didDetermineInteracting/I");
  fTree->Branch("didDetermineSignal" , &fDidDetermineSignal   ,"didDetermineSignal/I");
}


//----------------------------------------------------------------------------------------------------
// Fill truth info
void PionInelasticXS::fillTruthInfo(const sim::ParticleList& plist)
{
  // Loop over g4 particles 
  // Counter for the primary
  size_t iPrim(0);
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    auto mcParticle = plist.Particle(iG4);

    // If this is not a primary, skip it
    if (mcParticle->Process().find("primary") == std::string::npos) continue;
     
    // Store the track id 
    fG4PrimaryTrkId.push_back(mcParticle->TrackId());

    // Store the processes for this primary, only those that occur in TPC
    simb::MCTrajectory truetraj = mcParticle->Trajectory();
    auto thisTrjProcessMap = truetraj.TrajectoryProcesses();
  
    // If there's nothing, check the end of the track
    if (thisTrjProcessMap.size() == 0)
    {
      // Get the vertex
      if (!inActiveRegion(mcParticle->EndPosition().Vect())) continue;
      fG4PrimaryVtx.push_back( mcParticle->NumberTrajectoryPoints()-1 );

      auto tempMom = mcParticle->Momentum(mcParticle->NumberTrajectoryPoints()-2).Vect();
      fTrueInteractingKE = std::sqrt( tempMom*tempMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
 
      // If it doesn't have daughters, it's through going
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
        if (!inActiveRegion( truetraj.Position(interestingPoint).Vect()) ) continue;

        fG4PrimaryVtx.push_back( interestingPoint );
        fG4PrimaryProcesses.push_back( truetraj.KeyToProcess(couple.second) );

        auto tempMom = mcParticle->Momentum(interestingPoint-1).Vect();
        fTrueInteractingKE = std::sqrt( tempMom*tempMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
      }
    }//<-- End if map is not zero

    // Store the positions and momentum 
    fG4PrimaryPos0.push_back( mcParticle->Position().Vect() );
    fG4PrimaryPosf.push_back( mcParticle->EndPosition().Vect() );
    fG4PrimaryMom0.push_back( mcParticle->Momentum().Vect()*1000 ); // convert to MeV
    
    // ### Fill true length histo;
    hMCLength->Fill( (fG4PrimaryPosf[iPrim] - fG4PrimaryPos0[iPrim]).Mag() );
    
    // ### Project onto tpc 
    TVector3 thisPosProjVec = fG4PrimaryPos0[iPrim] - ( fG4PrimaryPos0[iPrim].Z()/fG4PrimaryMom0[iPrim].Z() )*fG4PrimaryMom0[iPrim];
    fG4PrimaryProjPos0.push_back(thisPosProjVec);

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
bool PionInelasticXS::filter(size_t& xsRecoTrkId, const std::vector<art::Ptr<recob::Track> >& tracklist)
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
    hRecoUpstreamZPos->Fill(tempMinZ);
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
// In active region
bool PionInelasticXS::inActiveRegion( const TVector3& thePos  )
{
  if ( FV_X_BOUND[0] < thePos.X() && thePos.X() < FV_X_BOUND[1] &&
       FV_Y_BOUND[0] < thePos.Y() && thePos.Y() < FV_Y_BOUND[1] &&
       FV_Z_BOUND[0] < thePos.Z() && thePos.Z() < FV_Z_BOUND[1] ) return true;

  return false;
}


//----------------------------------------------------------------------------------------------------
// Track matched track
void PionInelasticXS::trackMatchedTrack(const size_t& xsRecoTrkId, 
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
      if ( !inActiveRegion(theXYZ) ) continue;

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
      zPosVec.size()  != sc) throw cet::exception("PionInelasticXS") << "Calorimetry sizes are not the same!\n";

  std::cout << "Most downstream point from calorimetry is (" 
            << furthestInZCaloPoint.X() << ", "
            << furthestInZCaloPoint.Y() << ", "
            << furthestInZCaloPoint.Z() << ")\n\n";

  // ### Determine if the particle interacted 
  // ### This is simply if the furthest Z calo point is within active volume
  if (furthestInZCaloPointIndex >= 0)
  {
    // ### Fill histos
    hRecoFurthestInZCaloX->Fill(furthestInZCaloPoint.X());
    hRecoFurthestInZCaloY->Fill(furthestInZCaloPoint.Y());
    hRecoFurthestInZCaloZ->Fill(furthestInZCaloPoint.Z());
    if (inActiveRegion(furthestInZCaloPoint)) fDidDetermineInteracting = 1;
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
  if (fDidDetermineInteracting) determineInelasticity(xsRecoTrkId, furthestInZCaloPoint, tracklist);

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
  for (auto iKE : incidentKEVec) hRecoIncidentKE->Fill(iKE);

  fRecoInteractingKE = incidentKEVec.back();
  if (fDidDetermineSignal && incidentKEVec.size()) hRecoInteractingKE->Fill(incidentKEVec.back());
}

//----------------------------------------------------------------------------------------------------
// Determine inelasticity
void PionInelasticXS::determineInelasticity(const size_t& xsRecoTrkId,
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
    hRecoSecLength->Fill(length);
  
    std::cout << "Track " << iTrk << " reco info\n";
    std::cout << "Start Point "; PrintVec(dXYZstart);
    std::cout << "End Point   "; PrintVec(dXYZend);
    std::cout << "Length      " << length << std::endl;
    std::cout << "True info\n";
    for (const auto& v : fG4PrimaryVtx) PrintVec(fG4PrimaryTrTrjPos[0][v]);
    std::cout << std::endl;
 

    // ### Length cut
    if (length < fSecondaryLengthCut) continue;
    
    // ### Which endpoint is connected to vertex
    double diff1 = (furthestInZCaloPoint - dXYZstart).Mag();
    double diff2 = (furthestInZCaloPoint - dXYZend  ).Mag();

    std::cout << "Distances from vertex " << diff1 << "  " << diff2 << std::endl;
    double diffArr[2] = {diff1, diff2};    
    size_t minId    = diffArr[0] < diffArr[1] ? 0 : 1;

    // ### We mant the closest point
    if (diff1 < fVertexCut || diff2 < fVertexCut)
    {
      if (minId == 0) {theTracksLeaving.emplace(iTrk, false); hRecoVertexDiff->Fill(diffArr[0]);}
      if (minId == 1) {theTracksLeaving.emplace(iTrk, true); hRecoVertexDiff->Fill(diffArr[1]);}
    }

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

    std::cout << "Primary Start Point "; PrintVec(primPoints[0]);
    std::cout << "Primary End Point   "; PrintVec(primPoints[1]);
    std::cout << "Second. Start Point "; PrintVec(secPoints[0]);
    std::cout << "Second. End Point   "; PrintVec(secPoints[1]);
    std::cout << "THETA " << theta << std::endl;
  }//<-- End if 1 visible secondary

  // ### If there are at least 2 visible tracks leaving this vertex, yes
  if (theTracksLeaving.size() >= 2) fDidDetermineSignal = 1;

  std::cout << fDidDetermineSignal << std::endl;
  if (fDidDetermineSignal) std::cout << "Determined as inelastic!" << std::endl;
  else std::cout << "Determined NOT inelastic!\n";
}


DEFINE_ART_MODULE(PionInelasticXS)
}
