////////////////////////////////////////////////////////////////////////
// Class:       CalculatePionXS
// Module Type: analyzer
// File:        CalculatePionXS_module.cc
//
// Generated at Thu May 30 13:07:47 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


class CalculatePionXS;

class CalculatePionXS : public art::EDAnalyzer {
public:
  explicit CalculatePionXS(fhicl::ParameterSet const & p);

  CalculatePionXS(CalculatePionXS const &) = delete;
  CalculatePionXS(CalculatePionXS &&) = delete;

  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

// ====================================================================================================
// ====================================  CUTS AND CONSTANTS  ==========================================
// ====================================================================================================
// tpc boundaries
float TPC_X_BOUND[2] = {   0.0, 47.0 };
float TPC_Y_BOUND[2] = { -20.0, 20.0 };
float TPC_Z_BOUND[2] = {   0.0, 90.0 };

// fiducial volume definition
float FV_X_BOUND[2] = {   2.0, 45.0 };
float FV_Y_BOUND[2] = { -18.0, 18.0 };
float FV_Z_BOUND[2] = {   0.0, 88.0 };

// mass of pion in MeV
float PARTICLE_MASS(139.57); 

// number of centimeters in Z we require a track to have a spacepoint
float FIRST_SP_Z_POS(2.0);

// portion of upstream TPC which we will restrict the number of tracks
float UPPER_PART_OF_TPC(14.0);

// number of upper tpc tracks allowed
size_t N_UPPER_TPC_TRACKS(4);

// Delta X Between Wire Chamber Track and TPC Track 
float DELTA_X_BOUND[2] = { -2.0, 6.0 };

// Delta Y Between Wire Chamber Track and TPC Track 
float DELTA_Y_BOUND[2] = { -3.0, 6.0 };

// the assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(36); //MeV

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
float EVENT_WEIGHT(1.0);

// True  = Use the momentum based weighting  
// False = Don't weight events
bool USE_EVENT_WEIGHT(true);

// constants for cross section calculation
double RHO(1395);                 //kg/m^3
double MOLAR_MASS(39.948);          //g/mol
double G_PER_KG(1000);
double AVOGADRO(6.022140857e+23);        //number/mol
double NUMBER_DENSITY( (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO );
double SLAB_WIDTH(0.0045);        //in m
double PI(3.141592654);
double BARN_PER_M2(1e28);

// cut on angle between wc and tpc track (degrees)
float ALPHA_CUT(10.);

};

//----------------------------------------------------------------------------------------------------
// Contructor
CalculatePionXS::CalculatePionXS(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{}

//----------------------------------------------------------------------------------------------------
// Analyze
void CalculatePionXS::analyze(art::Event const & e)
{
  // #########################################
  // ### Increment our total event counter ###
  // #########################################
  nTotalEvents++; 
  if (nTotalEvents%1000 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl; 


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Looking at only MC truth/G4 information
 
  // =====================================================================================
  // ==============================  VARIABLES FOR G4 INFO  ==============================
  // =====================================================================================
  // count the number of primaries
  // this is to save us from a memory crash
  // which is not an easy bug to hunt down
  int nPrimary(0);                                       // number of primaries
  int maxTrTrjPoints(0);                                 // max number of truth traj points
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    if (process_primary[iG4] == 1) nPrimary++;
    if (process_primary[iG4] == 1 && maxTrTrjPoints < NTrTrajPts[iG4]) maxTrTrjPoints = NTrTrajPts[iG4];
  }
  std::vector<int> g4PrimaryTrkId;                       // track id 
  std::vector<std::vector<std::string>> g4PrimaryProcesses; // processes
  std::vector<TVector3> g4PrimaryPos0;                   // start pos
  std::vector<TVector3> g4PrimaryPosf;                   // final pos 
  std::vector<TVector3> g4PrimaryMom0;                   // momentum
  std::vector<TVector3> g4PrimaryProjPos0;               // projected position
  std::vector<std::vector<TVector3>> g4PrimaryTrTrjPos;  // trajectory sp
  std::vector<std::vector<TVector3>> g4PrimaryTrTrjMom;  // trajectory momentum 
  std::vector<std::vector<TVector3>> g4PrimaryVtx;       // primary interaction points
  
  // ##############################
  // ### Loop over g4 particles ###
  // ##############################
  {
  // counter for the primary
  size_t iPrim(0);
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    // #########################################
    // ### If this is not a primary, skip it ###
    // #########################################
    if (process_primary[iG4] == 0) continue;
     
    // store the track id 
    g4PrimaryTrkId.push_back(TrackId[iG4]);
    // store the processes for this primary, only those that occur in TPC
    std::vector<std::string> proc;
    std::vector<TVector3>    vtx;
    for (size_t iPt = 0; iPt < (*InteractionPoint).size(); iPt++)
    {
      // get the position of this interation point
      TVector3 thePos( MidPosX[iG4][(*InteractionPoint)[iPt]], 
                       MidPosY[iG4][(*InteractionPoint)[iPt]],
                       MidPosZ[iG4][(*InteractionPoint)[iPt]] );
      if (!InActiveRegion(thePos)) continue;
    
      proc.push_back( (*InteractionPointType)[iPt] );
      vtx.push_back(thePos); 
    }
    g4PrimaryProcesses.push_back(proc);
    g4PrimaryVtx.push_back(vtx);
    // store the positions and momentum 
    g4PrimaryPos0.push_back( TVector3(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]) );
    g4PrimaryPosf.push_back( TVector3(EndPointx[iG4],   EndPointy[iG4],   EndPointz[iG4]) );
    g4PrimaryMom0.push_back( TVector3(Px[iG4]*1000,     Py[iG4]*1000,     Pz[iG4]*1000) ); // convert to MeV
    // setting a global event weight 
    // setting Event weight 
    if(USE_EVENT_WEIGHT)
    {
      if(g4PrimaryMom0[iPrim].Z() > 0   && g4PrimaryMom0[iPrim].Z() < 100) EVENT_WEIGHT = 0.010;
      if(g4PrimaryMom0[iPrim].Z() > 100 && g4PrimaryMom0[iPrim].Z() < 200) EVENT_WEIGHT = 0.020;
      if(g4PrimaryMom0[iPrim].Z() > 200 && g4PrimaryMom0[iPrim].Z() < 300) EVENT_WEIGHT = 0.100;
      if(g4PrimaryMom0[iPrim].Z() > 300 && g4PrimaryMom0[iPrim].Z() < 400) EVENT_WEIGHT = 0.535;
      if(g4PrimaryMom0[iPrim].Z() > 400 && g4PrimaryMom0[iPrim].Z() < 500) EVENT_WEIGHT = 0.840;
      if(g4PrimaryMom0[iPrim].Z() > 500 && g4PrimaryMom0[iPrim].Z() < 600) EVENT_WEIGHT = 0.965;
      if(g4PrimaryMom0[iPrim].Z() > 600 && g4PrimaryMom0[iPrim].Z() < 700) EVENT_WEIGHT = 1.000;
      if(g4PrimaryMom0[iPrim].Z() > 700 && g4PrimaryMom0[iPrim].Z() < 800) EVENT_WEIGHT = 0.620;
      if(g4PrimaryMom0[iPrim].Z() > 800 && g4PrimaryMom0[iPrim].Z() < 900) EVENT_WEIGHT = 0.225;
      if(g4PrimaryMom0[iPrim].Z() > 900 && g4PrimaryMom0[iPrim].Z() <1000) EVENT_WEIGHT = 0.094;
      if(g4PrimaryMom0[iPrim].Z() >1000 && g4PrimaryMom0[iPrim].Z() <1100) EVENT_WEIGHT = 0.0275;
      if(g4PrimaryMom0[iPrim].Z() >1100)                             EVENT_WEIGHT = 0.010;
    }
    // fill momentum histos
    hMCPrimaryPx->Fill(g4PrimaryMom0[iPrim].X(),   EVENT_WEIGHT);
    hMCPrimaryPy->Fill(g4PrimaryMom0[iPrim].Y(),   EVENT_WEIGHT);
    hMCPrimaryPz->Fill(g4PrimaryMom0[iPrim].Z(),   EVENT_WEIGHT);
    hMCPrimaryP ->Fill(g4PrimaryMom0[iPrim].Mag(), EVENT_WEIGHT);
    // fill true length histo;
    hTrueLength->Fill( (g4PrimaryPosf[iPrim] - g4PrimaryPos0[iPrim]).Mag() );
    // project onto tpc 
    TVector3 thisPosProjVec = g4PrimaryPos0[iPrim] - ( g4PrimaryPos0[iPrim].Z()/g4PrimaryMom0[iPrim].Z() )*g4PrimaryMom0[iPrim];
    g4PrimaryProjPos0.push_back(thisPosProjVec);
    // fill the proj histos
    hMCPrimaryProjX0->Fill( g4PrimaryProjPos0[iPrim].X() );
    hMCPrimaryProjY0->Fill( g4PrimaryProjPos0[iPrim].Y() );
    hMCPrimaryProjZ0->Fill( g4PrimaryProjPos0[iPrim].Z() );
    // store the trajectory points and momentum
    std::vector<TVector3> temp;
    g4PrimaryTrTrjPos.push_back(temp);
    g4PrimaryTrTrjMom.push_back(temp);
    for (size_t iPoint = 0; iPoint < NTrTrajPts[iG4]; iPoint++)
    {
      g4PrimaryTrTrjPos[iPrim].push_back( TVector3(MidPosX[iG4][iPoint],    MidPosY[iG4][iPoint],    MidPosZ[iG4][iPoint]) );
      g4PrimaryTrTrjMom[iPrim].push_back( TVector3(MidPx[iG4][iPoint]*1000, MidPy[iG4][iPoint]*1000, MidPz[iG4][iPoint]*1000)   ); // convert to MeV
    }
    iPrim++;
  } 
  }

  // ==========================================================================
  // =================  LOOKING AT EVENTS THAT ENTER THE TPC  =================
  // ==========================================================================
  // ###########################
  // ### Loop over primaries ###
  // ###########################
  bool isGoodEvent(false);
  for (size_t iPrim = 0; iPrim < nPrimary; iPrim++)
  {
    // only look if a primary entered the tpc
    if ( g4PrimaryPosf[iPrim].Z() > 0 ) isGoodEvent = true;
    if ( !isGoodEvent ) 
    {
      hMCPrimaryMissedTpcX->Fill(g4PrimaryPosf[iPrim].X());
      hMCPrimaryMissedTpcY->Fill(g4PrimaryPosf[iPrim].Y());
      hMCPrimaryMissedTpcZ->Fill(g4PrimaryPosf[iPrim].Z());
      continue;
    }
    // ####################################
    // ### This primary entered the TPC ###
    // ####################################
    nPrimariesEntered++;
    // calculate energy loss
    float energyLoss(0);
    // loop over trj points for this primary
    for (size_t iPoint = 0; iPoint < g4PrimaryTrTrjPos[iPrim].size(); iPoint++)
    {
      // only look at the upstream portion
      if ( g4PrimaryTrTrjPos[iPrim][iPoint].Z() > 0 ) continue;
      // ignore last point
      if ( (iPoint+1) >= g4PrimaryTrTrjPos[iPrim].size() ) break;
      auto mom1Vec = g4PrimaryTrTrjMom[iPrim][iPoint];
      auto mom2Vec = g4PrimaryTrTrjMom[iPrim][iPoint+1];
      float energy1 = std::sqrt( mom1Vec.Mag()*mom1Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS );
      float energy2 = std::sqrt( mom2Vec.Mag()*mom2Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS ); 
      energyLoss += (energy1 - energy2);
    }//<--- End loop over true traj points
    hMCELossUpstream->Fill(energyLoss);
  }//<--- End loop over primaries
  if (!isGoodEvent) continue;
  nGoodMCEvents++;
  
 // End looking at only MC truth/G4 information
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
  for (size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
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
    // ###################################
    // ### Loop over trajectory points ###
    // ###################################
    for (size_t iTrjPoint = 0; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
    {
      // ######################
      // ### Fiducial check ###
      // ######################
      TVector3 theTrjPointVec( trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint] );
      TVector3 theTrjP0HatVec( pHat0_X[iTrk][iTrjPoint], pHat0_Y[iTrk][iTrjPoint], pHat0_Z[iTrk][iTrjPoint] );
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
    //hDataUpstreamZPos->Fill(tempMinZ);
  }//<-- End loop over reconstructed tracks
  // skip events that do not have a trk at the front
  if (frontFaceTrksId.size() == 0) continue;
  nEventsFrontTpcTrk++;
  // skip events that have too many tracks upstream        
  if(nUpstreamTpcTrks > N_UPPER_TPC_TRACKS || nUpstreamTpcTrks == 0) continue;
  nEventsUpperTpcTrkCount++;
  // ========================================================================
  // =================  MATCHING MC TO RECONSTRUCTED TRACK  =================
  // ========================================================================
  size_t nMcTpcMatch(0);
  // ###########################################
  // ### Loop over all the front face Tracks ###
  // ###########################################
  for(size_t iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
  {
    auto deltaVec = frontFaceTrksPos[iFrFaTrk] - g4PrimaryProjPos0[0];
    //hDeltaX->Fill(deltaVec.X());
    //hDeltaY->Fill(deltaVec.Y());
    //hDeltaZ->Fill(deltaVec.Z());
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
  TVector3 mcP0Hat = g4PrimaryMom0[0]; 
  mcP0Hat = mcP0Hat.Unit();
  float    mcPhi(0), mcTheta(0);

  // compute the angles for primary
  ComputeAngles( mcPhi, mcTheta, mcP0Hat );
  // ##############################################################
  // ### Calculating the angles for the front face tracks (TPC) ###
  // ##############################################################
  for(int iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
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
  size_t xsRecoTrkId;
  for(int iFrFaTrk = 0; iFrFaTrk < frontFaceTrksId.size(); iFrFaTrk++)
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
    //hAlpha->Fill(alpha);
    // checking if this made our cut
    if (alpha < ALPHA_CUT)
    {
      isAlphaMatch = true;
      xsRecoTrkId = frontFaceTrksId[iFrFaTrk];
    }
  }
  // did we get a match?
  if (nMcTpcMatch > 0) nEventsDeltaMatch++;
  // force there to be one match
  if (nMcTpcMatch != 1) continue;
  nEventsWcTpcUniqueMatch++;
  // force this unique match to pass our alpha cut
  if (!isAlphaMatch) continue;
  nEventsWcTpcUniqueMatchAlpha++;
  
  
  // ############################################################
  // ### We are know tracking a uniquely WC/TPC matched track ###
  // ###         The ID is in xsRecoTrkId variable            ###
  // ############################################################
  
  // ==========================================================================
  // =================  BEGIN TRACKING OUR MATCHED TPC TRACK  =================
  // ==========================================================================
  // define the containers of the important information about the track for the XS
  std::vector<double> pitchVec;
  std::vector<double> dEdXVec;
  std::vector<double> eDepVec; // Energy deposition at each slice. Simply dEdX*Pitch
  std::vector<double> resRanVec;
  std::vector<double> zPosVec;
  std::vector<double> incidentKEVec;

  // ### loop over the reconstructed tracks
  size_t furtherstInZCaloPointIndex  = 0;
  double furtherstInZCaloPointZ      = 0.;
  bool   isInteracting          = false;
  bool   isInelastic            = false;
  for(size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
  {
    // only look at the one that passed the WC_TPC match and cuts
    if(iTrk != xsRecoTrkId) continue;
    // we need the calorimetry information for this track
    size_t limit = sizeof(trkxyz[xsRecoTrkId][1])/sizeof(*trkxyz[xsRecoTrkId][1]);
    for (size_t k = 0; k < limit; k++)
    {
      // make sure we can actually see this
      TVector3 theXYZ(trkxyz[xsRecoTrkId][1][k][0], trkxyz[xsRecoTrkId][1][k][1], trkxyz[xsRecoTrkId][1][k][2]); 
      if ( !InActiveRegion(theXYZ) ) continue;
      // we're determining the most downstream point
      if (theXYZ.Z() > furtherstInZCaloPointZ) 
      {
        furtherstInZCaloPointIndex = k;
        furtherstInZCaloPointZ = theXYZ.Z();
      }
      pitchVec. push_back( trkpitchhit[xsRecoTrkId][1][k] );
      dEdXVec.  push_back( trkdedx[xsRecoTrkId][1][k]     );
      eDepVec.  push_back( trkdedx[xsRecoTrkId][1][k]*trkpitchhit[xsRecoTrkId][1][k] );
      resRanVec.push_back( trkrr[xsRecoTrkId][1][k] );
      zPosVec.  push_back( theXYZ.Z() );
    }//<-- End loop over track points 
  }//<--- End loop over tracks
  // determine if the particle interacted 
  if (furtherstInZCaloPointIndex)
  {
    TVector3 theXYZ( trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][0], 
                     trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][1], 
                     trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][2] ); 
    // fill histos
    //hFurtherstInZCaloX->Fill(theXYZ.X());
    //hFurtherstInZCaloY->Fill(theXYZ.Y());
    //hFurtherstInZCaloZ->Fill(theXYZ.Z());
    if (InActiveRegion(theXYZ)) isInteracting = true;
  }
  // determine if the interaction was inelastic
  if (isInteracting)
  {
    // i think for now, we will just look for any daughters around this 
    // furthest in Z calo point
    TVector3 theXYZ( trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][0], 
                     trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][1], 
                     trkxyz[xsRecoTrkId][1][furtherstInZCaloPointIndex][2] );
    // loop over tracks to find a starting point close to this
    size_t nTracksLeaving(0);
    for(size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // skip our primary track
      //if (iTrk == xsRecoTrkId) continue;
      // what determines the start and ending point here??
      TVector3 dXYZstart( trkvtxx[iTrk], trkvtxy[iTrk], trkvtxz[iTrk] );
      TVector3 dXYZend  ( trkendx[iTrk], trkendy[iTrk], trkendz[iTrk] );
      double diff = (theXYZ - dXYZstart).Mag();
      if ( diff < 0.01 ) nTracksLeaving++;
      for (const auto p : g4PrimaryProcesses[0]) std::cout << p << std::endl;
      std::cout << "Reco info\n";
      std::cout << "Calo Point  "; theXYZ.Print();
      std::cout << "Start Point "; dXYZstart.Print();
      std::cout << "End Point   "; dXYZend.Print();
      std::cout << "True info\n";
      for (const auto& v : g4PrimaryVtx[0]) v.Print();
      std::cout << std::endl;
    }
    // there must be at least 2 visible tracks leaving this vertex
    if (nTracksLeaving >= 2) isInelastic = true;
    std::cout << "\n" << nTracksLeaving << " " << ntracks_reco << std::endl;
  }
  // make sure we have stuff here to work with
  if (pitchVec.size() == 0) continue;
  // get the kinetic energy at "WC 4" 
  // (for single particle gun just use KE at z=0)
  TVector3 theWCMom(0,0,0);
  for (size_t iPt = 0; iPt < maxTrTrjPoints; iPt++)
  {
    if (g4PrimaryTrTrjPos[0][iPt].Z() > 0) break;
    theWCMom = g4PrimaryTrTrjMom[0][iPt];
  }
  double theWCKE = std::sqrt( theWCMom.Mag()*theWCMom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
  incidentKEVec.push_back(theWCKE);
  // fill the energy depositions
  double totalEnDep(0);
  for (size_t iDep = 0; iDep < eDepVec.size(); iDep++)
  {
    totalEnDep = totalEnDep + eDepVec[iDep];
    incidentKEVec.push_back(theWCKE - totalEnDep);
  }
  // fill the Incident and interacting histograms
  //for (auto iKE : incidentKEVec) hRecoMCIncidentKE->Fill(iKE);
  //if (isInelastic && incidentKEVec.size() != 0) hRecoMCInteractingKE->Fill(incidentKEVec.back());
}

//----------------------------------------------------------------------------------------------------
// Begin job
void CalculatePionXS::beginJob()
{
  // Implementation of optional member function here.
}

//----------------------------------------------------------------------------------------------------
// End job
void CalculatePionXS::endJob()
{
  // Implementation of optional member function here.
}

//----------------------------------------------------------------------------------------------------
// Reconfigure
void CalculatePionXS::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CalculatePionXS)
