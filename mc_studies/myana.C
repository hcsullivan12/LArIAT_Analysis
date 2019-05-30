#define myana_cxx
#include "myana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>

#include <assert.h>

// ============================================================================================
// ====================================  PROTOTYPES  ==========================================
// ============================================================================================
void MakePlots();
void ComputeAngles(float& phi, float& theta, const TVector3& p0Hat);
std::string ConvertProcessToString(const int& i);
void CalculateG4Xs();

// ==================================================================================================
// ====================================  USEFUL VARIABLES  ==========================================
// ==================================================================================================
// Counters
size_t nTotalEvents(0),      nPrimariesEntered(0),       nGoodMCEvents(0),
       nEventsInelastic(0),  nEventsFrontTpcTrk(0),      nEventsUpperTpcTrkCount(0),
       nEventsDeltaMatch(0), nEventsWcTpcUniqueMatch(0), nEventsWcTpcUniqueMatchAlpha(0),
       nXsG4Signal(0);

// slabs
std::vector<float> slabsInZ;


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


// ============================================================================================
// ====================================  HISTOGRAMS  ==========================================
// ============================================================================================
TH1D* hMCELossUpstream = new TH1D("hMCELossUpstream", "MC Energy Loss Upstream", 1000, 0, 1000);

TH1D* hMCPrimaryMissedTpcX = new TH1D("hMCPrimaryMissedTpcX", "MC Primary Missed TPC X", 200, -50, 50);
TH1D* hMCPrimaryMissedTpcY = new TH1D("hMCPrimaryMissedTpcY", "MC Primary Missed TPC Y", 200, -50, 50);
TH1D* hMCPrimaryMissedTpcZ = new TH1D("hMCPrimaryMissedTpcZ", "MC Primary Missed TPC Z", 200, -110, 10);

TH1D *hMCPrimaryProjX0 = new TH1D("hMCPrimaryProjX0", "Primary Particle X_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjY0 = new TH1D("hMCPrimaryProjY0", "Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjZ0 = new TH1D("hMCPrimaryProjZ0", "Primary Particle Z_{0}", 100, -5 , 5);

TH1D *hDeltaX = new TH1D("hDeltaX", "#Delta X_{0} of the most upstream Reco Track and the Projected Primary Particle X_{0}", 200, -50 , 50);
TH1D *hDeltaY = new TH1D("hDeltaY", "#Delta Y_{0} of the most upstream Reco Track and the Projected Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hDeltaZ = new TH1D("hDeltaZ", "#Delta Z_{0} of the most upstream Reco Track and the Projected Primary Particle Z_{0}", 200, -50 , 50);

TH1D *hMCPrimaryPx = new TH1D("hMCPrimaryPx", "Primary Particle P_{x}", 300, -150 , 150);
TH1D *hMCPrimaryPy = new TH1D("hMCPrimaryPy", "Primary Particle P_{y}", 300, -150 , 150);
TH1D *hMCPrimaryPz = new TH1D("hMCPrimaryPz", "Primary Particle P_{z}", 3000, -500 , 2000);
TH1D *hMCPrimaryP  = new TH1D("hMCPrimaryP", "Primary Particle P", 3000, -500 , 2000);

TH1D *hTrueLength = new TH1D("hTrueLength", "#True Length of the Primary Particle inside the TPC", 200, 0 , 100);

TH1D *hDataUpstreamZPos = new TH1D("hDataUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);

TH1D *hAlpha = new TH1D("hAlpha", "#alpha between MC Particle and TPC Track", 90, 0, 90);

// =======================================================================================
// ====================================  FILES  ==========================================
// =======================================================================================
TFile myRootFile("piMinusAna.root", "RECREATE");



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Main loop
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void myana::Loop(int inDebug)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    // ###########################################
    // ### If debug, only look at a sub sample ###
    // ###########################################
    if (inDebug == 1 && jentry%500 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

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
    std::vector<TVector3> g4PrimaryPos0;                   // start pos
    std::vector<TVector3> g4PrimaryPosf;                   // final pos 
    std::vector<TVector3> g4PrimaryMom0;                   // momentum
    std::vector<TVector3> g4PrimaryProjPos0;               // projected position
    std::vector<std::vector<TVector3>> g4PrimaryTrTrjPos;  // trajectory sp
    std::vector<std::vector<TVector3>> g4PrimaryTrTrjMom;  // trajectory momentum 
    
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
      g4PrimaryTrTrjPos.push_back(std::vector<TVector3>);
      g4PrimaryTrTrjMom.push_back(std::vector<TVector3>);

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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    


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
      hDataUpstreamZPos->Fill(tempMinZ);
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
      auto deltaVec = frontFaceTrkPos[iFrFaTrk] - g4PrimaryProjPos0[0];

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
      auto tpcP0Hat = frontFaceTrkP0Hat[iFrFaTrk];
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
      hAlpha->Fill(alpha);

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
    





    // #####################################################
    // ### 
    // #####################################################
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
    bool   isTrackInteracting          = false;
    bool   isTrackInelastic            = false;
    for(size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // only look at the one that passed the WC_TPC match and cuts
      if(iTrk != xsRecoTrkId) continue;

      // we need the calorimetry information for this track
      size_t limit = sizeof(trkxyz[xsREcoTrkId][1])/sizeof(*trkxyz[xsREcoTrkId][1]);
      for (size_t k = 0; k < limit; k++)
      {
        // make sure we can actually see this
        TVector3 theXYZ(trkxyz[xsRecoTrkId][1][k][0], trkxyz[xsRecoTrkId][1][k][1], trkxyz[xsRecoTrkId][1][k][2]); 
        if ( FV_X_BOUND[0] > theXYZ.X() || theXYZ.X() > FV_X_BOUND[1] ||
             FV_Y_BOUND[0] > theXYZ.Y() || theXYZ.Y() > FV_Y_BOUND[1] ||
             FV_Z_BOUND[0] > theXYZ.Z() || theXYZ.Z() > FV_Z_BOUND[1] ) continue;

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
      if ( FV_X_BOUND[0] > theXYZ.X() || theXYZ.X() > FV_X_BOUND[1] ||
           FV_Y_BOUND[0] > theXYZ.Y() || theXYZ.Y() > FV_Y_BOUND[1] ||
           FV_Z_BOUND[0] > theXYZ.Z() || theXYZ.Z() > FV_Z_BOUND[1] ) isTrackIntercting = true;
    }

    // determine if the interaction was inelastic
    if (isTrackInteracting)
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
        if (iTrk == xsRecoTrkId) continue;

        // what determines the start and ending point here??
        TVector dXYZstart( trkvtxx[iTrk], trkvtxy[iTrk], trkvtxz[iTrk] );
        double diff = (theXYZ - dXYZstart).Mag();
        if ( diff < 0.01 ) nTracksLeaving++;
      }
      // there must be at least 2 visible tracks leaving this vertex
      if (nTracksLeaving >= 2) isTrackInelastic = true;
    }

    // make sure we have stuff here to work with
    if (pitchVec.size() == 0) continue;

    // get the kinetic energy at "WC 4" 
    // (for single particle gun just use KE at z=0)
    TVector3 theWCMom(0,0,0);
    for (size_t iPt = 0; iPt < maxTrTrjPoints; iPt++)
    {
      if (g4PrimaryTrTrjZ[0][iPt] > 0) break;
      theWCMom = TVector3( g4PrimaryTrTrjPx[0][iPt], g4PrimaryTrTrjPy[0][iPt], g4PrimaryTrTrjPz[0][iPt]);
    }
    double theWCKE = std::sqrt( theWCMag()*theWCMom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
    incidentKEVec.push_back(theWCKE);

    // fill the energy depositions
    double totalEnDep(0);
    for (size_t iDep = 0; iDep < eDepVec.size(); iDep++)
    {
      totalEnDep = totalEnDep + eDepVec[iDep];
      incidentKEVec.push_back(theWCKE - totalEnDep);
    }

    // fill the Incident and interacting histograms
    for (auto iKE : incidentKEVec) hRecoMCIncidentKE->Fill(iKE);
    if (isInteracting && incideKEVec.size() != 0) hRecoMCInteractingKE->Fill(incidentKEVec.back());
  }//<---End loop over entries





  // ==========================================================
  // =================  EVENT REDUCTION TABLE =================
  // ==========================================================
  std::cout << endl
            << "Events simulated:                           " << nTotalEvents                 << endl
            << "Tracks entered TPC:                         " << nPrimariesEntered            << endl
            << "Events with > 0 entered tracks:             " << nGoodMCEvents                << endl
            << "Events with front TPC track:                " << nEventsFrontTpcTrk           << endl
            << "Events with < N upper TPC tracks:           " << nEventsUpperTpcTrkCount      << endl
            << "Events with WC_TPC match (nocut):           " << nEventsDeltaMatch            << endl
            << "Events with unique WC_TPC match (nocut):    " << nEventsWcTpcUniqueMatch      << endl
            << "Events with unique WC TPC match (alphacut): " << nEventsWcTpcUniqueMatchAlpha << endl
            << "Inelastic interactions:                     " << nEventsInelastic             << endl
            << endl;

  // ===============================================
  // =================  MAKE PLOTS =================
  // ===============================================
  MakePlots();

  gApplication->Terminate(0);
}//<---- End main loop





// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Make plots
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MakePlots()
{
  myRootFile.cd();

  hMCELossUpstream->Write();

  hMCPrimaryMissedTpcX->Write();  
  hMCPrimaryMissedTpcY->Write();  
  hMCPrimaryMissedTpcZ->Write();  

  hMCPrimaryProjX0->Write();
  hMCPrimaryProjY0->Write();
  hMCPrimaryProjZ0->Write();

  hDeltaX->Write();
  hDeltaY->Write();
  hDeltaZ->Write();

  hMCPrimaryPx->Write();
  hMCPrimaryPy->Write();
  hMCPrimaryPz->Write();
  hMCPrimaryP->Write();

  hTrueLength->Write();

  hDataUpstreamZPos->Write();

  hAlpha->Write();

  hXsG4IncidentKinEn->Write();
  hXsG4InteractingKinEn->Write();

  hXsG4->Write();

  hStepSize->Write();

  myRootFile.Close();
}




// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% compute angles
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ComputeAngles(float& phi, float& theta, const TVector3& p0Hat)
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




// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Convert process to string
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::string ConvertProcessToString(const int& i)
{
  std::string theProcess;

  switch(i)
  {
    case 0:  { theProcess = "None"; break; }
    case 1:  { theProcess = "Pi-Inelastic"; break; }
    case 2:  { theProcess = "Pi+Inelastic"; break; }
    case 3:  { theProcess = "NeutronInelastic"; break; }
    case 4:  { theProcess = "HadElastic"; break; }
    case 5:  { theProcess = "Ncapture"; break; }
    case 6:  { theProcess = "ChipsNuclearCaptureAtRest"; break; }
    case 7:  { theProcess = "Decay"; break; }
    case 8:  { theProcess = "Kaon0Linelastic"; break; }
    case 9:  { theProcess = "CoulombScat"; break; }
    case 10:  { theProcess = "MuMinusCaptureAtRest"; break; }
    case 11: { theProcess = "ProtonInelastic"; break; }
    case 12: { theProcess = "Kaon+Inelastic"; break; }
    case 13: { theProcess = "Kaon-Inelastic"; break; }
    case 14: { theProcess = "ProtonInelastic"; break; }
    default: { theProcess = "Unknown"; }

  }

  return theProcess;
}
