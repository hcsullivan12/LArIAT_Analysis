/////////////////////////////////////////////////////////////////
//  QUESTIONS
//
//  * Shouldn't we consider all species of particles that enter
//    the TPC instead of just the primary? Presumably this would 
//    give us a more realistic sample.
//  * Pion.C requires all primaries to enter the TPC
//
/////////////////////////////////////////////////////////////////


#define pi_minus_50000e_g4step_cxx
#include "pi_minus_50000e_g4step.h"
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
// option to compute xs using g4 info
bool IS_XS_USING_G4_INFO(true);

// option to compute xs using reco info
bool IS_XS_USING_RECO_INFO(false);

// the interaction channel we're interested in
std::string INTERACTION_CHANNEL("Pi-Inelastic");

// step size
float STEP_SIZE(0.0003); // m

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

TH1D *hXsG4IncidentKinEn = new TH1D("hXsG4IncidentKinEn", "MC Truth Incident Kinetic Energy", 20, 0, 1000);
TH1D *hXsG4InteractingKinEn = new TH1D("hXsG4InteractingKinEn", "MC Truth Interacting Kinetic Energy", 20, 0, 1000);

TH1D *hXsG4 = new TH1D("hXsG4", "XS", 20, 0, 1000);

TH1D *hStepSize = new TH1D("hStepSize", "G4 Step Size", 1000, 0, 0.0004);

// =======================================================================================
// ====================================  FILES  ==========================================
// =======================================================================================
TFile myRootFile("piMinusAna.root", "RECREATE");



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Main loop
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pi_minus_50000e_g4step::Loop(int inDebug)
{

  // ########################
  // ### Define our slabs ###
  // ########################
  if (IS_XS_USING_G4_INFO)
  {
    float iZ(TPC_Z_BOUND[0]);
    while (iZ <= TPC_Z_BOUND[1])
    {
      slabsInZ.push_back(iZ);
      iZ = iZ + STEP_SIZE*100; // convert to cm
    }
  }
 

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
    if (nTotalEvents%500 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl; 


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Looking at only MC truth/G4 information
 
    // =====================================================================================
    // ==============================  VARIABLES FOR G4 INFO  ==============================
    // =====================================================================================
    // count the number of primaries
    // this is to save us from a memory crash
    // which is not an easy bug to hunt down
    int nPrimary(0);                                                                                      // number of primaries
    int maxTrTrjPoints(0);                                                                                // max number of truth traj points
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (process_primary[iG4] == 1) nPrimary++;
      if (process_primary[iG4] == 1 && maxTrTrjPoints < NTrTrajPts[iG4]) maxTrTrjPoints = NTrTrajPts[iG4];
    }
    int nTrTrjPoints[nPrimary];                                                               // n truth trajectory points
    int g4PrimaryTrkId[nPrimary];                                                             // track id 
    int g4PrimaryProcess[nPrimary];                                                           // interaction process
    float g4PrimaryX0[nPrimary], g4PrimaryY0[nPrimary], g4PrimaryZ0[nPrimary];                // start pos
    float g4PrimaryXf[nPrimary], g4PrimaryYf[nPrimary], g4PrimaryZf[nPrimary];                // final pos 
    float g4PrimaryPx[nPrimary], g4PrimaryPy[nPrimary], g4PrimaryPz[nPrimary];                // momentum
    float g4PrimaryTrTrjX[nPrimary][maxTrTrjPoints],  g4PrimaryTrTrjY[nPrimary][maxTrTrjPoints],  g4PrimaryTrTrjZ[nPrimary][maxTrTrjPoints]; // trajectory sp
    float g4PrimaryTrTrjPx[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjPy[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjPz[nPrimary][maxTrTrjPoints]; // trajectory momentum 
    float g4PrimaryProjX0[nPrimary], g4PrimaryProjY0[nPrimary], g4PrimaryProjZ0[nPrimary];    // projected position

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
       
      // store the position and momentum 
      g4PrimaryX0[iPrim] = StartPointx[iG4];
      g4PrimaryY0[iPrim] = StartPointy[iG4];
      g4PrimaryZ0[iPrim] = StartPointz[iG4];

      g4PrimaryXf[iPrim] = EndPointx[iG4];
      g4PrimaryYf[iPrim] = EndPointy[iG4];
      g4PrimaryZf[iPrim] = EndPointz[iG4];

      g4PrimaryPx[iPrim] = Px[iG4] * 1000; // convert to MeV
      g4PrimaryPy[iPrim] = Py[iG4] * 1000; // convert to MeV
      g4PrimaryPz[iPrim] = Pz[iG4] * 1000; // convert to MeV

      // setting a global event weight 
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

      // setting Event weight 
      if(USE_EVENT_WEIGHT)
      {
        if(g4PrimaryPz[iPrim] > 0   && g4PrimaryPz[iPrim] < 100) EVENT_WEIGHT = 0.010;
        if(g4PrimaryPz[iPrim] > 100 && g4PrimaryPz[iPrim] < 200) EVENT_WEIGHT = 0.020;
        if(g4PrimaryPz[iPrim] > 200 && g4PrimaryPz[iPrim] < 300) EVENT_WEIGHT = 0.100;
        if(g4PrimaryPz[iPrim] > 300 && g4PrimaryPz[iPrim] < 400) EVENT_WEIGHT = 0.535;
        if(g4PrimaryPz[iPrim] > 400 && g4PrimaryPz[iPrim] < 500) EVENT_WEIGHT = 0.840;
        if(g4PrimaryPz[iPrim] > 500 && g4PrimaryPz[iPrim] < 600) EVENT_WEIGHT = 0.965;
        if(g4PrimaryPz[iPrim] > 600 && g4PrimaryPz[iPrim] < 700) EVENT_WEIGHT = 1.000;
        if(g4PrimaryPz[iPrim] > 700 && g4PrimaryPz[iPrim] < 800) EVENT_WEIGHT = 0.620;
        if(g4PrimaryPz[iPrim] > 800 && g4PrimaryPz[iPrim] < 900) EVENT_WEIGHT = 0.225;
        if(g4PrimaryPz[iPrim] > 900 && g4PrimaryPz[iPrim] <1000) EVENT_WEIGHT = 0.094;
        if(g4PrimaryPz[iPrim] >1000 && g4PrimaryPz[iPrim] <1100) EVENT_WEIGHT = 0.0275;
        if(g4PrimaryPz[iPrim] >1100)                             EVENT_WEIGHT = 0.010;
      }

      // fill momentum histos
      TVector3 thisMomVec( g4PrimaryPx[iPrim], g4PrimaryPy[iPrim], g4PrimaryPz[iPrim] );
      hMCPrimaryPx->Fill(thisMomVec.X(), EVENT_WEIGHT);
      hMCPrimaryPy->Fill(thisMomVec.Y(), EVENT_WEIGHT);
      hMCPrimaryPz->Fill(thisMomVec.Z(), EVENT_WEIGHT);
      hMCPrimaryP->Fill(thisMomVec.Mag(), EVENT_WEIGHT);

      // fill true length histo
      TVector3 thisPos0Vec( g4PrimaryX0[iPrim], g4PrimaryY0[iPrim], g4PrimaryZ0[iPrim] );
      TVector3 thisPos1Vec( g4PrimaryXf[iPrim], g4PrimaryYf[iPrim], g4PrimaryZf[iPrim] );

      float g4TrueLength = (thisPos1Vec - thisPos0Vec).Mag();
      hTrueLength->Fill(g4TrueLength);

      // project onto tpc 
      TVector3 thisPosProjVec = thisPos0Vec - ( thisPos0Vec.Z()/thisMomVec.Z() )*thisMomVec;

      g4PrimaryProjX0[iPrim] = thisPosProjVec.X();
      g4PrimaryProjY0[iPrim] = thisPosProjVec.Y();
      g4PrimaryProjZ0[iPrim] = thisPosProjVec.Z();

      // fill the proj histos
      hMCPrimaryProjX0->Fill(g4PrimaryProjX0[iPrim]);
      hMCPrimaryProjY0->Fill(g4PrimaryProjY0[iPrim]);
      hMCPrimaryProjZ0->Fill(g4PrimaryProjZ0[iPrim]);

      // store the track id 
      g4PrimaryTrkId[iPrim] = TrackId[iG4];

      // store the trajectory points and momentum
      nTrTrjPoints[iPrim] = NTrTrajPts[iG4];
      for (size_t iPoint = 0; iPoint < nTrTrjPoints[iPrim]; iPoint++)
      {
        g4PrimaryTrTrjX[iPrim][iPoint] = MidPosX[iG4][iPoint];
        g4PrimaryTrTrjY[iPrim][iPoint] = MidPosY[iG4][iPoint];
        g4PrimaryTrTrjZ[iPrim][iPoint] = MidPosZ[iG4][iPoint];

        g4PrimaryTrTrjPx[iPrim][iPoint] = MidPx[iG4][iPoint] * 1000; // convert to MeV
        g4PrimaryTrTrjPy[iPrim][iPoint] = MidPy[iG4][iPoint] * 1000; // convert to MeV
        g4PrimaryTrTrjPz[iPrim][iPoint] = MidPz[iG4][iPoint] * 1000; // convert to MeV
      }
      iPrim++;
    } 
    }

    // ###################################
    // ### Get the interaction process ###
    // ###################################
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      for (size_t iPrim = 0; iPrim < nPrimary; iPrim++)
      {
        if ( Mother[iG4] == g4PrimaryTrkId[iPrim] ) g4PrimaryProcess[iPrim] = Process[iG4];
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
      if ( g4PrimaryZf[iPrim] > 0 ) isGoodEvent = true;
      if ( !isGoodEvent ) 
      {
        hMCPrimaryMissedTpcX->Fill(g4PrimaryXf[iPrim]);
        hMCPrimaryMissedTpcY->Fill(g4PrimaryYf[iPrim]);
        hMCPrimaryMissedTpcZ->Fill(g4PrimaryZf[iPrim]);
        continue;
      }

      // ####################################
      // ### This primary entered the TPC ###
      // ####################################
      nPrimariesEntered++;

      // calculate energy loss
      float energyLoss(0);
      // loop over trj points for this primary
      for (size_t iPoint = 0; iPoint < nTrTrjPoints[iPrim]; iPoint++)
      {
        // only look at the upstream portion
        if ( g4PrimaryTrTrjZ[iPrim][iPoint] > 0 ) continue;
        // ignore last point
        if ( (iPoint+1) >= nTrTrjPoints[iPrim] ) break;
 
        TVector3 mom1Vec = TVector3( g4PrimaryTrTrjPx[iPrim][iPoint], 
                                     g4PrimaryTrTrjPy[iPrim][iPoint], 
                                     g4PrimaryTrTrjPz[iPrim][iPoint] );
        TVector3 mom2Vec = TVector3( g4PrimaryTrTrjPx[iPrim][iPoint+1], 
                                     g4PrimaryTrTrjPy[iPrim][iPoint+1], 
                                     g4PrimaryTrTrjPz[iPrim][iPoint+1] );
        float mom1    = mom1Vec.Mag();
        float energy1 = std::sqrt( mom1*mom1 + PARTICLE_MASS*PARTICLE_MASS ); 
        float mom2 = mom2Vec.Mag();
        float energy2 = std::sqrt( mom2*mom2 + PARTICLE_MASS*PARTICLE_MASS ); 

        energyLoss += (energy1 - energy2);
      }//<--- End loop over true traj points
      hMCELossUpstream->Fill(energyLoss);
    }//<--- End loop over primaries
    if (!isGoodEvent) continue;
    nGoodMCEvents++;





    // ================================================================
    // ================= MC TRUTH LEVEL CROSS SECTION =================
    // ================================================================
    if (IS_XS_USING_G4_INFO) 
    {
      // ###########################
      // ### Loop over primaries ###
      // ###########################
      for (size_t iPrim = 0; iPrim < nPrimary; iPrim++)
      {
        // ####################################################################
        // ### Loop over tr trj points to look for interaction in this slab ###
        // ####################################################################
        size_t nInter(0);
        double slabBoundZ[2] = { FV_Z_BOUND[0], FV_Z_BOUND[0] + STEP_SIZE };
        for (size_t iPoint = 1; iPoint < nTrTrjPoints[iPrim]; iPoint++)
        {
          TVector3 theTrTrjPointVec( g4PrimaryTrTrjX[iPrim][iPoint],  g4PrimaryTrTrjY[iPrim][iPoint],  g4PrimaryTrTrjZ[iPrim][iPoint] );
          TVector3 thePreTrTrjPointVec( g4PrimaryTrTrjX[iPrim][iPoint-1],  g4PrimaryTrTrjY[iPrim][iPoint-1],  g4PrimaryTrTrjZ[iPrim][iPoint-1] );
          double diff = (thePreTrTrjPointVec-theTrTrjPointVec).Mag();
          //if (diff < 1e-5) continue;
          TVector3 theTrTrjMomVec( g4PrimaryTrTrjPx[iPrim][iPoint-1], g4PrimaryTrTrjPy[iPrim][iPoint-1], g4PrimaryTrTrjPz[iPrim][iPoint-1] );
          double kineticEnergy(0);

          // compute the kinetic energy
          kineticEnergy = std::sqrt( theTrTrjMomVec.Mag()*theTrTrjMomVec.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
          // should we do this later?
          if (kineticEnergy < 0.001) break;

          // #######################################
          // ### Check if in the fiducial volume ###
          // #######################################
          if ( FV_X_BOUND[0] < theTrTrjPointVec.X() && theTrTrjPointVec.X() < FV_X_BOUND[1] &&
               FV_Y_BOUND[0] < theTrTrjPointVec.Y() && theTrTrjPointVec.Y() < FV_Y_BOUND[1] &&
               FV_Z_BOUND[0] < theTrTrjPointVec.Z() && theTrTrjPointVec.Z() < FV_Z_BOUND[1] )
          { 
            // ignore points not in this slab
            if ( theTrTrjPointVec.Z() < slabBoundZ[0] ) continue;
            if ( theTrTrjPointVec.Z() > slabBoundZ[1] )
            {
              // we've entered a new slab
              slabBoundZ[0] = slabBoundZ[1];
              slabBoundZ[1] = slabBoundZ[1] + STEP_SIZE;

              if (nInter <= 1)
              {
                // fill the incident
                hXsG4IncidentKinEn->Fill(kineticEnergy);
              }

              if (nInter == 1)
              {
                // we've got the one we're looking for
                // fill the interacting histogram
                nXsG4Signal++;
                hXsG4InteractingKinEn->Fill(kineticEnergy);
              }

              nInter = 0;
              //iPoint--;
              continue;
            }

            // we're in this slab!
          
            // did it interact here?
            auto it = std::find( InteractionPoint->begin(), InteractionPoint->end(), iPoint );
            if ( it != InteractionPoint->end() )
            {
              // we have an interaction at this point!
              size_t d = std::distance(InteractionPoint->begin(), it);
              std::string theProcess = ConvertProcessToString( (*InteractionPointType)[d] );
              // check if it's the process we're interested in
              if (theProcess == INTERACTION_CHANNEL) nInter++;
            }
          }
        }//<--- End loop over trj points  
      }//<--- End loop over primaries
    }//<--- End if xs using g4 


// End looking at only MC truth/G4 information
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





    // quit if we don't care about reconstructed variables
    if(!IS_XS_USING_RECO_INFO) continue;




    // ============================================================================
    // =================      APPLY FRONT FACE TPC TRACK CUT      =================
    // =================  APPLY CUT ON NUMBER OF UPSTREAM TRACKS  =================
    // =================          (Cuts on EM showers)            =================
    // ============================================================================
    std::vector<double> frontFaceTrkX,      frontFaceTrkY,      frontFaceTrkZ;
    std::vector<double> frontFaceTrkP0HatX, frontFaceTrkP0HatY, frontFaceTrkP0HatZ;
    std::vector<double> frontFaceTrkPhi,    frontFaceTrkTheta;
    std::vector<size_t> frontFaceTrkId;
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
      float dummyTrkX, dummyTrkY, dummyTrkZ;
      float dummyTrkP0HatX, dummyTrkP0HatY, dummyTrkP0HatZ;
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
            dummyTrkX = theTrjPointVec.X();
            dummyTrkY = theTrjPointVec.Y();
            dummyTrkZ = theTrjPointVec.Z();

            dummyTrkP0HatX = theTrjP0HatVec.X();
            dummyTrkP0HatY = theTrjP0HatVec.Y();
            dummyTrkP0HatZ = theTrjP0HatVec.Z();

            dummyTrkId = iTrk;
          }//<--- End if within first few cm of TPC

          // check if this point is in the upstream portion
          if (tempMinZ < UPPER_PART_OF_TPC) isUpstreamTpcTrk = true;

        }//<--- End if in FV
      }//<--- End loop over trj points 

      // we should have the observables attached to the most upstream trj point for this track
      // if it is at the front, append to our vectors
      if (isFrontTpcTrk)
      {
        frontFaceTrkX.push_back( dummyTrkX );
        frontFaceTrkY.push_back( dummyTrkY );
        frontFaceTrkZ.push_back( dummyTrkZ );

        frontFaceTrkP0HatX.push_back( dummyTrkP0HatX );
        frontFaceTrkP0HatY.push_back( dummyTrkP0HatY );        
        frontFaceTrkP0HatZ.push_back( dummyTrkP0HatZ );        

        frontFaceTrkId.push_back( dummyTrkId );
      }

      // increment if in upstream portion
      if (isUpstreamTpcTrk) nUpstreamTpcTrks++;

      // fill histos
      hDataUpstreamZPos->Fill(tempMinZ);
    }//<-- End loop over reconstructed tracks

    // skip events that do not have a trk at the front
    if (frontFaceTrkId.size() == 0) continue;
    nEventsFrontTpcTrk++;

    // skip events that have too many tracks upstream        
    if(nUpstreamTpcTrks > N_UPPER_TPC_TRACKS || nUpstreamTpcTrks == 0) continue;
    nEventsUpperTpcTrkCount++;






    // ========================================================================
    // =================  MATCHING MC TO RECONSTRUCTED TRACK  =================
    // ========================================================================
    if ( nPrimary > 1 ) std::cout << "\n\n\nTHERE IS MORE THAN ONE PRIMARY\n\n\n";
  
    size_t nMcTpcMatch(0);
    // ###########################################
    // ### Loop over all the front face Tracks ###
    // ###########################################
    for(size_t iFrFaTrk = 0; iFrFaTrk < frontFaceTrkId.size(); iFrFaTrk++)
    {
      float deltaX = frontFaceTrkX[iFrFaTrk] - g4PrimaryProjX0[0];
      float deltaY = frontFaceTrkY[iFrFaTrk] - g4PrimaryProjY0[0];
      float deltaZ = frontFaceTrkZ[iFrFaTrk] - g4PrimaryProjZ0[0];

      hDeltaX->Fill(deltaX);
      hDeltaY->Fill(deltaY);
      hDeltaZ->Fill(deltaZ);

      // matching in delta X and delta Y 
      if( DELTA_X_BOUND[0] < deltaX && deltaX < DELTA_X_BOUND[1] && 
          DELTA_Y_BOUND[0] < deltaY && deltaY < DELTA_Y_BOUND[1] )
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
    TVector3 mcP0Hat( g4PrimaryPx[0], g4PrimaryPy[0], g4PrimaryPz[0]);
    mcP0Hat = mcP0Hat.Unit();
    float    mcPhi(0), mcTheta(0);
  
    // compute the angles for primary
    ComputeAngles( mcPhi, mcTheta, mcP0Hat );

    // ##############################################################
    // ### Calculating the angles for the front face tracks (TPC) ###
    // ##############################################################
    for(int iFrFaTrk = 0; iFrFaTrk < frontFaceTrkId.size(); iFrFaTrk++)
    {
      // setting the TVector 
      TVector3 tpcP0Hat( frontFaceTrkP0HatX[iFrFaTrk], 
                         frontFaceTrkP0HatY[iFrFaTrk], 
                         frontFaceTrkP0HatZ[iFrFaTrk] );
      tpcP0Hat = tpcP0Hat.Unit();

      float thisPhi(0), thisTheta(0);
      ComputeAngles( thisPhi, thisTheta, tpcP0Hat );

      frontFaceTrkPhi.push_back(thisPhi);
      frontFaceTrkTheta.push_back(thisTheta);
    }//<--- End iFrFaTrk loop

    // sanity check
    assert( frontFaceTrkX.size()      == frontFaceTrkY.size()      == frontFaceTrkZ.size()      == 
            frontFaceTrkP0HatX.size() == frontFaceTrkP0HatY.size() == frontFaceTrkP0HatZ.size() ==
            frontFaceTrkPhi.size()    == frontFaceTrkTheta.size()  == frontFaceTrkId.size() );

    





    // ======================================================
    // =================  APPLY ANGLE CUTS  =================
    // ======================================================
    bool isAlphaMatch(false);
    size_t xsRecoTrkId;
    for(int iFrFaTrk = 0; iFrFaTrk < frontFaceTrkId.size(); iFrFaTrk++)
    {
      // compute the g4 primary momentum unit vector
      TVector3 theMCUnit( std::sin(mcTheta)*std::cos(mcPhi),
                          std::sin(mcTheta)*std::sin(mcPhi),
                          std::cos(mcTheta) );

      // compute the front face tpc trk momontum unit vector
      TVector3 theTpcUnit( std::sin(frontFaceTrkTheta[iFrFaTrk])*std::cos(frontFaceTrkPhi[iFrFaTrk]),
                           std::sin(frontFaceTrkTheta[iFrFaTrk])*std::sin(frontFaceTrkPhi[iFrFaTrk]),
                           std::cos(frontFaceTrkTheta[iFrFaTrk]) );
      
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
        xsRecoTrkId = frontFaceTrkId[iFrFaTrk];
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










 /*   // =================================================================
    // ================= ADDING TO CROSS SECTION PLOTS =================
    // =================================================================
   
    // ########################################################################
    // ### Variables for the track we are calculating the cross-section for ###
    // ########################################################################
    double xsdEdX[1000]    ={0.};
    double xsResRange[1000]={0.};
    double xsPitchHit[1000]={0.};
    int    xsSpts(0);
    double xsSumEnergy(0);
    float  xsTrackEndX(999), xsTrackEndY(999), xsTrackEndZ(999);
    bool   xsExitingTrack(false);
    float  xsKineticEnergy(0);
    float  xsMomentum(0);

    // calculating the momentum from the MC primary
    TVector3 thisMomVec( g4PrimaryPx[0], g4PrimaryPy[0], g4PrimaryPz[0] );
    xsMomentum = thisMomVec.Mag();

    // calculating the kinetic energy 
    // KE = E - mass 
    xsKineticEnergy = std::sqrt( xsMomentum*xsMomentum + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
    // subtract of the assumed energy loss
    xsKineticEnergy = xsKineticEnergy - ENTRY_TPC_ENERGY_LOSS;
    // fill the KE histo
    hMCInitialKE->Fill(xsKineticEnergy);

    // ################################
    // ### Loop over all TPC Tracks ###
    // ################################
    for(size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // only look at the one that passed the WC_TPC match and cuts
      if(iTrk != xsRecoTrackIndex) continue;

      // recording the end point of this track 
      xsTrackEndX = trkendx[iTrk];
      xsTrackEndY = trkendy[iTrk];
      xsTrackEndZ = trkendz[iTrk];

      hDataPionTrackEndX->Fill(xsTrackEndX);
      hDataPionTrackEndY->Fill(xsTrackEndY);
      hDataPionTrackEndZ->Fill(xsTrackEndZ);

      // recording the start-point of this track 
      hdataPionTrackStartX->Fill(trkvtxx[iTrk]);
      hdataPionTrackStartY->Fill(trkvtxy[iTrk]);
      hdataPionTrackStartZ->Fill(trkvtxz[iTrk]);

      xsRecoLength = trklength[iTrk];
      hRecoLength->Fill(xsRecoLength);

      // #################################################
      // ### If this tracks end point is at a boundary ###
      // ###   then tag this track as "through-going"  ###
      // #################################################
      if( xsTrackEndX < 1   || xsTrackEndX > 42.0 || xsTrackEndY > 19 ||
          xsTrackEndY < -19 || xsTrackEndZ > 89.0) xsExitingTrack = true;

      iPionSp = 0;
      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
      for(size_t iSpt = 0; iSpt < ntrkhits[iTrk]; iSpt++)
      {
        // Note: Format for this variable is:           
        //         [trk number][plane 0 = induction, 1 = collection][spts number] 
        xsdEdX[iPionSp] = trkdedx[iTrk][1][iSpt];  

        // putting in a fix in the case that the dE/dX is negative in this step 
        // then take the point before and the point after and average them
        if( xsdEdX[iPionSp] < 0 && iSpt < ntrkhits[iTrk] && iSpt > 0)
        {
          xsdEdX[iPionSp] =  ( trkdedx[iTrk][1][iSpt - 1] + trkdedx[iTrk][1][iSpt + 1] )*0.5;
        }

        // if this didn't fix it, then just put in a flat 2.4 MeV / cm fix 
        if(xsdEdX[iPionSp] < 0)
        {
          xsdEdX[iPionSp] = 2.4;
          continue;
        }

        xsResRange[iPionSp] = trkrr[iTrk][1][iSpt];
        xsPitchHit[iPionSp] = trkpitchhit[iTrk][1][iSpt];
        xsSumEnergy         = (xsdEdX[iPionSp] * xsPionPitchHit[iPionSp]) + xsPionSumEnergy;

        // recording the dE/dX
        hDataPiondEdX->Fill(xsdEdX[iPionSp], EVENT_WEIGHT);
        // recording the residual range 
        hDataPionRR->Fill(xsPionResRange[iPionSp]);
        // recording the Pitch 
        hDataPionTrkPitch->Fill(xsPionPitchHit[iPionSp], EVENT_WEIGHT);
        // filling 2d dE/dX vs RR 
        hDataPiondEdXvsRR->Fill(xsPionResRange[iPionSp], xsdEdX[iPionSp]);

        iPionSp++;
      }//<--- End iSpt loop
    }//<--- End iTrk loop
*/





  }//<---End loop over entries





  // ==========================================================
  // =================  EVENT REDUCTION TABLE =================
  // ==========================================================
  if (IS_XS_USING_RECO_INFO)
  {
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
  }
  if (IS_XS_USING_G4_INFO)
  {
    std::cout << endl
              << "Events simulated:                           " << nTotalEvents                 << endl
              << "Tracks entered TPC:                         " << nPrimariesEntered            << endl
              << "Events with > 0 entered tracks:             " << nGoodMCEvents                << endl
              << "Signal events:                              " << nXsG4Signal                  << endl
              << endl;
  }



  // ====================================================
  // =================  CALCULATE G4 XS =================
  // ====================================================
  if (IS_XS_USING_G4_INFO) CalculateG4Xs();



  // ===============================================
  // =================  MAKE PLOTS =================
  // ===============================================
  MakePlots();



  gApplication->Terminate(0);
}//<---- End main loop






// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Calculate the g4 cross section
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CalculateG4Xs()
{
  for (int bin = 1; bin <= hXsG4InteractingKinEn->GetXaxis()->GetNbins(); bin++)
  {
    auto num = hXsG4InteractingKinEn->GetBinContent(bin);
    auto den = hXsG4IncidentKinEn->GetBinContent(bin);
    auto err1 = std::sqrt(num);
    auto err2 = std::sqrt(den);

    double xs    = BARN_PER_M2*(1/NUMBER_DENSITY)*(1/STEP_SIZE)*(num/den);
    double error = xs*std::sqrt( (err1/num)*(err1/num) + (err2/den)*(err2/den)  );
    cout << xs << " " << error << endl;

    hXsG4->SetBinContent(bin, xs);
    hXsG4->SetBinError(bin, error);
  }
}





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
