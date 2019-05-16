/////////////////////////////////////////////////////////////////
//  QUESTIONS
//
//  * Shouldn't we consider all species of particles that enter
//    the TPC instead of just the primary? Presumably this would 
//    give us a more realistic sample.
//  * Pion.C requires all primaries to enter the TPC
//
/////////////////////////////////////////////////////////////////


#define pi_minus_100000e_cxx
#include "pi_minus_100000e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>

// ============================================================================================
// ====================================  PROTOTYPES  ==========================================
// ============================================================================================
void MakePlots();



// ==================================================================================================
// ====================================  USEFUL VARIABLES  ==========================================
// ==================================================================================================
// Counters
size_t nTotalEvents(0),      nPrimariesEntered(0),     nGoodMCEvents(0),
       nEventsInelastic(0),  nEventsFrontTpcTrk(0), nEventsUpperTpcTrkCount(0),
       nEventsDeltaMatch(0), nEventsWcTpcUniqueMatch(0);



// ====================================================================================================
// ====================================  CUTS AND CONSTANTS  ==========================================
// ====================================================================================================
// tpc dim 
TVector3 TPC_DIM_VEC( 47.0, 40.0, 90.0 );

// fiducial volume definition
float FV_X_BOUND[2] = [  2.0,                        TPC_DIM_VEC.X()-2.0 ];
float FV_Y_BOUND[2] = [ -0.5*TPC_DIM_VEC(1)+2.0, 0.5*TPC_DIM_VEC.Y()-2.0 ];
float FV_Z_BOUND[2] = [  0.0,                        TPC_DIM_VEC.Z()-2.0 ];

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
float ENTRY_TPC_ENERGY_LOSS = 36; //MeV

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
bool USE_EVENT_WEIGHT = true;

// constants for cross section calculation
float RHO           = 1400;                           //kg/m^3
float MOLAR_MASS     = 39.9;                           //g/mol
float G_PER_KG        = 1000;
float AVOGADRO      = 6.02e+23;                       //number/mol
float NUMBER_DENSITY = RHO*G_PER_KG/MOLAR_MASS*AVOGADRO;
float SLAB_WIDTH     = 0.0045;                         //in m




// ============================================================================================
// ====================================  HISTOGRAMS  ==========================================
// ============================================================================================
TH1F* hMCELossUpstream = new TH1F("hMCELossUpstream", "MC Energy Loss Upstream", 20, 0, 100);

TH1F* hMCPrimaryMissedTpcX = new TH1F("hMCPrimaryMissedTpcX", "MC Primary Missed TPC X", 200, -50, 50);
TH1F* hMCPrimaryMissedTpcY = new TH1F("hMCPrimaryMissedTpcY", "MC Primary Missed TPC Y", 200, -50, 50);
TH1F* hMCPrimaryMissedTpcZ = new TH1F("hMCPrimaryMissedTpcZ", "MC Primary Missed TPC Z", 200, -110, 10);

TH1D *hMCPrimaryProjX0 = new TH1D("hMCPrimaryProjX0", "Primary Particle X_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjY0 = new TH1D("hMCPrimaryProjY0", "Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjZ0 = new TH1D("hMCPrimaryProjZ0", "Primary Particle Z_{0}", 100, -5 , 5);

TH1D *hDeltaX = new TH1D("hDeltaX", "#Delta X_{0} of the most upstream Reco Track and the Projected Primary Particle X_{0}", 200, -50 , 50);
TH1D *hDeltaY = new TH1D("hDeltaY", "#Delta Y_{0} of the most upstream Reco Track and the Projected Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hDeltaZ = new TH1D("hDeltaZ", "#Delta Z_{0} of the most upstream Reco Track and the Projected Primary Particle Z_{0}", 200, -50 , 50);




// =======================================================================================
// ====================================  FILES  ==========================================
// =======================================================================================
TFile myRootFile("piMinusAna.root", "RECREATE");



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%% Main loop
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pi_minus_100000e::Loop(int inDebug)
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
    if (nTotalEvents%5000 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl;



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
    int nTrTrjPoints[nPrimary];                                                                           // n truth trajectory points
    int g4PrimaryTrkId[nPrimary];                                                                         // track id 
    int g4PrimaryProcess[nPrimary];                                                                       // interaction process
    float g4PrimaryX0[nPrimary], g4PrimaryY0[nPrimary], g4PrimaryZ0[nPrimary];                            // start pos
    float g4PrimaryXf[nPrimary], g4PrimaryYf[nPrimary], g4PrimaryZf[nPrimary];                            // final pos 
    float g4PrimaryPx[nPrimary], g4PrimaryPy[nPrimary], g4PrimaryPz[nPrimary];                            // momentum
    float g4PrimaryTrTrjX[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjY[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjZ[nPrimary][maxTrTrjPoints];    // trajectory sp
    float g4PrimaryTrTrjPx[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjPy[nPrimary][maxTrTrjPoints], g4PrimaryTrTrjPz[nPrimary][maxTrTrjPoints]; // trajectory momentum 
    float g4PrimaryProjX0[nPrimary], g4PrimaryProjY0[nPrimary], g4PrimaryProjZ0[nPrimary];  // projected position

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
        if(g4PrimaryPz[iPrim] >1100)                                  EVENT_WEIGHT = 0.010;
      }

      // fill momentum histos
      hMCPrimaryPx->Fill(g4PrimaryPx[iPrim], EVENT_WEIGHT);
      hMCPrimaryPy->Fill(g4PrimaryPy[iPrim], EVENT_WEIGHT);
      hMCPrimaryPz->Fill(g4PrimaryPz[iPrim], EVENT_WEIGHT);

      // project onto front of tpc
      // v = v0 + pHat*t
      // t = -1*(z0/pHat)
      TVector3 thisPos0Vec( g4PrimaryX0[iPrim], g4PrimaryY0[iPrim], g4PrimaryZ0[iPrim] );
      TVector3 thisPos1Vec( g4PrimaryXf[iPrim], g4PrimaryYf[iPrim], g4PrimaryZf[iPrim] );
      float g4TrueLength = (thisPos1Vec - thisPos0Vec).Mag();
      hTrueLength->Fill(g4TrueLength);

      TVector3 thisMomUnitVec( g4PrimaryPx[iPrim], g4PrimaryPy[iPrim], g4PrimaryPz[iPrim] );
      TVector3 thisPosProjVec = thisPos0Vec;
      
      thisMomUnitVec = thisMomUnitVec.Unit();
      float t        = -1*( g4PrimaryZ0[iPrim]/thisMomUnitVec.Z() );
      thisPosProjVec = thisPos0Vec + t*thisMomUnitVec;

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

// End looking at only MC truth/G4 information
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





    // ===================================================================
    // =================  APPLY FRONT FACE TPC TRACK CUT =================
    // ===================================================================
    bool isFrontTpcTrk(false);   
    std::vector<double> frontFaceTrkX,      frontFaceTrkY,      frontFaceTrkZ;
    std::vector<double> frontFaceTrkP0HatX, frontFaceTrkP0HatY, frontFaceTrkP0HatZ;
    std::vector<size_t> frontFaceTrkId;

    // ##########################################################
    // ### Only keeping events if there is a track within the ###
    // ###            the first x cm of the TPC               ###
    // ##########################################################
    for (size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // ###################################
      // ### Loop over trajectory points ###
      // ###################################
      // looking for most upstream point still in FV of TPC
      // most upstream point in FV
      float tempZMin(FV_Z_BOUND[1]);
      for (size_t iTrjPoint = 0; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
      {
        // ######################
        // ### Fiducial check ###
        // ######################
        TVector3 theTrjPointVec( trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint] );
        TVector3 theTrjP0HatVec( pHat0_X[iTrk][iTrjPoint], pHat0_Y[iTrk][iTrjPoint], pHat0_Z[iTrk][iTrjPoint] );

        if ( FV_X_BOUND[0] < theTrjPointVec.X() && theTrjPointVec.X() < FV_X_BOUND[1] &&
             FV_Y_BOUND[0] < theTrjPointVec.Y() && theTrjPointVec.Y() < FV_Y_BOUND[1] &&
             FV_Z_BOUND[0] < theTrjPointVec.Z() && theTrjPointVec.Z() < tempZMin        )
        {
          // this point is in the fiducial volume
          tempZMin = theTrjPointVec.Z();
          // we're storing all tracks that are in the front of the TPC
          if (tempZMin < FIRST_SP_Z_POS)
          {
            isFrontTpcTrk = true;
            
            // Store these points for later
            frontFaceTrkX.push_back( theTrjPointVec.X() );
            frontFaceTrkY.push_back( theTrjPointVec.Y() );
            frontFaceTrkZ.push_back( theTrjPointVec.Z() );

            frontFaceTrkP0HatX.push_back( theTrjP0HatVec.X() );
            frontFaceTrkP0HatX.push_back( theTrjP0HatVec.Y() );
            frontFaceTrkP0HatX.push_back( theTrjP0HatVec.Z() );

            frontFaceTrkId.push_back(iTrk);
          }//<--- End if within first few cm of TPC
        }//<--- End if in FV
      }//<--- End loop over trj points 
    }//<-- End loop over reconstructed tracks

    // skip events that do not pass this cut
    if (!isFrontTpcTrk) continue;
    nEventsFrontTpcTrk++;






    // =================================================================================
    // =================  APPLY CUT FOR THE NUMBER OF UPSTREAM TRACKS  =================
    // =================================================================================
    size_t nUpperTpcTrks(0);

    // #################################################################
    // ### Only keeping events if there is less than N tracks in the ###
    // ###    first y cm of the TPC (to help cut out EM Showers)     ###
    // #################################################################
    for (size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // ###################################
      // ### Loop over trajectory points ###
      // ###################################
      bool isUpperTpcTrk(false);
      
      // looking for most upstream point still in FV of TPC
      // most upstream point in FV
      float tempZMin(fvBoundVec[2](1));
      for (size_t iTrjPoint = 0; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
      {
        // ######################
        // ### Fiducial check ###
        // ######################
        TVector3 theTrjPointVec(trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint]);
        if ( fvBoundVec[0](0) < theTrjPointVec(0) && theTrjPointVec(0) < fvBoundVec[0](1) &&
             fvBoundVec[1](0) < theTrjPointVec(1) && theTrjPointVec(1) < fvBoundVec[1](1) &&
             fvBoundVec[2](0) < theTrjPointVec(2) && theTrjPointVec(2) < tempZMin        )
        {
          // this point is in the fiducial volume
          tempZMin = theTrjPointVec(2);
          if (tempZMin < UPPER_PART_OF_TPC) isUpperTpcTrk = true;
        }
      }
      if (isUpperTpcTrk) nUpperTpcTrks++; 
    }

    // Skipping the event if there are too many 
    // front TPC tracks in the event          
    if(nUpperTpcTrks > N_UPPER_TPC_TRACKS || nUpperTpcTrks == 0){continue;}
    nEventsUpperTpcTrkCount++;






   // ========================================================================
   // =================  MATCHING MC TO RECONSTRUCTED TRACK  =================
   // ========================================================================
  if ( nPrimary > 1 ) std::cout << "\n\n\nTHERE IS MORE THAN ONE PRIMARY\n\n\n";
  
  bool mcTpcMatch(false);
  size_t nMcTpcMatch(0);
  // ###################################
  // ### Loop over all the US Tracks ###
  // ###################################
  for(size_t iUpTrk = 0; iUpTrk < frontFaceTrkX.size(); iUpTrk++)
  {
    float deltaX = frontFaceTrkX[iUpTrk] - g4PrimaryProjX0[0];
    float deltaY = frontFaceTrkY[iUpTrk] - g4PrimaryProjY0[0];
    float deltaZ = frontFaceTrkZ[iUpTrk] - g4PrimaryProjZ0[0];

    hDeltaX->Fill(deltaX);
    hDeltaY->Fill(deltaY);
    hDeltaZ->Fill(deltaZ);

    // matching in delta X and delta Y 
    if( DELTA_X_BOUND[0] < deltaX && deltaX < DELTA_X_BOUND[1] && 
        DELTA_Y_BOUND[0] < deltaY && deltaY < DELTA_Y_BOUND[1] )
    {
      mcTpcMatch = true;
      nMcTpcMatch++;
    }
  }//<---End bb loop
  if (mcTpcMatch) nEventsDeltaMatch++;
  // force there to be one match
  if (nMcTpcMatch != 1) continue;
  nEventsWcTpcUniqueMatch++;

  // ######################################################
  // ### Calculating the angles for the upstream tracks ###
  // ######################################################
  /*TVector3 mcZHat(0,0,1);
  TVector3 mcP0Hat( g4PrimaryPx[0], g4PrimaryPy[0], g4PrimaryPz[0]);
  float    mcPhi(0), mcTheta(0), PI(3.141592654);

  // calculate theta 
  mcTheta = std::acos( mcZHat.Dot(mcP0Hat)/mcP0Hat.Mag() );

  // calculate phi
       if ( mcP0Hat.Y() > 0 && mcP0Hat.X() > 0 )   mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X());             
  else if ( mcP0Hat.Y() > 0 && mcP0Hat.X() < 0 )   mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+PI; 
  else if ( mcP0Hat.Y() < 0 && mcP0Hat.X() < 0 )   mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+PI; 
  else if ( mcP0Hat.Y() < 0 && mcP0Hat.X() > 0 )   mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+2*PI;     
  else if ( mcP0Hat.Y() == 0 && mcP0Hat.X() == 0 ) mcPhi = 0; //defined by convention
  else if ( mcP0Hat.Y() == 0 )
  {
    if ( p_hat_0_MC.X() > 0 ) mcPhi = 0; 
    else                      mcPhi = PI; 
  }
  else if ( p_hat_0_MC.X() == 0 )
  {
    if ( p_hat_0_MC.Y() > 0 ) mcPhi = PI/2; 
    else                      mcPhi = 3*PI/2; 
  }
  

  // ############################################################
  // ### Calculating the angles for the upstream tracks (TPC) ###
  // ############################################################
  TVector3 z_hat_TPC(0,0,1);
  TVector3 p_hat_0_TPC;
  for(int aa = 0; aa < nUpStreamTrk; aa++)
      {
      // ### Setting the TVector ###
      p_hat_0_TPC.SetX(dummyTrk_pHat0X[aa]);
      p_hat_0_TPC.SetY(dummyTrk_pHat0Y[aa]);
      p_hat_0_TPC.SetZ(dummyTrk_pHat0Z[aa]);


      // ### Calculating TPC track theta ###
      dummyTrk_Theta[aa] = acos(z_hat_TPC.Dot(p_hat_0_TPC)/p_hat_0_TPC.Mag());

      // ### Calculating TPC track phi ###
      //---------------------------------------------------------------------------------------------------------------------
      if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X()); }
      else if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+6.28318; }
      else if( p_hat_0_TPC.Y() == 0 && p_hat_0_TPC.X() == 0 ){ dummyTrk_Phi[aa] = 0; }//defined by convention
      else if( p_hat_0_TPC.Y() == 0 )
         {
         if( p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = 0; }

         else{ dummyTrk_Phi[aa] = 3.141592654; }

         }
      else if( p_hat_0_TPC.X() == 0 )
         {
         if( p_hat_0_TPC.Y() > 0 ){ dummyTrk_Phi[aa] = 3.141592654/2; }
         else{ dummyTrk_Phi[aa] = 3.141592654*3/2; }

         }
      //---------------------------------------------------------------------------------------------------------------------

      }//<---End aa loop

*/




    // =================================================================
    // ================= ADDING TO CROSS SECTION PLOTS =================
    // =================================================================
   
    // calculating the momentum from the MC primary
    TVector3 thisMomVec( g4PrimaryPx[0], g4PrimaryPy[0], g4PrimaryPz[0] );
    float momentum = thisMomVec.Mag();

    // calculating the kinetic energy 
    // KE = E - mass 
    float kineticEnergy = std::pow( momentum*momentum + PARTICLE_MASS*PARTICLE_MASS, 0.5 ) - PARTICLE_MASS;

    // subtract of the assumed energy loss
    kineticEnergy = kineticEnergy - ENTRY_TPC_ENERGY_LOSS;

    // fill the KE histo
    hMCInitialKE->Fill(kineticEnergy);

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

      // ### Recording the start-point of this track ###

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
          xsTrackEndY < -19 || xsTrackEndZ > 89.0) ExitingTrack = true;

      iPionSp = 0;
      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
      for(size_t iSp = 0; iSp < ntrkhits[iTrk]; iSp++)
      {
        // Note: Format for this variable is:           
        //         [trk number][plane 0 = induction, 1 = collection][spts number] 
        xsPion_dEdX[iPionSp] = trkdedx[iTrk][1][iSp];  

        // putting in a fix in the case that the dE/dX is negative in this step 
        // then take the point before and the point after and average them
        if( xsPion_dEdX[iPionSp] < 0 && iSp < ntrkhits[iTrk] && iSp > 0)
            {xsPion_dEdX[iPionSp] = ( (trkdedx[iTrk][1][iSp - 1] + trkdedx[iTrk][1][iSp + 1]) / 2);}

        // if this didn't fix it, then just put in a flat 2.4 MeV / cm fix 
        if(xsPion_dEdX[iPionSp] < 0)
        {
          xsPion_dEdX[iPionSp] = 2.4;
          continue;
        }

        xsPionResRange[iPionSp] = trkrr[iTrk][1][iSp];
        xsPionPitchHit[iPionSp] = trkpitchhit[iTrk][1][iSp];
        xsPionSumEnergy         = (xsPion_dEdX[iPionSp] * xsPionPitchHit[iPionSp]) + xsPionSumEnergy;

        // recording the dE/dX
        hDataPiondEdX->Fill(xsPion_dEdX[iPionSp], EVENT_WEIGHT);
        // recording the residual range 
        hDataPionRR->Fill(xsPionResRange[iPionSp]);
        // recording the Pitch 
        hDataPionTrkPitch->Fill(xsPionPitchHit[iPionSp], EVENT_WEIGHT);
        // filling 2d dE/dX vs RR 
        hDataPiondEdXvsRR->Fill(xsPionResRange[iPionSp], xsPion_dEdX[iPionSp]);

        iPionSp++;
      }//<--- End iSp loop
    }//<--- End iTrk loop

  }//<---End loop over entries





// ==========================================================
// =================  EVENT REDUCTION TABLE =================
// ==========================================================
  std::cout << endl
            << "Events simulated:                        " << nTotalEvents        << endl
            << "Tracks entered TPC:                      " << nPrimariesEntered      << endl
	    << "Good MC events:                          " << nGoodMCEvents       << endl
            << "Events with front TPC track:             " << nEventsFrontTpcTrk  << endl
            << "Events with < N upper TPC tracks:        " << nEventsUpperTpcTrkCount << endl
            << "Events with WC_TPC match (nocut):        " << nEventsDeltaMatch   << endl
            << "Events with unique WC_TPC match (nocut): " << nEventsWcTpcUniqueMatch << endl           
            << "Inelastic interactions:                  " << nEventsInelastic    << endl
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

  myRootFile.Close();
}
