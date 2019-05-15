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
bool InFV(const TVector3& );



// ==================================================================================================
// ====================================  USEFUL VARIABLES  ==========================================
// ==================================================================================================
// Counters
size_t nTotalEvents(0),      nTracksEntered(0),     nGoodMCEvents(0),
       nEventsInelastic(0),  nEventsFrontTPCTrk(0), nEventsNUpperTPCTrk(0),
       nEventsDeltaMatch(0), nEventsWC_TPCUniqueMatch(0);

// tpc dim and fiducial volume definition
TVector3 tpcDimVec( 47.0, 40.0, 90.0 );
using Vec2D = std::vector<TVector3>;
Vec2D fvBoundVec = { TVector3 (                   2.0,     tpcDimVec(0)-2.0, 0 ),
                     TVector3 ( -0.5*tpcDimVec(1)+2.0, 0.5*tpcDimVec(1)-2.0, 0 ),
                     TVector3 (                   0.0,     tpcDimVec(2)-2.0, 0 ) };

// Identifiers
std::set<std::string> processes;
std::set<int>         ufos;



// ====================================================================================================
// ====================================  CUTS AND CONSTANTS  ==========================================
// ====================================================================================================
// mass of pion in MeV
float PARTICLEMASS(139.57); 

// number of centimeters in Z we require a track to have a spacepoint
float FIRSTSPZPOS(2.0);

// portion of upstream TPC which we will restrict the number of tracks
float UPPERPARTOFTPC(14.0);

// number of upper tpc tracks allowed
size_t NUPPERTPCTRACKS(4);

// Delta X Between Wire Chamber Track and TPC Track 
float DELTAXBOUND[2] = { -2.0, 6.0 };

// Delta Y Between Wire Chamber Track and TPC Track 
float DELTAYBOUND[2] = { -3.0, 6.0 };

// the assumed energy loss between the cryostat and the TPC 
float ENTRYTPCENERGYLOSS = 36; //MeV




// ============================================================================================
// ====================================  HISTOGRAMS  ==========================================
// ============================================================================================
TH1F* hMCELossUpstream = new TH1F("hMCELossUpstream", "MC Energy Loss Upstream", 20, 0, 100);

TH1F* hMCPrimaryMissedTPCX = new TH1F("hMCPrimaryMissedTPCX", "MC Primary Missed TPC X", 200, -50, 50);
TH1F* hMCPrimaryMissedTPCY = new TH1F("hMCPrimaryMissedTPCY", "MC Primary Missed TPC Y", 200, -50, 50);
TH1F* hMCPrimaryMissedTPCZ = new TH1F("hMCPrimaryMissedTPCZ", "MC Primary Missed TPC Z", 200, -110, 10);

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
    // if debug, only look at a sub sample
    if (inDebug == 1 && jentry%500 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //#########################################
    //### Increment our total event counter ###
    //#########################################
    nTotalEvents++; 
    if (nTotalEvents%5000 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl;

// =====================================================================================
// ==============================  VARIABLES FOR G4 INFO  ==============================
// =====================================================================================
    int nPrimary(0);                                                                                      // number of primaries
    // count the number of primaries
    // this is to save us from a memory crash
    // which is not an easy bug to hunt down
    int maxTrjPoints(0);
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (process_primary[iG4] == 1) nPrimary++;
      if (process_primary[iG4] == 1 && maxTrjPoints < NTrTrajPts[iG4]) maxTrjPoints = NTrTrajPts[iG4];
    }
    int nTrjPoints[nPrimary];                                                                             // n trajectory sp
    int g4PrimaryTrkId[nPrimary];                                                                         // track id 
    int g4PrimaryProcess[nPrimary];                                                                       // interaction process
    float g4PrimaryX0[nPrimary], g4PrimaryY0[nPrimary], g4PrimaryZ0[nPrimary];                            // start pos
    float g4PrimaryXf[nPrimary], g4PrimaryYf[nPrimary], g4PrimaryZf[nPrimary];                            // final pos 
    float g4PrimaryPx[nPrimary], g4PrimaryPy[nPrimary], g4PrimaryPz[nPrimary];                            // momentum
    float g4PrimaryTrjX[nPrimary][maxTrjPoints], g4PrimaryTrjY[nPrimary][maxTrjPoints], g4PrimaryTrjZ[nPrimary][maxTrjPoints];    // trajectory sp
    float g4PrimaryTrjPx[nPrimary][maxTrjPoints], g4PrimaryTrjPy[nPrimary][maxTrjPoints], g4PrimaryTrjPz[nPrimary][maxTrjPoints]; // trajectory momentum 
    float g4PrimaryProjX0[nPrimary], g4PrimaryProjY0[nPrimary], g4PrimaryProjZ0[nPrimary];  // projected position

    // ##############################
    // ### Loop over g4 particles ###
    // ##############################
    {
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

      // project onto front of tpc
      // v = v0 + pHat*t
      // t = -1*(z0/pHat)
      TVector3 thisPos0( g4PrimaryX0[iPrim], g4PrimaryY0[iPrim], g4PrimaryZ0[iPrim] );
      TVector3 thisMomUnit( g4PrimaryPx[iPrim], g4PrimaryPy[iPrim], g4PrimaryPz[iPrim] );
      TVector3 thisPosProj = thisPos0;
      
      thisMomUnit = thisMomUnit.Unit();
      float t     = -1*( g4PrimaryZ0[iPrim]/thisMomUnit.Z() );
      thisPosProj = thisPos0 + t*thisMomUnit;

      g4PrimaryProjX0[iPrim] = thisPosProj.X();
      g4PrimaryProjY0[iPrim] = thisPosProj.Y();
      g4PrimaryProjZ0[iPrim] = thisPosProj.Z();

      // fill the proj histos
      hMCPrimaryProjX0->Fill(g4PrimaryProjX0[iPrim]);
      hMCPrimaryProjY0->Fill(g4PrimaryProjY0[iPrim]);
      hMCPrimaryProjZ0->Fill(g4PrimaryProjZ0[iPrim]);

      // store the track id 
      g4PrimaryTrkId[iPrim] = TrackId[iG4];

      // store the trajectory points and momentum
      nTrjPoints[iPrim] = NTrTrajPts[iG4];
      for (size_t iPoint = 0; iPoint < nTrjPoints[iPrim]; iPoint++)
      {
        g4PrimaryTrjX[iPrim][iPoint] = MidPosX[iG4][iPoint];
        g4PrimaryTrjY[iPrim][iPoint] = MidPosY[iG4][iPoint];
        g4PrimaryTrjZ[iPrim][iPoint] = MidPosZ[iG4][iPoint];

        g4PrimaryTrjPx[iPrim][iPoint] = MidPx[iG4][iPoint] * 1000; // convert to MeV
        g4PrimaryTrjPy[iPrim][iPoint] = MidPy[iG4][iPoint] * 1000; // convert to MeV
        g4PrimaryTrjPz[iPrim][iPoint] = MidPz[iG4][iPoint] * 1000; // convert to MeV
      }
      iPrim++;
    } 
    }

    // ###################################
    // ### Get the interaction process ###
    // ###################################
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if ( Mother[iG4] == g4PrimaryTrkId[nPrimary-1] ) g4PrimaryProcess[nPrimary-1] = Process[iG4];
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
        hMCPrimaryMissedTPCX->Fill(g4PrimaryXf[iPrim]);
        hMCPrimaryMissedTPCY->Fill(g4PrimaryYf[iPrim]);
        hMCPrimaryMissedTPCZ->Fill(g4PrimaryZf[iPrim]);
        continue;
      }

      // ####################################
      // ### This primary entered the TPC ###
      // ####################################
      nTracksEntered++;

      // calculate energy loss
      float energyLoss(0);
      // loop over trj points for this primary
      for (size_t iPoint = 0; iPoint < nTrjPoints[iPrim]; iPoint++)
      {
        // only look at the upstream portion
        if ( g4PrimaryTrjZ[iPrim][iPoint] > 0 ) continue;
        // ignore last point
        if ( (iPoint+1) >= nTrjPoints[iPrim] ) break;
 
        TVector3 momVec1 = TVector3( g4PrimaryTrjPx[iPrim][iPoint], 
                                     g4PrimaryTrjPy[iPrim][iPoint], 
                                     g4PrimaryTrjPz[iPrim][iPoint] );
        TVector3 momVec2 = TVector3( g4PrimaryTrjPx[iPrim][iPoint+1], 
                                     g4PrimaryTrjPy[iPrim][iPoint+1], 
                                     g4PrimaryTrjPz[iPrim][iPoint+1] );
        float mom1    = momVec1.Mag();
        float energy1 = std::sqrt( mom1*mom1 + PARTICLEMASS*PARTICLEMASS ); 
        float mom2 = momVec2.Mag();
        float energy2 = std::sqrt( mom2*mom2 + PARTICLEMASS*PARTICLEMASS ); 

        energyLoss += (energy1 - energy2);
      }
      hMCELossUpstream->Fill(energyLoss);
    }
    if (!isGoodEvent) continue;
    nGoodMCEvents++;

// ===================================================================
// =================  APPLY FRONT FACE TPC TRACK CUT =================
// ===================================================================
    bool isFrontTPCTrk(false);   
    std::vector<double> dummyTrkX,      dummyTrkY,      dummyTrkZ;
    std::vector<double> dummyTrkP0HatX, dummyTrkP0HatY, dummyTrkP0HatZ;

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
      float tempZMin(fvBoundVec[2](1));
      for (size_t iTrjPoint = 0; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
      {
        // ######################
        // ### Fiducial check ###
        // ######################
        TVector3 theTrjPoint( trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint] );
        TVector3 theTrjP0Hat( pHat0_X[iTrk][iTrjPoint], pHat0_Y[iTrk][iTrjPoint], pHat0_Z[iTrk][iTrjPoint] );
        if ( fvBoundVec[0](0) < theTrjPoint(0) && theTrjPoint(0) < fvBoundVec[0](1) &&
             fvBoundVec[1](0) < theTrjPoint(1) && theTrjPoint(1) < fvBoundVec[1](1) &&
             fvBoundVec[2](0) < theTrjPoint(2) && theTrjPoint(2) < tempZMin        )
        {
          // this point is in the fiducial volume
          tempZMin = theTrjPoint(2);
          // we're storing all tracks that are in the front of the TPC
          if (tempZMin < FIRSTSPZPOS)
          {
            isFrontTPCTrk = true;
            
            // Store these points for later
            dummyTrkX.push_back( theTrjPoint(0) );
            dummyTrkY.push_back( theTrjPoint(1) );
            dummyTrkZ.push_back( theTrjPoint(2) );
            dummyTrkP0HatX.push_back( theTrjP0Hat(0) );
            dummyTrkP0HatX.push_back( theTrjP0Hat(1) );
            dummyTrkP0HatX.push_back( theTrjP0Hat(2) );
          }
        }//<--- End if in FV
      }//<--- End loop over trj points 
    }//<-- End loop over reconstructed tracks

    // skip events that do not pass this cut
    if (!isFrontTPCTrk) continue;
    nEventsFrontTPCTrk++;

// =================================================================================
// =================  APPLY CUT FOR THE NUMBER OF UPSTREAM TRACKS  =================
// =================================================================================
    size_t nUpperTPCTrks(0);

    // #################################################################
    // ### Only keeping events if there is less than N tracks in the ###
    // ###    first y cm of the TPC (to help cut out EM Showers)     ###
    // #################################################################
    for (size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
    {
      // ###################################
      // ### Loop over trajectory points ###
      // ###################################
      bool isUpperTPCTrk(false);
      
      // looking for most upstream point still in FV of TPC
      // most upstream point in FV
      float tempZMin(fvBoundVec[2](1));
      for (size_t iTrjPoint = 0; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
      {
        // ######################
        // ### Fiducial check ###
        // ######################
        TVector3 theTrjPoint(trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint]);
        if ( fvBoundVec[0](0) < theTrjPoint(0) && theTrjPoint(0) < fvBoundVec[0](1) &&
             fvBoundVec[1](0) < theTrjPoint(1) && theTrjPoint(1) < fvBoundVec[1](1) &&
             fvBoundVec[2](0) < theTrjPoint(2) && theTrjPoint(2) < tempZMin        )
        {
          // this point is in the fiducial volume
          tempZMin = theTrjPoint(2);
          if (tempZMin < UPPERPARTOFTPC) isUpperTPCTrk = true;
        }
      }
      if (isUpperTPCTrk) nUpperTPCTrks++; 
    }

    // Skipping the event if there are too many 
    // fron TPC tracks in the event          
    if(nUpperTPCTrks > NUPPERTPCTRACKS || nUpperTPCTrks == 0){continue;}
    nEventsNUpperTPCTrk++;

// ========================================================================
// =================  MATCHING MC TO RECONSTRUCTED TRACK  =================
// ========================================================================
  if ( nPrimary > 1 ) std::cout << "\n\n\nTHERE IS MORE THAN ONE PRIMARY\n\n\n";
  
  bool MC_TPCMatch(false);
  size_t nMC_TPCMatch(0);
  // ###################################
  // ### Loop over all the US Tracks ###
  // ###################################
  for(size_t iUpTrk = 0; iUpTrk < dummyTrkX.size(); iUpTrk++)
  {
    float deltaX = dummyTrkX[iUpTrk] - g4PrimaryProjX0[0];
    float deltaY = dummyTrkY[iUpTrk] - g4PrimaryProjY0[0];
    float deltaZ = dummyTrkZ[iUpTrk] - g4PrimaryProjZ0[0];

    hDeltaX->Fill(deltaX);
    hDeltaY->Fill(deltaY);
    hDeltaZ->Fill(deltaZ);

    // matching in delta X and delta Y 
    if( DELTAXBOUND[0] < deltaX && deltaX < DELTAXBOUND[1] && 
        DELTAYBOUND[0] < deltaY && deltaY < DELTAYBOUND[1] )
    {
      MC_TPCMatch = true;
      nMC_TPCMatch++;
    }
  }//<---End bb loop
  if (MC_TPCMatch) nEventsDeltaMatch++;
  // force there to be one match
  if (nMC_TPCMatch != 1) continue;
  nEventsWC_TPCUniqueMatch++;

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


// ==============================================================
// =================  CALCULATING CROSS SECTION =================
// ==============================================================




  }//<---End loop over entries


// ==========================================================
// =================  EVENT REDUCTION TABLE =================
// ==========================================================
  std::cout << endl
            << "Events simulated:                        " << nTotalEvents        << endl
            << "Tracks entered TPC:                      " << nTracksEntered      << endl
	    << "Good MC events:                          " << nGoodMCEvents       << endl
            << "Events with front TPC track:             " << nEventsFrontTPCTrk  << endl
            << "Events with < N upper TPC tracks:        " << nEventsNUpperTPCTrk << endl
            << "Events with WC_TPC match (nocut):        " << nEventsDeltaMatch   << endl
            << "Events with unique WC_TPC match (nocut): " << nEventsWC_TPCUniqueMatch << endl           
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

  hMCPrimaryMissedTPCX->Write();  
  hMCPrimaryMissedTPCY->Write();  
  hMCPrimaryMissedTPCZ->Write();  

  hMCPrimaryProjX0->Write();
  hMCPrimaryProjY0->Write();
  hMCPrimaryProjZ0->Write();

  hDeltaX->Write();
  hDeltaY->Write();
  hDeltaZ->Write();

  myRootFile.Close();
}
