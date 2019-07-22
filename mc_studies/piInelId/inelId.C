/**
 * @file inelId.C
 * @brief MC studies for identifying pi inelastic
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#define inelId_cxx
#include "inelId.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <iostream>

// ============================================================================================
// ====================================  PROTOTYPES  ==========================================
// ============================================================================================
void MakePlots();
void ComputeAngles(float& phi, float& theta, const TVector3& p0Hat);
std::string ConvertProcessToString(const int& i);
void CalculateG4Xs();
bool InActiveRegion(const TVector3& thePos);
bool InTPCRegion(const TVector3& thePos);
void BeginJob();
inline void PrintVec(const TVector3& pos) {std::cout<<"\t("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")\n";};

// mass of pion in MeV
float PARTICLE_MASS(139.57); 
inline double ToKineticEnergy(const TVector3& mom){return std::sqrt( mom.Mag()*mom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; }


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

int IN_DEBUG(0);

// tpc boundaries
float TPC_X_BOUND[2] = {   0.0, 47.0 };
float TPC_Y_BOUND[2] = { -20.0, 20.0 };
float TPC_Z_BOUND[2] = {   0.0, 90.0 };

// fiducial volume definition
float FV_X_BOUND[2] = {   2.0, 45.0 };
float FV_Y_BOUND[2] = { -18.0, 18.0 };
float FV_Z_BOUND[2] = {   0.0, 88.0 };

// number of centimeters in Z we require a track to have a spacepoint
float FIRST_SP_Z_POS(4.0);

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
bool USE_EVENT_WEIGHT(false);

// constants for cross section calculation
double RHO(1395);                 //kg/m^3
double MOLAR_MASS(39.948);          //g/mol
double G_PER_KG(1000);
double AVOGADRO(6.022140857e+23);        //number/mol
double NUMBER_DENSITY( (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO );
double SLAB_WIDTH(0.0047);        //in m
double PI(3.141592654);
double M2_PER_BARN(1e-28);

// cut on angle between wc and tpc track (degrees)
float ALPHA_CUT(10.);

// threshold for hadron dEdx
double HIT_DEDX_THRESHOLD(40.);

// cut on length of secondary
double SECONDARY_LENGTH_CUT(3.0);

// cut on distance from vertex
double VERTEX_CUT(4.0);

// cut on angle between secondary and primary
double SECONDARY_ANGLE_CUT(10.);

// cut for tagging throughgoing tracks
double THROUGHGOING_Z_CUT(86.);

// ============================================================================================
// ====================================  HISTOGRAMS  ==========================================
// ============================================================================================
// WC to TPC matching/beamline
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
TH1D *hMCPrimaryPx = new TH1D("hMCPrimaryPx", "Primary Particle P_{x}", 22, 0, 1100);
TH1D *hMCPrimaryPy = new TH1D("hMCPrimaryPy", "Primary Particle P_{y}", 22, 0, 1100);
TH1D *hMCPrimaryPz = new TH1D("hMCPrimaryPz", "Primary Particle P_{z}", 22, 0 , 1100);
TH1D *hMCPrimaryP  = new TH1D("hMCPrimaryP", "Primary Particle P", 22, 0, 1100);
TH1D *hTrueLength = new TH1D("hTrueLength", "#True Length of the Primary Particle inside the TPC", 200, 0 , 100);
TH1D *hMostUpstreamZPos = new TH1D("hMostUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);
TH1D *hAlpha = new TH1D("hAlpha", "#alpha between MC Particle and TPC Track", 90, 0, 90);

// Reco
TH1D *hRecoIncidentKE = new TH1D("hRecoIncidentKE", "Incident for Reconstructed MC", 22, 0, 1100);
TH1D *hRecoInteractingKE = new TH1D("hRecoInteractingKE", "Interacting for Reconstructed MC", 22, 0, 1100);
TH1D *hRecoXSKE          = new TH1D("hRecoXSKE", "Reconstruced XS", 22, 0, 1100);
TH1D *hRecoFurthestInZCaloX = new TH1D("hRecoFurthestInZCaloX", "Most downstream in TPC Calorimetry X", 1300, TPC_X_BOUND[0]-20, TPC_X_BOUND[1]+20);
TH1D *hRecoFurthestInZCaloY = new TH1D("hRecoFurthestInZCaloY", "Most downstream in TPC Calorimetry Y", 1300, TPC_Y_BOUND[0]-20, TPC_Y_BOUND[1]+20);
TH1D *hRecoFurthestInZCaloZ = new TH1D("hRecoFurthestInZCaloZ", "Most downstream in TPC Calorimetry Z", 1300, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hRecoSecondaryLength    = new TH1D("hRecoSecondaryLength",     "Length of Secondary Tracks", 100, 0, 100);
TH1D* hRecoTrackVertexDiff    = new TH1D("hRecoTrackVertexDiff",          "Distance of Candidate Endpoints", 100, 0, 50); 
TH1D* hRecoVertexDiff        = new TH1D("hRecoVertexDiff",         "Vertex Reco - MC", 100, -50, 50); 
TH1D* hRecoSecondaries  = new TH1D("hRecoSecondaries", "Reconstructed number of tracks leaving vertex", 10, 0, 10);
TH1D* hRecoOneSecondaryTheta = new TH1D("hRecoOneSecondaryTheta", "Angle between primary and single secondary", 180, 0, 180);
TH2D* hRecoSmearingMatrix = new TH2D("hRecoSmearingMatrix", "Kinetic Energy smearing matrix", 22, 0, 1100, 22, 0, 1100);
TH1D* hRecoIntTypeBkg = new TH1D("hRecoIntTypeBkg", "Interacting Type Background", 22, 0, 1100);
TH1D* hRecoIntVertexBkg = new TH1D("hRecoIntVertexBkg", "Interacting Vertex Background", 22, 0, 1100);
TH1D* hRecoIntSignalBkg = new TH1D("hRecoIntSignalBkg", "Interacting Signal", 22, 0, 1100);
TH1D* hRecoFirstInTpcPointX = new TH1D("hRecoFirstInTpcPointX", "Reco First Point in TPC X", 400, TPC_X_BOUND[0]-20, 20);
TH1D* hRecoFirstInTpcPointY = new TH1D("hRecoFirstInTpcPointY", "Reco First Point in TPC Y", 400, TPC_Y_BOUND[0]-20, 20);
TH1D* hRecoFirstInTpcPointZ = new TH1D("hRecoFirstInTpcPointZ", "Reco First Point in TPC Z", 400, TPC_Z_BOUND[0]-20, 20);
TH1D* hRecoTrackPitch = new TH1D("hRecoTrackPitch", "Reco Track Pitch", 1000, TPC_Z_BOUND[0], 10);
TH1D* hRecoIntBkgDecay = new TH1D("hRecoIntBkgDecay", "Interaction Type Background from Decay", 22, 0, 1100);
TH1D* hRecoIntBkgElastic = new TH1D("hRecoIntBkgElastic", "Interaction Type Background from Elastic", 22, 0, 1100);
TH1D* hRecoIntBkgCapture = new TH1D("hRecoIntBkgCapture", "Interaction Type Background from Capture", 22, 0, 1100);

// mc
TH1D* hMCPreWCtoTPCInteractingKE = new TH1D("hMCPreWCtoTPCInteractingKE", "Pre WC to TPC Match True Interacting", 22, 0, 1100);
TH1D* hMCPreWCtoTPCIncidentKE = new TH1D("hMCPreWCtoTPCIncidentKE", "Pre WC to TPC Match True Incident", 22, 0, 1100);
TH1D *hMCPreWCtoTPCXSKE         = new TH1D("hMCPreWCtoTPCXSKE", "True XS", 22, 0, 1100);
TH1D* hMCLastPosFirstIntPosDiff = new TH1D("hMCLastPosFirstIntPosDiff", "MC Last Position - First Interacting Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCLastPosZ = new TH1D("hMCLastPosZ", "MC Last Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCUniformZ = new TH1D("hMCUniformZ", "MC Uniform Z", 1000, TPC_Z_BOUND[0], 10);
TH1D* hMCEnDep = new TH1D("hMCEnDep", "MC Energy Deposits", 400, -100, 100);
TH1D* hMCdEdX = new TH1D("hMCdEdX", "MC dEdX", 500, 0, 50);
TH1D* hMCInteractingKE = new TH1D("hMCInteractingKE", "True Interacting with cuts", 22, 0, 1100);
TH1D* hMCIncidentKE = new TH1D("hMCIncidentKE", "True Incident with cuts", 22, 0, 1100);
TH1D *hMCNoCutXSKE       = new TH1D("hMCNoCutXSKE", "True XS", 22, 0, 1100);
TH1D* hMCElasticAngle       = new TH1D("hMCElasticAngle", "Elastic Scattering Angle of Primary", 180, 0, 180);
TH1D* hMCInelasticAngle     = new TH1D("hMCInelasticAngle", "Angle Between Secondaries and Primary for Inelastic", 180, 0, 180);
TH1D* hMCInelasticOneVisDAngle = new TH1D("hMCInelasticOneVisDAngle", "Angle Between Single Visible Secondary and Primary for Inelastic", 180, 0, 180);
TH1D* hMCSecondaries      = new TH1D("hMCSecondaries",     "True number of tracks leaving vertex", 10, 0, 10);
TH1D* hMCNoCutIncidentKE = new TH1D("hMCNoCutIncidentKE", "True Incident", 22, 0, 1100);
TH1D* hMCNoCutInteractingKE = new TH1D("hMCNoCutInteractingKE", "True Interacting", 22, 0, 1100);
TH1D* hMCFirstInTpcPointX = new TH1D("hMCFirstInTpcPointX", "MC First Point in TPC X", 400, TPC_X_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointY = new TH1D("hMCFirstInTpcPointY", "MC First Point in TPC Y", 400, TPC_Y_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointZ = new TH1D("hMCFirstInTpcPointZ", "MC First Point in TPC Z", 400, TPC_Z_BOUND[0]-20, 20);
TH1D* hMCLastInTpcPointX = new TH1D("hMCLastInTpcPointX", "MC Last Point in TPC X", 1300, TPC_X_BOUND[0]-20, TPC_X_BOUND[1]+20);
TH1D* hMCLastInTpcPointY = new TH1D("hMCLastInTpcPointY", "MC Last Point in TPC Y", 1300, TPC_Y_BOUND[0]-20, TPC_Y_BOUND[1]+20);
TH1D* hMCLastInTpcPointZ = new TH1D("hMCLastInTpcPointZ", "MC Last Point in TPC Z", 1300, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);




// =======================================================================================
// ====================================  FILES  ==========================================
// =======================================================================================
TFile myRootFile("piMinusAna.root", "RECREATE");


// =====================================================================================
// ==============================  VARIABLES FOR G4 INFO  ==============================
// =====================================================================================
std::vector<size_t> g4PrimaryTrkId;                       // track id 
std::vector<TVector3> g4PrimaryPos0;                   // start pos
std::vector<TVector3> g4PrimaryPosf;                   // final pos 
std::vector<TVector3> g4PrimaryMom0;                   // momentum
std::vector<TVector3> g4PrimaryMomf;                   // final momentum
std::vector<TVector3> g4PrimaryProjPos0;               // projected position
std::vector<std::vector<TVector3>> g4PrimaryTrTrjPos;  // trajectory sp
std::vector<std::vector<TVector3>> g4PrimaryTrTrjMom;  // trajectory momentum 
std::vector<std::map<size_t, std::string>> g4PrimaryInteractions;
std::vector<bool> g4IsTrackInteracting, g4IsTrackSignal; 
std::vector<double> g4PrimaryKEFF;
std::vector<std::string> g4PrimarySubProcess;


/**
 * @brief Main loop
 * 
 * @param inDebug Option to run over only a few events
 */
void myana::Loop(int inDebug)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0/*47000*/; jentry<100000/*nentries*/;jentry++) 
  {
    // ###########################################
    // ### If debug, only look at a sub sample ###
    // ###########################################
    if (inDebug == 1 && jentry%1000 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << jentry << endl;
    IN_DEBUG = inDebug;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // #########################################
    // ### Increment our total event counter ###
    // #########################################
    nTotalEvents++; 
    if (nTotalEvents%1000 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl; 

    // TEMP
    //if (event == 1) cout << jentry << endl;
    if(IN_DEBUG) std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    //std::cout << "Event " << event << "\n\n";
    //if (nTotalEvents > 500) break;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Looking at only MC truth/G4 information

    // ### CLEAR G4 CONTAINERS
    // =====================================================================================
    // ==============================  VARIABLES FOR G4 INFO  ==============================
    // =====================================================================================
    g4PrimaryTrkId.clear();                     
    g4PrimaryInteractions.clear(); 
    g4PrimaryPos0.clear();                   
    g4PrimaryPosf.clear();                   
    g4PrimaryMom0.clear();                   
    g4PrimaryMomf.clear();                   
    g4PrimaryProjPos0.clear();               
    g4PrimaryTrTrjPos.clear();  
    g4PrimaryTrTrjMom.clear();      
    g4IsTrackInteracting.clear();
    g4IsTrackSignal.clear();
    g4PrimaryKEFF.clear();
    g4PrimarySubProcess.clear();

    // ##############################
    // ### Loop over g4 particles ###
    // ##############################
    {
    // ### Counter for the primary
    size_t iPrim(0);
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // ### If this is not a primary, skip it
      if (process_primary[iG4] == 0) continue;

      // ### We've got a new primary
      g4PrimaryInteractions.push_back(std::map<size_t,std::string>());
      g4PrimaryTrTrjPos.  push_back(std::vector<TVector3>());
      g4PrimaryTrTrjMom.  push_back(std::vector<TVector3>());
      g4PrimarySubProcess.push_back("none");

      // ### Store the track id 
      g4PrimaryTrkId.push_back(TrackId[iG4]); 

      // ### Store the interaction points
      for (size_t iPt = 0; iPt < (*InteractionPoint).size(); iPt++)
      {
        auto p     = (*InteractionPoint)[iPt];
        auto ptype = (*InteractionPointType)[iPt]; 
        TVector3 pos(MidPosX[iPrim][p], MidPosY[iPrim][p], MidPosZ[iPrim][p]);

        g4PrimaryInteractions[iPrim].emplace(p, ptype);
      }
      if (g4PrimaryInteractions[iPrim].size()) g4PrimarySubProcess[iPrim] = (*InteractionPointSubType)[iPrim];

      // ### Store the positions and momentum 
      g4PrimaryPos0.push_back( TVector3( StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]) );
      g4PrimaryPosf.push_back( TVector3( EndPointx[iG4],   EndPointy[iG4],   EndPointz[iG4]) );
      g4PrimaryMom0.push_back( 1000*TVector3( Px[iG4],     Py[iG4],          Pz[iG4]) );
      g4PrimaryMomf.push_back( 1000*TVector3( EndPx[iG4],  EndPy[iG4],       EndPz[iG4]) );

      // ### Setting a global event weight 
      // ### Setting Event weight 
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

      // ### Fill momentum histos
      hMCPrimaryPx->Fill(g4PrimaryMom0[iPrim].X(),   EVENT_WEIGHT);
      hMCPrimaryPy->Fill(g4PrimaryMom0[iPrim].Y(),   EVENT_WEIGHT);
      hMCPrimaryPz->Fill(g4PrimaryMom0[iPrim].Z(),   EVENT_WEIGHT);
      hMCPrimaryP ->Fill(g4PrimaryMom0[iPrim].Mag(), EVENT_WEIGHT);
      
      // ### Fill true length histo;
      hTrueLength->Fill( (g4PrimaryPosf[iPrim] - g4PrimaryPos0[iPrim]).Mag() );
      
      // ### Project onto tpc 
      TVector3 thisPosProjVec = g4PrimaryPos0[iPrim] - ( g4PrimaryPos0[iPrim].Z()/g4PrimaryMom0[iPrim].Z() )*g4PrimaryMom0[iPrim];
      g4PrimaryProjPos0.push_back(thisPosProjVec);

      // ### Fill the proj histos
      hMCPrimaryProjX0->Fill( g4PrimaryProjPos0[iPrim].X() );
      hMCPrimaryProjY0->Fill( g4PrimaryProjPos0[iPrim].Y() );
      hMCPrimaryProjZ0->Fill( g4PrimaryProjPos0[iPrim].Z() );
  
      // ### Store the trajectory points and momentum
      std::map<float, float> keMap;
      for (size_t iPoint = 0; iPoint < NTrTrajPts[iG4]; iPoint++)
      {
        TVector3 pos(     MidPosX[iPrim][iPoint],    MidPosY[iPrim][iPoint],    MidPosZ[iPrim][iPoint] );
        TVector3 mom(1000*MidPx[iPrim][iPoint], 1000*MidPy[iPrim][iPoint], 1000*MidPz[iPrim][iPoint] );
        if (pos.Z() < 0) keMap.emplace(pos.Z(), ToKineticEnergy(mom));
        
        g4PrimaryTrTrjPos[iPrim].push_back( pos );
        g4PrimaryTrTrjMom[iPrim].push_back( mom );   
      }

      // ### Set the KE at z = 0
      if(!keMap.empty()) g4PrimaryKEFF.push_back( (--keMap.end())->second );

      iPrim++;
    }//<--- End loop over G4 particles 
    }

    // ### Store if interacting/signal
    for (size_t iPrim = 0; iPrim < g4PrimaryInteractions.size(); iPrim++)
    {
      // ### Set the interacting variables
      bool isInterInTpc = false;
      bool isSignal     = true;
      for (const auto& p : g4PrimaryInteractions[iPrim])
      { 
        if (InActiveRegion(g4PrimaryTrTrjPos[iPrim][p.first])) isInterInTpc = true;
        isSignal = p.second.find("pi") != std::string::npos && p.second.find("Inelastic") != std::string::npos;
      }
      g4IsTrackInteracting.push_back( isInterInTpc );
      g4IsTrackSignal.push_back( isSignal );
    }

    // ### Sanity check...
    // ( asserts won't work with the interpreter :( )
    size_t sc = g4PrimaryTrkId.size();
    if (g4PrimaryTrkId.size()       != sc || g4PrimaryInteractions.size() != sc ||
        g4PrimaryPos0.size()        != sc || g4PrimaryPosf.size()         != sc ||   
        g4PrimaryMom0.size()        != sc || g4PrimaryProjPos0.size()     != sc ||
        g4PrimaryTrTrjPos.size()    != sc || g4PrimaryTrTrjMom.size()     != sc ||
        g4IsTrackSignal.size()      != sc || 
        g4IsTrackInteracting.size() != sc || g4PrimaryKEFF.size()         != sc ||
        g4PrimaryMomf.size()        != sc || g4PrimarySubProcess.size()   != sc) {cout << "CHECK PRIMARY INFO SIZES!\n"; return;};
    
// End looking at only MC truth/G4 information
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin apply "WC to TPC" cuts

    // ==========================================================================
    // =================  LOOKING AT EVENTS THAT ENTER THE TPC  =================
    // ==========================================================================
    // ### Loop over primaries 
    bool isGoodEvent(false);
    for (size_t iPrim = 0; iPrim < g4PrimaryTrkId.size(); iPrim++)
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
        if ( g4PrimaryTrTrjPos[iPrim][iPoint].Z() > 0 ) break;
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
    for (int iTrk = 0; iTrk < ntracks_reco; iTrk++)
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
      for (int iTrjPoint = 1; iTrjPoint < nTrajPoint[iTrk]; iTrjPoint++)
      {
        // ######################
        // ### Fiducial check ###
        // ######################
        TVector3 theTrjPointVec( trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint] );
        TVector3 tempVec       ( trjPt_X[iTrk][iTrjPoint], trjPt_Y[iTrk][iTrjPoint], trjPt_Z[iTrk][iTrjPoint] );
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
      hMostUpstreamZPos->Fill(tempMinZ);
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
      auto tpcP0Hat = frontFaceTrksP0Hat[iFrFaTrk];
      tpcP0Hat = tpcP0Hat.Unit();

      float thisPhi(0), thisTheta(0);
      ComputeAngles( thisPhi, thisTheta, tpcP0Hat );

      frontFaceTrksPhi.push_back(thisPhi);
      frontFaceTrksTheta.push_back(thisTheta);
    }//<--- End iFrFaTrk loop

    // sanity check
    sc = frontFaceTrksPos.size();
    if ( frontFaceTrksPos.size() != sc || frontFaceTrksP0Hat.size() != sc || 
         frontFaceTrksPhi.size() != sc || frontFaceTrksTheta.size() != sc ||
         frontFaceTrksId.size()  != sc ) {cout << "CHECK FRONT FACE CONTAINER SIZES\n"; return;};


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

// End applying "WC to TPC" cuts
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // ############################################################
    // ### We are know tracking a uniquely WC/TPC matched track ###
    // ###         The ID is in xsRecoTrkId variable            ###
    // ############################################################
    



  }//<---End loop over entries


  // ==========================================================
  // =================  EVENT REDUCTION TABLE =================
  // ==========================================================
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
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



/**
 * @brief Method to determine if track underwent inelastic interaction
 *
 */
bool myana::DetermineInelasticity(const size_t& xsRecoTrkId,
                                  const TVector3& furthestInZCaloPoint)
{
  // ### Look for any daughters around this furthest in Z calo point
  // ### Loop over tracks to find a starting point close to this
  if(IN_DEBUG) std::cout << "There are " << ntracks_reco << " reconstructed tracks...\n" << std::endl;
  std::map<size_t, bool> theTracksLeaving; // trk id, is tracking inverted
  bool isSignal = false;

  for(size_t iTrk = 0; iTrk < ntracks_reco; iTrk++)
  {
    // ### Skip our primary track
    if (iTrk == xsRecoTrkId) continue;

    // ### What determines the start and ending point here??
    TVector3 dXYZstart( trkvtxx[iTrk],  trkvtxy[iTrk], trkvtxz[iTrk] );
    TVector3 dXYZend  ( trkendx[iTrk],  trkendy[iTrk], trkendz[iTrk] );
    double length = (dXYZstart-dXYZend).Mag();
    
    // ### Which endpoint is connected to vertex
    double diff1 = (furthestInZCaloPoint - dXYZstart).Mag();
    double diff2 = (furthestInZCaloPoint - dXYZend  ).Mag();

    double diffArr[2] = {diff1, diff2};
    size_t minId    = diffArr[0] < diffArr[1] ? 0 : 1;
    hRecoTrackVertexDiff->Fill(diffArr[minId]);

    if(IN_DEBUG)
    {
      std::cout << "Track " << iTrk << " reco info\n";
      std::cout << "\tStart Point "; PrintVec(dXYZstart);
      std::cout << "\tEnd Point   "; PrintVec(dXYZend);
      std::cout << "\tLength " << length << std::endl;
      std::cout << "\tDistances from vertex " << diff1 << "  " << diff2 << std::endl;
      std::cout << "\tMinId is " << minId << std::endl;
      std::cout << "True info\n";
      for (size_t iPrim = 0; iPrim < g4PrimaryInteractions.size(); iPrim++)
      {
        for (const auto& v : g4PrimaryInteractions[iPrim]) PrintVec(g4PrimaryTrTrjPos[iPrim][v.first]);
      }
      std::cout << std::endl;
    }

    // ### Vertex cut and length cut
    if (diffArr[minId] < VERTEX_CUT) hRecoSecondaryLength->Fill(length);
    if (diffArr[minId] < VERTEX_CUT && length > SECONDARY_LENGTH_CUT)
    {
      if (minId == 0) theTracksLeaving.emplace(iTrk, false);
      if (minId == 1) theTracksLeaving.emplace(iTrk, true);
    }
  }//<-- End loop over trks

  if(IN_DEBUG) std::cout << "Tracks leaving = " << theTracksLeaving.size() << std::endl;
  hRecoSecondaries->Fill(theTracksLeaving.size());

  // ### If there are zero visible tracks leaving this vertex, 
  if (theTracksLeaving.size() == 0) 
  {
    // ### Maybe the primary scattered but tracking didn't break up the track
    // ### Check the scattering angle
   
  }

  // ### If there is just 1
  if (theTracksLeaving.size() == 1)
  {
    // ### I'm still interested in supplementing this with 
    // ### small energy deposits around vertex
    auto primTrkId  = xsRecoTrkId;
    auto secTrkId   = theTracksLeaving.begin()->first;
    bool isInverted = theTracksLeaving.begin()->second;
    
    // ### To compute the angle between tracks, I'll only use the first n
    // ### points, in case there's something funky going on in tracking

    // ### Grab the end points and the nth point from them
    std::vector<TVector3> primPoints, secPoints;
    int primNpointsBack = 10; int primNpointsTot = nTrajPoint[primTrkId];
    int secNpointsBack  = 10; int secNpointsTot  = nTrajPoint[secTrkId];

    primNpointsBack = primNpointsBack < primNpointsTot ? primNpointsBack : primNpointsTot;
    secNpointsBack  = secNpointsBack  < secNpointsTot  ? secNpointsBack  : secNpointsTot;

    // ### For primary, just in case, grab the most downstream
    if (trkvtxz[primTrkId] > trkendz[primTrkId])
    {
      // ### It's inverted, grab the first few points
      primPoints.push_back( TVector3(trjPt_X[primTrkId][primNpointsBack-1], 
                                     trjPt_Y[primTrkId][primNpointsBack-1], 
                                     trjPt_Z[primTrkId][primNpointsBack-1]) );      
      primPoints.push_back( TVector3(trjPt_X[primTrkId][0], 
                                     trjPt_Y[primTrkId][0], 
                                     trjPt_Z[primTrkId][0]) );
    }
    else
    {
      // ### It's not inverted, grab the last few points
      primPoints.push_back( TVector3(trjPt_X[primTrkId][primNpointsTot-primNpointsBack], 
                                     trjPt_Y[primTrkId][primNpointsTot-primNpointsBack], 
                                     trjPt_Z[primTrkId][primNpointsTot-primNpointsBack]) );
      primPoints.push_back( TVector3(trjPt_X[primTrkId][primNpointsTot-1], 
                                     trjPt_Y[primTrkId][primNpointsTot-1], 
                                     trjPt_Z[primTrkId][primNpointsTot-1]) );                                     
    }

    // ### Now handle secondary
    if (!isInverted) 
    {
      // ### It's the start point
      secPoints.push_back( TVector3(trjPt_X[secTrkId][0], 
                                    trjPt_Y[secTrkId][0], 
                                    trjPt_Z[secTrkId][0]) );
      secPoints.push_back( TVector3(trjPt_X[secTrkId][secNpointsBack-1], 
                                    trjPt_Y[secTrkId][secNpointsBack-1], 
                                    trjPt_Z[secTrkId][secNpointsBack-1]) );      
    }
    else 
    {
      // ### It's the end point
      secPoints.push_back( TVector3(trjPt_X[secTrkId][secNpointsTot-1], 
                                    trjPt_Y[secTrkId][secNpointsTot-1], 
                                    trjPt_Z[secTrkId][secNpointsTot-1]) );
      secPoints.push_back( TVector3(trjPt_X[secTrkId][secNpointsTot-secNpointsBack], 
                                    trjPt_Y[secTrkId][secNpointsTot-secNpointsBack], 
                                    trjPt_Z[secTrkId][secNpointsTot-secNpointsBack]) );
    }

    // ### Now form our direction vectors
    TVector3 primDir(primPoints[1] - primPoints[0]);
    TVector3 secDir (secPoints[1]  - secPoints[0]);

    double theta = (180/TMath::Pi())*std::acos( primDir.Unit().Dot(secDir.Unit()) );
    hRecoOneSecondaryTheta->Fill(theta);
    if (theta > SECONDARY_ANGLE_CUT) isSignal = true;

    if(IN_DEBUG)
    {
      std::cout << "Primary Start Point "; PrintVec(primPoints[0]);
      std::cout << "Primary End Point   "; PrintVec(primPoints[1]);
      std::cout << "Second. Start Point "; PrintVec(secPoints[0]);
      std::cout << "Second. End Point   "; PrintVec(secPoints[1]);
      std::cout << "THETA " << theta << std::endl;
    }
  }//<-- End if 1 visible secondary

  // ### If there are at least 2 visible tracks leaving this vertex, yes
  if (theTracksLeaving.size() >= 2) isSignal = true;

  if(IN_DEBUG)
  {
    if (isSignal) std::cout << "Determined as inelastic!" << std::endl;
    else std::cout << "Determined NOT inelastic!\n";
  }
  
  return isSignal;
}//<-- End DertermineInelasticity


/**
 * @brief Make plots and save to output file
 * 
 */
void MakePlots()
{
  // ### make the pre wc to TPC xs plot
  for (int iBin = 1; iBin <= hMCPreWCtoTPCIncidentKE->GetNbinsX(); iBin++)
  {
    if (hMCPreWCtoTPCIncidentKE->GetBinContent(iBin) == 0) continue;

    // ### our cross section
    float tempXS = (hMCPreWCtoTPCInteractingKE->GetBinContent(iBin)/hMCPreWCtoTPCIncidentKE->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hMCPreWCtoTPCXSKE->SetBinContent(iBin, tempXS);

    // ### incident taken as poissonian
    float denomError = std::sqrt(hMCPreWCtoTPCIncidentKE->GetBinContent(iBin));
    float denom      = hMCPreWCtoTPCIncidentKE->GetBinContent(iBin);
    if (denom == 0) continue;
    float term2 = denomError/denom;

    auto intCounts = hMCPreWCtoTPCInteractingKE->GetBinContent(iBin);
    auto incCounts = hMCPreWCtoTPCIncidentKE->GetBinContent(iBin);
    float var      = intCounts*( 1 - intCounts/incCounts );
    float numError = std::sqrt(var);
    float num      = intCounts;

    if (num)
    {
      float term1 = numError/num;
      float xs    = hMCPreWCtoTPCXSKE->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 );

      hMCPreWCtoTPCXSKE->SetBinError(iBin,totalError);
    }
  }

  // ### Make the mc no cuts XS plot
  for (int iBin = 1; iBin <= hMCNoCutIncidentKE->GetNbinsX(); iBin++)
  {
    if (hMCNoCutIncidentKE->GetBinContent(iBin) == 0) continue;

    // ### our cross section
    float tempXS = (hMCNoCutInteractingKE->GetBinContent(iBin)/hMCNoCutIncidentKE->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hMCNoCutXSKE->SetBinContent(iBin, tempXS);

    // ### incident taken as poissonian
    float denomError = std::sqrt(hMCNoCutIncidentKE->GetBinContent(iBin));
    float denom      = hMCNoCutIncidentKE->GetBinContent(iBin);
    if (denom == 0) continue;
    float term2 = denomError/denom;

    auto intCounts = hMCNoCutInteractingKE->GetBinContent(iBin);
    auto incCounts = hMCNoCutIncidentKE->GetBinContent(iBin);
    float var      = intCounts*( 1 - intCounts/incCounts );
    float numError = std::sqrt(var);
    float num      = intCounts;

    if (num)
    {
      float term1 = numError/num;
      float xs    = hMCNoCutXSKE->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 );

      hMCNoCutXSKE->SetBinError(iBin,totalError);
    }
  }

  // ### NEED TO MAKE MC WITH CUTS

  // ### Make the reco XS plot
  for (int iBin = 1; iBin <= hRecoIncidentKE->GetNbinsX(); iBin++)
  {
    if (hRecoIncidentKE->GetBinContent(iBin) == 0) continue;

    // ### our cross section
    float tempXS = (hRecoInteractingKE->GetBinContent(iBin)/hRecoIncidentKE->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hRecoXSKE->SetBinContent(iBin, tempXS);

    // ### incident taken as poissonian
    float denomError = std::sqrt(hRecoIncidentKE->GetBinContent(iBin));
    float denom      = hRecoIncidentKE->GetBinContent(iBin);
    if (denom == 0) continue;
    float term2 = denomError/denom;

    auto intCounts = hRecoInteractingKE->GetBinContent(iBin);
    auto incCounts = hRecoIncidentKE->GetBinContent(iBin);
    float var      = intCounts*( 1 - intCounts/incCounts );
    float numError = std::sqrt(var);
    float num      = intCounts;

    if (num)
    {
      float term1 = numError/num;
      float xs    = hRecoXSKE->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 );

      hRecoXSKE->SetBinError(iBin,totalError);
    }
  }

  // ### Normalize smearing matrix
  for (size_t iBinX = 1; iBinX <= hRecoSmearingMatrix->GetXaxis()->GetNbins(); iBinX++)
  {
    double sum = 0;
    for(size_t iBinY = 1; iBinY <= hRecoSmearingMatrix->GetYaxis()->GetNbins(); iBinY++)
    {
      sum = sum + hRecoSmearingMatrix->GetBinContent(iBinX, iBinY);
    }
    for(size_t iBinY = 1; iBinY <= hRecoSmearingMatrix->GetYaxis()->GetNbins(); iBinY++)
    {
      auto content = hRecoSmearingMatrix->GetBinContent(iBinX, iBinY);
      hRecoSmearingMatrix->SetBinContent(iBinX, iBinY, content/sum);
    }
  }

  myRootFile.cd();
  myRootFile.mkdir("wcToTpc/");
  myRootFile.cd("wcToTpc/");
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
  hMostUpstreamZPos->Write();
  hAlpha->Write();

  myRootFile.cd();
  myRootFile.mkdir("reco/");
  myRootFile.cd("reco/");
  hRecoIncidentKE->Write();
  hRecoInteractingKE->Write();
  hRecoXSKE->Write();
  hRecoFurthestInZCaloX->Write();
  hRecoFurthestInZCaloY->Write();
  hRecoFurthestInZCaloZ->Write();
  hRecoSecondaryLength->Write();
  hRecoTrackVertexDiff->Write();
  hRecoSecondaries->Write();
  hRecoOneSecondaryTheta->Write();
  hRecoSmearingMatrix->Write();
  hRecoIntTypeBkg  ->Write();
  hRecoIntVertexBkg->Write();
  hRecoIntSignalBkg->Write();
  hRecoIntBkgDecay->Write();
  hRecoIntBkgElastic->Write();
  hRecoIntBkgCapture->Write();
  hRecoVertexDiff->Write();
  hRecoFirstInTpcPointX->Write();
  hRecoFirstInTpcPointY->Write();
  hRecoFirstInTpcPointZ->Write();
  hRecoTrackPitch->Write();

  myRootFile.cd();
  myRootFile.mkdir("mc/");
  myRootFile.cd("mc/");
  hMCElasticAngle->Write();
  hMCInelasticAngle->Write();
  hMCInelasticOneVisDAngle->Write();
  hMCSecondaries->Write();
  hMCNoCutIncidentKE->Write();
  hMCNoCutInteractingKE->Write();
  hMCNoCutXSKE->Write();
  hMCFirstInTpcPointX->Write();
  hMCFirstInTpcPointY->Write();
  hMCFirstInTpcPointZ->Write();
  hMCLastInTpcPointX->Write();
  hMCLastInTpcPointY->Write();
  hMCLastInTpcPointZ->Write();
  hMCPreWCtoTPCIncidentKE->Write();
  hMCPreWCtoTPCInteractingKE->Write();
  hMCPreWCtoTPCXSKE->Write();
  hMCLastPosFirstIntPosDiff->Write();
  hMCLastPosZ->Write();
  hMCUniformZ->Write();
  hMCEnDep->Write();
  hMCdEdX->Write();
  hMCIncidentKE->Write();
  hMCInteractingKE->Write();

  myRootFile.Close();
}



/**
 * @brief Compute angle for beamline matching
 * 
 */
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


/**
 * @brief Check if in active region
 *
 */
bool InActiveRegion( const TVector3& thePos  )
{
  if ( FV_X_BOUND[0] < thePos.X() && thePos.X() < FV_X_BOUND[1] &&
       FV_Y_BOUND[0] < thePos.Y() && thePos.Y() < FV_Y_BOUND[1] &&
       FV_Z_BOUND[0] < thePos.Z() && thePos.Z() < FV_Z_BOUND[1] ) return true;

  return false;
}



/**
 * @brief Check if position in TPC
 * 
 */
bool InTPCRegion( const TVector3& thePos  )
{
  if ( TPC_X_BOUND[0] < thePos.X() && thePos.X() < TPC_X_BOUND[1] &&
       TPC_Y_BOUND[0] < thePos.Y() && thePos.Y() < TPC_Y_BOUND[1] &&
       TPC_Z_BOUND[0] < thePos.Z() && thePos.Z() < TPC_Z_BOUND[1] ) return true;

  return false;
}



