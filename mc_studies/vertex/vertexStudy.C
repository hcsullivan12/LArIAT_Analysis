/**
 * @file vertexStudy.C
 * @brief Vertex study
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#define vertexStudy_cxx
#include "vertexStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <iostream>
#include <TVector3.h>

/// @name Prototypes
/// @{
void MakePlots();
bool InActiveRegion(const TVector3& thePos);
bool InTPCRegion(const TVector3& thePos);
inline void PrintVec(const TVector3& pos) {std::cout<<"\t("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")\n";};
/// Mass of pion in MeV
float PARTICLE_MASS(139.57); 
inline double ToKineticEnergy(const TVector3& mom){return std::sqrt( mom.Mag()*mom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; }
/// @}

/// @name Useful variables
/// @{
/// Counters
size_t nTotalEvents(0),  nPrimariesEntered(0), 
       nGoodMCEvents(0), nEventsInelastic(0),
       nElasticLikeEvents(0);

/// Primary track Id
int primTrkId(-1);

/// Vector of g4 Ids for visible secondaries
std::vector<size_t> visSec;

/// Container for pdg codes for charged particles (cuz I can't remember) 
struct PdgCodes_t
{
  int kElectron = 11;
  int kMuon     = 13;
  int kPion     = 211;
  int kKaon     = 321;
  int kProton   = 2212;
};
PdgCodes_t pdgcodes;
/// Method to check if charged
inline bool IsCharged(int const& pdg) 
{ 
  return (pdg == pdgcodes.kElectron || pdg == pdgcodes.kMuon || 
          pdg == pdgcodes.kPion     || pdg == pdgcodes.kKaon || 
          pdg == pdgcodes.kProton);
}

/// vector of TrackIdes 
std::vector<TrackIde_t> ides;

/// vector of vertices
std::vector<Vertex_t> vertices;
/// @}

/// @name Cuts and constants
/// @{
int IN_DEBUG(0);

/// TPC boundaries
float TPC_X_BOUND[2] = {   0.0, 47.0 };
float TPC_Y_BOUND[2] = { -20.0, 20.0 };
float TPC_Z_BOUND[2] = {   0.0, 90.0 };

/// Fiducial volume definition
float FV_X_BOUND[2] = {   1.0, 46.0 };
float FV_Y_BOUND[2] = { -18.0, 18.0 };
float FV_Z_BOUND[2] = {   0.0, 88.0 };

/// Track length cut for secondaries
float SECONDARY_LENGTH_CUT(2.0);

/// The assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(36); //MeV
/// @}

/// @name Histograms
/// @{
/// WC to TPC matching/beamline
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

/// mc
TH1D* hMCLastPosFirstIntPosDiff = new TH1D("hMCLastPosFirstIntPosDiff", "MC Last Position - First Interacting Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCLastPosZ = new TH1D("hMCLastPosZ", "MC Last Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCElasticAngle       = new TH1D("hMCElasticAngle", "Elastic Scattering Angle of Primary", 180, 0, 180);
TH1D* hMCInelasticAngle     = new TH1D("hMCInelasticAngle", "Angle Between Secondaries and Primary for Inelastic", 180, 0, 180);
TH1D* hMCInelasticOneVisDAngle = new TH1D("hMCInelasticOneVisDAngle", "Angle Between Single Visible Secondary and Primary for Inelastic", 180, 0, 180);
TH1D* hMCElasticConeAngle = new TH1D("hMCElasticConeAngle", "Cone Angle Between Single Visible Secondary and Primary for Elastic", 90, 0, 90);
TH1D* hMCInelasticConeAngle = new TH1D("hMCInelasticConeAngle", "Cone Angle Between Single Visible Secondary and Primary for Inelastic", 90, 0, 90);
TH1D* hMCSecondaries      = new TH1D("hMCSecondaries",     "True number of tracks leaving vertex", 10, 0, 10);
TH1D* hMCFirstInTpcPointX = new TH1D("hMCFirstInTpcPointX", "MC First Point in TPC X", 400, TPC_X_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointY = new TH1D("hMCFirstInTpcPointY", "MC First Point in TPC Y", 400, TPC_Y_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointZ = new TH1D("hMCFirstInTpcPointZ", "MC First Point in TPC Z", 400, TPC_Z_BOUND[0]-20, 20);
TH1D* hMCLastInTpcPointX = new TH1D("hMCLastInTpcPointX", "MC Last Point in TPC X", 1300, TPC_X_BOUND[0]-20, TPC_X_BOUND[1]+20);
TH1D* hMCLastInTpcPointY = new TH1D("hMCLastInTpcPointY", "MC Last Point in TPC Y", 1300, TPC_Y_BOUND[0]-20, TPC_Y_BOUND[1]+20);
TH1D* hMCLastInTpcPointZ = new TH1D("hMCLastInTpcPointZ", "MC Last Point in TPC Z", 1300, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCSecondaryTrkLength = new TH1D("hMCSecondaryTrkLength", "Secondary Track Lengths", 100, 0, 50);
TH2D* hMCIdeVertexElastic = new TH2D("hMCIdeVertexElastic", "IDEs versus bubble radius elastic", 20, 0, 10, 200, 0, 200);
TH2D* hMCIdeVertexInelastic = new TH2D("hMCIdeVertexInelastic", "IDEs versus bubble radius inelastic", 20, 0, 10, 200, 0, 200);
TH2D* hMCIdesXvsZ = new TH2D("hMCIdesXvsZ", "IDEs for X vs Z ", 200, 0, 100, 100, 0, 50);
TH2D* hMCIdesSubXvsZ = new TH2D("hMCIdesSubXvsZ", "IDEs subtracted for X vs Z ", 200, 0, 100, 100, 0, 50);
TH1I* hMCNumInteractions = new TH1I("hMCNumInteractions", "Number of interactions in TPC for elastic like", 5, 0, 5);
/// @}

/// Output root file
TFile myRootFile("piMinusAna.root", "RECREATE");

Long64_t jentry(0);

/// Fiducial volume definition
float FV_X_BOUND[2] = {   1.0, 46.0 };
float FV_Y_BOUND[2] = { -18.0, 18.0 };
float FV_Z_BOUND[2] = {   0.0, 88.0 };

/**
 * @brief Main loop
 * 
 * @param inDebug Option to run over only a few events
 */
void vertexStudy::Loop(int inDebug)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (jentry=2; jentry<=2/*nentries*/;jentry++) 
  {
    // If debug, only look at a sub sample 
    //if (inDebug == 1 && jentry%1000 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << jentry << endl;
    IN_DEBUG = inDebug;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // Increment our total event counter 
    nTotalEvents++; 
    if (nTotalEvents%1000 == 0) std::cout << "EVENT = " << nTotalEvents << std::endl; 

    // TEMP
    //if (event == 1) cout << jentry << endl;
    if(IN_DEBUG) std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    //std::cout << "Event " << event << "\n\n";
    //if (nTotalEvents > 500) break;


    // Filling primary information 
    primTrkId = -1;
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // If this is not a primary, skip it
      if (process_primary[iG4] == 0) continue;
      primTrkId = TrackId[iG4];

      // ### Store the positions and momentum 
      TVector3 g4PrimaryPos0( StartPointx[iG4], StartPointy[iG4], StartPointz[iG4] );
      TVector3 g4PrimaryPosf( EndPointx[iG4],   EndPointy[iG4],   EndPointz[iG4] );
      TVector3 g4PrimaryMom0(1000*Px[iG4],    1000*Py[iG4],    1000*Pz[iG4]);
      TVector3 g4PrimaryMomf(1000*EndPx[iG4], 1000*EndPy[iG4], 1000*EndPz[iG4]);

      // Fill momentum histos
      hMCPrimaryPx->Fill(g4PrimaryMom0.X());
      hMCPrimaryPy->Fill(g4PrimaryMom0.Y());
      hMCPrimaryPz->Fill(g4PrimaryMom0.Z());
      hMCPrimaryP ->Fill(g4PrimaryMom0.Mag());
      
      // Fill true length histo;
      hTrueLength->Fill( (g4PrimaryPosf - g4PrimaryPos0).Mag() );
      
      // Project onto tpc 
      TVector3 g4PrimaryProjPos0 = g4PrimaryPos0 - ( g4PrimaryPos0.Z()/g4PrimaryMom0.Z() )*g4PrimaryMom0;

      // Fill the proj histos
      hMCPrimaryProjX0->Fill( g4PrimaryProjPos0.X() );
      hMCPrimaryProjY0->Fill( g4PrimaryProjPos0.Y() );
      hMCPrimaryProjZ0->Fill( g4PrimaryProjPos0.Z() );
    }//<--- End loop over G4 particles
    if (primTrkId != 1) std::cout << "WARNING! Check the primary track ID!\n";
    
    // Check if entered TPC
    bool isGoodEvent(false);
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // If this is not a primary, skip it
      if (process_primary[iG4] == 0) continue;
      
      // only look if a primary entered the tpc
      if ( EndPointz[iG4] > 0 ) isGoodEvent = true;
      if ( !isGoodEvent ) 
      {
        hMCPrimaryMissedTpcX->Fill(EndPointx[iG4]);
        hMCPrimaryMissedTpcY->Fill(EndPointy[iG4]);
        hMCPrimaryMissedTpcZ->Fill(EndPointz[iG4]);
        continue;
      }

      nPrimariesEntered++;
      // calculate energy loss
      float energyLoss(0);
      // loop over trj points for this primary
      for (size_t iPoint = 0; iPoint < NTrTrajPts[iG4]; iPoint++)
      {
        // only look at the upstream portion
        if ( MidPosZ[iG4][iPoint] > 0 ) break;

        // ignore last point
        if ( (iPoint+1) >= NTrTrajPts[iG4] ) break;
        TVector3 mom1Vec(MidPx[iG4][iPoint], MidPy[iG4][iPoint], MidPz[iG4][iPoint]);
        TVector3 mom2Vec(MidPx[iG4][iPoint+1], MidPy[iG4][iPoint+1], MidPz[iG4][iPoint+1]);

        float energy1 = std::sqrt( mom1Vec.Mag()*mom1Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS );
        float energy2 = std::sqrt( mom2Vec.Mag()*mom2Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS ); 
        energyLoss += (energy1 - energy2);
      }//<--- End loop over true traj points
      hMCELossUpstream->Fill(energyLoss);
    }//<--- End loop over primaries
    if (!isGoodEvent) continue;
    nGoodMCEvents++;

    if (nTotalEvents == 1)
    {
      for (int iIDE = 0; iIDE < IDETrackId->size(); iIDE++) 
      {
        auto xBin = hMCIdesXvsZ->GetXaxis()->FindBin(IDEPos->at(iIDE)[2]);
        auto yBin = hMCIdesXvsZ->GetYaxis()->FindBin(IDEPos->at(iIDE)[0]);
        hMCIdesXvsZ->SetBinContent(xBin, yBin, IDEEnergy->at(iIDE));
      }
    }
  

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin study

    // Get first interaction in TPC
    if (IN_DEBUG) std::cout << "Getting first interaction in TPC...\n";
    Vertex_t firstVertexInTpc(-1, "none", TVector3(0,0,-100));
    vertices.clear();
    GetFirstInteractionInTPC(firstVertexInTpc);
    //firstVertexInTpc.position.Print();
    //for (const auto& v : vertices) cout << v.position.Z() << " " << v.process << std::endl;
    if (vertex.process.find("none") != std::string::npos) continue;

    // Identify visible secondaries
    visSec.clear();
    IdentifyVisibleSecondaries(firstVertexInTpc);

    // Angle study
    if (IN_DEBUG) std::cout << "Starting angle study...\n";
    bool isElasticLike(false);
    AngleStudy(firstVertexInTpc, isElasticLike);

    // Option to only do events for which there is only one elastic like interaction in TPC
    // Is we have something other than elastic, CoulombScat, or pi-Inelastic, skip
    bool doSkip(false);
    {
      bool isElastic(false);
      bool isInelastic(false);
      bool isOther(false);

      // Check the vertices
      //cout << "\n";
      for (const auto& v : vertices)
      {
        if      (v.process.find("hadElastic") != std::string::npos || v.process.find("CoulombScat") != std::string::npos) isElastic = true;
        else if (v.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;
        else isOther = true;
        //std::cout << v.process << "\n";// v.position.Print();
      }
      //// Check the interactionxxx vectors too in case we missed something
      //for (const auto& v : *InteractionPointType)
      //{
      //       if (v.find("hadElastic") != std::string::npos || v.find("CoulombScat") != std::string::npos) isElastic = true;
      //  else if (v.find("pi-Inelastic") != std::string::npos) isInelastic = true;
      //  else isOther = true;
      //}
      //if (isElastic && isInelastic) doSkip = true;
      if (!isElastic && !isInelastic) doSkip = true;
      if (isOther) doSkip = true;

      if (IN_DEBUG) 
      {
        std::cout << "isElastic = " << isElastic << "  isInelastic = " << isInelastic << "  isOther = " << isOther << std::endl;
        std::cout << "Vertices:\n";
        for (const auto& v : vertices) std::cout << v.process << " " << v.position.Z() << std::endl;
        std::cout << "Interactionxxx:\n";
        for (const auto& v : *InteractionPointType) std::cout << v << std::endl;
      }
    }
    if (doSkip) continue;

    if (isElasticLike) nElasticLikeEvents++;
    if (IN_DEBUG) 
    {
      std::cout << "\njEntry                 = " << jentry 
                << "\nFirst int point in TPC = " << MidPosZ[0][firstVertexInTpc.p]
                << "\nFirst process in TPC   = " << firstVertexInTpc.process
                << "\nIs elastic like        = ";
      if (isElasticLike) std::cout << "yes" << std::endl;
      else               std::cout << "no"  << std::endl;
      std::cout << "\nContains both          = ";
      if (doSkip)        std::cout << "yes" << std::endl;
      else               std::cout << "no"  << std::endl;
      for (const auto& p : *InteractionPointType) std::cout << p << endl;
    }

    // Vertex study for elastic like events
    if (isElasticLike) 
    {
      hMCNumInteractions->Fill(InteractionPoint->size());

      ides.clear();
      ides.reserve(IDETrackId->size());
      for (int iIDE = 0; iIDE < IDETrackId->size(); iIDE++) 
      {
        TrackIde_t tide( (*IDETrackId)[iIDE],
                         DetermineG4Id(tide.trackId),
                         (*IDEEnergy)[iIDE], 
                         TVector3( (*IDEPos)[iIDE][0], (*IDEPos)[iIDE][1], (*IDEPos)[iIDE][2] ) );
        ides.push_back(tide);
      }
     
      if (IN_DEBUG) std::cout << "Starting vertex study...\n";
      VertexStudy(firstVertexInTpc);
    }

// End study
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  }//<---End loop over entries


  // Event reduction table
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events simulated:                           " << nTotalEvents                 << endl
            << "Tracks entered TPC:                         " << nPrimariesEntered            << endl
            << "Events with > 0 entered tracks:             " << nGoodMCEvents                << endl
            << "Elastic like events:                        " << nElasticLikeEvents           << endl
            << "Inelastic interactions:                     " << nEventsInelastic             << endl
            << endl;

  // Make plots
  MakePlots();

  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Method to get the first interaction in tpc
 * 
 * @note Something is wrong with the InteractionPoint vectors.
 *       I have seen events for which the Interactionxxx vectors
 *       contain only a hadElastic but see a piInelastic when 
 *       checking the daughters. So I will check the Process
 *       labels for each daughter of primary.
 * 
 * @note When an elastic shows signs of > 0 ide, it's most likely
 *       the case that there is a subsequent inelastic interaction.
 * 
 * @param point The primary point
 * @param process The process type
 */
void vertexStudy::GetFirstInteractionInTPC(Vertex_t& point)
{
  // Check Interactionxxx just in case
  for (int iPt = 0; iPt < InteractionPoint->size(); iPt++)
  {
    Vertex_t thisVertex( InteractionPoint->at(iPt),
                         InteractionPointType->at(iPt),
                         TVector3(MidPosX[0][InteractionPoint->at(iPt)], MidPosY[0][InteractionPoint->at(iPt)], MidPosZ[0][InteractionPoint->at(iPt)]) );
    // Check if in FV
    if (!InActiveRegion(thisVertex.position)) continue;

    bool weGotIt(false);
    for (const auto& v : vertices)
    {
      float diff = (v.position-thisVertex.position).Mag();
      if (  diff < 0.01 && v.process == thisVertex.process) weGotIt = true;
    }
    if (!weGotIt) vertices.push_back(thisVertex);
  }
  std::sort(vertices.begin(), vertices.end(), [](const auto& l, const auto& r) {return l.position.Z() < r.position.Z();});

  if (vertices.size()) point = vertices[0];
}

/**
 * @brief Method to get G4 Id from track id
 * 
 * @param tid Track id
 * @return int The particle's G4 Id
 */
int vertexStudy::DetermineG4Id(const int& tid)
{
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    if (TrackId[iG4] == tid) return iG4;
  }
  return -999;
} 

/**
 * @brief Method to study angles for elastic like events
 * 
 * @param vertex The first vertex in the tpc
 * @param isElasticLike Is this process elastic like (one visible secondary)
 */
void vertexStudy::AngleStudy(const Vertex_t& vertex, bool& isElasticLike)
{
  // Handle the case of just the primary
  bool isElastic(false);
  bool isInelastic(false);
  if (vertex.process.find("hadElastic")   != std::string::npos || vertex.process.find("CoulombScat") != std::string::npos) isElastic = true;
  if (vertex.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;

  // Only looking at elastic or inelastic
  if (!isElastic && !isInelastic) return;

  // Get incident momentum
  TVector3 primIncMom( MidPx[0][vertex.point-1], MidPy[0][vertex.point-1], MidPz[0][vertex.point-1] );

  // If the primary doesn't end here
  if ( (vertex.point+1) < NTrTrajPts[0] ) 
  {
    TVector3 trailMom( MidPx[0][vertex.point], MidPy[0][vertex.point], MidPz[0][vertex.point] );
    double theta = (180/TMath::Pi())*std::acos( primIncMom.Unit().Dot(trailMom.Unit()) );
    if (isElastic)   
    {
      isElasticLike = true;
      hMCElasticAngle->Fill(theta);
    }
    else if (isInelastic) hMCInelasticAngle->Fill(theta); 
  }
  
  // Fill inelastic scattering angles for all relevant daughters
  if (isInelastic)
  {
    for (const auto& iG4 : visSec)
    {
      if (iG4 == DetermineG4Id(primTrkId)) continue;

      // This is a charged daughter
      TVector3 dMom0(Px[iG4], Py[iG4], Pz[iG4]);

      // Compute angle between daughter and incoming primary
      double theta = (180/TMath::Pi())*std::acos( primIncMom.Unit().Dot(dMom0.Unit()) );
      hMCInelasticAngle->Fill(theta);
    }//<-- End loop over visible particles
  }//<--End if inelastic

  // Do it again, except for cases in which inelastic but one visible daughter
  if (isInelastic)
  {
    // If visSec.size() == 1 this looks elastic
    hMCSecondaries->Fill(visSec.size());
    if (visSec.size() == 1)
    {
      size_t iG4 = visSec[0];
      isElasticLike = true;
      TVector3 trailMom = primIncMom;

      // If this visible secondary is the primary
      if (iG4 == DetermineG4Id(primTrkId)) trailMom = TVector3( MidPx[0][vertex.point], MidPy[0][vertex.point], MidPz[0][vertex.point] );
      else                                 trailMom = TVector3( Px[iG4], Py[iG4], Pz[iG4] );
    
      double theta = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(trailMom.Unit()));
      hMCInelasticOneVisDAngle->Fill(theta);
    }//<-- End if one visible daughter
  }//<-- End if is inelastic

  // Look at the cone angle 
  if (visSec.size() == 1)
  {
    size_t iG4 = visSec[0];

    TVector3 trailMom = primIncMom;

    // If this visible secondary is the primary
    if (iG4 == DetermineG4Id(primTrkId)) trailMom = TVector3( MidPx[0][vertex.point], MidPy[0][vertex.point], MidPz[0][vertex.point] );
    else                                 trailMom = TVector3( Px[iG4], Py[iG4], Pz[iG4] );
    
    auto resultant = primIncMom + trailMom;
    double theta = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(resultant.Unit()));

         if (isElastic)   hMCElasticConeAngle->Fill(theta);
    else if (isInelastic) hMCInelasticConeAngle->Fill(theta);
  }//<-- End if one visible daughter

  // Sanity checks
  if (isElastic && visSec.size() != 1) 
  {
    std::cout << "WARNING: Is elastic but visible secondaries does not make sense!\n";
    std::cout << visSec.size() << " " << jentry << std::endl;
  }
}//<-- End truth studies

/**
 * @brief Method to identify visible secondaries attached to the interaction vertex
 * 
 * @param vertex The interaction vertex
 */
void vertexStudy::IdentifyVisibleSecondaries(const Vertex_t& vertex)
{
  visSec.clear();

  // If the primary doesn't end here
  size_t primiG4 = DetermineG4Id(primTrkId);
  if ((vertex.p+1) < NTrTrajPts[0]) visSec.push_back(primiG4);

  // Fill inelastic scattering angles for all relevant daughters
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    // @todo Take absolute value here?
    if (Mother[iG4] != primTrkId) continue;

    // Make sure she is charged 
    if ( !IsCharged(std::abs(pdg[iG4])) ) continue;

    // This is a charged daughter
    TVector3 dPos0(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]);
    TVector3 dPosf(EndPointx[iG4], EndPointy[iG4], EndPointz[iG4]);

    // Check track length
    // @todo Do we need this cut?
    if ((dPosf-dPos0).Mag() < SECONDARY_LENGTH_CUT) continue;

    // This particle needs to be attached to the vertex
    if ((vertex.position-dPos0).Mag() < 0.01) visSec.push_back(iG4);
  }//<-- End loop over G4 particles
}


/**
 * @brief Area to study IDEs near vertex.
 * 
 */
void vertexStudy::VertexStudy(const Vertex_t& vertex)
{
  bool isElastic(false);
  bool isInelastic(false);
  if (vertex.process.find("hadElastic") != std::string::npos || vertex.process.find("CoulombScat") != std::string::npos) isElastic = true;
  if (vertex.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;
  if (!isElastic && !isInelastic) return;
  if (isElastic && isInelastic) std::cout << "WARNING: Check the process types!\n";

  // Get the trj point here
  TVector3 posHere = vertex.position;

  // First subtract out ides which belong to primary or single secondary
  SubtractIdes(vertex.p);

  if (nTotalEvents == 1)
  {
    for (const auto& ide : ides) 
    {
      auto xBin = hMCIdesSubXvsZ->GetXaxis()->FindBin(ide.pos.Z());
      auto yBin = hMCIdesSubXvsZ->GetYaxis()->FindBin(ide.pos.X());
      hMCIdesSubXvsZ->SetBinContent(xBin, yBin, ide.energy);
    }
  }
  
  // Total energy within some radius
  for (int iBin = 1; iBin <= hMCIdeVertexElastic->GetNbinsX(); iBin++)
  {
    float r = hMCIdeVertexElastic->GetXaxis()->GetBinCenter(iBin);
    float totEnergy = GetIDENearVertex(posHere, r);

    if (isElastic)   hMCIdeVertexElastic  ->Fill(r, totEnergy);
    if (isInelastic) hMCIdeVertexInelastic->Fill(r, totEnergy);

    //if (isElastic) 
    //{
    //  cout << "ENTRY = " << jentry << endl;
    //  cout << MidPosZ[0][vertex.p] << endl;
    //  for (size_t p = 0; p < InteractionPointType->size(); p++) cout << MidPosZ[0][InteractionPoint->at(p)] << " " << InteractionPointType->at(p) << endl;
    //}
  }
}

/**
 * @brief Method to subtract out ides which belong to primary or single secondary.
 * 
 */
void vertexStudy::SubtractIdes(const int& p)
{
  if (visSec.size() != 1) cout << "\nWHAT ARE YOU DOING HERE? " << visSec.size() << "\n";

  // Loop over ides 
  // @note Negative tid means EM activity, -1*tid is mother
  for (auto& tide : ides)
  {
    // If this is the primary or the single secondary set the energy to zero
    if (TrackId[tide.g4Id] == primTrkId || tide.g4Id == visSec[0]) tide.energy = 0;
   
    // Check vicinity to tracks
    CheckVicinity(tide, p);

    // @todo Don't quite understand this, only looking for neutrons or gammas
    //if (tide.trackId < 1 || pdg[tide.g4Id] > 10000) tide.energy = 0;
  }//<-- End loop over ides
}

/**
 * @brief Check vicinity to primary and secondary. 
 * 
 * @param tide The track ide
 * @param p The interaction point
 */
void vertexStudy::CheckVicinity(TrackIde_t& tide, const int& p)
{
  // @note using argoneut numbers
  float CONE_DISTANCE_CUT(2.4);  // cm
  float CONE_ANGLE_CUT(120);     // degrees
  float CYLINDER_DISTANCE_CUT(3); // argoneut used 5cm

  // Get direction vectors
  TVector3 point = tide.pos;
  TVector3 vertex(MidPosX[0][p], MidPosY[0][p], MidPosZ[0][p]);
  TVector3 primDir(MidPosX[0][p-1], MidPosY[0][p-1], MidPosZ[0][p-1]);
  primDir = primDir - vertex;
  primDir = primDir.Unit();
  TVector3 secDir(EndPointx[visSec[0]], EndPointy[visSec[0]], EndPointz[visSec[0]]);
  secDir = secDir - vertex;
  secDir = secDir.Unit();

  // First check for primary
  auto u = point - vertex;
  auto uDotPrim  = u.Dot(primDir);
  auto uCosPrim  = uDotPrim/u.Mag();
  auto uSinPrim  = std::sqrt(1 - uCosPrim*uCosPrim);
  auto uDistPrim = u.Mag() * uSinPrim;

  // check if in cone first
  if (uDotPrim < CONE_DISTANCE_CUT)
  {
    if (10 <= point.Z() && point.Z() < 16) 
    {
      //cout << "A " << uDotPrim << " " << CONE_DISTANCE_CUT << " " << uDistPrim << endl;
      //u.Print();
      //point.Print();
      //vertex.Print();
      //primDir.Print();
    }
    float coneDist = uDotPrim * std::tan( 0.5*CONE_ANGLE_CUT*TMath::Pi()/180 );
    if (uDistPrim < coneDist) 
    {
      tide.energy = 0;
    }
  }
  // check if in cylinder
  else if (uDistPrim < CYLINDER_DISTANCE_CUT) { tide.energy = 0;}

  // Now check for secondary
  auto uDotSec   = u.Dot(secDir);
  auto uCosSec   = uDotSec/u.Mag();
  auto uSinSec   = std::sqrt(1 - uCosSec*uCosSec);
  auto uDistSec  = u.Mag() * uSinSec;
  
  // check if in cone first
  if (uDotSec < CONE_DISTANCE_CUT)
  {
    float coneDist = uDotSec * std::tan( 0.5*CONE_ANGLE_CUT*TMath::Pi()/180 );
    if (uDistSec < coneDist) tide.energy = 0;
  }
  // check if in cylinder
  else if (uDistSec < CYLINDER_DISTANCE_CUT) { tide.energy = 0; }
  //if (tide.energy > 0 && uDistPrim < 5) cout << uDistPrim << " " << tide.energy << " " << tide.pos.X() << " " << tide.pos.Y() << " " << tide.pos.Z() << endl;
}

/**
 * @brief Method to get the energy deposited inside bubble near vertex.
 * 
 * Subtracts out all ides from primary particle and single secondary. The remaining
 * energy should be from neutrons and/or dexcitation gammas. 
 * 
 * @note We should've set other ides equal to zero
 * @todo Should we consider a track length cut? 
 * 
 * @param vertex Interaction vertex
 * @param radius Bubble radius 
 * @return float Total energy deposited within bubble radius of interaction vertex
 */
float vertexStudy::GetIDENearVertex(const TVector3& vertex, const float& radius)
{
  float totalEnergy(0);
  for (const auto& tide : ides) 
  {
    if ( (tide.pos-vertex).Mag() < radius ) totalEnergy += tide.energy;
  }
  return totalEnergy;
}


/**
 * @brief Make plots and save to output file
 * 
 */
void MakePlots()
{
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

  myRootFile.cd();
  myRootFile.mkdir("mc/");
  myRootFile.cd("mc/");
  hMCElasticAngle->Write();
  hMCInelasticAngle->Write();
  hMCInelasticOneVisDAngle->Write();
  hMCSecondaries->Write();
  hMCFirstInTpcPointX->Write();
  hMCFirstInTpcPointY->Write();
  hMCFirstInTpcPointZ->Write();
  hMCLastInTpcPointX->Write();
  hMCLastInTpcPointY->Write();
  hMCLastInTpcPointZ->Write();
  hMCLastPosFirstIntPosDiff->Write();
  hMCLastPosZ->Write();
  hMCSecondaryTrkLength->Write();
  hMCIdeVertexElastic->Write();
  hMCIdeVertexInelastic->Write();
  hMCIdesXvsZ->Write();
  hMCIdesSubXvsZ->Write();
  hMCNumInteractions->Write();
  hMCInelasticConeAngle->Write();
  hMCElasticConeAngle->Write();

  myRootFile.Close();
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

/**
 * @brief Convert position to primary point
 * 
 * @param position The position
 * @return int The point id of primary
 */
int vertexStudy::ConvertToPrimaryPoint(const TVector3& position)
{
  for (int iPt = 0; iPt < NTrTrajPts[0]; iPt++)
  {
    TVector3 testPosition(MidPosX[0][iPt], MidPosY[0][iPt], MidPosZ[0][iPt]);

    if ( (testPosition-position).Mag() < 0.01 ) return iPt;
  }
  return 0;
}



