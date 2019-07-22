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
struct PdgCodes
{
  int kElectron = 11;
  int kMuon     = 13;
  int kPion     = 211;
  int kKaon     = 321;
  int kProton   = 2212;
};
PdgCodes pdgcodes;
/// Method to check if charged
inline bool IsCharged(int const& pdg) 
{ 
  return (pdg == pdgcodes.kElectron || pdg == pdgcodes.kMuon || 
          pdg == pdgcodes.kPion     || pdg == pdgcodes.kKaon || 
          pdg == pdgcodes.kProton);
}

/// Stucture for IDE
struct TrackIde
{
  int trackId;
  int g4Id;
  float energy;
  TVector3 pos;
};
/// vector of TrackIdes 
std::vector<TrackIde> ides;
/// @}

/// @name Cuts and constants
/// @{
int IN_DEBUG(0);

/// TPC boundaries
float TPC_X_BOUND[2] = {   0.0, 47.0 };
float TPC_Y_BOUND[2] = { -20.0, 20.0 };
float TPC_Z_BOUND[2] = {   0.0, 90.0 };

/// Fiducial volume definition
float FV_X_BOUND[2] = {   2.0, 45.0 };
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
TH1D* hMCSecondaries      = new TH1D("hMCSecondaries",     "True number of tracks leaving vertex", 10, 0, 10);
TH1D* hMCFirstInTpcPointX = new TH1D("hMCFirstInTpcPointX", "MC First Point in TPC X", 400, TPC_X_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointY = new TH1D("hMCFirstInTpcPointY", "MC First Point in TPC Y", 400, TPC_Y_BOUND[0]-20, 20);
TH1D* hMCFirstInTpcPointZ = new TH1D("hMCFirstInTpcPointZ", "MC First Point in TPC Z", 400, TPC_Z_BOUND[0]-20, 20);
TH1D* hMCLastInTpcPointX = new TH1D("hMCLastInTpcPointX", "MC Last Point in TPC X", 1300, TPC_X_BOUND[0]-20, TPC_X_BOUND[1]+20);
TH1D* hMCLastInTpcPointY = new TH1D("hMCLastInTpcPointY", "MC Last Point in TPC Y", 1300, TPC_Y_BOUND[0]-20, TPC_Y_BOUND[1]+20);
TH1D* hMCLastInTpcPointZ = new TH1D("hMCLastInTpcPointZ", "MC Last Point in TPC Z", 1300, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCSecondaryTrkLength = new TH1D("hMCSecondaryTrkLength", "Secondary Track Lengths", 100, 0, 50);
TH1D* hMCIdeVertexElastic = new TH1D("hMCIdeVertexElastic", "IDEs versus bubble radius elastic", 100, 0, 10);
TH1D* hMCIdeVertexInelastic = new TH1D("hMCIdeVertexInelastic", "IDEs versus bubble radius inelastic", 100, 0, 10);
/// @}

/// Output root file
TFile myRootFile("piMinusAna.root", "RECREATE");


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
  for (Long64_t jentry=0/*47000*/; jentry<100000/*nentries*/;jentry++) 
  {
    // If debug, only look at a sub sample 
    if (inDebug == 1 && jentry%1000 != 0) continue;
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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin study

    // Make distribution of secondary track lengths
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // Only look at daughters of primary
      if (Mother[iG4] != primTrkId) continue;
      // Make sure she's charged
      if ( !IsCharged(std::abs(pdg[iG4])) ) continue;

      TVector3 startPoint(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]);
      TVector3 endPoint(EndPointx[iG4], EndPointy[iG4], EndPointz[iG4]);
      hMCSecondaryTrkLength->Fill((startPoint-endPoint).Mag());
    }

    // Get first interaction in TPC
    int firstProcessPointInTpc(NTrTrajPts[0]);
    std::string firstProcessInTpc("none");
    GetFirstInteractionInTPC(firstProcessPointInTpc, firstProcessInTpc);

    // Angle study
    bool isElasticLike(false);
    visSec.clear();
    AngleStudy(firstProcessPointInTpc, firstProcessInTpc, isElasticLike);
    if (isElasticLike) nElasticLikeEvents++;

    // Vertex study for elastic like events
    if (isElasticLike) 
    {
      ides.clear();
      ides.reserve(maxTrackIDE);
      for (int iIDE = 0; iIDE < maxTrackIDE; iIDE++) 
      {
        TrackIde tide;
        tide.trackId = IDETrackId[iIDE];
        tide.energy  = IDEEnergy[iIDE];
        tide.pos     = TVector3( IDEPos[iIDE][0], IDEPos[iIDE][1], IDEPos[iIDE][2] );
        tide.g4Id    = DetermineG4Id(tide.trackId);
        ides.push_back(tide);
      }
      VertexStudy();
    }

// End study
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  }//<---End loop over entries


  // Event reduction table
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events simulated:                           " << nTotalEvents                 << endl
            << "Tracks entered TPC:                         " << nPrimariesEntered            << endl
            << "Events with > 0 entered tracks:             " << nGoodMCEvents                << endl
            << "Elastic like events:                        " << nElasticLikeEvents++         << endl
            << "Inelastic interactions:                     " << nEventsInelastic             << endl
            << endl;

  // Make plots
  MakePlots();

  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Method to get the first interaction in tpc
 * 
 * @param point The primary point
 * @param process The process type
 */
void vertexStudy::GetFirstInteractionInTPC(int& point, std::string& process)
{
  for (size_t iPr = 0; iPr < (*InteractionPoint).size(); iPr++)
  {
    // The point and process
    auto p    = (*InteractionPoint)[iPr];
    auto proc = (*InteractionPointType)[iPr];

    // Get the trj point here and before here
    TVector3 posBefore( MidPosX[0][p-1], MidPosY[0][p-1], MidPosZ[0][p-1] );
    TVector3 posHere  ( MidPosX[0][p],   MidPosY[0][p],   MidPosZ[0][p] );
    auto primIncDir = (posHere-posBefore).Unit();
    if (!InActiveRegion(posHere)) continue;

    point   = p;
    process = proc;
    return;
  }
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
} 

/**
 * @brief Method to study angles for elastic like events
 * 
 * @param p Interaction point
 * @param proc Interaction process type
 * @param isElasticLike Is this process elastic like (one visible secondary)
 */
void vertexStudy::AngleStudy(const int& p, const std::string& proc, bool& isElasticLike)
{
  // Make the distribution of scattering angles
  // We need direction vectors before interaction points
  size_t primiG4(0);
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++) {if (process_primary[iG4]==1) primiG4=iG4;}

  // Handle the case of just the primary
  bool isElastic(false);
  bool isInelastic(false);
  if (proc.find("hadElastic") != std::string::npos || proc.find("CoulombScat") != std::string::npos) isElastic = true;
  if (proc.find("pi-Inelastic") != std::string::npos) isInelastic = true;

  // Only looking at elastic or inelastic
  if (!isElastic && !isInelastic) return;
  if (isElastic && isInelastic) std::cout << "WARNING: Check the process types!\n";

  // Get the trj point here and before here
  TVector3 posBefore( MidPosX[0][p-1], MidPosY[0][p-1], MidPosZ[0][p-1] );
  TVector3 posHere  ( MidPosX[0][p],   MidPosY[0][p],   MidPosZ[0][p] );
  auto primIncDir = (posHere-posBefore).Unit();

  // After
  if ( (p+1) < NTrTrajPts[0]) 
  {
    TVector3 posAfter( MidPosX[0][p+1], MidPosY[0][p+1], MidPosZ[0][p+1] );
    auto primTrailDir = posAfter-posHere;
    double theta = (180/TMath::Pi())*std::acos( primIncDir.Unit().Dot(primTrailDir.Unit()) );
    
    if (isElastic)   
    {
      isElasticLike = true;
      hMCElasticAngle->Fill(theta);
    }

    // In the inelastic case, add this to the visible daughters list
    // since the primary still lives
    // @todo Add this case to the track length distribution
    if (isInelastic) 
    {
      hMCInelasticAngle->Fill(theta); 
      visSec.push_back(primiG4);
    }
  }
  
  // Handle inelastic type with all relevant daughters
  if (isInelastic)
  {
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (Mother[iG4] != primTrkId) continue;
      // Make sure she is charged 
      if ( !IsCharged(std::abs(pdg[iG4])) ) continue;

      // This is a charged daughter
      TVector3 dPos0(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]);
      TVector3 dPosf(EndPointx[iG4],   EndPointy[iG4],   EndPointz[iG4]);
      TVector3 dMom0(Px[iG4], Py[iG4], Pz[iG4]);

      // Check track length
      if ((dPosf-dPos0).Mag() < SECONDARY_LENGTH_CUT) continue;

      // This needs to be attached to this vertex
      if ((posHere-dPos0).Mag() < 0.01)
      {
        visSec.push_back(iG4);
        // Compute angle between daughter and incoming primary
        double theta = (180/TMath::Pi())*std::acos( primIncDir.Unit().Dot(dMom0.Unit()) );
        hMCInelasticAngle->Fill(theta);
      }
    }//<-- End loop over G4 particles

    // Do it again, except for cases in which inelastic but one visible daughter
    // If visSec.size() == 1 this looks elastic
    hMCSecondaries->Fill(visSec.size());
    if (visSec.size() == 1)
    {
      size_t iG4 = visSec[0];
      isElasticLike = true;

      // If this visible secondary is the primary
      if (iG4 == primiG4)
      {
        TVector3 posAfter(MidPosX[0][p + 1], MidPosY[0][p + 1], MidPosZ[0][p + 1]);
        auto primTrailDir = posAfter - posHere;
        double theta = (180 / TMath::Pi()) * std::acos(primIncDir.Unit().Dot(primTrailDir.Unit()));
        hMCInelasticOneVisDAngle->Fill(theta);
      }
      else
      {
        TVector3 dPos0(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]);
        TVector3 dPosf(EndPointx[iG4], EndPointy[iG4], EndPointz[iG4]);
        TVector3 dMom0(Px[iG4], Py[iG4], Pz[iG4]);
      
        // Compute angle between daughter and incoming primary
        double theta = (180 / TMath::Pi()) * std::acos(primIncDir.Unit().Dot(dMom0.Unit()));
        hMCInelasticOneVisDAngle->Fill(theta);
      }
    }//<-- End if one visible daughter
  }//<-- End if is inelastic

}//<-- End truth studies


/**
 * @brief Area to study IDEs near vertex.
 * 
 * Uses the first vertex in the TPC for the study. Currently looking 
 * at IDEs within some radius of vertex.
 * 
 */
void vertexStudy::VertexStudy()
{
  for (size_t iPr = 0; iPr < (*InteractionPoint).size(); iPr++)
  {
    // The point and process
    auto p    = (*InteractionPoint)[iPr];
    auto proc = (*InteractionPointType)[iPr];

    bool isElastic(false);
    bool isInelastic(false);
    if (proc.find("hadElastic") != std::string::npos || proc.find("CoulombScat") != std::string::npos) isElastic = true;
    if (proc.find("pi-Inelastic") != std::string::npos) isInelastic = true;
    if (!isElastic && !isInelastic) continue;
    if (isElastic && isInelastic) std::cout << "WARNING: Check the process types!\n";

    // Get the trj point here
    TVector3 posHere( MidPosX[0][p], MidPosY[0][p], MidPosZ[0][p] );
    if (!InActiveRegion(posHere)) continue;

    // Total energy within some radius
    for (int iBin = 1; iBin <= hMCIdeVertexElastic->GetNbinsX(); iBin++)
    {
      float r = hMCIdeVertexElastic->GetBinCenter(iBin);
      float totEnergy = GetIDENearVertex(posHere, r);
      if (isElastic)   hMCIdeVertexElastic->SetBinContent(iBin, totEnergy);
      if (isInelastic) hMCIdeVertexInelastic->SetBinContent(iBin, totEnergy);
    }
    break;
  }
}

/**
 * @brief Method to get the energy deposited inside bubble near vertex.
 * 
 * Subtracts out all ides from visible particle (e.g. charged). The remaining
 * energy is from neutrons and/or dexcitation gammas. 
 * 
 * @todo Should we consider a track length cut? 
 * 
 * @param vertex 
 * @param radius 
 * @return float 
 */
float vertexStudy::GetIDENearVertex(const TVector3& vertex, const float& radius)
{
  // Loop over ides and subtract out those which belong to all visible particles
  for (auto& tide : ides)
  {
    if ( !IsCharged( pdg[tide.g4Id] ) ) continue;


  }

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



