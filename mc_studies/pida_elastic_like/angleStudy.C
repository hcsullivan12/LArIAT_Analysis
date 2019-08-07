/**
 * @file angleStudy.C
 * @brief Angle study
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#define angleStudy_cxx
#include "angleStudy.h"
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
std::string IntProcessToString(const int& p);
std::string GetParticle(int pdg);

/// Mass of pion in MeV
float PARTICLE_MASS(139.57); 
inline double ToKineticEnergy(const TVector3& mom){return std::sqrt( mom.Mag()*mom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; }
/// @}

/// @name Useful variables
/// @{
/// Counters
size_t _nTotalEvents(0),  _nPrimariesEntered(0), 
       _nGoodMCEvents(0), _nEventsInelastic(0),
       _nElasticLikeEvents(0), _nRangedOut(0),
       _nElasticRangedOut(0), _nInelasticRangedOut(0);

/// Primary track Id
int _primTrkId(-1);

/// Vector of g4 Ids for visible secondaries
std::vector<size_t> _visSec;

/// Container for pdg codes for charged particles (cuz I can't remember) 
struct PdgCodes_t
{
  int kElectron = 11;
  int kMuon     = 13;
  int kPion     = 211;
  int kKaon     = 321;
  int kProton   = 2212;
};
PdgCodes_t _pdgcodes;
/// Method to check if charged
inline bool IsCharged(int const& p) 
{ 
  int pdg = std::abs(p);
  return (pdg == _pdgcodes.kElectron || pdg == _pdgcodes.kMuon || 
          pdg == _pdgcodes.kPion     || pdg == _pdgcodes.kKaon || 
          pdg == _pdgcodes.kProton);
}

/// unknown processes
std::set<std::string> _processes;

/// vector of _vertices
std::vector<Vertex_t> _vertices;
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
TH1D* hMCPrimaryMissedTpcX = new TH1D("hMCPrimaryMissedTpcX", "MC Primary Missed TPC X", 200, -50, 50);
TH1D* hMCPrimaryMissedTpcY = new TH1D("hMCPrimaryMissedTpcY", "MC Primary Missed TPC Y", 200, -50, 50);
TH1D* hMCPrimaryMissedTpcZ = new TH1D("hMCPrimaryMissedTpcZ", "MC Primary Missed TPC Z", 200, -110, 10);
TH1D *hMCPrimaryProjX0 = new TH1D("hMCPrimaryProjX0", "Primary Particle X_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjY0 = new TH1D("hMCPrimaryProjY0", "Primary Particle Y_{0}", 200, -50 , 50);
TH1D *hMCPrimaryProjZ0 = new TH1D("hMCPrimaryProjZ0", "Primary Particle Z_{0}", 100, -5 , 5);
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
TH1I* hMCNumInteractions = new TH1I("hMCNumInteractions", "Number of interactions in TPC for elastic like", 5, 0, 5);
TH1D* hMCFirstInteractionInTpcZ = new TH1D("hMCFirstInteractionInTpcZ", "First interaction Z point in TPC", 100, 0, 100);

/// @}

/// Output root file
TFile myRootFile("piMinusAna.root", "RECREATE");

Long64_t jentry(0);

/**
 * @brief Main loop
 * 
 * @param inDebug Option to run over only a few events
 */
void angleStudy::Loop(int inDebug)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (jentry=0; jentry<=nentries;jentry++) 
  {
    // If debug, only look at a sub sample 
    //if (inDebug == 1 && jentry%1000 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << jentry << endl;
    IN_DEBUG = inDebug;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // Increment our total event counter 
    _nTotalEvents++; 
    if (_nTotalEvents%1000 == 0) std::cout << "EVENT = " << _nTotalEvents << std::endl; 

    // TEMP
    //if (event == 1) cout << jentry << endl;
    if(IN_DEBUG) std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    //std::cout << "Event " << event << "\n\n";
    //if (_nTotalEvents > 500) break;


    // Filling primary information 
    _primTrkId = -1;
    int nPrimaries = 0;
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (TrackId->at(iG4) == 0) std::cout << "TRACK ID = 0\n";
      // If this is not a primary, skip it
      if (process_primary->at(iG4) == 0) continue;
      _primTrkId = TrackId->at(iG4);
      nPrimaries++;

      // Store the positions and momentum 
      TVector3 g4PrimaryPos0( StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4) );
      TVector3 g4PrimaryPosf( EndPointx  ->at(iG4), EndPointy  ->at(iG4), EndPointz  ->at(iG4) );
      TVector3 g4PrimaryMom0(1000*StartPx->at(iG4), 1000*StartPy->at(iG4), 1000*StartPz->at(iG4));
      TVector3 g4PrimaryMomf(1000*EndPx->at(iG4), 1000*EndPy->at(iG4), 1000*EndPz->at(iG4));

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
    if (_primTrkId != 1) std::cout << "WARNING! Check the primary track ID!\n";
    if (nPrimaries > 1) std::cout << "WARNING! Check the number of primaries!\n";
    
    // Check if entered TPC
    bool isGoodEvent(false);
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // If this is not a primary, skip it
      if (process_primary->at(iG4) == 0) continue;
      
      // only look if a primary entered the tpc
      if ( EndPointz->at(iG4) > 0 ) isGoodEvent = true;
      if ( !isGoodEvent ) 
      {
        hMCPrimaryMissedTpcX->Fill( EndPointx->at(iG4) );
        hMCPrimaryMissedTpcY->Fill( EndPointy->at(iG4) );
        hMCPrimaryMissedTpcZ->Fill( EndPointz->at(iG4) );
        continue;
      }

      _nPrimariesEntered++;
    }//<--- End loop over primaries
    if (!isGoodEvent) continue;
    _nGoodMCEvents++;
  

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin study

    if (IN_DEBUG) std::cout << "\nFilling track lengths...";
    // Make distribution of secondary track lengths
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // Only look at daughters of primary
      if (Mother->at(iG4) != _primTrkId) continue;
  
      // Make sure she's charged
      if ( !IsCharged(PDG->at(iG4)) ) continue;
      TVector3 startPoint(StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4)); 
      TVector3 endPoint(EndPointx->at(iG4), EndPointy->at(iG4), EndPointz->at(iG4));
      hMCSecondaryTrkLength->Fill((startPoint-endPoint).Mag());
    }

    // Get interactions in TPC
    if (IN_DEBUG) std::cout << "Getting interactions in TPC...\n";
    _vertices.clear();
    GetInteractionsInTpc();

    // Get first interaction
    Vertex_t firstVertexInTpc(-1, "none", TVector3(0,0,200));
    if (_vertices.size()) firstVertexInTpc = _vertices[0];
    if (firstVertexInTpc.process.find("none") != std::string::npos) continue;
    hMCFirstInteractionInTpcZ->Fill(firstVertexInTpc.position.Z());

    // Identify visible secondaries
    _visSec.clear();
    IdentifyVisibleSecondaries(firstVertexInTpc);

    // Angle study
    if (IN_DEBUG) std::cout << "Starting angle study...\n";
    bool isElasticLike(false);
    AngleStudy(firstVertexInTpc, isElasticLike);

    if (isElasticLike) _nElasticLikeEvents++;
    if (IN_DEBUG) 
    {
      std::cout << "\nEvent                 = " << event 
                << "\nFirst int point in TPC = " << MidPosZ->at(0)[firstVertexInTpc.point]
                << "\nFirst process in TPC   = " << firstVertexInTpc.process
                << "\nVisible secondaries    = " << _visSec.size() 
                << "\nIs elastic like        = ";
      if (isElasticLike) std::cout << "yes" << std::endl;
      else               std::cout << "no"  << std::endl;
      for (const auto& p : *InteractionPointType) std::cout << p << endl;
    }

    // If this wasn't elastic like, we're done
    if (!isElasticLike)continue;

    // Do pid study
    PidStudy(firstVertexInTpc);

// End study
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  }//<---End loop over entries

  // Event reduction table
  std::cout << "\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"
            << "Events simulated:                           " << _nTotalEvents                 << endl
            << "Tracks entered TPC:                         " << _nPrimariesEntered            << endl
            << "Events with > 0 entered tracks:             " << _nGoodMCEvents                << endl
            << "Elastic like events:                        " << _nElasticLikeEvents           << endl
            << "Ranged out:                                 " << _nRangedOut                   << endl
            << "Elastic ranged out:                         " << _nElasticRangedOut            << endl
            << "Inelastic ranged out:                       " << _nInelasticRangedOut          << endl
            << endl;

  // Make plots
  MakePlots();
  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Method to get the interactions in tpc
 * 
 * @note When an elastic shows signs of > 0 ide, it's most likely
 *       the case that there is a subsequent inelastic interaction.
 */
void angleStudy::GetInteractionsInTpc()
{
  // If the interaction container is empty, check the daughters
  if (!InteractionPoint->size())
  {
    if (IN_DEBUG) std::cout << "\nChecking Daughters!!!\n";
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (Mother->at(iG4) != _primTrkId) continue;
      // Check the process
      std::string proc = Process->at(iG4);
      // We skip elastic and inelastic
      if (proc.find("hadElastic") != std::string::npos ||
          proc.find("pi-Inelastic") != std::string::npos ||
          proc.find("hIon") != std::string::npos ||
          proc.find("CoulombScat") != std::string::npos) continue;

      // Should be decay or capture
      _processes.insert(proc);

      TVector3 pos(StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4));
      Vertex_t thisVertex( ConvertToPrimaryPoint(pos),
                           proc,
                           pos );
      // Check if in FV
      if (!InActiveRegion(thisVertex.position)) continue;
      
      // This should be it!
      _vertices.push_back(thisVertex);
      break;
    }
  }
  // G4 gave us something
  else
  {
    for (int iPt = 0; iPt < InteractionPoint->size(); iPt++)
    {
      Vertex_t thisVertex(InteractionPoint->at(iPt),
                          InteractionPointType->at(iPt),
                          TVector3(MidPosX->at(0)[InteractionPoint->at(iPt)], MidPosY->at(0)[InteractionPoint->at(iPt)], MidPosZ->at(0)[InteractionPoint->at(iPt)]));
      // Check if in FV
      if (!InActiveRegion(thisVertex.position))continue;

      bool weGotIt(false);
      for (const auto &v : _vertices)
      {
        float diff = (v.position - thisVertex.position).Mag();
        if (diff < 0.01 && v.process == thisVertex.process)weGotIt = true;
      }
      if (!weGotIt)_vertices.push_back(thisVertex);
    }
  }
  std::sort(_vertices.begin(), _vertices.end(), [](const auto& l, const auto& r) {return l.position.Z() < r.position.Z();});
}

/**
 * @brief Method to get G4 Id from track id
 * 
 * @param tid Track id
 * @return int The particle's G4 Id
 */
int angleStudy::DetermineG4Id(const int& tid)
{
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    if (TrackId->at(iG4 == tid)) return iG4;
  }
  return -999;
} 

/**
 * @brief Method to study angles for elastic like events
 * 
 * @param vertex The first vertex in the tpc
 * @param isElasticLike Is this process elastic like (one visible secondary)
 */
void angleStudy::AngleStudy(const Vertex_t& vertex, bool& isElasticLike)
{
  // Handle the case of just the primary
  bool isElastic(false);
  bool isInelastic(false);
  if (vertex.process.find("hadElastic")   != std::string::npos || vertex.process.find("CoulombScat") != std::string::npos) isElastic = true;
  if (vertex.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;

  // Only looking at elastic or inelastic
  if (!isElastic && !isInelastic) return;

  // Get incident momentum
  TVector3 primIncMom( MidPx->at(0)[vertex.point-1], MidPy->at(0)[vertex.point-1], MidPz->at(0)[vertex.point-1] );

  // If the primary doesn't end here
  if ( (vertex.point+1) < (NTrTrajPts->at(0)-1) ) 
  {
    TVector3 trailMom( MidPx->at(0)[vertex.point], MidPy->at(0)[vertex.point], MidPz->at(0)[vertex.point] );
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
    for (const auto& iG4 : _visSec)
    {
      if (iG4 == DetermineG4Id(_primTrkId)) continue;

      // This is a charged daughter
      TVector3 dMom0(StartPx->at(iG4), StartPy->at(iG4), StartPz->at(iG4));

      // Compute angle between daughter and incoming primary
      double theta = (180/TMath::Pi())*std::acos( primIncMom.Unit().Dot(dMom0.Unit()) );
      hMCInelasticAngle->Fill(theta);
    }//<-- End loop over visible particles
  }//<--End if inelastic

  // Do it again, except for cases in which inelastic but one visible daughter
  if (isInelastic)
  {
    // If _visSec.size() == 1 this looks elastic
    hMCSecondaries->Fill(_visSec.size());
    if (_visSec.size() == 1)
    {
      size_t iG4 = _visSec[0];
      isElasticLike = true;
      TVector3 trailMom = primIncMom;

      // If this visible secondary is the primary
      if (iG4 == DetermineG4Id(_primTrkId)) trailMom = TVector3( MidPx->at(0)[vertex.point], MidPy->at(0)[vertex.point], MidPz->at(0)[vertex.point] );
      else                                 trailMom = TVector3( StartPx->at(iG4), StartPy->at(iG4), StartPz->at(iG4) );
    
      double theta = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(trailMom.Unit()));
      hMCInelasticOneVisDAngle->Fill(theta);
    }//<-- End if one visible daughter
  }//<-- End if is inelastic

  // Look at the cone angle 
  if (_visSec.size() == 1)
  {
    size_t iG4 = _visSec[0];

    TVector3 trailMom = primIncMom;

    // If this visible secondary is the primary
    if (iG4 == DetermineG4Id(_primTrkId)) trailMom = TVector3( MidPx->at(0)[vertex.point], MidPy->at(0)[vertex.point], MidPz->at(0)[vertex.point] );
    else                                 trailMom = TVector3( StartPx->at(iG4), StartPy->at(iG4), StartPz->at(iG4) );
    
    auto resultant = primIncMom + trailMom;
    double theta = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(resultant.Unit()));

         if (isElastic)   hMCElasticConeAngle->Fill(theta);
    else if (isInelastic) hMCInelasticConeAngle->Fill(theta);
  }//<-- End if one visible daughter

  // Sanity checks
  if (isElastic && _visSec.size() != 1) 
  {
    std::cout << "WARNING: Is elastic but visible secondaries does not make sense!\n";
    std::cout << _visSec.size() << " " << jentry << std::endl;
  }
}//<-- End truth studies


/**
 * @brief PID study of visible secondary
 * 
 * @param vertex The interaction vertex
 */
void angleStudy::PidStudy(const Vertex_t& vertex)
{
  if (_visSec.size() != 1) std::cout << "WHY ARE YOU HERE!\n";
  int secG4id = _visSec[0];

  std::cout << "\nEvent " << event << std::endl;
  std::cout << "Single daughter is " << abs(PDG->at(secG4id)) << " " << GetParticle(std::abs(PDG->at(secG4id))) << std::endl;
  // Let's see how many of these range out in detector
  auto endPos = TVector3(EndPointx->at(secG4id), EndPointy->at(secG4id), EndPointz->at(secG4id));
  
  // Make sure this isn't a pion
  if (InActiveRegion(endPos) && std::abs(PDG->at(secG4Id)) != kPion)
  {
    std::cout << "RANGED OUT!\n";
    _nRangedOut++;
    if (vertex.process == "hadElastic") _nElasticRangedOut++;
    if (vertex.process == "pi-Inelastic") _nInelasticRangedOut++;
  }
  if (InActiveRegion(endPos) && std::abs(PDG->at(secG4Id)) == kPion)
  {
    // Make sure the pion doesn't undergo another interaction in tpc
  }
}

/**
 * @brief Convert pdg to particle name
 * 
 * @param pdg 
 * @return std::string 
 */
std::string GetParticle(int pdg)
{
  pdg = std::abs(pdg);
  std::string type = "not charged";
  switch (pdg)
  {
    case 11:   {type = "electron";break;}
    case 13:   {type = "muon"    ;break;}
    case 211:  {type = "pion"    ;break;}
    case 321:  {type = "kaon"    ;break;}
    case 2212: {type = "proton"  ;break;}
  }
  return type;
}

/**
 * @brief Method to identify visible secondaries attached to the interaction vertex
 * 
 * @param vertex The interaction vertex
 */
void angleStudy::IdentifyVisibleSecondaries(const Vertex_t& vertex)
{
  _visSec.clear();

  // If the primary doesn't end here
  size_t primiG4 = DetermineG4Id(_primTrkId);
  if ( (vertex.point+1) < (NTrTrajPts->at(0)-1) ) _visSec.push_back(primiG4);

  // Fill inelastic scattering angles for all relevant daughters
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    // @todo Take absolute value here?
    if (Mother->at(iG4) != _primTrkId) continue;

    // Make sure she is charged 
    if ( !IsCharged(PDG->at(iG4)) ) continue;

    // This is a charged daughter
    TVector3 dPos0(StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4));
    TVector3 dPosf(EndPointx  ->at(iG4), EndPointy  ->at(iG4), EndPointz  ->at(iG4));
    double length = (dPosf-dPos0).Mag();
    if (IN_DEBUG) std::cout << "Particle " << GetParticle(PDG->at(iG4)) << "  at " << dPos0.X() << " " << dPos0.Y() << " " << dPos0.Z() << "  Daughter length = " << length << std::endl;

    // Check track length
    // @todo Do we need this cut?
    if (length < SECONDARY_LENGTH_CUT) continue;

    // This particle needs to be attached to the vertex
    if ((vertex.position-dPos0).Mag() < 0.01) _visSec.push_back(iG4);
  }//<-- End loop over G4 particles
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
  hMCPrimaryMissedTpcX->Write();  
  hMCPrimaryMissedTpcY->Write();  
  hMCPrimaryMissedTpcZ->Write();  
  hMCPrimaryProjX0->Write();
  hMCPrimaryProjY0->Write();
  hMCPrimaryProjZ0->Write();
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
int angleStudy::ConvertToPrimaryPoint(const TVector3& position)
{
  for (int iPt = 0; iPt < NTrTrajPts->at(0); iPt++)
  {
    TVector3 testPosition(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);

    if ( (testPosition-position).Mag() < 0.01 ) return iPt;
  }
  return 0;
}



