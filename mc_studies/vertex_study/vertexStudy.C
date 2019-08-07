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
void makePlots();
bool inActiveRegion(const TVector3& thePos);
bool inTPCRegion(const TVector3& thePos);
inline void printVec(const TVector3& pos) {std::cout<<"\t("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")\n";};
void writeIdes(int stage);

/// Mass of pion in MeV
float PARTICLE_MASS(139.57); 
inline double toKineticEnergy(const TVector3& mom){return std::sqrt( mom.Mag()*mom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; }
/// @}

/// @name Useful variables
/// @{
/// Counters
size_t _nTotalEvents(0),  _nPrimariesEntered(0), 
       _nGoodMCEvents(0), _nElasticLikeEvents(0);

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
PdgCodes_t pdgcodes;
/// Method to check if charged
inline bool isCharged(int const& p) 
{ 
  int pdg = std::abs(p);
  return (pdg == pdgcodes.kElectron || pdg == pdgcodes.kMuon || 
          pdg == pdgcodes.kPion     || pdg == pdgcodes.kKaon || 
          pdg == pdgcodes.kProton);
}

/// vector of TrackIdes 
std::vector<TrackIde_t> _ides;

/// vector of _vertices
std::vector<Vertex_t> _vertices;

/// set of unknown processes
std::set<std::string> _processes;
/// @}

/// @name Cuts and constants
/// @{
int IN_DEBUG(0);

int EVENT_FOR_EVT_DSP = 999;

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
TH1D *hMCPrimaryPx = new TH1D("hMCPrimaryPx", "Primary Particle P_{x}", 22, 0, 1100);
TH1D *hMCPrimaryPy = new TH1D("hMCPrimaryPy", "Primary Particle P_{y}", 22, 0, 1100);
TH1D *hMCPrimaryPz = new TH1D("hMCPrimaryPz", "Primary Particle P_{z}", 22, 0 , 1100);
TH1D *hMCPrimaryP  = new TH1D("hMCPrimaryP", "Primary Particle P", 22, 0, 1100);
TH1D *hTrueLength = new TH1D("hTrueLength", "#True Length of the Primary Particle inside the TPC", 200, 0 , 100);
TH1D *hMostUpstreamZPos = new TH1D("hMostUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);

/// mc
TH1D* hMCLastPosFirstIntPosDiff = new TH1D("hMCLastPosFirstIntPosDiff", "MC Last Position - First Interacting Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
TH1D* hMCLastPosZ = new TH1D("hMCLastPosZ", "MC Last Position", 1000, TPC_Z_BOUND[0]-20, TPC_Z_BOUND[1]+20);
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
TH2D* hMCXvsZ = new TH2D("hMCXvsZ", "X vs Z ", 200, 0, 100, 100, 0, 50);
TH1D* hMCFirstInteractionInTpcZ = new TH1D("hMCFirstInteractionInTpcZ", "First interaction Z point in TPC", 100, 0, 100);

/// @}

/// Output root file
TFile myRootFile("results.root", "RECREATE");

Long64_t _jentry(0);

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
  for (_jentry=0; _jentry<=nentries;_jentry++) 
  {
    // If debug, only look at a sub sample 
    //if (inDebug == 1 && _jentry%1000 != 0) continue;
    if (inDebug == 1) cout << "InDubug: " << _jentry << endl;
    IN_DEBUG = inDebug;

    Long64_t ientry = LoadTree(_jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(_jentry);   nbytes += nb;
    // Increment our total event counter 
    _nTotalEvents++; 
    if (_nTotalEvents%1000 == 0) std::cout << "EVENT = " << _nTotalEvents << std::endl; 

    // TEMP
    //if (event == 1) cout << _jentry << endl;
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

      // ### Store the positions and momentum 
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
      // calculate energy loss
      float energyLoss(0);
      // loop over trj points for this primary
      for (size_t iPoint = 0; iPoint < NTrTrajPts->at(iG4); iPoint++)
      {
        // only look at the upstream portion
        if ( MidPosZ->at(iG4)[iPoint] > 0 ) break;

        // ignore last point
        if ( (iPoint+1) >= NTrTrajPts->at(iG4) ) break;
        TVector3 mom1Vec( MidPx->at(iG4)[iPoint],   MidPy->at(iG4)[iPoint],   MidPz->at(iG4)[iPoint]  );
        TVector3 mom2Vec( MidPx->at(iG4)[iPoint+1], MidPy->at(iG4)[iPoint+1], MidPz->at(iG4)[iPoint+1]);

        float energy1 = std::sqrt( mom1Vec.Mag()*mom1Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS );
        float energy2 = std::sqrt( mom2Vec.Mag()*mom2Vec.Mag() + PARTICLE_MASS*PARTICLE_MASS ); 
        energyLoss += (energy1 - energy2);
      }//<--- End loop over true traj points
      hMCELossUpstream->Fill(energyLoss);
    }//<--- End loop over primaries
    if (!isGoodEvent) continue;
    _nGoodMCEvents++;

    //if (_nTotalEvents == 1)
    //{
    //  for (int iIDE = 0; iIDE < IDETrackId->size(); iIDE++) 
    //  {
    //    auto xBin = hMCIdesXvsZ->GetXaxis()->FindBin(IDEPos->at(iIDE)[2]);
    //    auto yBin = hMCIdesXvsZ->GetYaxis()->FindBin(IDEPos->at(iIDE)[0]);
    //    hMCIdesXvsZ->SetBinContent(xBin, yBin, IDEEnergy->at(iIDE));
    //  }
    //}

    if (_nTotalEvents == 1)
    {
      for (int iPt = 0; iPt < NTrTrajPts->at(0); iPt++)
      {
        auto xBin = hMCXvsZ->GetXaxis()->FindBin(MidPosZ->at(0)[iPt]);
        auto yBin = hMCXvsZ->GetYaxis()->FindBin(MidPosX->at(0)[iPt]);
        hMCXvsZ->SetBinContent(xBin, yBin, 10);
      }
      for (size_t iDtr = 0; iDtr < NDTrTrajPts->size(); iDtr++)
      {
        if (!isCharged(DPdgCode->at(iDtr))) continue;
        for (int iPt = 0; iPt < NTrTrajPts->at(iDtr); iPt++)
        {
          auto xBin = hMCXvsZ->GetXaxis()->FindBin(DMidPosZ->at(iDtr)[iPt]);
          auto yBin = hMCXvsZ->GetYaxis()->FindBin(DMidPosX->at(iDtr)[iPt]);
          hMCXvsZ->SetBinContent(xBin, yBin, 10);
        }
      }
    }
  

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Begin study

    if (IN_DEBUG) std::cout << "\nFilling track lengths...";
    // Make distribution of secondary track lengths
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      // Only look at daughters of primary
      if (Mother->at(iG4) != _primTrkId) continue;
  
      // Make sure she's charged
      if ( !isCharged(PDG->at(iG4)) ) continue;
      TVector3 startPoint(StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4)); 
      TVector3 endPoint(EndPointx->at(iG4), EndPointy->at(iG4), EndPointz->at(iG4));
      hMCSecondaryTrkLength->Fill((startPoint-endPoint).Mag());
    }

    // Get first interaction in TPC
    if (IN_DEBUG) std::cout << "\nGetting first interaction in TPC...\n";
    _vertices.clear();
    getInteractionsInTpc();
    Vertex_t firstVertexInTpc(-1, "none", TVector3(0,0,-100));
    if (_vertices.size()) firstVertexInTpc = _vertices[0];
    hMCNumInteractions->Fill(_vertices.size());
    if (firstVertexInTpc.process.find("none") != std::string::npos) continue;

    // Identify visible secondaries
    _visSec.clear();
    identifyVisibleSecondaries(firstVertexInTpc);
    if (_visSec.size() != 1) continue;
    _nElasticLikeEvents++;

    if (IN_DEBUG) 
    {
      std::cout << "\nEvent                  = " << event 
                << "\nFirst int point in TPC = " << MidPosZ->at(0)[firstVertexInTpc.point]
                << "\nFirst process in TPC   = " << firstVertexInTpc.process
                << "\nProcesses:\n";
      for (const auto& p : *InteractionPointType) std::cout << p << endl;
    }
    
    // Prepare for vertex study
    if (IN_DEBUG) std::cout << "Initializing ides...\n";
    _ides.clear();
    _ides.reserve(TrackIdes_x->size());
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      for (size_t iIde = 0; iIde < TrackIdes_e->at(iG4).size(); iIde++)
      {
        //cout << "  " << iIde << endl;
        TrackIde_t tide( TrackId->at(iG4),
                         iG4,
                         TrackIdes_e->at(iG4)[iIde], 
                         TVector3( TrackIdes_x->at(iG4)[iIde], TrackIdes_y->at(iG4)[iIde], TrackIdes_z->at(iG4)[iIde] )
                       );
        _ides.push_back(tide);
      }
    }
    
    if (IN_DEBUG) std::cout << "Starting vertex study...\n";
    doVertexStudy(firstVertexInTpc);

    // Make the plot
    for (int rbin = 1; rbin <= hMCIdeVertexElastic->GetNbinsX(); rbin++)
    {
      float rcut = hMCIdeVertexElastic->GetXaxis()->GetBinCenter(rbin);
      float en = getIdeNearVertex(firstVertexInTpc.position, rcut);
      if (firstVertexInTpc.process == "hadElastic") 
      {
        int ebin = hMCIdeVertexElastic->GetYaxis()->FindBin(en);
        float content = hMCIdeVertexElastic->GetBinContent(rbin, ebin);
        hMCIdeVertexElastic->SetBinContent(rbin, ebin, content+en);
      }
      if (firstVertexInTpc.process == "pi-Inelastic")
      {
        int ebin = hMCIdeVertexInelastic->GetYaxis()->FindBin(en);
        float content = hMCIdeVertexInelastic->GetBinContent(rbin, ebin);
        hMCIdeVertexInelastic->SetBinContent(rbin, ebin, content+en);
      }
    }

// End study
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  }//<---End loop over entries


  // Event reduction table
  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events simulated:                           " << _nTotalEvents                 << endl
            << "Tracks entered TPC:                         " << _nPrimariesEntered            << endl
            << "Events with > 0 entered tracks:             " << _nGoodMCEvents                << endl
            << "Elastic like events:                        " << _nElasticLikeEvents           << endl
            << "Unknown processes:\n";
  for (const auto& p : _processes) cout << p << endl;

  // Make plots
  makePlots();
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
 */
void vertexStudy::getInteractionsInTpc()
{
  if (InteractionPoint->size())
  {
    for (int iPt = 0; iPt < InteractionPoint->size(); iPt++)
    {
      Vertex_t thisVertex( InteractionPoint->at(iPt),
                           InteractionPointType->at(iPt),
                           TVector3(MidPosX->at(0)[InteractionPoint->at(iPt)], MidPosY->at(0)[InteractionPoint->at(iPt)], MidPosZ->at(0)[InteractionPoint->at(iPt)]) );
      // Check if in FV
      if (!inActiveRegion(thisVertex.position)) continue;
      
      bool weGotIt(false);
      for (const auto& v : _vertices)
      {
        float diff = (v.position-thisVertex.position).Mag();
        if (  diff < 0.01 && v.process == thisVertex.process) weGotIt = true;
      }
      if (!weGotIt) _vertices.push_back(thisVertex);
    }
  }
  else
  {
    if (IN_DEBUG) std::cout << "Checking daughter!\n";
    // Check the daughters
    std::map<int, std::string> temp;
    for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
      if (Mother->at(iG4) != _primTrkId)continue;

      auto proc = Process->at(iG4);
      temp.emplace(iG4, proc);
    }

    // Look for decay and capture, this should only leave inelastic
    bool isDecay(false);
    bool isCapture(false);
    bool isInelastic(false);
    int tempG4 = -1;
    for (const auto& m : temp)
    {
      if (m.second.find("Decay") != std::string::npos) isDecay = true;
      if (m.second.find("Capture") != std::string::npos) isCapture = true;
      if (m.second.find("pi-Inelastic") != std::string::npos) {tempG4=m.first; isInelastic = true;}
    }
    if (isDecay && isCapture){std::cerr << "\n\nFAILING!\n"; exit(1);}
    if (!isDecay && !isCapture)
    {
      if (isInelastic) 
      {
        auto proc = Process->at(tempG4);
        TVector3 pos(StartPointx->at(tempG4), StartPointy->at(tempG4), StartPointz->at(tempG4));
        Vertex_t thisVertex( convertToPrimaryPoint(pos),
                             proc,
                             pos );
        // Check if in FV
        if (inActiveRegion(thisVertex.position)) _vertices.push_back(thisVertex);
      }
      else exit(1);
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
int vertexStudy::determineG4Id(const int& tid)
{
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    if (TrackId->at(iG4 == tid)) return iG4;
  }
  return -999;
} 

/**
 * @brief Method to identify visible secondaries attached to the interaction vertex
 * 
 * @param vertex The interaction vertex
 */
void vertexStudy::identifyVisibleSecondaries(const Vertex_t& vertex)
{
  // If the primary doesn't end here
  size_t primiG4 = determineG4Id(_primTrkId);
  if ( (vertex.point+1) < (NTrTrajPts->at(0)-1) ) _visSec.push_back(primiG4);

  // Fill inelastic scattering angles for all relevant daughters
  for (size_t iG4 = 0; iG4 < geant_list_size; iG4++)
  {
    // @todo Take absolute value here?
    if (Mother->at(iG4) != _primTrkId) continue;

    // Make sure she is charged 
    if ( !isCharged(PDG->at(iG4)) ) continue;

    // This is a charged daughter
    TVector3 dPos0(StartPointx->at(iG4), StartPointy->at(iG4), StartPointz->at(iG4));
    TVector3 dPosf(EndPointx  ->at(iG4), EndPointy  ->at(iG4), EndPointz  ->at(iG4));

    // Check track length
    // @todo Do we need this cut?
    if ((dPosf-dPos0).Mag() < SECONDARY_LENGTH_CUT) continue;

    // This particle needs to be attached to the vertex
    if ((vertex.position-dPos0).Mag() < 0.01) _visSec.push_back(iG4);
  }//<-- End loop over G4 particles
}


/**
 * @brief Area to study IDEs near vertex.
 * 
 */
void vertexStudy::doVertexStudy(const Vertex_t& vertex)
{
  bool isElastic(false);
  bool isInelastic(false);
  if (vertex.process.find("hadElastic") != std::string::npos || vertex.process.find("CoulombScat") != std::string::npos) isElastic = true;
  if (vertex.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;
  if (!isElastic && !isInelastic) return;
  if (isElastic && isInelastic) std::cout << "WARNING: Check the process types!\n";

  // Get the trj point here
  TVector3 posHere = vertex.position;

  // Option to save ides 
  if (event == EVENT_FOR_EVT_DSP) writeIdes(0);

  // First subtract out _ides which belong to primary or single secondary
  subtractIdes(vertex.point);

  // Option to save ides
  if (event == EVENT_FOR_EVT_DSP) writeIdes(1);

  // Total energy within some radius
  for (int iBin = 1; iBin <= hMCIdeVertexElastic->GetNbinsX(); iBin++)
  {
    float r = hMCIdeVertexElastic->GetXaxis()->GetBinCenter(iBin);
    float totEnergy = getIdeNearVertex(posHere, r);

    if (isElastic)   hMCIdeVertexElastic  ->Fill(r, totEnergy);
    if (isInelastic) hMCIdeVertexInelastic->Fill(r, totEnergy);
  }
}

/**
 * @brief Method to subtract out _ides which belong to primary or single secondary.
 * 
 */
void vertexStudy::subtractIdes(const int& p)
{
  if (_visSec.size() != 1) cout << "\nWHAT ARE YOU DOING HERE? " << _visSec.size() << "\n";

  // Loop over _ides 
  // @note Negative tid means EM activity, -1*tid is mother
  for (auto& tide : _ides)
  {
    // If this is the primary or the single secondary set the energy to zero
    if (std::abs(TrackId->at(tide.g4Id)) == _primTrkId || tide.g4Id == _visSec[0]) tide.energy = 0;
   
    // Check vicinity to tracks
    checkVicinity(tide, p);
  }//<-- End loop over _ides
}

/**
 * @brief Writing ides to text file for viewing later
 * 
 * @param stage 
 */
void writeIdes(int stage)
{
  std::string name = "ides_for_viewing";
  if (stage) name = name+"_sub.txt";
  else name = name+".txt";

  std::ofstream file(name.c_str());
  for (const auto& ide : _ides)
  {
    file << ide.pos.X() << " " << ide.pos.Y() << " " << ide.pos.Z() << " " << ide.energy << std::endl;
  }
  file.close();
}

/**
 * @brief Check vicinity to primary and secondary. 
 * 
 * @param tide The track ide
 * @param p The interaction point
 */
void vertexStudy::checkVicinity(TrackIde_t& tide, const int& p)
{
  // @note using argoneut numbers
  float CONE_DISTANCE_CUT(2.4);  // cm
  float CONE_ANGLE_CUT(120);     // degrees
  float CYLINDER_DISTANCE_CUT(3); // argoneut used 5cm

  // Get direction vectors
  TVector3 point = tide.pos;
  TVector3 vertex( MidPosX->at(0)[p],   MidPosY->at(0)[p],   MidPosZ->at(0)[p]);
  TVector3 primDir(MidPosX->at(0)[p-1], MidPosY->at(0)[p-1], MidPosZ->at(0)[p-1]);
  primDir = primDir - vertex;
  primDir = primDir.Unit();
  TVector3 secDir(EndPointx->at(_visSec[0]), EndPointy->at(_visSec[0]), EndPointz->at(_visSec[0]));
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
 * Subtracts out all _ides from primary particle and single secondary. The remaining
 * energy should be from neutrons and/or dexcitation gammas. 
 * 
 * @note We should've set other _ides equal to zero
 * @todo Should we consider a track length cut? 
 * 
 * @param vertex Interaction vertex
 * @param radius Bubble radius 
 * @return float Total energy deposited within bubble radius of interaction vertex
 */
float vertexStudy::getIdeNearVertex(const TVector3& vertex, const float& radius)
{
  float totalEnergy(0);
  for (const auto& tide : _ides) 
  {
    if ( (tide.pos-vertex).Mag() < radius ) totalEnergy += tide.energy;
  }
  return totalEnergy;
}


/**
 * @brief Make plots and save to output file
 * 
 */
void makePlots()
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
  hMCPrimaryPx->Write();
  hMCPrimaryPy->Write();
  hMCPrimaryPz->Write();
  hMCPrimaryP->Write();
  hTrueLength->Write();
  hMostUpstreamZPos->Write();

  myRootFile.cd();
  myRootFile.mkdir("mc/");
  myRootFile.cd("mc/");
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
  hMCNumInteractions->Write();
  hMCXvsZ->Write();

  myRootFile.Close();
}



/**
 * @brief Check if in active region
 *
 */
bool inActiveRegion( const TVector3& thePos  )
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
bool inTPCRegion( const TVector3& thePos  )
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
int vertexStudy::convertToPrimaryPoint(const TVector3& position)
{
  for (int iPt = 0; iPt < NTrTrajPts->at(0); iPt++)
  {
    TVector3 testPosition(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);

    if ( (testPosition-position).Mag() < 0.01 ) return iPt;
  }
  return 0;
}



