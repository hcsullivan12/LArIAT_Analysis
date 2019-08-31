/**
 * @file xs_ana_truth.C
 * @brief XS making script for truth information 
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#define xs_ana_truth_cxx
#include "xs_ana_truth.h"
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
bool InTpcRegion(const TVector3& thePos);
inline void PrintVec(const TVector3& pos) {std::cout<<"\t("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")\n";};
std::string IntProcessToString(const int& p);

/// Mass of pion in MeV
float PARTICLE_MASS(139.57); 
inline double toKineticEnergy(const TVector3& mom){return std::sqrt( mom.Mag()*mom.Mag() + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS; }
/// @}

/// @name Useful variables
/// @{
/// Container for pdg codes for charged particles (cuz I can't remember) 
enum PdgCodes_t
{
  kElectron = 11,
  kMuon     = 13,
  kPion     = 211,
  kKaon     = 321,
  kProton   = 2212
};
/// Method to check if charged
inline bool isCharged(int const& p) 
{ 
  int pdg = std::abs(p);
  return (pdg == PdgCodes_t::kElectron || pdg == PdgCodes_t::kMuon || 
          pdg == PdgCodes_t::kPion     || pdg == PdgCodes_t::kKaon || 
          pdg == PdgCodes_t::kProton);
}
/// Method to get particle name
inline std::string getParticle(int pdg)
{
  pdg = std::abs(pdg);
  std::string name = "unknown";
  switch(pdg)
  {
    case PdgCodes_t::kElectron:{name="electron";break;}
    case PdgCodes_t::kMuon:    {name="muon";break;}
    case PdgCodes_t::kPion:    {name="pion";break;}
    case PdgCodes_t::kKaon:    {name="kaon";break;}
    case PdgCodes_t::kProton:  {name="proton";break;}
  }
  return name;
}
/// Counters
int _nTotalEvents(0), _nEventsInel(0),
    _nEventsUniqueMatch(0), _nEventsMatch(0),
    _nEventsInelSec(0), _nEventsInelZeroSec(0),
    _nEventsInelOneSec(0), _nMuonBkg(0),
    _nProtonBkg(0), _nOtherBkg(0),
    _nElectronBkg(0), _nPionBkg(0),
    _nGoodEvents(0);

/// vertices in tpc
std::vector<Vertex_t> _vertices;
/// Background in matching
std::set<std::string> _bkgCandidates;

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

/// max on dedx
float DEDX_MAX(40);

/// The assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(36); //MeV
/// @}

/// @name Histograms
/// @{
TH1D* hTrueInteractingKe = new TH1D("hTrueInteractingKe",  "True Interacting",   24, 0, 1200);
TH1D* hTrueIncidentKe    = new TH1D("hTrueIncidentKe",     "True Incident",      24, 0, 1200);
TH1D* hTruePitch         = new TH1D("hTruePitch", "hTruePitch", 1000, 0, 100);
TH1D* hTrueEnDep         = new TH1D("hTrueEnDep", "hTrueEnDep", 1000, 0, 100);

TH1D* hTruFirstPosInTpcX = new TH1D("hTruFirstPosInTpcX", "First Position In Tpc X", 120, -5, 55);
TH1D* hTruFirstPosInTpcY = new TH1D("hTruFirstPosInTpcY", "First Position In Tpc Y", 100, -25, 25);
TH1D* hTruFirstPosInTpcZ = new TH1D("hTruFirstPosInTpcZ", "First Position In Tpc Z", 240, -10, 110);
TH1D* hTruIntPosInTpcX = new TH1D("hTruIntPosInTpcX", "Last Position In Tpc X", 120, -5, 55);
TH1D* hTruIntPosInTpcY = new TH1D("hTruIntPosInTpcY", "Last Position In Tpc Y", 100, -25, 25);
TH1D* hTruIntPosInTpcZ = new TH1D("hTruIntPosInTpcZ", "Last Position In Tpc Z", 240, -10, 110);

TH1D* hTruWc4KeCapture = new TH1D("hTruWc4KeCapture", "True Wc4 Ke Capture", 24, 0, 1200);
TH1D* hTruWc4KeDecay   = new TH1D("hTruWc4KeDecay",   "True Wc4 Ke Decay",   24, 0, 1200);

TH1D* hRecoFirstPosInTpcX = new TH1D("hRecoFirstPosInTpcX", "Reco First Pos In Tpc X", 120, -5, 55);
TH1D* hRecoFirstPosInTpcY = new TH1D("hRecoFirstPosInTpcY", "Reco First Pos In Tpc Y", 100, -25, 25);
TH1D* hRecoFirstPosInTpcZ = new TH1D("hRecoFirstPosInTpcZ", "Reco First Pos In Tpc Z", 240, -10, 110);
TH1D* hRecoIntPosInTpcX = new TH1D ("hRecoIntPosInTpcX",  "Reco Last Pos In Tpc X",  120, -5, 55);
TH1D* hRecoIntPosInTpcY = new TH1D ("hRecoIntPosInTpcY",  "Reco Last Pos In Tpc Y",  100, -25, 25);
TH1D* hRecoIntPosInTpcZ = new TH1D ("hRecoIntPosInTpcZ",  "Reco Last Pos In Tpc Z",  240, -10, 110);

TH1D* hTruLength  = new TH1D("hTruLength",  "Tru Length",  200, 0, 100);
TH1D* hRecoLength = new TH1D("hRecoLength", "Reco Length", 200, 0, 100);
TH2D* hTruVsRecoLength = new TH2D("hTruVsRecoLength", "Tru Vs Reco Length", 200, 0, 100, 200, 0, 100);

TH1D* hDiffLength           = new TH1D("hDiffLength", "hDiffLength", 400, -100, 100);
TH1D* hDiffFirstPosInTpcX   = new TH1D("hDiffFirstPosInTpcX",   "hDiffFirstPosInTpcX",   80, -20, 20);
TH1D* hDiffFirstPosInTpcY   = new TH1D("hDiffFirstPosInTpcY",   "hDiffFirstPosInTpcY",   80, -20, 20);
TH1D* hDiffFirstPosInTpcZ   = new TH1D("hDiffFirstPosInTpcZ",   "hDiffFirstPosInTpcZ",   40, -10, 10);
TH1D* hDiffFirstPosInTpcMag = new TH1D("hDiffFirstPosInTpcMag", "hDiffFirstPosInTpcMag", 200, 0, 100);
TH1D* hDiffIntPosInTpcX    = new TH1D("hDiffIntPosInTpcX",    "hDiffIntPosInTpcX",    80, -20, 20);
TH1D* hDiffIntPosInTpcY    = new TH1D("hDiffIntPosInTpcY",    "hDiffIntPosInTpcY",    80, -20, 20);
TH1D* hDiffIntPosInTpcZ    = new TH1D("hDiffIntPosInTpcZ",    "hDiffIntPosInTpcZ",    400, -100, 100);
TH1D* hDiffIntPosInTpcMag  = new TH1D("hDiffIntPosInTpcMag",  "hDiffIntPosInTpcMag",  200, 0, 100);
/// @}

/// Output root file
TFile myRootFile("XS_ANA_TRUTH.root", "RECREATE");

/**
 * @brief Main loop
 * 
 * @param inDebug For debugging
 */
void xs_ana_truth::Loop(int inDebug)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (int jentry=0; jentry<=nentries;jentry++) 
  {
    IN_DEBUG = inDebug;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // Increment our total event counter 
    _nTotalEvents++;
    if(_nTotalEvents%500==0)cout<<_nTotalEvents<<" / "<<nentries<<endl;
    if(IN_DEBUG)cout<<"\nProcessing run #"<<run<<" event#"<<event<<"\n";

    // get the interactions in the tpc
    _vertices.clear();
    getInteractionsInTpc();
    if(IN_DEBUG){for(const auto& v : _vertices)v.Dump();}

    // Measure true cross section
    //doTrueXs();

    // loop over tracks
    size_t nTrksMatched = 0;
    int    idMatchedTrk = -1;
    for (int iTrk=0; iTrk<ntracks_reco; iTrk++)
    {
      if (track_WC2TPC_match->at(iTrk))
      {
        idMatchedTrk=iTrk;
        nTrksMatched++;
      }
    }
    // leave if we didn't get any tracks
    if (nTrksMatched>0) _nEventsMatch++;
    if (nTrksMatched!=1)continue;
    _nEventsUniqueMatch++;

    // apply filters
    //doFilter(idMatchedTrk);

    // look at primary truth information
    // get first point in tpc
    TVector3 firstPosInTpc(0,0,-100);
    for (int iPt=0; iPt<MidPosX->at(0).size(); iPt++)
    {
      TVector3 pos(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);
      if(!InTpcRegion(pos))continue;
      firstPosInTpc = pos;
      break;
    }
    if(!InTpcRegion(firstPosInTpc))
    {
      // What did then?
      // make sure candidates are charged, and proj onto front face
      if(IN_DEBUG)cout<<"Warning: Run #"<<run<<" Event #"<<event<<" primary did not enter tpc but a track was matched\n";
      if(IN_DEBUG)cout<<"Primary ended at "<<EndPointx->at(0)<<" " <<EndPointy->at(0)<<" " <<EndPointz->at(0)<<"\n";
      if(IN_DEBUG)cout<<"Candidates...\n";
      std::set<std::string> cand;
      for (int ig4=0; ig4<geant_list_size; ig4++)
      {
        if(process_primary->at(ig4)==1)continue;
        if(!isCharged(PDG->at(ig4)))continue;
        TVector3 sp(StartPointx->at(ig4), StartPointy->at(ig4), StartPointz->at(ig4)); 
        TVector3 mom(StartPx->at(ig4), StartPy->at(ig4), StartPz->at(ig4)); 
        TVector3 proj = sp - (sp.Z()/mom.Z())*mom;
        if(TPC_X_BOUND[0]>proj.X() || proj.X()>TPC_X_BOUND[1])continue;
        if(TPC_Y_BOUND[0]>proj.Y() || proj.Y()>TPC_Y_BOUND[1])continue;
        cand.insert(getParticle(PDG->at(ig4)));
        _bkgCandidates.insert(getParticle(PDG->at(ig4)));
      }
      if(IN_DEBUG)for(const auto& c : cand){cout<<"\t"<<c<<endl;}

      // some particles take precedence
      if(cand.find("pion")!=cand.end()){_nPionBkg++;}
      else if(cand.find("muon")!=cand.end()){_nMuonBkg++;}
      else if(cand.find("proton")!=cand.end()){_nProtonBkg++;}
      else if(cand.find("electron")!=cand.end()){_nElectronBkg++;}
      else {_nOtherBkg++;}

      // if there are no candidates, wth did we do wrong?
      if(!cand.size()){cerr<<"Error. No candidates for matched track!\n";exit(1);}
      // we're done here, we can't track this particle
      continue;
    }
    // get last position in tpc
    int prim_pts = MidPosX->at(0).size();
    TVector3 lastPosInTpc(MidPosX->at(0)[prim_pts-1],MidPosY->at(0)[prim_pts-1],MidPosZ->at(0)[prim_pts-1]);
    for (int iPt=MidPosX->at(0).size()-1; iPt>=0; iPt--)
    {
      TVector3 pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
      if(!InTpcRegion(pos))continue;
      lastPosInTpc=pos;
      break;
    }

    // cases of decay and capture
    for(const auto& v : _vertices)
    {
      // plotting kinetic energy at WC 4 for decay and capture
      TVector3 mom(1000*MidPx->at(0)[0],1000*MidPy->at(0)[0],1000*MidPz->at(0)[0]);
      float kinEn = std::sqrt(mom*mom + PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;
      if(v.process=="CaptureAtRest")hTruWc4KeCapture->Fill(kinEn);
      if(v.process=="Decay")hTruWc4KeDecay->Fill(kinEn);
    }

    // looking at reco
    // check if track is inverted
    TVector3 trk_startpos(track_start_x->at(idMatchedTrk),track_start_y->at(idMatchedTrk),track_start_z->at(idMatchedTrk));
    TVector3 trk_endpos(track_end_x->at(idMatchedTrk),track_end_y->at(idMatchedTrk),track_end_z->at(idMatchedTrk));
    if(trk_endpos.Z()<trk_startpos.Z())
    {
      auto temp=trk_startpos;
      trk_startpos=trk_endpos;
      trk_endpos=temp;
    }

    // get closest interaction point in tpc, default is last tpc position
    Vertex_t vtx(getPrimaryPoint(lastPosInTpc),"none",lastPosInTpc);
    float minDist = std::numeric_limits<float>::max();
    for(const auto& v : _vertices)
    {
      float dist = (v.position-trk_endpos).Mag();
      //if(v.process=="hadElastic")continue;
      if(dist<minDist){minDist=dist;vtx=v;}
    }
    if(vtx.position.Z()<firstPosInTpc.Z())cout<<event<<"Warning: Interaction point is before first point. "<<vtx.position.Z()<<" "<<firstPosInTpc.Z()<<"\n";

    //cout<<col_track_z->at(idMatchedTrk)[0]<<" " << col_track_z->at(idMatchedTrk)[col_track_z->at(idMatchedTrk).size()-1]<<endl;

    // should be properly matched pions
    _nGoodEvents++;
    double trueLength = (firstPosInTpc-vtx.position).Mag();
    double recoLength = (trk_startpos-trk_endpos).Mag();

    hTruFirstPosInTpcX->Fill(firstPosInTpc.X());
    hTruFirstPosInTpcY->Fill(firstPosInTpc.Y());
    hTruFirstPosInTpcZ->Fill(firstPosInTpc.Z());
    hTruIntPosInTpcX->Fill(vtx.position.X());
    hTruIntPosInTpcY->Fill(vtx.position.Y());
    hTruIntPosInTpcZ->Fill(vtx.position.Z());
    hTruLength->Fill(trueLength);

    hRecoFirstPosInTpcX->Fill(trk_startpos.X());
    hRecoFirstPosInTpcY->Fill(trk_startpos.Y());
    hRecoFirstPosInTpcZ->Fill(trk_startpos.Z());
    hRecoIntPosInTpcX->Fill(trk_endpos.X());
    hRecoIntPosInTpcY->Fill(trk_endpos.Y());
    hRecoIntPosInTpcZ->Fill(trk_endpos.Z());
    hRecoLength->Fill(recoLength);
    
    hTruVsRecoLength->Fill( recoLength, trueLength );
    hDiffLength->Fill( recoLength - trueLength );
    hDiffFirstPosInTpcX->Fill(   (trk_startpos-firstPosInTpc).X() );
    hDiffFirstPosInTpcY->Fill(   (trk_startpos-firstPosInTpc).Y() );
    hDiffFirstPosInTpcZ->Fill(   (trk_startpos-firstPosInTpc).Z() );
    hDiffFirstPosInTpcMag->Fill( (trk_startpos-firstPosInTpc).Mag() );
    hDiffIntPosInTpcX->Fill(   (trk_endpos-vtx.position).X() );
    hDiffIntPosInTpcY->Fill(   (trk_endpos-vtx.position).Y() );
    hDiffIntPosInTpcZ->Fill(   (trk_endpos-vtx.position).Z() );
    hDiffIntPosInTpcMag->Fill( (trk_endpos-vtx.position).Mag() );

    //if((trk_endpos-vtx.position).Mag()>20)cout<<event<<" "<<trk_endpos.Z()<<" "<<vtx.position.Z() <<endl;

  }//<---End loop over entries


  // Event reduction table
  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events:                      " << _nTotalEvents       << endl
            << "Events with WC match:        " << _nEventsMatch       << endl
            << "Events with unique WC match: " << _nEventsUniqueMatch << endl
            << "Muon background:             " << _nMuonBkg           << endl
            << "Proton background:           " << _nProtonBkg         << endl
            << "Electron background:         " << _nElectronBkg       << endl
            << "Other background:            " << _nOtherBkg          << endl
            << "Good events:                 " << _nGoodEvents        << endl;
            
            
  std::cout<<"Candidates:"<<std::endl;
  for(const auto& c : _bkgCandidates)cout<<"\t"<<c<<endl;

  // Make plots
  MakePlots();
  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Getting interactions in tpc
 * 
 */
void xs_ana_truth::getInteractionsInTpc()
{
  // @note We should've only filled the interactions with 
  //       what g4 spat out. We handle the zero case here.
  // @todo Do we need to check the end point still?
  if (InteractionPoint->size())
  {
    for (int iInt=0; iInt<InteractionPoint->size(); iInt++)
    {
      int prim_pt = InteractionPoint->at(iInt);
      std::string proc = InteractionPointType->at(iInt);
      TVector3 pos(MidPosX->at(0)[prim_pt], MidPosY->at(0)[prim_pt], MidPosZ->at(0)[prim_pt]);
      if(!InActiveRegion(pos))continue;
      Vertex_t vtx(prim_pt, proc, pos);
      _vertices.push_back(vtx);
    }
  }
  else
  {
    if(IN_DEBUG)cout<<event<<"Checking end process and daughters...\n";
    // this has to be something catastrophic
    std::string proc_maybe = "none";
    int primTrkId = -999999;
    int primG4Id  = -999999;
    for (int ig4=0; ig4<geant_list_size; ig4++)
    {
      if(process_primary->at(ig4)!=1)continue;

      primTrkId=TrackId->at(ig4);
      primG4Id=ig4;

      // make sure this is in the tpc, if not, she is throughgoing
      TVector3 pos(EndPointx->at(ig4), EndPointy->at(ig4), EndPointz->at(ig4));
      if(!InActiveRegion(pos))break;
      proc_maybe = G4FinalProcess->at(ig4);
      break;
    }
    if(primTrkId==-999999){cerr<<"Error: Did not get primary trk id,\n";exit(1);}

    // @note For some reason, we must check the end process for Decay.
    //       But others show up too, for the others, check the daughter's processes.
    //       Don't know why this happens, but we have to hack something up to circumvent this.
    if(proc_maybe!="Decay" && proc_maybe!="none" && proc_maybe!="LArVoxelReadoutScoringProcess")
      {cerr<<"\nSomething happened at EVENT="<<event<<" PROCESS="<<proc_maybe<<"\n"<<endl;exit(1);}

    // hopefully geant4 got this part correct...
    if(proc_maybe=="none")return;

    // if we found decay, get out of here
    if (proc_maybe=="Decay")
    {
      int prim_pt = MidPosX->at(0).size()-1;
      TVector3 pos(EndPointx->at(primG4Id), EndPointy->at(primG4Id), EndPointz->at(primG4Id));
      Vertex_t vtx(prim_pt, proc_maybe, pos);
      _vertices.push_back(vtx);
      return;
    }

    // check again... 
    if(proc_maybe!="LArVoxelReadoutScoringProcess")
      {cerr<<"\nSomething happened at EVENT="<<event<<" PROCESS="<<proc_maybe<<"\n"<<endl;exit(1);}

    // @note Event displays suggest that LArVoxelReadoutScoringProcess is actually CaptureAtRest.
    //       There is a Bragg peak at the end of tracks. We will tag these as capture.
    int prim_pt = MidPosX->at(0).size()-1;
    TVector3 pos(EndPointx->at(primG4Id), EndPointy->at(primG4Id), EndPointz->at(primG4Id));
    Vertex_t vtx(prim_pt, "CaptureAtRest", pos);
    _vertices.push_back(vtx);
  }
}

/**
 * @brief Convert position to primary point
 * 
 */
int xs_ana_truth::getPrimaryPoint(const TVector3& pos)
{
  for (int iPt=0; iPt<NTrTrajPts->at(0); iPt++)
  {
    TVector3 testPosition(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);
    if ( (testPosition-pos).Mag() < 0.01 )return iPt;
  }
  cout<<"Warning: Return 0\n";
  return 0;
}

  /**
 * @brief Area to apply thin slab to true energy deposits
 *
 */
void xs_ana_truth::doTrueXs()
{
  // get first position in tpc
  TVector3 firstTpcPos(0,0,-100);
  TVector3 initialMom(0,0,0);
  for (int iPt=0; iPt < MidPosX->at(0).size(); iPt++)
  {
    TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
    if(!InActiveRegion(true_pos))continue;
    firstTpcPos=true_pos;
    initialMom = TVector3(1000*MidPx->at(0)[iPt], 1000*MidPy->at(0)[iPt], 1000*MidPz->at(0)[iPt]);
    break;
  }
  if(!InActiveRegion(firstTpcPos))return;
  // get last position closest to end position
  TVector3 lastTpcPos(0,0,200);
  for (int iPt=MidPosX->at(0).size()-1; iPt>=0; iPt--)
  {
    TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
    if(!InActiveRegion(true_pos))continue;
    lastTpcPos=true_pos;
    break;
  }
  if(!InActiveRegion(lastTpcPos))return;
  if(firstTpcPos.Z()>lastTpcPos.Z())cout<<"Warning: Check the tpc positions.\n";

  // get unifomly spaced points to the nearest inelastic vertex
  std::vector<TVector3> orderedPoints;
  bool isInel = false;
  if(_vertices.size())
  {
    TVector3 current_pos = firstTpcPos;
    for(const auto& v : _vertices)
    {
      orderedPoints.push_back(current_pos);
      // adding in range (current, v.pos)
      addUniformPoints(current_pos, v.position, orderedPoints);
      current_pos = v.position;
      if(v.process.find("pi-Inelastic")!=std::string::npos){isInel=true;break;}
    }
  }
  else 
  {
    orderedPoints.push_back(firstTpcPos);
    addUniformPoints(firstTpcPos, lastTpcPos, orderedPoints);
    orderedPoints.push_back(lastTpcPos);
  }

  // ides need to be ordered
  std::map<double, std::vector<double>> orderedIdes;
  for(int iIde=0; iIde<TrackIdes_x->at(0).size(); iIde++)
  {
    std::vector<double> temp = {TrackIdes_x->at(0)[iIde], TrackIdes_y->at(0)[iIde], TrackIdes_z->at(0)[iIde], TrackIdes_e->at(0)[iIde]};
    orderedIdes.emplace(TrackIdes_z->at(0)[iIde], temp);
  }

  // now we can fill
  double kinEn = std::sqrt(initialMom.Mag()*initialMom.Mag()+PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;

  int oldPt = 0;
  for (int iPt=1; iPt<orderedPoints.size(); iPt++, oldPt++)
  {
    auto oldPos     = orderedPoints[oldPt];
    auto currentPos = orderedPoints[iPt];

    auto itOldIde = orderedIdes.begin();
    double eDep = 0;
    for(auto itIde = orderedIdes.begin(); itIde!=orderedIdes.end(); itIde++,itOldIde++)
    {
      auto thisIde = itIde->second;
      if(thisIde[2]<oldPos.Z())continue;
      if(thisIde[2]>currentPos.Z())continue;
      eDep += thisIde[3];
    }
    kinEn -= eDep;

    hTruePitch->Fill((oldPos-currentPos).Mag());
    hTrueEnDep->Fill(eDep);
    hTrueIncidentKe->Fill(kinEn);
  }
  if(isInel)hTrueInteractingKe->Fill(kinEn);
}

/**
 * @brief Add ordered points between two positions
 *
 */
void xs_ana_truth::addUniformPoints(const TVector3& x0, const TVector3& xf, std::vector<TVector3>& orderedPoints)
{
  // we're not including the end points
  float length = (x0-xf).Mag();
  float trackPitch = 0.47;
  int nPts = (int)(length/trackPitch);
  for(int iPt=1;iPt<=nPts;iPt++)
  {
    auto newPoint = x0 + iPt*(trackPitch/length) * (xf - x0);
    orderedPoints.push_back(newPoint);
  }
  // remove second to last point if it's too close to the last point
  auto secondToLastPoint = orderedPoints.end()-2;
  int distance = std::distance(orderedPoints.begin(), secondToLastPoint);
  if( (*secondToLastPoint-xf).Mag() < trackPitch*0.5 )orderedPoints.erase(secondToLastPoint);
}



/**
 * @brief Make plots and save to output file
 * 
 */
void MakePlots()
{
  myRootFile.cd();
  hTrueInteractingKe->Write();
  hTrueIncidentKe->Write();
  hTruePitch->Write();
  hTrueEnDep->Write();
  hTruFirstPosInTpcX->Write();
  hTruFirstPosInTpcY->Write();
  hTruFirstPosInTpcZ->Write();
  hTruIntPosInTpcX->Write();
  hTruIntPosInTpcY->Write();
  hTruIntPosInTpcZ->Write();
  hTruWc4KeCapture->Write();
  hTruWc4KeDecay->Write();
  hRecoFirstPosInTpcX->Write();
  hRecoFirstPosInTpcY->Write();
  hRecoFirstPosInTpcZ->Write();
  hRecoIntPosInTpcX->Write();
  hRecoIntPosInTpcY->Write();
  hRecoIntPosInTpcZ->Write();
  hTruLength->Write();
  hRecoLength->Write();
  hTruVsRecoLength->Write();
  hDiffLength->Write();
  hDiffFirstPosInTpcX->Write();
  hDiffFirstPosInTpcY->Write();
  hDiffFirstPosInTpcZ->Write();
  hDiffFirstPosInTpcMag->Write();
  hDiffIntPosInTpcX->Write();
  hDiffIntPosInTpcY->Write();
  hDiffIntPosInTpcZ->Write();
  hDiffIntPosInTpcMag->Write();

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
bool InTpcRegion( const TVector3& thePos  )
{
  if ( TPC_X_BOUND[0] < thePos.X() && thePos.X() < TPC_X_BOUND[1] &&
       TPC_Y_BOUND[0] < thePos.Y() && thePos.Y() < TPC_Y_BOUND[1] &&
       TPC_Z_BOUND[0] < thePos.Z() && thePos.Z() < TPC_Z_BOUND[1] ) return true;

  return false;
}




