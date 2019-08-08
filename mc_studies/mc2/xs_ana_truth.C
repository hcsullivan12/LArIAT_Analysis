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
bool InTPCRegion(const TVector3& thePos);
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

/// Track length cut for secondaries
float SECONDARY_LENGTH_CUT(2.0);

/// Vertex cut
float VERTEX_CUT(3.0);

/// angle cut
float ANGLE_CUT(20);

/// downstream cut
float DOWNSTREAM_Z_CUT(88);

/// max on dedx
float DEDX_MAX(40);

/// The assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(36); //MeV
/// @}

/// @name Histograms
/// @{
TH1D* hInteractingKe = new TH1D("hInteractingKe",  "Interacting",   24, 0, 1200);
TH1D* hIncidentKe    = new TH1D("hIncidentKe",     "Incident",      24, 0, 1200);         
TH1D* hTrkZ          = new TH1D("hTrkZ",           "Z pos in tpc",  50, 0, 100);
TH1D* hDeDx          = new TH1D("hDeDx",           "dEdX",          200, 0, 50);
TH1D* hPitch         = new TH1D("hPitch",          "Track pitch",   100, 0, 5);
TH1D* hTruFirstPosInTpcX = new TH1D("hTruFirstPosInTpcX", "First Position In Tpc X", 60, -5, 55);
TH1D* hTruFirstPosInTpcY = new TH1D("hTruFirstPosInTpcY", "First Position In Tpc Y", 50, -25, 25);
TH1D* hTruFirstPosInTpcZ = new TH1D("hTruFirstPosInTpcZ", "First Position In Tpc Z", 120, -10, 110);
/// @}

/// Output root file
TFile myRootFile("XS_ANA.root", "RECREATE");

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

    // this is where reconstruction and truth split
    // look at primary truth information
    // get first point in tpc
    TVector3 firstPosInTpc(0,0,-100);
    for (int iPt=0; iPt<=MidPosX->at(0).size(); iPt++)
    {
      TVector3 pos(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);
      if(!InActiveRegion(pos))continue;
      firstPosInTpc = pos;
      break;
    }
    if(!InActiveRegion(firstPosInTpc))
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

    // should be properly matched pions
    _nGoodEvents++;
    hTruFirstPosInTpcX->Fill(firstPosInTpc.X());
    hTruFirstPosInTpcY->Fill(firstPosInTpc.Y());
    hTruFirstPosInTpcZ->Fill(firstPosInTpc.Z());

    // get the interactions in the tpc
    _vertices.clear();
    getInteractionsInTpc();
    if(IN_DEBUG){for(const auto& v : _vertices)v.Dump();}

    // now what?

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
 * @brief Area to determine inelasticity
 *
 */
bool xs_ana_truth::isInelastic(int const& idMatchedTrk)
{
  // three cases to consider
  // 1) If more than one track attached to vertex, yes
  // 2) If one, check angle
  // 3) If none, where is the end point? 

  // get the start and enpoint
  int mtrk_n_points = col_track_x->at(idMatchedTrk).size();
  TVector3 mtrk_start_point( col_track_x->at(idMatchedTrk)[0],
                             col_track_y->at(idMatchedTrk)[0],
                             col_track_z->at(idMatchedTrk)[0] );
  TVector3 mtrk_end_point( col_track_x->at(idMatchedTrk)[mtrk_n_points-1],
                           col_track_y->at(idMatchedTrk)[mtrk_n_points-1],
                           col_track_z->at(idMatchedTrk)[mtrk_n_points-1] );
    
  // check if inverted
  if(mtrk_start_point.Z()>mtrk_end_point.Z())
  {
    if(IN_DEBUG)cout<<"Track inverted: "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<endl;
    auto temp = mtrk_start_point;
    mtrk_start_point = mtrk_end_point;
    mtrk_end_point = temp;
  }

  // loop over tracks to see if any are attached to vertex
  std::vector<int> sec_ids;
  for (int iTrk=0; iTrk<ntracks_reco; iTrk++)
  {
    // skip the matched trk
    if(iTrk==idMatchedTrk)continue;

    // make sure this track has points
    int n_points = col_track_x->at(iTrk).size();
    if(n_points<2)continue;
      
    TVector3 start_point( col_track_x->at(iTrk)[0],
                          col_track_y->at(iTrk)[0],
                          col_track_z->at(iTrk)[0] );
    TVector3 end_point( col_track_x->at(iTrk)[n_points-1],
                        col_track_y->at(iTrk)[n_points-1],
                        col_track_z->at(iTrk)[n_points-1] );
    float dist1 = (mtrk_end_point - start_point).Mag();
    float dist2 = (mtrk_end_point - end_point).Mag();
    float dist_min = dist1 < dist2 ? dist1 : dist2;

    // checl inversion
    if(dist2<dist1)
    {
      if(IN_DEBUG)cout<<"Track inverted: "<<start_point.Z()<<" "<<end_point.Z()<<endl;
      auto temp = start_point;
      start_point = end_point;
      end_point = temp;
    }
        
    // we're done with this track if it's not close enough
    if(dist_min>VERTEX_CUT)continue;

    // we got one!
    // @todo can we not use vertex information from reconstruction?
    sec_ids.push_back(iTrk);
  }
  // we should have all secondaries connected to vertex
  // Case 1)
  if(sec_ids.size() > 1){_nEventsInelSec++; return true;}
  else if(sec_ids.size()==1)
  {
    if(idMatchedTrk==sec_ids[0])cerr<<"Error: Sec id = matched\n";
    float angle = getAngle(idMatchedTrk, sec_ids[0]);
    cout<<angle<<endl;
    if(angle>ANGLE_CUT){_nEventsInelOneSec++; return true;}
    // @todo anything else here?
    else return false;
  }    

  //
  // @todo Think about what to do if sec = 0 but the track endpoint 
  //       is in the middle of the tpc. It's not obvious to me right 
  //       what's better to do: assume no interaction, assume inelastic
  //       'or some other third thing' - SS
  //       For now, assume inelastic.
  //
  if(sec_ids.size())cerr<<"Error. Sec ids > 0\n";
  
  if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT){_nEventsInelZeroSec++;return true;}
  return false;
}

/**
 * @brief Getting interactions in tpc
 * 
 */
void xs_ana_truth::getInteractionsInTpc()
{
  // @note We should've only filled the interactions with 
  //       what g4 spat out. We handle the zero case here.
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
    if(IN_DEBUG)cout<<"Checking end process and daughters...\n";
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
    if (proc_maybe!="Decay" && proc_maybe!="none" && proc_maybe!="LArVoxelReadoutScoringProcess")
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
 * @brief Get angle between to tracks
 * 
 * @todo This currently assumes the tracks have no kinks,
 *       namely, it uses start and end points.
 * 
 */
float xs_ana_truth::getAngle(int const& idm, int const& ids)
{
  TVector3 trk1_sp( track_start_x->at(idm),
                    track_start_y->at(idm),
                    track_start_z->at(idm) );
  TVector3 trk1_ep( track_end_x->at(idm),
                    track_end_y->at(idm),
                    track_end_z->at(idm) );
  // check inversion
  if(trk1_sp.Z()>trk1_ep.Z())
  {
    if(IN_DEBUG)cout<<"Track inverted: "<<trk1_sp.Z()<<" "<<trk1_ep.Z()<<endl;
    auto temp = trk1_sp;
    trk1_sp = trk1_ep;
    trk1_ep = temp;
  }
  auto dir1 = (trk1_ep-trk1_sp).Unit();

  TVector3 trk2_sp( track_start_x->at(ids),
                    track_start_y->at(ids),
                    track_start_z->at(ids) );
  TVector3 trk2_ep( track_end_x->at(ids),
                    track_end_y->at(ids),
                    track_end_z->at(ids) );
  // check inversion
  float dist1 = (trk1_ep - trk2_sp).Mag();
  float dist2 = (trk1_ep - trk2_ep).Mag();
  if(dist2<dist1)
  {
    if(IN_DEBUG)cout<<"Track inverted: "<<trk2_sp.Z()<<" "<<trk2_ep.Z()<<endl;
    auto temp = trk2_sp;
    trk2_sp = trk2_ep;
    trk2_ep = temp;
  }

  auto dir2 = (trk2_ep-trk2_sp).Unit();
  return std::acos(dir1.Dot(dir2))*180./M_PI;
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
 * @brief Make plots and save to output file
 * 
 */
void MakePlots()
{
  myRootFile.cd();
  hInteractingKe->Write();
  hIncidentKe->Write();
  hTrkZ->Write();
  hDeDx->Write();
  hPitch->Write();

  hTruFirstPosInTpcX->Write();
  hTruFirstPosInTpcY->Write();
  hTruFirstPosInTpcZ->Write();

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




