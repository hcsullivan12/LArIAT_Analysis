/**
 * @file xs_ana.C
 * @brief XS making script
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#define xs_ana_cxx
#include "xs_ana.h"
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
/// Counters
int _nTotalEvents(0), _nEventsInel(0),
    _nEventsUniqueMatch(0), _nEventsMatch(0),
    _nEventsInelSec(0), _nEventsInelZeroSec(0),
    _nEventsInelOneSec(0);

/// Vector of g4 Ids for visible secondaries
std::vector<size_t> visSec;
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
/// @}

/// Output root file
TFile myRootFile("XS_ANA.root", "RECREATE");

/**
 * @brief Main loop
 * 
 * @param inDebug For debugging
 */
void xs_ana::Loop(int inDebug)
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
    //if(_nTotalEvents==10)break;

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

    // check if inelastic
    bool isInel = isInelastic(idMatchedTrk);
    if(isInel)_nEventsInel++;

    // ke at front face
    // @todo what do I do about the WC momentum?
    float wcMom = wctrk_momentum[0];
    float kinEn = std::sqrt( wcMom*wcMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
    kinEn -= ENTRY_TPC_ENERGY_LOSS;

    // fill incident
    // @note The equals sign here <=
    int slabs = col_track_x->at(idMatchedTrk).size();
    for (int iSb = 0; iSb<=slabs; iSb++)
    {
      float dedx  = col_track_dedx->at(idMatchedTrk)[iSb];
      float pitch = col_track_pitch_hit->at(idMatchedTrk)[iSb];
      // protection against large dedx
      if(dedx>DEDX_MAX)continue;

      hIncidentKe->Fill(kinEn);
      hDeDx->Fill(dedx);
      hPitch->Fill(pitch);
      kinEn -= dedx * pitch;

      // check for < 0
      if(kinEn<=0){cout<<"Warning: Filling zero\n";break;}
    }
    // fill interacting
    if(isInel)hInteractingKe->Fill(kinEn);
  }//<---End loop over entries


  // Event reduction table
  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events:                      " << _nTotalEvents       << endl
            << "Events with WC match:        " << _nEventsMatch       << endl
            << "Events with unique WC match: " << _nEventsUniqueMatch << endl
            << "Inelastic interactions:      " << _nEventsInel        << endl
            << "Inelastic by > 2 sec:        " << _nEventsInelSec     << endl
            << "Inelastic by 1 sec:          " << _nEventsInelOneSec  << endl
            << "Inelastic by 0 sec:          " << _nEventsInelZeroSec << endl
            << endl;

  // Make plots
  MakePlots();
  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Area to determine inelasticity
 *
 */
bool xs_ana::isInelastic(int const& idMatchedTrk)
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
 * @brief Get angle between to tracks
 * 
 * @todo This currently assumes the tracks have no kinks,
 *       namely, it uses start and end points.
 * 
 */
float xs_ana::getAngle(int const& idm, int const& ids)
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




