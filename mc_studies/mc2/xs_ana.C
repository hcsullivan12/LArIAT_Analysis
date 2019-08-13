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
void makePlots();
bool inActiveRegion(const TVector3& thePos);
bool inTpcRegion(const TVector3& thePos);
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
    _nEventsInelOneSec(0), _nCorrect_alg1(0), 
    _nIncorrect_alg1(0), _nMissed_alg1(0),
    _nCorrect_alg2(0), _nIncorrect_alg2(0), 
    _nMissed_alg2(0), _n_correct(0),
    _n_incorrect(0), _n_missed(0), 
    _nMuonBkg(0),
    _nProtonBkg(0), _nOtherBkg(0),
    _nElectronBkg(0), _nPionBkg(0);

/// Vector of g4 Ids for visible secondaries
std::vector<size_t> _visSec;

/// interaction vertices
std::vector<Vertex_t> _vertices;

/// the reason for tagging
std::string _reason_alg1 = "";
std::string _reason_alg2 = "";

/// true particle in the tpc
std::string _trueTpcParticle = "pion";

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
float SECONDARY_LENGTH_CUT(1.0);
float SECONDARY_LENGTH_CUT2(3.0);

/// Vertex cut
float VERTEX_CUT(3.0);

/// angle cut
float ANGLE_CUT(10);

/// downstream cut
float DOWNSTREAM_Z_CUT(80);

/// max on dedx
float DEDX_MAX(20);

/// max on pitch
float PITCH_MAX(2.0);

/// The assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(42); //MeV
/// @}

/// @name Histograms
/// @{
TH1D* hInteractingKe = new TH1D("hInteractingKe",  "Interacting",   24, 0, 1200);
TH1D* hIncidentKe    = new TH1D("hIncidentKe",     "Incident",      24, 0, 1200);         
TH1D* hTrkZ          = new TH1D("hTrkZ",           "Z pos in tpc",  50, 0, 100);
TH1D* hDeDx          = new TH1D("hDeDx",           "dEdX",          200, 0, 50);
TH1D* hPitch         = new TH1D("hPitch",          "Track pitch",   100, 0, 5);
TH1D* hVertexDistSec = new TH1D("hVertexDistSec", "Vertex Dist Sec", 100, 0, 10);

TH1D* hDiffInitKe = new TH1D("hDiffInitKe", "hDiffInitKe", 1000, -500, 500);
TH1D* hEnLossUpstream = new TH1D("hEnLossUpstream", "hEnLossUpstream", 1000, -500, 500);

TH1D* hDiffFirstPosInTpcX   = new TH1D("hDiffFirstPosInTpcX",   "hDiffFirstPosInTpcX",   80, -20, 20);
TH1D* hDiffFirstPosInTpcY   = new TH1D("hDiffFirstPosInTpcY",   "hDiffFirstPosInTpcY",   80, -20, 20);
TH1D* hDiffFirstPosInTpcZ   = new TH1D("hDiffFirstPosInTpcZ",   "hDiffFirstPosInTpcZ",   40, -10, 10);
TH1D* hDiffFirstPosInTpcMag = new TH1D("hDiffFirstPosInTpcMag", "hDiffFirstPosInTpcMag", 200, 0, 100);
TH1D* hDiffIntPosInTpcX   = new TH1D("hDiffIntPosInTpcX",   "hDiffIntPosInTpcX",   80, -20, 20);
TH1D* hDiffIntPosInTpcY   = new TH1D("hDiffIntPosInTpcY",   "hDiffIntPosInTpcY",   80, -20, 20);
TH1D* hDiffIntPosInTpcZ   = new TH1D("hDiffIntPosInTpcZ",   "hDiffIntPosInTpcZ",   400, -100, 100);
TH1D* hDiffIntPosInTpcMag = new TH1D("hDiffIntPosInTpcMag", "hDiffIntPosInTpcMag", 200, 0, 100);

TH2D* hTruVsRecoLength = new TH2D("hTruVsRecoLength", "Tru Vs Reco Length", 200, 0, 100, 200, 0, 100);

TH1D* hDiffLength           = new TH1D("hDiffLength", "hDiffLength", 400, -100, 100);
/// @}

/// Output root file
TFile myRootFile("XS_ANA.root", "RECREATE");

/**
 * @brief Main loop
 * 
 * @param inDebug For debugging
 */
void xs_ana::Loop(int inDebug, int isMc)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (int jentry=223; jentry<=223/*nentries*/;jentry++) 
  {
    IN_DEBUG = inDebug;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Increment our total event counter 
    _nTotalEvents++;
    if(_nTotalEvents%500==0)cout<<_nTotalEvents<<" / "<<nentries<<endl;
    if(_nTotalEvents%10==0)continue;

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

    // checking the particle in the tpc
    if(isMc)
    {
      TVector3 firstPosInTpc(0,0,-100);
      for (int iPt=0; iPt<MidPosX->at(0).size(); iPt++)
      {
        TVector3 pos(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);
        if(!inTpcRegion(pos))continue;
        firstPosInTpc = pos;
        break;
      }
      if(!inTpcRegion(firstPosInTpc))
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
        if(cand.find("pion")!=cand.end()){_nPionBkg++;_trueTpcParticle="sec_pion";}
        else if(cand.find("muon")!=cand.end()){_nMuonBkg++;_trueTpcParticle="muon";}
        else if(cand.find("proton")!=cand.end()){_nProtonBkg++;_trueTpcParticle="proton";}
        else if(cand.find("electron")!=cand.end()){_nElectronBkg++;_trueTpcParticle="electron";}
        else {_nOtherBkg++;}

        // if there are no candidates, wth did we do wrong?
        if(!cand.size()){cerr<<"Error. No candidates for matched track!\n";exit(1);}
        if(_trueTpcParticle!="pion")continue;
      }
    }

    // Fill these temporary containers 
    std::vector<double> t_trk_x, t_trk_y, t_trk_z, t_trk_p, t_trk_dedx;
    int slabs = col_track_x->at(idMatchedTrk).size();
    for (int iSb=0; iSb<slabs; iSb++)
    {
      t_trk_x.push_back(col_track_x->at(idMatchedTrk)[iSb]);
      t_trk_y.push_back(col_track_y->at(idMatchedTrk)[iSb]);
      t_trk_z.push_back(col_track_z->at(idMatchedTrk)[iSb]);
      t_trk_p.push_back(col_track_pitch_hit->at(idMatchedTrk)[iSb]);
      t_trk_dedx.push_back(col_track_dedx->at(idMatchedTrk)[iSb]);
    }

    // check for inversion
    if(t_trk_z[0] > t_trk_z.back())
    {
      std::reverse(t_trk_x.begin(),    t_trk_x.end());
      std::reverse(t_trk_y.begin(),    t_trk_y.end());
      std::reverse(t_trk_z.begin(),    t_trk_z.end());
      std::reverse(t_trk_p.begin(),    t_trk_p.end());
      std::reverse(t_trk_dedx.begin(), t_trk_dedx.end());
    }

    // check if inelastic
    cout<<"\nEvent #"<<event<<endl;
    cout<<"Kinks "<<track_kink_x->at(idMatchedTrk).size()<<endl;
    TVector3 end_position(t_trk_x.back(), t_trk_y.back(), t_trk_z.back());
    TVector3 recobVertex(t_trk_x.back(), t_trk_y.back(), t_trk_z.back());
    bool isInel_1 = isInelastic_Alg1(idMatchedTrk);
    bool isInel_2 = isInelastic_Alg2(idMatchedTrk, recobVertex);
    if(isInel_2 && recobVertex.Z()>t_trk_z.front()) end_position = recobVertex;

    TVector3 kink_position;
    bool passKink = checkKink(idMatchedTrk, kink_position);
    if(!isInel_1 && passKink && kink_position.Z()>t_trk_z.back()) {isInel_1=true;end_position=kink_position;}
    if(!isInel_2 && passKink && kink_position.Z()>t_trk_z.back()) {isInel_2=true;end_position=kink_position;}

    // now we have to decide which one to take and the end point
    bool isInel = isInel_1 || isInel_2;

    // Fill these containers 
    std::vector<double> trk_x, trk_y, trk_z, trk_p, trk_dedx;
    for (int iSb=0; iSb<slabs; iSb++)
    {
      if(t_trk_z[iSb] > end_position.Z())break;
      trk_x.push_back(t_trk_x[iSb]);
      trk_y.push_back(t_trk_y[iSb]);
      trk_z.push_back(t_trk_z[iSb]);
      trk_p.push_back(t_trk_p[iSb]);
      trk_dedx.push_back(t_trk_dedx[iSb]);
    }

    // ke at front face
    // @todo what do I do about the WC momentum?
    float wcMom = wctrk_momentum[0];
    float kinEn = std::sqrt( wcMom*wcMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
    kinEn -= ENTRY_TPC_ENERGY_LOSS;
    if(isMc)
    {
      // get vertices
      _vertices.clear();
      getInteractionsInTpc();
      bool foundInel(false);
      bool foundOther(false);
      // @note for now, skip decay and capture at rest
      for(const auto& v : _vertices)
      {
        v.Dump();
        if(v.process=="pi-Inelastic")foundInel=true;
        if(v.process=="Decay"||v.process=="CaptureAtRest")foundOther=true;
      }
      if(foundOther)continue;
      
      cout<<"First alg...\n";
      cout<<"\tThe reason "<<_reason_alg1<<endl;
      if(isInel_1 && foundInel){cout<<"\tCorrect1!\n";_nCorrect_alg1++;}
      if(isInel_1 && !foundInel){cout<<"\tDetermined but not inelastic!\n";_nIncorrect_alg1++;}
      if(!isInel_1 && foundInel){cout<<event<<"\tMissed!\n";_nMissed_alg1++;}
      if(!isInel_1 && !foundInel){cout<<"\tCorrect2!\n";_nCorrect_alg1++;}
      cout<<"Second alg...\n";
      cout<<"\tThe reason "<<_reason_alg2<<endl;
      if(isInel_2 && foundInel){cout<<"\tCorrect1!\n";_nCorrect_alg2++;}
      if(isInel_2 && !foundInel){cout<<"\tDetermined but not inelastic!\n";_nIncorrect_alg2++;}
      if(!isInel_2 && foundInel){cout<<event<<"\tMissed!\n";_nMissed_alg2++;}
      if(!isInel_2 && !foundInel){cout<<"\tCorrect2!\n";_nCorrect_alg2++;}
      
      cout<<"Ntracks = "<<ntracks_reco<<endl;
      cout<<"Start point ("<<trk_x[0]<<","<<trk_y[0]<<","<<trk_z[0]<<")\n";
      cout<<"End point   ("<<trk_x[trk_x.size()-1]<<","<<trk_y[trk_x.size()-1]<<","<<trk_z[trk_x.size()-1]<<")\n";
      if(isInel && foundInel)_n_correct++;
      if(isInel_1 && !foundInel)_n_incorrect++;
      if(!isInel_1 && foundInel)_n_missed++;
      if(!isInel_1 && !foundInel)_n_correct++;

      // find nearest pi-inelastic vertex
      TVector3 trk_startpos(trk_x.front(), trk_y.front(), trk_z.front());
      TVector3 trk_endpos(trk_x.back(), trk_y.back(), trk_z.back());
      float minDist = std::numeric_limits<float>::max();
      Vertex_t closest_vtx(-1, "none", trk_endpos);
      for(const auto& v : _vertices)
      {
        if(v.process.find("pi-Inelastic")==std::string::npos)continue;
        float dist = (trk_endpos-v.position).Mag();
        if(dist<minDist){minDist=dist; closest_vtx=v;}
      }
      
      // get first position closest to reco start position
      TVector3 firstTruePos(0,0,-100);
      TVector3 lastTruePos(0,0,200);
      TVector3 initialMom(0,0,0);
      minDist = std::numeric_limits<float>::max();
      for (int iPt=0; iPt < MidPosX->at(0).size(); iPt++)
      {
        TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
        TVector3 true_mom(1000*MidPx->at(0)[iPt],1000*MidPy->at(0)[iPt],1000*MidPz->at(0)[iPt]);
        float dist = (true_pos-trk_startpos).Mag();
        if(dist<minDist){minDist=dist;firstTruePos=true_pos;initialMom=true_mom;}
      }
      // get last position closest to reco end position
      minDist = std::numeric_limits<float>::max();
      for (int iPt=MidPosX->at(0).size()-1; iPt>=0; iPt--)
      {
        TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
        float dist = (true_pos-trk_endpos).Mag();
        if(dist<minDist){minDist=dist;lastTruePos=true_pos;}
      }
      if(firstTruePos.Z()<0)cout<<"Warning: Check true first position\n";
      if(lastTruePos.Z()>90)cout<<"Warning: Check true last position\n";
      if(firstTruePos.Z()>lastTruePos.Z())cout<<"Warning: Check true positions\n";

      hDiffFirstPosInTpcX->Fill(   (trk_startpos-firstTruePos).X() );
      hDiffFirstPosInTpcY->Fill(   (trk_startpos-firstTruePos).Y() );
      hDiffFirstPosInTpcZ->Fill(   (trk_startpos-firstTruePos).Z() );
      hDiffFirstPosInTpcMag->Fill( (trk_startpos-firstTruePos).Mag() );

      auto last_interesting_point = closest_vtx.position;
      if(closest_vtx.process=="none") last_interesting_point = lastTruePos;

      cout<<firstTruePos.Z()<<" "<<last_interesting_point.Z()<<endl;
      if((firstTruePos-last_interesting_point).Mag()<1)cout<<"\nHEYYY "<<event<<endl;

      hDiffIntPosInTpcX->Fill(   (trk_endpos-last_interesting_point).X() );
      hDiffIntPosInTpcY->Fill(   (trk_endpos-last_interesting_point).Y() );
      hDiffIntPosInTpcZ->Fill(   (trk_endpos-last_interesting_point).Z() );
      hDiffIntPosInTpcMag->Fill( (trk_endpos-last_interesting_point).Mag() );
      hDiffLength->Fill(      (trk_startpos-trk_endpos).Mag() - (firstTruePos-last_interesting_point).Mag());
      hTruVsRecoLength->Fill( (trk_startpos-trk_endpos).Mag(),  (firstTruePos-last_interesting_point).Mag());
      //if((trk_endpos-lastTruePos).Z()>60)cout<<"EVENT "<<event<<endl;

      // compare the kinetic energies
      TVector3 startingMom(1000*MidPx->at(0)[0], 1000*MidPy->at(0)[0], 1000*MidPz->at(0)[0]);
      float startKe   = std::sqrt(startingMom*startingMom + PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;
      float initialKe = std::sqrt(initialMom*initialMom + PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;
      hEnLossUpstream->Fill(startKe-initialKe);
      hDiffInitKe->Fill(initialKe-kinEn);
    }

    // interpolate
    //float last_dedx = trk_dedx.front();
    //float enLoss(0);
    //for (int iSb = 1; iSb<(trk_x.size()-1); iSb++)
    //{
    //  float dedx = trk_dedx[iSb];
    //  float pitch = trk_p[iSb];
    //  enLoss += dedx*pitch;
    //  float next_dedx = trk_dedx[iSb+1];
    //  //if(dedx<DEDX_MAX)continue;
    //}

    // fill incident
    for (int iSb = 0; iSb<trk_x.size(); iSb++)
    {
      float dedx  = trk_dedx[iSb];
      float pitch = trk_p[iSb];
      // protection against large dedx and pitch
      //if(pitch>PITCH_MAX){cout<<"HEYYY "<<pitch<<" "<<iSb<<"/"<<trk_x.size()<<endl;continue;}
      if(dedx>DEDX_MAX)continue;

      hIncidentKe->Fill(kinEn);
      hDeDx->Fill(dedx);
      hPitch->Fill(pitch);
      if(iSb!=(trk_x.size()-1))kinEn -= dedx * pitch;

      // check for < 0
      if(kinEn<=0){cout<<"Warning: "<<event<<" Filling zero\n";break;}
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
            << "Correct_1:                   " << _nCorrect_alg1      << endl
            << "Incorrect_1:                 " << _nIncorrect_alg1    << endl
            << "Missed_1:                    " << _nMissed_alg1       << endl
            << "Correct_2:                   " << _nCorrect_alg2      << endl
            << "Incorrect_2:                 " << _nIncorrect_alg2    << endl
            << "Missed_2:                    " << _nMissed_alg2       << endl
            << "Correct:                     " << _n_correct          << endl
            << "Incorrect:                   " << _n_incorrect        << endl
            << "Missed:                      " << _n_missed           << endl
            << endl;

  // Make plots
  makePlots();
  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Area to determine inelasticity using crude method
 *
 */
bool xs_ana::isInelastic_Alg1(int const& idMatchedTrk)
{
  // A few cases to consider
  // 1) If more than one track attached to vertex, yes
  // 2) If one, check angle
  // 3) If none, where is the end point? 

  // get the start and endpoint
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

  // loop over tracks to look for candidate secondaries
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
    float length = (start_point-end_point).Mag();
    if(length<SECONDARY_LENGTH_CUT)continue;
    float dist1 = (mtrk_end_point - start_point).Mag();
    float dist2 = (mtrk_end_point - end_point).Mag();
    float dist_min = dist1 < dist2 ? dist1 : dist2;

    // check inversion
    if(dist2<dist1)
    {
      if(IN_DEBUG)cout<<"Track inverted: "<<start_point.Z()<<" "<<end_point.Z()<<endl;
      auto temp = start_point;
      start_point = end_point;
      end_point = temp;
    }
    hVertexDistSec->Fill(dist_min);

    // begin checking for connection to vertex
    if(dist_min>VERTEX_CUT)continue;
    // we got one!
    sec_ids.push_back(iTrk);
  }

  // we should have all secondaries connected to vertex
  // Case 1)
  if(sec_ids.size() > 1){_nEventsInelSec++;_reason_alg1=">2 secondaries";return true;}
  else if(sec_ids.size()==1)
  {
    if(idMatchedTrk==sec_ids[0])cerr<<"Error: Sec id = matched\n";
    float angle = getAngle(idMatchedTrk, sec_ids[0]);
    _reason_alg1="angle="+std::to_string(angle);
    if(angle>ANGLE_CUT){_nEventsInelOneSec++;return true;}
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
  _reason_alg1="track endpoint";
  if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT && inActiveRegion(mtrk_end_point)){_nEventsInelZeroSec++;return true;}
  return false;
}

/**
 * @brief Area to determine inelasticity based on reco info
 * 
 */
bool xs_ana::isInelastic_Alg2(int const& idMatchedTrk, TVector3& recobVertex)
{
  // get the start and endpoint
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

  // check the vertex information
  struct recob_vertex_t 
  {
    std::vector<int> tids;
    double x,y,z;
    void dump()const{cout<<"\tVertex at ("<<x<<","<<y<<","<<z<<") with "<<tids.size()<<" tracks\n";};
  };
  std::vector<recob_vertex_t> vtxs;
  for (int iVtx=0; iVtx< vertex_x->size(); iVtx++)
  {
    recob_vertex_t vtx;
    vtx.x = vertex_x->at(iVtx);
    vtx.y = vertex_y->at(iVtx);
    vtx.z = vertex_z->at(iVtx);
    //cout<<vtx.z<<" "<<vertex_track_ids->at(iVtx).size()<<endl;
    // we only care about those which contain the matched track
    bool pass(false);
    for (int iTrk = 0; iTrk<vertex_track_ids->at(iVtx).size(); iTrk++)
    {
      if(vertex_track_ids->at(iVtx)[iTrk]==idMatchedTrk)pass=true;
      // check track length
      int tid = vertex_track_ids->at(iVtx)[iTrk];
      if(track_length->at(tid) < SECONDARY_LENGTH_CUT2)continue;
      vtx.tids.push_back(vertex_track_ids->at(iVtx)[iTrk]);
    }
    if(pass&&vtx.tids.size())vtxs.push_back(vtx);
    
  }
  if(!vtxs.size())
  {
    _reason_alg2 = "vtx = 0";
    if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT && inActiveRegion(mtrk_end_point)) return true;
    return false;
  }
  //for(const auto& v : vtxs)v.dump();

  // sort in increasing z
  std::sort(vtxs.begin(), vtxs.end(), [](const auto& l, const auto& r){return l.z<r.z;});
  if(vtxs.back().tids.size()>1){_reason_alg2="> 1";recobVertex = TVector3(vtxs.back().x, vtxs.back().y, vtxs.back().z);return true;}

  _reason_alg2 = "vtx = 1";
  if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT  && inActiveRegion(mtrk_end_point))return true;
  return false;
}

/**
 * @brief Check for kinks in primary track 
 *  
 */
bool xs_ana::checkKink(int const& idMatchedTrk, TVector3& kink_position)
{
  // get the start and endpoint
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

  // get the kink information
  int nKinks = track_kink_x->at(idMatchedTrk).size();
  if(!nKinks)return false;
  int kinkCount = 1;
  while (kinkCount <= nKinks)
  {
    // get points before and after kink
    TVector3 beforeKinkPos, afterKinkPos;
    TVector3 kinkPos(track_kink_x->at(idMatchedTrk)[kinkCount-1], track_kink_y->at(idMatchedTrk)[kinkCount-1], track_kink_z->at(idMatchedTrk)[kinkCount-1]);
    getPosNearKink(idMatchedTrk, beforeKinkPos, kinkPos, afterKinkPos);

    // compare 
    TVector3 incDir = (kinkPos-beforeKinkPos).Unit();
    TVector3 outDir = (afterKinkPos-kinkPos).Unit();
    float angle = std::acos(incDir.Dot(outDir))*180./M_PI;

    if(angle > ANGLE_CUT)
    {
      _reason_alg1=_reason_alg1+" angle="+std::to_string(angle);
      _reason_alg2=_reason_alg2+" angle="+std::to_string(angle);
      kink_position = kinkPos;
      return true;
    }
    kinkCount++;
  }
  return false;
}

/**
 * @brief Area to get points near kink position
 * 
 */
void xs_ana::getPosNearKink(int const& idMatchedTrk, TVector3& beforeKinkPos, const TVector3& kinkPos, TVector3& afterKinkPos)
{
  std::vector<double> t_trk_x, t_trk_y, t_trk_z;
  for (int iSb=0; iSb<col_track_x->at(idMatchedTrk).size(); iSb++)
  {
    t_trk_x.push_back(col_track_x->at(idMatchedTrk)[iSb]);
    t_trk_y.push_back(col_track_y->at(idMatchedTrk)[iSb]);
    t_trk_z.push_back(col_track_z->at(idMatchedTrk)[iSb]);
  }
  // check for inversion
  if(t_trk_z[0] > t_trk_z.back())
  {
    std::reverse(t_trk_x.begin(),t_trk_x.end());
    std::reverse(t_trk_y.begin(),t_trk_y.end());
    std::reverse(t_trk_z.begin(),t_trk_z.end());
  }

  // get pos before kink
  for(int iPt=0; iPt<t_trk_x.size(); iPt++)
  {
    if(t_trk_z[iPt] >= kinkPos.Z())break;
    beforeKinkPos = TVector3(t_trk_x[iPt], t_trk_y[iPt], t_trk_z[iPt]);
  }

  // get pos after kink
  for(int iPt=0; iPt<t_trk_x.size(); iPt++)
  {
    if(t_trk_z[iPt] <= kinkPos.Z())continue;
    afterKinkPos = TVector3(t_trk_x[iPt], t_trk_y[iPt], t_trk_z[iPt]);
    break;
  }
  if(beforeKinkPos.Z()>afterKinkPos.Z())
  {cout<<"Warning: "<<event<<" Check kink positions!\n";}
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
 * @brief Getting interactions in tpc
 * 
 */
void xs_ana::getInteractionsInTpc()
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
      if(!inActiveRegion(pos))continue;
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
      if(!inActiveRegion(pos))break;
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
 * @brief Make plots and save to output file
 * 
 */
void makePlots()
{
  myRootFile.cd();
  hInteractingKe->Write();
  hIncidentKe->Write();
  hTrkZ->Write();
  hDeDx->Write();
  hPitch->Write();
  hVertexDistSec->Write();
  hDiffInitKe->Write();
  hEnLossUpstream->Write();
  hDiffFirstPosInTpcX->Write();
  hDiffFirstPosInTpcY->Write();
  hDiffFirstPosInTpcZ->Write();
  hDiffFirstPosInTpcMag->Write();
  hDiffIntPosInTpcX->Write();
  hDiffIntPosInTpcY->Write();
  hDiffIntPosInTpcZ->Write();
  hDiffIntPosInTpcMag->Write();
  hTruVsRecoLength->Write();
  hDiffLength->Write();

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
bool inTpcRegion( const TVector3& thePos  )
{
  if ( TPC_X_BOUND[0] < thePos.X() && thePos.X() < TPC_X_BOUND[1] &&
       TPC_Y_BOUND[0] < thePos.Y() && thePos.Y() < TPC_Y_BOUND[1] &&
       TPC_Z_BOUND[0] < thePos.Z() && thePos.Z() < TPC_Z_BOUND[1] ) return true;

  return false;
}




