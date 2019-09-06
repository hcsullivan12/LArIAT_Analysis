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
void reset();
void dumpHeader(const int& runId, const int& evtId, const int& jent);
void fillXsPlots(float                      recoKinEn, 
                 const int&                 evtId,
                 const std::vector<double>& trk_x, 
                 const std::vector<double>& trk_y, 
                 const std::vector<double>& trk_z,
                 const std::vector<double>& trk_p, 
                 const std::vector<double>& trk_dedx, 
                 const bool& isInel);

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
  kProton   = 2212,
  kNeutron  = 2112,
  kPi0      = 111
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
int _nTotalEvents(0),       _nEventsInel(0),
    _nEventsUniqueMatch(0), _nEventsMatch(0),
    _nEventsInelSec(0),     _nEventsInelZeroSec(0),
    _nEventsInelOneSec(0),  _nCorrect_alg1(0), 
    _nIncorrect_alg1(0),    _nMissed_alg1(0),
    _nCorrect_alg2(0),      _nIncorrect_alg2(0), 
    _nMissed_alg2(0),       _n_correct(0),
    _n_incorrect(0),        _n_missed(0), 
    _nMuonBkg(0),           _nProtonBkg(0), 
    _nOtherBkg(0),          _nElectronBkg(0), 
    _nPionBkg(0),           _nEventsPrimPion(0);

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

/// Output root file
TFile myRootFile("XS_ANA.root", "RECREATE");

/// entry in tree
int _jentry(0);

/// flag to fill elastic background
bool _fillElasticBkg = true;
/// @}

/// @name Cuts and constants
/// @{
/// verbosity schema:
///     >=1 Warnings
///     >=2 Tracking 
///     >=3 Identification
int VERBOSE=(0);

/// TPC boundaries
float TPC_X_BOUND[2] = {   0.0, 47.0 };
float TPC_Y_BOUND[2] = { -20.0, 20.0 };
float TPC_Z_BOUND[2] = {   0.0, 90.0 };

/// Fiducial volume definition
float FV_X_BOUND[2] = {   1.0, 46.0 };
float FV_Y_BOUND[2] = { -18.0, 18.0 };
float FV_Z_BOUND[2] = {   0.0, 88.0 };

/// Track length cut for secondaries
float SECONDARY_LENGTH_CUT(1.0); // for secondaries at the end point
float SECONDARY_LENGTH_CUT2(3.0); // for secondaries attached to a recob::Vertex

/// Vertex cut
float VERTEX_CUT(3.0);

/// angle cut
float ANGLE_CUT(10);

/// downstream cut
float DOWNSTREAM_Z_CUT(80);

/// max on dedx
float DEDX_MAX(40);

/// max on pitch
float PITCH_MAX(2.0);

/// The assumed energy loss between the cryostat and the TPC 
float ENTRY_TPC_ENERGY_LOSS(42); //MeV
/// @}

/// @name Histograms
/// @{
// plots for xs
TH1D* hInteractingKe         = new TH1D("hInteractingKe",  "Interacting",   24, 0, 1200);
TH1D* hIncidentKe            = new TH1D("hIncidentKe",     "Incident",      24, 0, 1200);         
TH1D* hWellRecoInteractingKe = new TH1D("hWellRecoInteractingKe",  "WellRecoInteracting",   24, 0, 1200);
TH1D* hWellRecoIncidentKe    = new TH1D("hWellRecoIncidentKe",     "WellRecoIncident",      24, 0, 1200);         
TH1D* hTrueIncidentKe            = new TH1D("hTrueIncidentKe", "hTrueIncidentKe", 24, 0, 1200);
TH1D* hTrueInteractingKe         = new TH1D("hTrueInteractingKe", "hTrueInteractingKe", 24, 0, 1200);
TH1D* hWellRecoTrueIncidentKe    = new TH1D("hWellRecoTrueIncidentKe", "hWellRecoTrueIncidentKe", 24, 0, 1200);
TH1D* hWellRecoTrueInteractingKe = new TH1D("hWellRecoTrueInteractingKe", "hWellRecoTrueInteractingKe", 24, 0, 1200);
TH1D* hIncidentKeElasticBkg = new TH1D("hIncidentKeElasticBkg", "hIncidentKeElasticBkg", 24, 0, 1200);
TH1D* hInteractingKeElasticBkg = new TH1D("hInteractingKeElasticBkg", "hInteractingKeElasticBkg", 24, 0, 1200);
// reco
TH1D* hTrkZ          = new TH1D("hTrkZ",           "Z pos in tpc",  50, 0, 100);
TH1D* hDeDx          = new TH1D("hDeDx",           "dEdX",          200, 0, 50);
TH1D* hPitch         = new TH1D("hPitch",          "Track pitch",   100, 0, 5);
TH1D* hRecoLength      = new TH1D("hRecoLength", "hRecoLength", 110, 0, 110);
TH2D* hRecoAngleVsRecoKe  = new TH2D("hRecoAngleVsRecoKe", "hRecoAngleVsRecoKe", 24, 0, 1200, 180, 0, 180);
// well reco stuff
TH1D* hWellRecoLength     = new TH1D("hWellRecoLength", "hWellRecoLength", 110, 0, 110);
TH1D* hWellRecoTrueLength = new TH1D("hWellRecoTrueLength", "hWellRecoTrueLength", 110, 0, 110);
TH1D* hWellRecoDiffLength = new TH1D("hWellRecoDiffLength", "hWellRecoDiffLength", 400, -100, 100);
TH2D* hWellRecoTruVsRecoLength = new TH2D("hWellRecoTruVsRecoLength", "Well Reco Tru Vs Reco Length", 200, 0, 100, 200, 0, 100);
// truth
TH1D* hEnLossUpstream = new TH1D("hEnLossUpstream", "hEnLossUpstream", 1000, -500, 500);
TH1D* hVertexDistSec = new TH1D("hVertexDistSec", "Vertex Dist Sec", 100, 0, 10);
TH1D* hTrueLength      = new TH1D("hTrueLength", "hTrueLength", 110, 0, 110);
TH1D* hKeWc4Decay   = new TH1D("hKeWc4Decay", "hKeWc4Decay", 24, 0, 1200);
TH1D* hKeWc4Capture = new TH1D("hKeWc4Capture", "hKeWc4Capture", 24, 0, 1200);
TH1D* hTruePitch      = new TH1D("hTruePitch",      "hTruePitch",      100, 0, 5);
TH1D* hTrueEnDep      = new TH1D("hTrueEnDep",      "hTrueEnDep",      1000, 0, 1000);
TH1D* hWellRecoTruePitch       = new TH1D("hWellRecoTruePitch", "hWellRecoTruePitch", 100, 0, 5);
TH1D* hWellRecoTrueEnDep       = new TH1D("hWellRecoTrueEnDep", "hWellRecoTrueEnDep", 1000, 0, 1000);
TH2D* hTrueAngleVsTrueKe       = new TH2D("hTrueAngleVsTrueKe", "hTrueAngleVsTrueKe", 24, 0, 1200, 180, 0, 180);
// diffs or reco/truth comparison
TH1D* hDiffInitKe = new TH1D("hDiffInitKe", "hDiffInitKe", 1000, -500, 500);
TH1D* hDiffFirstPosInTpcX   = new TH1D("hDiffFirstPosInTpcX",   "hDiffFirstPosInTpcX",   80, -20, 20);
TH1D* hDiffFirstPosInTpcY   = new TH1D("hDiffFirstPosInTpcY",   "hDiffFirstPosInTpcY",   80, -20, 20);
TH1D* hDiffFirstPosInTpcZ   = new TH1D("hDiffFirstPosInTpcZ",   "hDiffFirstPosInTpcZ",   40, -10, 10);
TH1D* hDiffFirstPosInTpcMag = new TH1D("hDiffFirstPosInTpcMag", "hDiffFirstPosInTpcMag", 200, 0, 100);
TH1D* hDiffIntPosInTpcX   = new TH1D("hDiffIntPosInTpcX",   "hDiffIntPosInTpcX",   80, -20, 20);
TH1D* hDiffIntPosInTpcY   = new TH1D("hDiffIntPosInTpcY",   "hDiffIntPosInTpcY",   80, -20, 20);
TH1D* hDiffIntPosInTpcZ   = new TH1D("hDiffIntPosInTpcZ",   "hDiffIntPosInTpcZ",   400, -100, 100);
TH1D* hDiffIntPosInTpcMag = new TH1D("hDiffIntPosInTpcMag", "hDiffIntPosInTpcMag", 200, 0, 100);
TH1D* hDiffLength      = new TH1D("hDiffLength", "hDiffLength", 400, -100, 100);
TH2D* hTruVsRecoLength = new TH2D("hTruVsRecoLength", "Tru Vs Reco Length", 200, 0, 100, 200, 0, 100);
TH2D* hTrueVsRecoKinEn = new TH2D("hTrueVsRecoKinEn", "hTrueVsRecoKinEn", 24, 0, 1200, 24, 0, 1200);
/// @}

/**
 * @brief Main loop
 * @todo Plot z position of mis identified interactions.
 * @todo Event 27211 shows that we need to check end point interactions in all cases. 
 * @todo Angle algs showing 90.00000, something must be wrong...
 * 
 * @param v Verbosity
 */
void xs_ana::Loop(int v, int isMc)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (int jentry=20000; jentry<=nentries;jentry++) 
  {
    VERBOSE = v;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Reset global variables 
    reset();
    if(VERBOSE)dumpHeader(run,event,jentry);
    _jentry = jentry;

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
    if (isMc) 
    {
      bool isPrimPion = checkParticleInTpc();
      if(!isPrimPion)
      {
        if(VERBOSE>=1){cout<<"Skipping due to non entering primary pion...\n\n";}
        continue;
      }
    }
    _nEventsPrimPion++;

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

    // Begin checking if inelastic
    // We have seperate algorithms here...
    //   Alg1) Checks end point using crude proximity method
    //   Alg2) Check the reconstruction variables (recob::Vertex)
    //   Kink) Check for kinks found in pmtrack
    // check if inelastic
    bool isInel_1 = isInelastic_Alg1(idMatchedTrk);
    TVector3 end_position(t_trk_x.back(), t_trk_y.back(), t_trk_z.back());
    TVector3 recobVertex = end_position;
    bool isInel_2 = isInelastic_Alg2(idMatchedTrk, recobVertex);
    TVector3 kink_position;
    bool passKink = checkKink(idMatchedTrk, kink_position);
    if(VERBOSE>=2 && passKink){cout<<"Kink found at";PrintVec(kink_position);}

    // Now we decide...
    // We must select the candidate vertex most upstream,
    // otherise this could bias the XS.
    // The only exception would be if developed a method to further discriminate 
    // against elastic scatters, (e.g. vertex activity)
    // trust the vertex reconstruction
    if(isInel_2 && recobVertex.Z()>t_trk_z.front()) end_position = recobVertex;

    // check to see if kink is before current prediction if it passed the inelasticity test 
    if(isInel_1 && passKink && kink_position.Z()>t_trk_z.front() && kink_position.Z()<end_position.Z())end_position = kink_position;
    if(isInel_2 && passKink && kink_position.Z()>t_trk_z.front() && kink_position.Z()<end_position.Z())end_position = kink_position;

    // this is case if we didn't tag as inelastic
    if(!isInel_1 && passKink && kink_position.Z()>t_trk_z.front()) {isInel_1=true;end_position=kink_position;}
    if(!isInel_2 && passKink && kink_position.Z()>t_trk_z.front()) {isInel_2=true;end_position=kink_position;}
    

    // Final decision
    bool isInel = isInel_1 || isInel_2;

    // Fill these containers 
    std::vector<double> trk_x, trk_y, trk_z, trk_p, trk_dedx;
    // Where do we stop?
    float min_dist =std::numeric_limits<float>::max();
    int stop_sb = -1;
    for (int iSb=0; iSb<slabs; iSb++)
    {
      TVector3 pos(t_trk_x[iSb], t_trk_y[iSb], t_trk_z[iSb]);
      float dist = (pos-end_position).Mag();
      if(dist<min_dist){min_dist=dist;stop_sb=iSb;}
    }
    for (int iSb=0; iSb<=stop_sb; iSb++)
    {
      //if(t_trk_z[iSb] > end_position.Z())break;
      trk_x.push_back(t_trk_x[iSb]);
      trk_y.push_back(t_trk_y[iSb]);
      trk_z.push_back(t_trk_z[iSb]);
      trk_p.push_back(t_trk_p[iSb]);
      trk_dedx.push_back(t_trk_dedx[iSb]);
    }

    // look for bragg peak, only if we determined it to be inelastic 
    //if(isInel && isBraggPeak(trk_x, trk_y, trk_z, trk_p, trk_dedx))isInel = false;

    // ke at front face
    float wcMom = wctrk_momentum[0];
    float recoKinEn = std::sqrt( wcMom*wcMom + PARTICLE_MASS*PARTICLE_MASS ) - PARTICLE_MASS;
    recoKinEn -= ENTRY_TPC_ENERGY_LOSS;
    if(isMc) 
    {
      // get vertices
      _vertices.clear();
      getInteractionsInTpc();
      // look for decay and capture
      bool isDecayOrCap(false);
      for(const auto& v : _vertices)
      {
        //v.Dump();
        if(v.process=="Decay"){hKeWc4Decay->Fill(wcMom);isDecayOrCap=true;}
        else if(v.process=="CaptureAtRest"){hKeWc4Capture->Fill(wcMom);isDecayOrCap=true;}
      }
      if(VERBOSE>=3)dumpEventInfo();

      if(isDecayOrCap)
      {
        if(VERBOSE>=1){cout<<"Skipping due to decay or capture...\n\n";}
        continue;
      }
      studyMc(wcMom, recoKinEn, isInel, isInel_1, isInel_2, trk_x, trk_y, trk_z, trk_p, trk_dedx);
    }

    // fill incident interacting histos
    fillXsPlots(recoKinEn, event, trk_x, trk_y, trk_z, trk_p, trk_dedx, isInel);
  }//<---End loop over entries


  // Event reduction table
  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
            << "Events:                      " << _nTotalEvents       << endl
            << "Events with WC match:        " << _nEventsMatch       << endl
            << "Events with unique WC match: " << _nEventsUniqueMatch << endl
            << "Events primary pion:         " << _nEventsPrimPion    << endl
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
            << "Muon Bkg:                    " << _nMuonBkg           << endl
            << "Proton Bkg:                  " << _nProtonBkg         << endl
            << "Electron Bkg:                " << _nElectronBkg       << endl
            << "Sec pion Bkg:                " << _nPionBkg           << endl
            << "Other Bkg:                   " << _nOtherBkg          << endl
            << endl;


  // Make plots
  makePlots();
  gApplication->Terminate(0);
}//<---- End main loop

/**
 * @brief Reset global variables (playing with fire)
 * 
 */
void reset()
{
  _visSec.clear();
  _vertices.clear();
  _reason_alg1 = "";
  _reason_alg2 = "";
  _trueTpcParticle = "pion";
  _bkgCandidates.clear();
  _fillElasticBkg = false;
}

/**
 * @brief Make it pretty 
 *
 */
void dumpHeader(const int& runId, const int& evtId, const int& jent)
{
  cout<<"\n#########################################################\n";
  cout<<"RUN #"<<runId<<", EVENT #"<<evtId<<", JENTRY #"<<jent<<endl;
}

/**
 * @brief Area to determine inelasticity using crude method
 * @todo Anything else for case 1?
 * @todo Think about what to do if sec = 0 but the track endpoint 
 *       is in the middle of the tpc. It's not obvious to me right 
 *       what's better to do: assume no interaction, assume inelastic
 *       'or some other third thing' - SS
 *       For now, assume inelastic.
 */
bool xs_ana::isInelastic_Alg1(int const& idMatchedTrk)
{
  // A few cases to consider
  // 1) If more than one track attached to vertex, yes
  // 2) If one, check angle
  // 3) If none, where is the end point?
  if(VERBOSE>=3)cout<<"Running alg1...\n"; 

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
    if(VERBOSE>=2)cout<<"\tTrack inverted: "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<" ";
    auto temp = mtrk_start_point;
    mtrk_start_point = mtrk_end_point;
    mtrk_end_point = temp;
    if(VERBOSE>=2)cout<<" ==> "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<endl;
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
    if(VERBOSE>=3)cout<<"\tSec length = "<<length<<endl;
    if(length<SECONDARY_LENGTH_CUT)continue;
    float dist1 = (mtrk_end_point - start_point).Mag();
    float dist2 = (mtrk_end_point - end_point).Mag();
    float dist_min = dist1 < dist2 ? dist1 : dist2;
    if(VERBOSE>=3)cout<<"\tdist1 = "<<dist1<<"  dist2 = "<<dist2<<endl;

    // check inversion
    if(dist2<dist1)
    {
      if(VERBOSE>=2)cout<<"\tSec track inverted: "<<start_point.Z()<<" "<<end_point.Z()<<" ";
      auto temp = start_point;
      start_point = end_point;
      end_point = temp;
      if(VERBOSE>=2)cout<<" ==> "<<start_point.Z()<<" "<<end_point.Z()<<endl;
    }
    hVertexDistSec->Fill(dist_min);

    // begin checking for connection to vertex
    if(VERBOSE>=3)cout<<"\tdist_min = "<<dist_min<<endl;
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
    float angle = getAngle(idMatchedTrk, sec_ids[0], mtrk_end_point);
    _reason_alg1="angle="+std::to_string(angle);
    if(angle>ANGLE_CUT){_nEventsInelOneSec++;return true;}
    else return false;
  }    

  if(sec_ids.size())cerr<<"Error. Sec ids > 0\n";
  _reason_alg1="track endpoint";
  if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT && inActiveRegion(mtrk_end_point)){_nEventsInelZeroSec++;return true;}
  return false;
}

/**
 * @brief Area to determine inelasticity based on reco info
 * @todo Check matched track end point for case of tids == 2.
 *       Perhaps need to handle the angle better. 
 * @todo The case of vtx = 1 needs to be readdressed. 
 * 
 */
bool xs_ana::isInelastic_Alg2(int const& idMatchedTrk, TVector3& recobVertex)
{
  if(VERBOSE>=3)cout<<"Running alg2...\n"; 
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
    if(VERBOSE>=2)cout<<"\tTrack inverted: "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<" ";
    auto temp = mtrk_start_point;
    mtrk_start_point = mtrk_end_point;
    mtrk_end_point = temp;
    if(VERBOSE>=2)cout<<" ==> "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<endl;
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
  if(VERBOSE>=3){cout<<"\tRecob::Vertices...\n";for(const auto& v : vtxs)v.dump();}

  // sort in increasing z
  std::sort(vtxs.begin(), vtxs.end(), [](const auto& l, const auto& r){return l.z<r.z;});

  // we take the first vertex with > 1 tids attached
  for(const auto& v : vtxs)
  {
    if(v.tids.size()>1)
    {
      if(v.tids.size()==2)
      {
        //cout<<track_length->at(v.tids[0])<<" "<<track_length->at(v.tids[1])<<endl;
        // this is a scatter, check the angle
        int id1 = v.tids[0] == idMatchedTrk ? v.tids[0] : v.tids[1];
        int id2 = v.tids[1] == idMatchedTrk ? v.tids[0] : v.tids[1];
        if(id1==id2 || id1!=idMatchedTrk){cerr<<"ERROR. Check the ids.\n"; exit(1);}
        float angle = getAngle(id1, id2, TVector3(v.x,v.y,v.z));
        if(angle>ANGLE_CUT)
        {
          _reason_alg2="angle = "+std::to_string(angle);
          recobVertex = TVector3(v.x, v.y, v.z);
          return true;
        }
      }
      _reason_alg2="> 2";
      recobVertex = TVector3(v.x, v.y, v.z);
      return true;
    }
  }

  // This was either a mistake or we missed something
  _reason_alg2 = "vtx = 1";
  //recobVertex = TVector3(vtxs[0].x, vtxs[0].y, vtxs[0].z);
  if(mtrk_end_point.Z() < DOWNSTREAM_Z_CUT  && inActiveRegion(mtrk_end_point))return true;
  return false;
}

/**
 * @brief Check that the primary entered in the tpc, if it didn't, what did?
 * 
 * @return true Is our primary pion
 * @return false Is not our primary pion
 */
bool xs_ana::checkParticleInTpc()
{
  TVector3 firstPosInTpc(0,0,-100);
  bool enteredTpc(false);
  for (int iPt=0; iPt<MidPosX->at(0).size(); iPt++)
  {
    TVector3 pos(MidPosX->at(0)[iPt], MidPosY->at(0)[iPt], MidPosZ->at(0)[iPt]);
    if(!inTpcRegion(pos))continue;
    firstPosInTpc = pos;
    enteredTpc=true;
    break;
  }

  if(!enteredTpc)
  {
    // What did then?
    // make sure candidates are charged, and proj onto front face
    if(VERBOSE>=1)cout<<"Warning: Run #"<<run<<" Event #"<<event<<" primary did not enter tpc but a track was matched\n";
    if(VERBOSE>=1)cout<<"Primary ended at ("<<EndPointx->at(0)<<", " <<EndPointy->at(0)<<", " <<EndPointz->at(0)<<")\n";
    if(VERBOSE>=1)cout<<"Candidates...\n";
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
    if(VERBOSE>=1)for(const auto& c : cand){cout<<"\t"<<c<<endl;}
    
    // some particles take precedence
    if(cand.find("pion")!=cand.end()){_nPionBkg++;_trueTpcParticle="sec_pion";}
    else if(cand.find("muon")!=cand.end()){_nMuonBkg++;_trueTpcParticle="muon";}
    else if(cand.find("proton")!=cand.end()){_nProtonBkg++;_trueTpcParticle="proton";}
    else if(cand.find("electron")!=cand.end()){_nElectronBkg++;_trueTpcParticle="electron";}
    else {_nOtherBkg++;_trueTpcParticle="other";}
    if(VERBOSE>=1)cout<<"Candidate decided to be "<<_trueTpcParticle<<endl;
    
    // if there are no candidates, wth did we do wrong?
    if(!cand.size()){cerr<<"Error. No candidates for matched track!\n";exit(1);}
  }
  return enteredTpc;
}

/**
 * @brief This is where we fill our interacting and incident plots using 
 *        reconstructed information.
 * 
 */
void fillXsPlots(float                      recoKinEn,
                 const int&                 evtId, 
                 const std::vector<double>& trk_x, 
                 const std::vector<double>& trk_y, 
                 const std::vector<double>& trk_z,
                 const std::vector<double>& trk_p, 
                 const std::vector<double>& trk_dedx, 
                 const bool& isInel)
{
  float lastPosKinEn = recoKinEn;
  for (int iSb = 0; iSb<trk_x.size(); iSb++)
  {
    float dedx  = trk_dedx[iSb];
    float pitch = trk_p[iSb];
    // protection against large dedx and pitch
    //if(pitch>PITCH_MAX){cout<<"HEYYY "<<pitch<<" "<<iSb<<"/"<<trk_x.size()<<endl;continue;}
    if(dedx>DEDX_MAX)dedx=2.1; // mean value from dedx curve
    
    hIncidentKe->Fill(recoKinEn);
    if(_fillElasticBkg)hIncidentKeElasticBkg->Fill(recoKinEn);
    hDeDx->Fill(dedx);
    hPitch->Fill(pitch);
    hTrkZ->Fill(trk_z[iSb]);
    if(iSb!=(trk_x.size()-1))recoKinEn -= dedx * pitch;

    //keeping track of the last > 0 ke
    if(recoKinEn>0)lastPosKinEn=recoKinEn;
    
    // check for < 0
    if(recoKinEn<=0 && VERBOSE>=1){cout<<"Warning: "<<evtId<<" Filling zero\n";break;}
  }
  
  // fill interacting
  if(isInel)hInteractingKe->Fill(lastPosKinEn);
  if(_fillElasticBkg)hInteractingKeElasticBkg->Fill(lastPosKinEn);
}

/**
 * @brief Check energy reconstruction
 * 
 */
void xs_ana::makeSmearingPlot(const std::vector<double>& trk_x,
                              const std::vector<double>& trk_y, 
                              const std::vector<double>& trk_z, 
                              const std::vector<double>& trk_p, 
                              const std::vector<double>& trk_dedx, 
                              float recoKinEn, 
                              float trueKinEn)
{
  // for each trk point, find the nearest true point
  int iOldPt = 0;
  for(int iPt=1; iPt<trk_x.size(); iPt++, iOldPt++)
  {
    float dedx  = trk_dedx[iPt];
    float pitch = trk_p[iPt];
    if(dedx>DEDX_MAX)dedx = 2.1; // mean of dedx curve
    recoKinEn -= dedx * pitch;

    TVector3 oldPos(trk_x[iOldPt], trk_y[iOldPt], trk_z[iOldPt]);
    TVector3 currentPos(trk_x[iPt], trk_y[iPt], trk_z[iPt]);
    for(int iTrPt=0; iTrPt<TrackIdes_x->at(0).size(); iTrPt++)
    {
      TVector3 true_pos(TrackIdes_x->at(0)[iTrPt],TrackIdes_y->at(0)[iTrPt],TrackIdes_z->at(0)[iTrPt]);
      if(true_pos.Z() < oldPos.Z())continue;
      if(true_pos.Z() > currentPos.Z())continue;
      trueKinEn -= TrackIdes_e->at(0)[iTrPt];
    }
    hTrueVsRecoKinEn->Fill(trueKinEn, recoKinEn);

    float diff = std::abs(trueKinEn-recoKinEn);
    if(VERBOSE>=1 && diff > 100)cout<<"\nWarning: KE difference > 100"<<event<<endl;
    //if(event==223)cout<<oldPos.Z()<<" True = "<<trueKinEn<<" Reco = "<<recoKinEn<<" Diff = "<<diff<<" dedx = "<<dedx<<" pitch = "<<pitch<<endl;
  }
}

/**
 * @brief A place to make study mc information
 * 
 * @return true Skip this event
 * 
 */
void xs_ana::studyMc(const float& wcMom, const float& recoKinEn, const bool& isInel,
                     const bool& isInel_1, const bool& isInel_2,
                     const std::vector<double>& trk_x, const std::vector<double>& trk_y, const std::vector<double>& trk_z,
                     const std::vector<double>& trk_p, const std::vector<double>& trk_dedx)
{
  // find nearest pi-inelastic vertex
  TVector3 trk_startpos(trk_x.front(), trk_y.front(), trk_z.front());
  TVector3 trk_endpos(trk_x.back(), trk_y.back(), trk_z.back());
  float minDist = std::numeric_limits<float>::max();
  Vertex_t closest_vtx(-1, "none", trk_endpos);
  int closest_vtx_id(-1);

  if(VERBOSE>=3)cout<<"Inelastic interactions in tpc:\n";
  bool foundInel(false);
  for(int iVtx=0; iVtx<_vertices.size(); iVtx++)
  {
    if(_vertices[iVtx].process.find("pi-Inelastic")==std::string::npos)continue;
    foundInel = true;
    float dist = (trk_endpos-_vertices[iVtx].position).Mag();
    if(dist<minDist){minDist=dist; closest_vtx=_vertices[iVtx]; closest_vtx_id=iVtx;}
    if(VERBOSE>=3){cout<<"\t";_vertices[iVtx].Dump();}
  }
  if(VERBOSE>=3){cout<<"Closest (inelastic) vertex: "; closest_vtx.Dump();}

  // get first position in tpc (3cm from front face hit removal)
  TVector3 firstTruePos(0,0,-100);
  TVector3 lastTruePos(0,0,200);
  TVector3 initialMom(0,0,0);
  int hit_removal_low(3);
  int hit_removal_high(87);
  for (int iPt=0; iPt < MidPosX->at(0).size(); iPt++)
  {
    TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
    TVector3 true_mom(1000*MidPx->at(0)[iPt],1000*MidPy->at(0)[iPt],1000*MidPz->at(0)[iPt]);
    if(true_pos.Z()<0)continue;
    firstTruePos=true_pos;
    if(true_pos.Z()<hit_removal_low)continue;
    firstTruePos=true_pos;
    initialMom=true_mom;
    break;
  }
  // get last position in tpc (87cm from front face hit removal)
  minDist = std::numeric_limits<float>::max();
  for (int iPt=MidPosX->at(0).size()-1; iPt>=0; iPt--)
  {
    TVector3 true_pos(MidPosX->at(0)[iPt],MidPosY->at(0)[iPt],MidPosZ->at(0)[iPt]);
    if(!inTpcRegion(true_pos))continue;
    if(true_pos.Z()>hit_removal_high)continue;
    lastTruePos = true_pos;
    break;
  }
  if(firstTruePos.Z()<0 && VERBOSE>=1){cout<<event<<"Warning: Check true first position";PrintVec(firstTruePos);}
  if(lastTruePos.Z()>hit_removal_high && VERBOSE>=1){cout<<"Warning: Check true last position";PrintVec(lastTruePos);}
  if(firstTruePos.Z()>lastTruePos.Z() && VERBOSE>=1)cout<<"Warning: Check true positions\n";

  hDiffFirstPosInTpcX->Fill(   (trk_startpos-firstTruePos).X() );
  hDiffFirstPosInTpcY->Fill(   (trk_startpos-firstTruePos).Y() );
  hDiffFirstPosInTpcZ->Fill(   (trk_startpos-firstTruePos).Z() );
  hDiffFirstPosInTpcMag->Fill( (trk_startpos-firstTruePos).Mag() );

  auto last_interesting_point = closest_vtx.position;
  if(closest_vtx.process=="none") last_interesting_point = lastTruePos;

  float reco_length = (trk_endpos-trk_startpos).Mag();
  float true_length = (firstTruePos-last_interesting_point).Mag();
  hRecoLength->Fill(reco_length);
  hTrueLength->Fill(true_length);
  hDiffIntPosInTpcX->Fill(   (trk_endpos-last_interesting_point).X() );
  hDiffIntPosInTpcY->Fill(   (trk_endpos-last_interesting_point).Y() );
  hDiffIntPosInTpcZ->Fill(   (trk_endpos-last_interesting_point).Z() );
  hDiffIntPosInTpcMag->Fill( (trk_endpos-last_interesting_point).Mag() );
  hDiffLength->Fill( reco_length - true_length );
  hTruVsRecoLength->Fill( reco_length,  true_length );
  //if(40 < reco_length && reco_length<50 && true_length > 80)cout<<"\n\nsadfsadfsdafHEYYYYY"<<event<<endl;

  // compare the kinetic energies
  TVector3 startingMom(1000*MidPx->at(0)[0], 1000*MidPy->at(0)[0], 1000*MidPz->at(0)[0]);
  float startKe   = std::sqrt(startingMom*startingMom + PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;
  float initialKe = std::sqrt(initialMom*initialMom + PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;
  hEnLossUpstream->Fill(startKe-initialKe);
  hDiffInitKe->Fill(initialKe-recoKinEn);

  float merit = (trk_endpos-last_interesting_point).Mag();
  if(VERBOSE>=3){cout<<"First alg report...\n";}
  if(VERBOSE>=3){cout<<"\tThe reason "<<_reason_alg1<<endl;}
  if(VERBOSE>=3){if(isInel_1  && foundInel){if(merit < 5.0){cout<<"\tCorrect1!\n";}
                                            else{cout<<"\tMissed1 by merit\n";}}}
  if(VERBOSE>=3){if(isInel_1  && !foundInel){cout<<"\tDetermined but not inelastic!\n";}}
  if(VERBOSE>=3){if(!isInel_1 && foundInel ){cout<<"\tMissed!\n";}}
  if(VERBOSE>=3){if(!isInel_1 && !foundInel){cout<<"\tCorrect2!\n";}}
  if(VERBOSE>=3){cout<<"Second alg report...\n";}
  if(VERBOSE>=3){cout<<"\tThe reason "<<_reason_alg2<<endl;}
  if(VERBOSE>=3){if(isInel_2  && foundInel){if(merit < 5.0){cout<<"\tCorrect1!\n";}
                                            else{cout<<"\tMissed2 by merit\n";}}}
  if(VERBOSE>=3){if(isInel_2  && !foundInel){cout<<"\tDetermined but not inelastic!\n";}}
  if(VERBOSE>=3){if(!isInel_2 && foundInel ){cout<<"\tMissed!\n";}}
  if(VERBOSE>=3){if(!isInel_2 && !foundInel){cout<<"\tCorrect2!\n";}}
  if(VERBOSE>=3){cout<<"Ntracks = "<<ntracks_reco<<endl;}
  if(VERBOSE>=3){cout<<"Start point ("<<trk_x[0]<<","<<trk_y[0]<<","<<trk_z[0]<<")\n";}
  if(VERBOSE>=3){cout<<"End point   ("<<trk_x.back()<<","<<trk_y.back()<<","<<trk_z.back()<<")\n";}

  if(isInel_1 && foundInel){if(merit < 5.0){_nCorrect_alg1++;}
                            else{_nMissed_alg1++;}}
  if(isInel_1  && !foundInel){_nIncorrect_alg1++;}
  if(!isInel_1 && foundInel) {_nMissed_alg1++;}
  if(!isInel_1 && !foundInel){_nCorrect_alg1++;}
  if(isInel_2  && foundInel){if(merit > 5.0){_nCorrect_alg2++;}
                             else{_nMissed_alg2++;}}
  if(isInel_2  && !foundInel){_nIncorrect_alg2++;}
  if(!isInel_2 && foundInel ){_nMissed_alg2++;}
  if(!isInel_2 && !foundInel){_nCorrect_alg2++;}

  if(isInel    && foundInel){if(merit < 5.0){_n_correct++;}
                            else{_n_missed++;}}
  if(isInel    && !foundInel){_n_incorrect++;}
  if(!isInel   && foundInel ){_n_missed++;}
  if(!isInel   && !foundInel){_n_correct++;}

  //####################################
  // Area to check events 
  //if( true_length < 80 && reco_length > 80 ){cout<<"entry: "<<_jentry<<" event: "<<event<<endl;dumpEventInfo();}
  //####################################

  // Tag elastic background
  if(isInel && !foundInel)
  {
    if(VERBOSE>=1)cout<<"Warning: Tagged as inelastic but closest is "<<closest_vtx.process<<endl;
    _fillElasticBkg = true;
  }


  // make the true xs plots
  doTrueXs(initialMom, firstTruePos, lastTruePos);

  // how about the well reco tracks?
  bool wellReco = abs(reco_length-true_length) < 2.0 && abs(recoKinEn-initialKe)<10;
  if(wellReco)
  {
    hWellRecoLength->Fill(reco_length);
    hWellRecoTrueLength->Fill(true_length);
    hWellRecoDiffLength->Fill(reco_length-true_length);
    hWellRecoTruVsRecoLength->Fill(reco_length, true_length);

    doTrueXsWellReco(initialMom, firstTruePos, last_interesting_point, closest_vtx_id); 
    makeSmearingPlot(trk_x, trk_y, trk_z, trk_p, trk_dedx, recoKinEn, initialKe);
    
    // fill incident
    float recoKinEn_copy = recoKinEn;
    float lastPosKinEn = recoKinEn;
    for (int iSb = 0; iSb<trk_x.size(); iSb++)
    {
      float dedx  = trk_dedx[iSb];
      float pitch = trk_p[iSb];
      //cout<<"recoKinEn = " <<recoKinEn_copy<<" dedx = "<<dedx<<" pitch = "<<pitch<<endl;
      // protection against large dedx and pitch
      //if(pitch>PITCH_MAX){cout<<"HEYYY "<<pitch<<" "<<iSb<<"/"<<trk_x.size()<<endl;continue;}
      if(dedx>DEDX_MAX)dedx=2.1; // mean value from dedx curve
      
      hWellRecoIncidentKe->Fill(recoKinEn_copy);
      hDeDx->Fill(dedx);
      hPitch->Fill(pitch);
      if(iSb!=(trk_x.size()-1))recoKinEn_copy -= dedx * pitch;

      //keeping track of the last > 0 ke
      if(recoKinEn_copy>0)lastPosKinEn=recoKinEn_copy;
      
      // check for < 0
      if(recoKinEn_copy<=0 && VERBOSE>=1){cout<<"Warning: "<<event<<" Got zero at "<<iSb<<" / "<<trk_x.size()-1<<"\n";break;}
    }

    // fill interacting
    if(isInel)hWellRecoInteractingKe->Fill(lastPosKinEn);
  }
  return false;
}

/**
 * @brief Check for kinks in primary track 
 *  
 */
bool xs_ana::checkKink(int const& idMatchedTrk, TVector3& kink_position)
{
  if(VERBOSE>=2)cout<<"Checking the "<<track_kink_x->at(idMatchedTrk).size()<<" kinks...\n";
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
    if(VERBOSE>=2)cout<<"\tTrack inverted: "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<" ";
    auto temp = mtrk_start_point;
    mtrk_start_point = mtrk_end_point;
    mtrk_end_point = temp;
    if(VERBOSE>=2)cout<<" ==> "<<mtrk_start_point.Z()<<" "<<mtrk_end_point.Z()<<endl;
  }

  // get the kink information
  int nKinks = track_kink_x->at(idMatchedTrk).size();
  if(!nKinks)return false;

  // (old, I think this biases the result) check all kinks and use the maximum????
  // better to take the first kink that passes kink cut
  int kinkCount = 1;
  float max_angle(0);
  bool didPass(false);
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

    if(angle > ANGLE_CUT && angle > max_angle)
    {
      if(VERBOSE>=2)cout<<"\tKink with angle = "<<angle<<endl;
      max_angle = angle;
      kink_position = kinkPos;
      didPass = true;
      break;
    }
    kinkCount++;
  }
  if(didPass)
  {
    _reason_alg1=_reason_alg1+" kinkangle="+std::to_string(max_angle);
    _reason_alg2=_reason_alg2+" kinkangle="+std::to_string(max_angle);
    return true;
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
  if(beforeKinkPos.Z()>afterKinkPos.Z() && VERBOSE>=1)
  {cout<<"Warning: Check kink positions!\n";}
}

/**
 * @brief For debug
 * 
 */
void xs_ana::dumpEventInfo()
{
  cout<<"Interactions...\n";
  for(int iPt=0; iPt<InteractionPoint->size(); iPt++)cout<<"\tz = "<<MidPosZ->at(0)[InteractionPoint->at(iPt)]<<"  process = "<<InteractionPointType->at(iPt)<<endl;

  //cout<<"Daughters...\n";
  //for(int ig4 = 0; ig4 < geant_list_size; ig4++)
  //{
  //  if(std::abs(Mother->at(ig4))!=1)continue;
  //  cout<<"PDG = "<<PDG->at(ig4)<<" startz = "<<StartPointz->at(ig4)<<" process = "<<Process->at(ig4)<<endl;
  //}

}

/**
 * @brief Get angle between two tracks
 * @note This assumes the first id is the matched track.
 *       This also assumes we haven't shaved off points from the primary track.
 */
float xs_ana::getAngle(int const& idm, int const& ids, const TVector3& vtx)
{
  if(VERBOSE>=3){cout<<"\tgetting angle...\n";}
  int nTrkPts = col_track_x->at(idm).size();
  TVector3 trk1_sp( col_track_x->at(idm)[0],
                    col_track_y->at(idm)[0],
                    col_track_z->at(idm)[0] );
  TVector3 trk1_ep( col_track_x->at(idm)[nTrkPts-1],
                    col_track_y->at(idm)[nTrkPts-1],
                    col_track_z->at(idm)[nTrkPts-1] );

  // check inversion
  bool isInverted = false;
  if(trk1_sp.Z()>trk1_ep.Z())
  {
    if(VERBOSE>=2)cout<<"\tTrack inverted: "<<trk1_sp.Z()<<" "<<trk1_ep.Z()<<" ";
    auto temp = trk1_sp;
    trk1_sp = trk1_ep;
    trk1_ep = temp;
    if(VERBOSE>=2)cout<<" ==> "<<trk1_sp.Z()<<" "<<trk1_ep.Z()<<endl;
    isInverted=true;
  }
  // Redo... Using points near vertex
  int iPt=0;
  for(; iPt<nTrkPts; iPt++)
  {
    TVector3 pos(col_track_x->at(idm)[iPt], col_track_y->at(idm)[iPt], col_track_z->at(idm)[iPt]);
    if(isInverted)pos = TVector3(col_track_x->at(idm)[nTrkPts-1 - iPt], col_track_y->at(idm)[nTrkPts-1 - iPt], col_track_z->at(idm)[nTrkPts-1 - iPt]);
    if(pos.Z() > vtx.Z())break;
    iPt++;
  }
  // this is somewhat aribitrary
  // We want to use points "far" from the vertex in case of error in reconstruction
  // But close enough to eliminate MCS effects
  int npts = 10;
  int iSp = iPt-npts > 0 ? iPt-npts : 0;
  int iEp = iSp+3;
  trk1_sp = TVector3(col_track_x->at(idm)[iSp], col_track_y->at(idm)[iSp], col_track_z->at(idm)[iSp]);
  trk1_ep = TVector3(col_track_x->at(idm)[iEp], col_track_y->at(idm)[iEp], col_track_z->at(idm)[iEp]);
  if(isInverted)
  {
    trk1_sp = TVector3(col_track_x->at(idm)[nTrkPts-1 - iSp], col_track_y->at(idm)[nTrkPts-1 - iSp], col_track_z->at(idm)[nTrkPts-1 - iSp]);
    trk1_ep = TVector3(col_track_x->at(idm)[nTrkPts-1 - iEp], col_track_y->at(idm)[nTrkPts-1 - iEp], col_track_z->at(idm)[nTrkPts-1 - iEp]);
  }
  if(VERBOSE>=3){cout<<"\tsp = "<<iSp<<" ";PrintVec(trk1_sp);}
  if(VERBOSE>=3){cout<<"\tep = "<<iEp<<" ";PrintVec(trk1_ep);}


  auto dir1 = (trk1_ep-trk1_sp).Unit();

  int nSecTrkPts = col_track_x->at(ids).size();
  TVector3 trk2_sp( col_track_x->at(ids)[0],
                    col_track_y->at(ids)[0],
                    col_track_z->at(ids)[0] );
  TVector3 trk2_ep( col_track_x->at(ids)[nSecTrkPts-1],
                    col_track_y->at(ids)[nSecTrkPts-1],
                    col_track_z->at(ids)[nSecTrkPts-1] );
  // check inversion
  float dist1 = (trk1_ep - trk2_sp).Mag();
  float dist2 = (trk1_ep - trk2_ep).Mag();
  isInverted = false;
  if(dist2<dist1)
  {
    if(VERBOSE>=2)cout<<"\tSec track inverted: "<<trk2_sp.Z()<<" "<<trk2_ep.Z()<<" ";
    auto temp = trk2_sp;
    trk2_sp = trk2_ep;
    trk2_ep = temp;
    if(VERBOSE>=2)cout<<" ==> "<<trk2_sp.Z()<<" "<<trk2_ep.Z()<<endl;
    isInverted=true;
  }
  int extend = 3 < nSecTrkPts ? 3 : 2;
  trk2_ep = TVector3(col_track_x->at(ids)[extend], col_track_y->at(ids)[extend], col_track_z->at(ids)[extend]);
  if(isInverted)
  {
    trk2_ep = TVector3(col_track_x->at(ids)[nSecTrkPts-1 - extend], col_track_y->at(ids)[nSecTrkPts-1 - extend], col_track_z->at(ids)[nSecTrkPts-1 - extend]);
  }
  auto dir2 = (trk2_ep-trk2_sp).Unit();
  return std::acos(dir1.Dot(dir2))*180./M_PI;
}

/**
 * @brief Getting interactions in tpc
 * @note We should've only filled the interactions with 
 *       what g4 spat out. We handle the zero case here.
 * @todo Do we need to check the end point still?
 */
void xs_ana::getInteractionsInTpc()
{
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
    if(VERBOSE>=2)cout<<"Checking end process and daughters...\n";
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
  bool foundInel(false);
  for(const auto& v :_vertices)
  {
    if(v.process.find("pi-Inelastic") != std::string::npos)foundInel=true;
  }
  if(foundInel)tagInelasticChannel();
}

/**
 * @brief Tagging the inelastic channel: quasi, abs, CX
 * 
 */
void xs_ana::tagInelasticChannel()
{
  //
  // We define sub channels for inelastic processes:
  //      1) Absorption (no pions in final state)
  //      2) Quasi (no neutral pions, 1 same charged pion in final state)
  //      3) Charge exchange (1 nuetral pion, no charged pions in final state)
  //      4) Pion production (what's left?)
  //

  bool isPrimaryChargedPion(false);
  TVector3 primaryFinalPosition;
  int primTrkId(-1);

  // Loop over the g4 particles
  for (size_t ig4 = 0; ig4 < geant_list_size; ig4++)
  {
    if (!process_primary->at(ig4)) continue;
    primTrkId = TrackId->at(ig4);
    // if charged pion
    if (std::abs(PDG->at(ig4)) == PdgCodes_t::kPion)
    {
      // change flag
      isPrimaryChargedPion = true;
      // get last point
      primaryFinalPosition = TVector3(MidPosX->at(ig4).back(), MidPosY->at(ig4).back(), MidPosZ->at(ig4).back());
    } 
  }
  if (!isPrimaryChargedPion)return;

  // Loop over the particles again checking daughters
  for (size_t ig4 = 0; ig4 < geant_list_size; ig4++)
  {
    // we only care about particles that start in the tpc
    TVector3 dtrStartPos(StartPointx->at(ig4), StartPointy->at(ig4), StartPointz->at(ig4)); 
    if (!inTpcRegion(dtrStartPos)) continue;
    //
    //// skip if particle is not a child of the primary
    if (Mother->at(ig4) != primTrkId) continue;
    // get pdg code
    int pdgCode = std::abs(PDG->at(ig4));

    // which vertex is she attached to
    float minDist = std::numeric_limits<float>::max();
    int index = -1;
    for(int iVtx = 0; iVtx < _vertices.size(); iVtx++)
    {
      if((_vertices[iVtx].position-dtrStartPos).Mag()<0.01)
      {
        index=iVtx;
        break;
      }
    }
    if(index>=0)
    {
      _vertices[index].secondaries.push_back(pdgCode);
    }
  }

  for(int iVtx = 0; iVtx < _vertices.size(); iVtx++)
  {
    if(_vertices[iVtx].process.find("pi-Inelastic") == std::string::npos)continue;
    bool isTheEnd = (primaryFinalPosition-_vertices[iVtx].position).Mag() < 0.01 ? true : false;
    _vertices[iVtx].subtype = "error";
    int nChargedPionDaughters(0), nNeutralPionDaughters(0),
        nProtonDaughters(0),      nNeutronDaughters(0),
        nMuonDaughters(0),        nKaonDaughters(0);

    // add up all the secondaries at this vertex
    for(const auto& sec : _vertices[iVtx].secondaries)
    {
      // we will classify the interaction based on pions, kaons, muons, protons, and neutrons
      switch(sec)
      {
        case PdgCodes_t::kPi0:{ nNeutralPionDaughters++; break; }
        case PdgCodes_t::kPion:{ nChargedPionDaughters++; break; }
        case PdgCodes_t::kProton:{ nProtonDaughters++; break; }
        case PdgCodes_t::kNeutron:{ nNeutronDaughters++; break; }
        case PdgCodes_t::kKaon:{ nKaonDaughters++; break; }
        case PdgCodes_t::kMuon:{ nMuonDaughters++; break; }
      }
    }
    if(VERBOSE>=3)
    {
      std::cout << "\nNumber of protons:           " << nProtonDaughters
                << "\nNumber of neutrons:          " << nNeutronDaughters
                << "\nNumber of charged pions:     " << nChargedPionDaughters
                << "\nNumber of neutral pions:     " << nNeutralPionDaughters
                << "\nNumber of charged kaons:     " << nKaonDaughters
                << "\nNumber of muons:             " << nMuonDaughters
                << "\nPrimary end vertex: (" 
                << primaryFinalPosition.X() << ", " 
                << primaryFinalPosition.Y() << ", "
                << primaryFinalPosition.Z() << ")\n";
    }
    // Handle absorption and inelastic
    if (nChargedPionDaughters == 0 && nNeutralPionDaughters == 0)_vertices[iVtx].subtype = "PionAbsorption";
    // Handle charge exchange
    else if (nNeutralPionDaughters == 1 && nChargedPionDaughters == 0)_vertices[iVtx].subtype = "ChargeExchange";
    // Special case
    else if (nChargedPionDaughters == 1 && nNeutralPionDaughters == 0)_vertices[iVtx].subtype = "QuasiElastic";
    // Everything else should be pion production
    else {_vertices[iVtx].subtype = "PionProduction";}
  }
}

/**
 * @brief Fill true incident and interacting histograms for well reco events
 * @todo True pitch histogram shows > 1000 enrtries
 *       with < 0.47. 
 *
 */
void xs_ana::doTrueXsWellReco(const TVector3& initialMom, 
                              const TVector3& firstPos, 
                              const TVector3& lastPos, 
                              const int& closest_vtx_id)
{
  // for now we will assume this is a straight track
  // get unifomly spaced points to the last pos
  std::vector<TVector3> orderedPoints;
  bool isInel = false;
  if(closest_vtx_id>=0)
  {
    TVector3 current_pos = firstPos;
    for(int iVtx=0; iVtx<=closest_vtx_id; iVtx++)
    {
      if(_vertices[iVtx].position.Z()>lastPos.Z())break;
      orderedPoints.push_back(current_pos);
      // adding in range (current, v.pos)
      addUniformPoints(current_pos, _vertices[iVtx].position, orderedPoints);
      current_pos = _vertices[iVtx].position;
    }
    // add last point
    orderedPoints.push_back(current_pos);
    if(_vertices[closest_vtx_id].process.find("pi-Inelastic")!=std::string::npos)isInel=true;
  }
  else 
  {
    orderedPoints.push_back(firstPos);
    addUniformPoints(firstPos, lastPos, orderedPoints);
    orderedPoints.push_back(lastPos);
  }
  // remove second to last point if it's too close to the last point
  auto secondToLastPoint = orderedPoints.end()-2;
  float trackPitch = 0.47;
  int distance = std::distance(orderedPoints.begin(), secondToLastPoint);
  if( (*secondToLastPoint-orderedPoints.back()).Mag() < trackPitch*0.5 )orderedPoints.erase(secondToLastPoint);

  // now we can fill
  double kinEn = std::sqrt(initialMom.Mag()*initialMom.Mag()+PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;

  int oldPt = 0;
  for (int iPt=1; iPt<orderedPoints.size(); iPt++, oldPt++)
  {
    auto oldPos     = orderedPoints[oldPt];
    auto currentPos = orderedPoints[iPt];

    double eDep = 0;
    for(int iIde=0; iIde<TrackIdes_x->at(0).size(); iIde++)
    {
      if(TrackIdes_z->at(0)[iIde]<oldPos.Z())continue;
      if(TrackIdes_z->at(0)[iIde]>currentPos.Z())continue;
      eDep += TrackIdes_e->at(0)[iIde];
    }
    kinEn -= eDep;

    hWellRecoTruePitch->Fill((oldPos-currentPos).Mag());
    hWellRecoTrueEnDep->Fill(eDep);
    hWellRecoTrueIncidentKe->Fill(kinEn);
  }
  if(isInel)hWellRecoTrueInteractingKe->Fill(kinEn);
}

/**
 * @brief Fill true incident and interacting histograms
 *
 */
void xs_ana::doTrueXs(const TVector3& initialMom, 
                      const TVector3& firstPos, 
                      const TVector3& lastPos)
{

  // get unifomly spaced points to the last pos
  std::vector<TVector3> orderedPoints;
  bool isInel = false;
  if(_vertices.size())
  {
    TVector3 current_pos = firstPos;
    for(int iVtx=0; iVtx<_vertices.size(); iVtx++)
    {
      orderedPoints.push_back(current_pos);
      // adding in range (current, v.pos)
      addUniformPoints(current_pos, _vertices[iVtx].position, orderedPoints);
      current_pos = _vertices[iVtx].position;

      // we're finished if this was an inelastic process
      if(_vertices[iVtx].process.find("pi-Inelastic") != std::string::npos){isInel=true;break;}
    }
    // add last point
    orderedPoints.push_back(current_pos);
  }
  else 
  {
    orderedPoints.push_back(firstPos);
    addUniformPoints(firstPos, lastPos, orderedPoints);
    orderedPoints.push_back(lastPos);
  }
  // remove second to last point if it's too close to the last point
  auto secondToLastPoint = orderedPoints.end()-2;
  float trackPitch = 0.47;
  int distance = std::distance(orderedPoints.begin(), secondToLastPoint);
  if( (*secondToLastPoint-orderedPoints.back()).Mag() < trackPitch*0.5 )orderedPoints.erase(secondToLastPoint);

  // now we can fill
  double kinEn = std::sqrt(initialMom.Mag()*initialMom.Mag()+PARTICLE_MASS*PARTICLE_MASS)-PARTICLE_MASS;

  int oldPt = 0;
  for (int iPt=1; iPt<orderedPoints.size(); iPt++, oldPt++)
  {
    auto oldPos     = orderedPoints[oldPt];
    auto currentPos = orderedPoints[iPt];

    double eDep = 0;
    for(int iIde=0; iIde<TrackIdes_x->at(0).size(); iIde++)
    {
      if(TrackIdes_z->at(0)[iIde]<oldPos.Z())continue;
      if(TrackIdes_z->at(0)[iIde]>currentPos.Z())continue;
      eDep += TrackIdes_e->at(0)[iIde];
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
void xs_ana::addUniformPoints(const TVector3& x0, const TVector3& xf, std::vector<TVector3>& orderedPoints)
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
  hWellRecoInteractingKe->Write();
  hWellRecoIncidentKe->Write();
  hTrueIncidentKe->Write();
  hTrueInteractingKe->Write();
  hWellRecoTrueIncidentKe->Write();
  hWellRecoTrueInteractingKe->Write();
  hIncidentKeElasticBkg->Write();
  hInteractingKeElasticBkg->Write();

  hTrkZ->Write();
  hDeDx->Write();
  hPitch->Write();
  hRecoLength->Write();

  hWellRecoLength->Write();
  hWellRecoTrueLength->Write();
  hWellRecoDiffLength->Write();
  hWellRecoTruVsRecoLength->Write();

  hEnLossUpstream->Write();
  hVertexDistSec->Write();
  hTrueLength->Write();
  hKeWc4Decay->Write();
  hKeWc4Capture->Write();
  hTruePitch->Write();
  hTrueEnDep->Write();
  hWellRecoTruePitch->Write();
  hWellRecoTrueEnDep->Write();

  hDiffInitKe->Write();
  hDiffFirstPosInTpcX->Write();
  hDiffFirstPosInTpcY->Write();
  hDiffFirstPosInTpcZ->Write();
  hDiffFirstPosInTpcMag->Write();
  hDiffIntPosInTpcX->Write();
  hDiffIntPosInTpcY->Write();
  hDiffIntPosInTpcZ->Write();
  hDiffIntPosInTpcMag->Write();
  hDiffLength->Write();
  hTruVsRecoLength->Write();
  hTrueVsRecoKinEn->Write();

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




