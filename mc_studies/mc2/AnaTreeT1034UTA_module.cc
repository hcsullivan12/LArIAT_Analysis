// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
//#include "larcorealg/Geometry/CryostatGeo.h"
//#include "larcorealg/Geometry/TPCGeo.h"
//#include "larcorealg/Geometry/PlaneGeo.h"
//#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
//#include "lardata/RecoBaseArt/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTOF        = 100;      // :/
const int kMaxWCTracks   = 1000;
const int kMaxTrack      = 1000;     //maximum number of tracks
const int kMaxHits       = 20000;    //maximum number of hits
const int kMaxTrackHits  = 1000;     //maximum number of space points
const int kMaxHitIDs     = 100;      //maximum number of space points ids
const int kMaxTrajHits   = 1000;     //maximum number of trajectory points
const int kMaxCluster    = 1000;     //maximum number of clusters
const int kMaxPrimaryPart = 10;	     //maximum number of true primary particles
const int kMaxDaughterPart = 100;     //maximum number of true daughter particles
const int kMaxPrimaries  = 20000;    //maximum number of true particles tracked
const int kMaxTruePrimaryPts = 5000; //maximum number of points in the true primary trajectory 
const int kMaxTrueDaughterPts = 5000; //maximum number of points in the true daughter trajectory 
const int kMaxIDE = 5000; //maximum number of points in the true primary trajectory 

namespace lariat 
{
  class AnaTreeT1034UTA;
}

class lariat::AnaTreeT1034UTA : public art::EDAnalyzer 
{
public:
  explicit AnaTreeT1034UTA(fhicl::ParameterSet const & p);
  virtual ~AnaTreeT1034UTA();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void reconfigure(fhicl::ParameterSet const & p);


private:

  // === Function used to reset all the variables  ===
  void ResetVars();
  
  // === Storing information into TTree ====
  TTree* fTree;
  // === Histos! ===
  TH2D* E_vs_e;
   
  //=== Storing Run Information ===
  int run;			//<---Run Number
  int subrun;			//<---SubRun Number
  int event;			//<---Event Number
  double evttime;		//<---Event Time Stamp
  double efield[3];		//<---Electric Field 
  int t0;
   
  // === Storing Track Information ===
  int ntracks_reco;		//<---Number of reconstructed tracks
  std::vector<int>    track_WC2TPC_match;
  std::vector<double> track_start_x;
  std::vector<double> track_start_y;
  std::vector<double> track_start_z;
  std::vector<double> track_end_x;
  std::vector<double> track_end_y;
  std::vector<double> track_end_z;
  
  std::vector<double> track_length;
  std::vector<int> track_primary;	//<---Is this particle primary (primary = 1, non-primary = 0)

  // === Storing the tracks SpacePoints (Individual 3D points)
  std::vector<int> ntrack_hits;
  std::vector< std::vector<double> > track_xpos;
  std::vector< std::vector<double> > track_ypos;
  std::vector< std::vector<double> > track_zpos;
  std::vector< std::vector<int> > nhit_ids;

  // === Storing track kinks
  std::vector<std::vector<double>> track_kink_x;
  std::vector<std::vector<double>> track_kink_y;
  std::vector<std::vector<double>> track_kink_z;

  // === Storing vertices
  std::vector<std::vector<int>> vertex_track_ids;
  std::vector<double> vertex_x;
  std::vector<double> vertex_y;
  std::vector<double> vertex_z;

  // === Storing the tracks Calorimetry Information
  std::vector<int> ind_track_hits;
  std::vector<double> ind_track_ke;
  std::vector< std::vector<double> > ind_track_wire;
  std::vector< std::vector<double> > ind_track_dedx;
  std::vector< std::vector<double> > ind_track_dqdx;
  std::vector< std::vector<double> > ind_track_rr;
  std::vector< std::vector<double> > ind_track_pitch_hit;
  std::vector<int> col_track_hits;
  std::vector<double> col_track_ke;
  std::vector< std::vector<double> > col_track_x;
  std::vector< std::vector<double> > col_track_y;
  std::vector< std::vector<double> > col_track_z;
  std::vector< std::vector<double> > col_track_wire;
  std::vector< std::vector<double> > col_track_dedx;
  std::vector< std::vector<double> > col_track_dqdx;
  std::vector< std::vector<double> > col_track_rr;
  std::vector< std::vector<double> > col_track_pitch_hit;


  // === hit info ===
  int nhits;
  std::vector<double> hit_time;
  std::vector<double> hit_wire;
  std::vector<double> hit_view;
  std::vector<double> hit_amp;
  std::vector<double> hit_charge;

  std::vector<int>         InteractionPoint;         //<---Geant 4 Primary Trj Point Corresponding to the Interaction
  std::vector<std::string> InteractionPointType;     //<---Geant 4 Primary Interaction Type


  // === Storing Geant4 MC Truth Information ===
  int no_primaries;				//<---Number of primary Geant4 particles in the event
  int geant_list_size;			//<---Number of Geant4 particles tracked
  double primary_p;				//<---Primary particle momentum
  std::vector<int> PDG;
  std::vector<double> StartPointx; 
  std::vector<double> StartPointy;
  std::vector<double> StartPointz;
  std::vector<double> StartEnergy;
  std::vector<double> StartKE;
  std::vector<double> LastKE;
  std::vector<double> StartPx;
  std::vector<double> StartPy;
  std::vector<double> StartPz;

  std::vector<double> EndPointx; 
  std::vector<double> EndPointy;
  std::vector<double> EndPointz;
  std::vector<double> EndEnergy;
  std::vector<double> EndPx;
  std::vector<double> EndPy;
  std::vector<double> EndPz;

  std::vector<std::vector<double>> TrackIdes_x;
  std::vector<std::vector<double>> TrackIdes_y;
  std::vector<std::vector<double>> TrackIdes_z;
  std::vector<std::vector<double>> TrackIdes_e;
  
  std::vector<std::string> Process;

  std::vector<int> NumberDaughters;
  std::vector<int> TrackId;
  std::vector<int> Mother;
  std::vector<int> process_primary;	//<---Is this particle primary (primary = 1, non-primary = 0)
  std::vector<std::string> G4Process;         //<---The process which created this particle
  std::vector<std::string> G4FinalProcess;    //<---The last process which this particle went under

  // === Storing additional Geant4 MC Truth Information for the primary track only ===	   
  std::vector<int> NTrTrajPts;
  std::vector< std::vector<double> > MidPosX;
  std::vector< std::vector<double> > MidPosY;
  std::vector< std::vector<double> > MidPosZ;
  std::vector< std::vector<double> > MidPx;
  std::vector< std::vector<double> > MidPy;
  std::vector< std::vector<double> > MidPz;


  std::vector<double> NDTrTrajPts;
  int NProtonDaughters = 0;
  int NNeutronDaughters = 0;
  std::vector<int> DTrackId;
  std::vector<int> DPdgCode;
  std::vector<double> DStartKE;
  std::vector<double> DStartEnergy;
  std::vector<double> DStartP;
  std::vector< std::vector<double> > DMidPosX;
  std::vector< std::vector<double> > DMidPosY;
  std::vector< std::vector<double> > DMidPosZ;
  // === Storing additional Geant4 MC Truth Information for the daughter tracks ===  
   

  // === beamline info ===
  int num_tof_objects;
  double tofObject[kMaxTOF];

  int num_wctracks;
  double wctrk_momentum[kMaxWCTracks];
  double wctrk_XFace[kMaxWCTracks];
  double wctrk_YFace[kMaxWCTracks];
  double wctrk_theta[kMaxWCTracks];
  double wctrk_phi[kMaxWCTracks];
  int wctrk_missed[kMaxWCTracks];
  int wctrk_picky[kMaxWCTracks];
  int wctrk_quality[kMaxWCTracks];
  double wctrk_x_proj_3cm[kMaxWCTracks];
  double wctrk_y_proj_3cm[kMaxWCTracks];
  double wctrk_z_proj_3cm[kMaxWCTracks];
  double wctrk_wc1_x[kMaxWCTracks];
  double wctrk_wc1_y[kMaxWCTracks];
  double wctrk_wc1_z[kMaxWCTracks];
  double wctrk_wc2_x[kMaxWCTracks];
  double wctrk_wc2_y[kMaxWCTracks];
  double wctrk_wc2_z[kMaxWCTracks];
  double wctrk_wc3_x[kMaxWCTracks];
  double wctrk_wc3_y[kMaxWCTracks];
  double wctrk_wc3_z[kMaxWCTracks];
  double wctrk_wc4_x[kMaxWCTracks];
  double wctrk_wc4_y[kMaxWCTracks];
  double wctrk_wc4_z[kMaxWCTracks];

  double electron_lifetime;

  std::string fTreeName;
  std::string fClusterModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fG4ModuleLabel;
  std::string fTOFModuleLabel;
  std::string fWCTrackLabel;
  std::string fWC2TPCModuleLabel;
  std::string fWCQualityProducerLabel;

  calo::CalorimetryAlg fCalorimetryAlg;

};


lariat::AnaTreeT1034UTA::AnaTreeT1034UTA(fhicl::ParameterSet const & pset) 
  : EDAnalyzer(pset)
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(pset);
}

lariat::AnaTreeT1034UTA::~AnaTreeT1034UTA()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034UTA::reconfigure(fhicl::ParameterSet const & pset)
{
  fTreeName                 = pset.get< std::string >("TreeName", "anatreeuc");

  fHitsModuleLabel          = pset.get< std::string >("HitsModuleLabel");
  fTrackModuleLabel         = pset.get< std::string >("TrackModuleLabel");
  fCalorimetryModuleLabel   = pset.get< std::string >("CalorimetryModuleLabel");
  fParticleIDModuleLabel    = pset.get< std::string >("ParticleIDModuleLabel");
  fClusterModuleLabel       = pset.get< std::string >("ClusterModuleLabel");
  fG4ModuleLabel            = pset.get< std::string >("G4ModuleLabel");

  fTOFModuleLabel           = pset.get< std::string >("TOFModuleLabel");
  fWCTrackLabel             = pset.get< std::string >("WCTrackLabel");
  fWC2TPCModuleLabel        = pset.get< std::string >("WC2TPCModuleLabel", "WC2TPCtrk");
  fWCQualityProducerLabel   = pset.get< std::string >("WCQualityProducerLabel", "wcquality");

  return;
}

void lariat::AnaTreeT1034UTA::analyze(art::Event const & evt)
{
  // #############################################
  // ### Reset variables before we get started ###
  // #############################################
  ResetVars();
  //std::cout<<"Check1"<<std::endl;

  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
  // === Liquid Argon Properties Services ===
  //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  // === Detector properties service ===
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // === BackTracker service ===
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
   
   
  // === Run Number ===
  run = evt.run();
  // === Sub-Run Number ===
  subrun = evt.subRun();
  // === Event Number ===
  event = evt.id().event();
   
  std::cout<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout << "Run = "         << run 
            << ", SubRun = "    << subrun 
            << ", Evt = "       << event    << std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<std::endl;
   
  // === Event Time ===
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();

  // === Electron lifetime ===
  electron_lifetime = detprop->ElectronLifetime();

  // === Electric Field ===
  // Note: LArProperties::Efield() has moved to DetectorProperties/DetectorPropertiesService
  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);
   
  // === Trigger Offset ====
  t0 = detprop->TriggerOffset();
   
  // #####################################
  // ### Getting the Track Information ###
  // #####################################
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
   
  // === Filling the tracklist from the tracklistHandle ===
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}

  // ###################################
  // ### Getting kink information ###
  // ################################### 
  art::Handle< std::vector<recob::Vertex> > kinkListHandle;
  std::vector<art::Ptr<recob::Vertex> > kinklist; 
  if (evt.getByLabel(fTrackModuleLabel,"kink",kinkListHandle))
    {art::fill_ptr_vector(kinklist, kinkListHandle);}

  // ###################################
  // ### Getting vertex information ###
  // ################################### 
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  std::vector<art::Ptr<recob::Vertex> > vertexlist; 
  if (evt.getByLabel(fTrackModuleLabel,vertexListHandle))
    {art::fill_ptr_vector(vertexlist, vertexListHandle);}

  // ###################################
  // ### Getting the Hit Information ###
  // ###################################
  art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::hits objects
  std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define hitklist as a pointer to recob::hit
   
  // === Filling the hitlist from the hitlistHandle ===
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}


  // #############################
  // ### beam line information ###
  // #############################
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  std::vector< art::Ptr<ldp::TOF> > tof;
  if(evt.getByLabel(fTOFModuleLabel, TOFColHandle))
    {art::fill_ptr_vector(tof, TOFColHandle);}

  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {art::fill_ptr_vector(wctrack, wctrackHandle);}

  // art::ValidHandle< std::vector< recob::PFParticle > >
  auto wcquality_handle = evt.getValidHandle< std::vector< recob::PFParticle > >
      (fWCQualityProducerLabel);

   
  // ##########################################################
  // ### Grabbing associations for use later in the AnaTool ###
  // ##########################################################
  //std::cout<<"Check2"<<std::endl;
  // === Associations between hits and raw digits ===
  art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
  // === Association between SpacePoints and Tracks ===
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  // === Association between Tracks and 2d Hits ===
  art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
  // === Association between Calorimetry objects and Tracks ===
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  // === Association between Particle ID objects (PID) and Tracks ===
  art::FindManyP<anab::ParticleID>  fmpid(trackListHandle, evt, fParticleIDModuleLabel);
  // ==== Association between Tracks and Hits
  art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
  // === Association between vertices and Tracks ===
  art::FindManyP<recob::Track> fmvtx(vertexListHandle, evt, fTrackModuleLabel);
 
  // Fill vertex information
  if(fmvtx.isValid())
  {
    for(int iVtx=0; iVtx<(int)vertexListHandle->size(); iVtx++)
    {
      std::cout<<"\nNEW VERTEX"<<std::endl;
      std::cout<<vertexlist.at(iVtx)->position().Z()<<std::endl;
      auto trks = fmvtx.at(iVtx);
      std::vector<int> trk_ids;
      for(int iTrk=0; iTrk<(int)trks.size(); iTrk++)trk_ids.push_back(trks.at(iTrk)->ID());
      
      vertex_track_ids.push_back(trk_ids);
      vertex_x.push_back(vertexlist.at(iVtx)->position().X());
      vertex_y.push_back(vertexlist.at(iVtx)->position().Y());
      vertex_z.push_back(vertexlist.at(iVtx)->position().Z());
      std::cout<<trk_ids.size()<<std::endl;
    }
  }



   
  // ###################################################################
  // ### Setting a boolian to only output MC info if this is MC-info ###
  // ###################################################################
  bool isdata = false;
  if (evt.isRealData())
    {isdata = true;}
	
  else isdata = false;
   

   
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE MCTruth Geant4 INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isdata)
    {
    // ### Making a vector of MCParticles ###
    std::vector<const simb::MCParticle* > geant_part;
    for(size_t p = 0; p < plist.size(); ++p) 
	  {geant_part.push_back(plist.Particle(p));}
    
    // ### Setting a string for primary ###
    std::string pri("primary");
    
    int primary=0;
    int geant_particle=0;
    int numTrue_inTPC = 0;
    
    // ### Determine the number of primary particles from geant ###
    for( unsigned int i = 0; i < geant_part.size(); ++i )
	  {
	  geant_particle++;
	  if(geant_part[i]->Process()==pri)
	    {primary++;}
	  }//<---End i loop
	
    // ### Saving the number of primary particles ###
    no_primaries=primary;
    // ### Saving the number of Geant4 particles ###
    geant_list_size=geant_particle;
     
    // ### Looping over all the Geant4 particles ###
    int iPrim = 0;
    //int iDaught = 0;
    std::cout<<"geant par size: "<<geant_part.size()<<std::endl;
    for(unsigned int i = 0; i < geant_part.size(); ++i)
      {
        // Fill track ides from backtracker
        auto simIDE = bt_serv->TrackIdToSimIDEs_Ps(geant_part[i]->TrackId(), geom->View(0));
        std::vector<double> _ides_x;
        std::vector<double> _ides_y;
        std::vector<double> _ides_z;
        std::vector<double> _ides_e;
        for (const auto& ide : simIDE) 
        {
          _ides_x.push_back(ide->x); 
          _ides_y.push_back(ide->y); 
          _ides_z.push_back(ide->z); 
          _ides_e.push_back(ide->energy); 
        }
        TrackIdes_x.push_back(_ides_x);
        TrackIdes_y.push_back(_ides_y);
        TrackIdes_z.push_back(_ides_z);
        TrackIdes_e.push_back(_ides_e);


	  if(geant_part[i]->Process()==pri){process_primary.push_back(1);}
	  else                             {process_primary.push_back(0);}
	   
	  Process.push_back(geant_part[i]->Process());
	  
      // ### Saving other info ###
      PDG.push_back(geant_part[i]->PdgCode());
      Mother.push_back(geant_part[i]->Mother());
      TrackId.push_back(geant_part[i]->TrackId());
      StartEnergy.push_back(geant_part[i]->E());
	  EndEnergy.push_back(geant_part[i]->EndE());
      double startp = sqrt( pow(geant_part[i]->Px(),2) + pow(geant_part[i]->Py(),2) + pow(geant_part[i]->Pz(),2));
      double mass = geant_part[i]->Mass(); 
      double startke = sqrt( pow(startp,2) + pow(mass,2) ) - mass;
      StartKE.push_back(startke);
      double lastke = 0;
	  simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
      int counter_ii = 0;
	  for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj) {
        if(truetraj.Z(counter_ii) > 0   && truetraj.Z(counter_ii) < 90   &&
          truetraj.X(counter_ii) > 0   && truetraj.X(counter_ii) < 47.5 &&
          truetraj.Y(counter_ii) > -20 && truetraj.Y(counter_ii) < 20) {
          double thisp;
          if(counter_ii) {
            thisp = sqrt( pow(truetraj.Px(counter_ii-1),2) + 
                          pow(truetraj.Py(counter_ii-1),2) + 
                          pow(truetraj.Pz(counter_ii-1),2));
          }
          else {
            thisp = sqrt( pow(truetraj.Px(counter_ii),2) + 
                          pow(truetraj.Py(counter_ii),2) + 
                          pow(truetraj.Pz(counter_ii),2));
          }
          lastke = sqrt( pow(thisp,2) + pow(mass,2) ) - mass;
        }//<--End if in tpc
        counter_ii++;
	  }//<--End loop on true trajectory points
      LastKE.push_back(lastke);


	   
	  // ### Saving the start and end Px, Py, Pz info ###
	  StartPx.push_back(geant_part[i]->Px());
	  StartPy.push_back(geant_part[i]->Py());
	  StartPz.push_back(geant_part[i]->Pz());
	  EndPx.push_back(geant_part[i]->EndPx());
	  EndPy.push_back(geant_part[i]->EndPy());
	  EndPz.push_back(geant_part[i]->EndPz());
	   
      //std::cout << "    p(x,y,z): (" << geant_part[i]->Px() << ", " << geant_part[i]->Py() << ", " << geant_part[i]->Px() << ")     " 
      //          << "    start E: " << geant_part[i]->E() << std::endl;
	  // ### Saving the Start and End Point for this particle ###
	  StartPointx.push_back(geant_part[i]->Vx());
	  StartPointy.push_back(geant_part[i]->Vy());
	  StartPointz.push_back(geant_part[i]->Vz());
	  EndPointx.push_back(geant_part[i]->EndPosition()[0]);
	  EndPointy.push_back(geant_part[i]->EndPosition()[1]);
	  EndPointz.push_back(geant_part[i]->EndPosition()[2]);

	  // ### Saving the processes for this particle ###
	  G4Process.push_back( geant_part[i]->Process() );
	  G4FinalProcess.push_back( geant_part[i]->EndProcess() );
 	   
	  // ### Saving the number of Daughters for this particle ###
	  NumberDaughters.push_back(geant_part[i]->NumberDaughters());

	  // ### Save intermediary information for the primary track
	  if(geant_part[i]->Process() == pri)
        {
        //std::cout << "    p(x,y,z): (" << geant_part[i]->Px() << ", " << geant_part[i]->Py() << ", " << geant_part[i]->Px() << ")" << std::endl; 
        //std::cout << "    sqrt(p2):  " << sqrt( pow(geant_part[i]->Px(), 2) + pow(geant_part[i]->Py(), 2) + pow(geant_part[i]->Pz(), 2) ) << std::endl;
        //std::cout << "    p       :  " << geant_part[i]->P() << std::endl;
        primary_p = geant_part[i]->P();
	    NTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
	    //simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
        std::vector<double> midx;  std::vector<double> midy;  std::vector<double> midz;
        std::vector<double> midpx; std::vector<double> midpy; std::vector<double> midpz;
	    int iPrimPt = 0;	
        double true_total_distance = 0;
        double previous_mid_xpt = geant_part[i]->Vx();    
        double previous_mid_ypt = geant_part[i]->Vy();    
        double previous_mid_zpt = geant_part[i]->Vz();    
	    for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	      {
          // ### true pt vars ###
          double mid_xpt = truetraj.X(iPrimPt);
          double mid_ypt = truetraj.Y(iPrimPt);
          double mid_zpt = truetraj.Z(iPrimPt);
          double true_distance = sqrt( pow(mid_xpt - previous_mid_xpt, 2) 
                                     + pow(mid_ypt - previous_mid_ypt, 2) 
                                     + pow(mid_zpt - previous_mid_zpt, 2));
          if(mid_zpt > 0   && mid_zpt < 90   &&
             mid_xpt > 0   && mid_xpt < 47.5 &&
             mid_ypt > -20 && mid_ypt < 20){numTrue_inTPC++;}
          true_total_distance += true_distance; 
                                 
                                 
          previous_mid_xpt = mid_xpt;
          previous_mid_ypt = mid_ypt;
          previous_mid_zpt = mid_zpt;
          // ### pushing back vars ###
          midx.push_back(truetraj.X(iPrimPt));
          midy.push_back(truetraj.Y(iPrimPt));
          midz.push_back(truetraj.Z(iPrimPt));
          midpx.push_back(truetraj.Px(iPrimPt));
          midpy.push_back(truetraj.Py(iPrimPt));
          midpz.push_back(truetraj.Pz(iPrimPt));
	      iPrimPt++;
	      }//<--End loop on true trajectory points

        MidPosX.push_back(midx); MidPosY.push_back(midy); MidPosZ.push_back(midz);
        MidPx.push_back(midpx);  MidPy.push_back(midpy);  MidPz.push_back(midpz);
        midx.clear(); midy.clear(); midz.clear();
        midpx.clear(); midpy.clear(); midpz.clear();
	    
	    auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
	    // Ok, let's only store the processes that were saved in g4
      // We will handle the special case of 0 elsewhere
          if (thisTracjectoryProcessMap.size())
	        {
	        // The map is not zero: somthing interesting might happen in the middle of the track!!
	        for(auto const& couple: thisTracjectoryProcessMap) 
	          {
	          int interestingPoint = (int) couple.first;
	          InteractionPoint.push_back(interestingPoint);         	   
	          InteractionPointType.push_back((truetraj.KeyToProcess(couple.second)));           
	          }
	        }

          std::cout << "\nInteractions..." << std::endl;
          for (int iPt = 0; iPt < (int)InteractionPoint.size(); iPt++)
          {
            std::cout << "Interaction at " << geant_part[i]->Position(InteractionPoint[iPt]).Vect().Z() << " was " << InteractionPointType[iPt] << std::endl;
          }
   
	    iPrim++;
	    }//<--End if primary

	  if(geant_part[i]->Process() != pri)
        {
        if(geant_part[i]->Mother() == 1) 
          {
          //std::cout << "    DPdgCode: " << geant_part[i]->PdgCode() << std::endl;
          //if(geant_part[i]->PdgCode() < 10000){DPdgCode.push_back(geant_part[i]->PdgCode());}
          //DPdgCode.push_back(geant_part[i]->PdgCode());
          //DStartEnergy.push_back(geant_part[i]->E());
	      //NDTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
          if(geant_part[i]->PdgCode() == 2112)
            {NNeutronDaughters++;}
          if(geant_part[i]->PdgCode() == 2212)
            {NProtonDaughters++;}
          if(geant_part[i]->PdgCode() < 10000)
            {
	        NDTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
            DTrackId.push_back(geant_part[i]->TrackId());
            DPdgCode.push_back(geant_part[i]->PdgCode());
            DStartEnergy.push_back(geant_part[i]->E());
            DStartP.push_back(geant_part[i]->P());

            std::cout<<"daughter info: "<<std::endl;
            std::cout<<"\tpdg: "<<geant_part[i]->PdgCode()<<std::endl;
            std::cout<<"\tenergy: "<<geant_part[i]->E()<<std::endl;
            std::cout<<"\tmomentum: "<<geant_part[i]->P()<<std::endl;

	        simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
	        int jDaughtPt = 0;	
            std::vector<double> midx; std::vector<double> midy; std::vector<double> midz;
	        for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	          {
              midx.push_back(truetraj.X(jDaughtPt));
              midy.push_back(truetraj.Y(jDaughtPt));
              midz.push_back(truetraj.Z(jDaughtPt));
	          jDaughtPt++;
	          }//<--End loop on true trajectory points
            DMidPosX.push_back(midx); DMidPosY.push_back(midy); DMidPosZ.push_back(midz);
            midx.clear(); midy.clear(); midz.clear();
            }//<-- End if proton
          }//<-- End if primary daughter
        }//<-- End if not primary
      }//<--End loop on geant particles
    }//<---End checking if this is MC 

// ------------------------------------------------------   
//                      tof stuff 
// ------------------------------------------------------
    num_tof_objects = tof.size();
    size_t tof_counter = 0; // book-keeping
    for(size_t i = 0; i < tof.size(); i++){
      size_t number_tof = tof[i]->NTOF();
      for(size_t tof_idx = 0; tof_idx < number_tof; ++tof_idx){
        tofObject[tof_counter] =  tof[i]->SingleTOF(tof_idx);
        ++tof_counter;
      } // loop over TOF
    }//<---End tof_count loop

// --------------------------------------------------------
//                      wc track stuff
// --------------------------------------------------------
    num_wctracks = wctrack.size();
    for(size_t wct_count = 0; wct_count < wctrack.size(); wct_count++){
      wctrk_momentum[wct_count] = wctrack[wct_count]->Momentum();
      wctrk_XFace[wct_count] = wctrack[wct_count]->XYFace(0);
      wctrk_YFace[wct_count] = wctrack[wct_count]->XYFace(1);
      wctrk_theta[wct_count] = wctrack[wct_count]->Theta();
      wctrk_phi[wct_count] = wctrack[wct_count]->Phi();
      wctrk_missed[wct_count] = wctrack[wct_count]->WCMissed();
      wctrk_picky[wct_count] = static_cast< int > (wctrack[wct_count]->IsPicky());
      wctrk_quality[wct_count] = wcquality_handle->size();

      wctrk_x_proj_3cm[wct_count] = wctrack[wct_count]->ProjectionAtZ(3, false).X();
      wctrk_y_proj_3cm[wct_count] = wctrack[wct_count]->ProjectionAtZ(3, false).Y();
      wctrk_z_proj_3cm[wct_count] = wctrack[wct_count]->ProjectionAtZ(3, false).Z();

      wctrk_wc1_x[wct_count] = wctrack[wct_count]->HitPosition(0, 0);
      wctrk_wc1_y[wct_count] = wctrack[wct_count]->HitPosition(0, 1);
      wctrk_wc1_z[wct_count] = wctrack[wct_count]->HitPosition(0, 2);
      wctrk_wc2_x[wct_count] = wctrack[wct_count]->HitPosition(1, 0);
      wctrk_wc2_y[wct_count] = wctrack[wct_count]->HitPosition(1, 1);
      wctrk_wc2_z[wct_count] = wctrack[wct_count]->HitPosition(1, 2);
      wctrk_wc3_x[wct_count] = wctrack[wct_count]->HitPosition(2, 0);
      wctrk_wc3_y[wct_count] = wctrack[wct_count]->HitPosition(2, 1);
      wctrk_wc3_z[wct_count] = wctrack[wct_count]->HitPosition(2, 2);
      wctrk_wc4_x[wct_count] = wctrack[wct_count]->HitPosition(3, 0);
      wctrk_wc4_y[wct_count] = wctrack[wct_count]->HitPosition(3, 1);
      wctrk_wc4_z[wct_count] = wctrack[wct_count]->HitPosition(3, 2);
    }


  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE 3-D TRACK INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
   
  // ### Calling the track momentum calculator ###
  trkf::TrackMomentumCalculator trkm;
  trkm.SetMinLength(10); //change the minimal track length requirement to 10 cm

  // === Saving the number of tracks per event ===
  ntracks_reco=tracklist.size();
  //TVector3 larStart;
  //TVector3 larEnd;
  
  // === Association between WC Tracks and TPC Tracks ===
  int TempTrackMatchedID = -1;
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {
    std::cout<<"\t\t\tGOT LABEL\n";
    art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
    if(fWC2TPC.isValid())
      {
      std::cout<<"\t\t\t\t\tmodule is valid!\n";
      // === Loop on all the Assn WC-TPC tracks === 
      for(unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn )
        {
        // =========================                                                                                       
        // === Get the TPC track ===
        // =========================                                                                      
        cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));

        if(!trackWC2TPC) continue;
        recob::Track const& aTrack(trackWC2TPC.ref());
        TempTrackMatchedID = aTrack.ID();
        std::cout<<"\t\t\t\t\t\t\tTEMPTRACKID: "<<TempTrackMatchedID<<std::endl;

        }//<----End indexAssn loop          
      }//<---End checking that the WC2TPC   
    }



  // ### Looping over tracks ###
  //double maxtrackenergy = -1;



  //bool is_primary = false;
  //std::cout << "Reco info! trying to track down calo info..." << std::endl;
  for(int iTrack = 0; iTrack < ntracks_reco; iTrack++)
    {
    bool is_primary = false;

    // ### Storing an integer for the match of a WC to TPC track ###
    if(TempTrackMatchedID == tracklist[iTrack]->ID() )
      {track_WC2TPC_match.push_back(1);std::cout<<"\n\n\t\tMATCHED\n";}//<---End match
    else{track_WC2TPC_match.push_back(0);std::cout<<"\t\tNOTMATCHED\n";}


    
    // ### Setting the track information into memory ###
    auto trackStartEnd = tracklist[iTrack]->Extent();
    //larStart = tracklist[iTrack]->VertexDirection();
    //larEnd = tracklist[iTrack]->EndDirection();
    
    // ### track start and end ###
    track_start_x.push_back(trackStartEnd.first.X()); 
    track_start_y.push_back(trackStartEnd.first.Y()); 
    track_start_z.push_back(trackStartEnd.first.Z()); 
    track_end_x.push_back(trackStartEnd.second.X()); 
    track_end_y.push_back(trackStartEnd.second.Y()); 
    track_end_z.push_back(trackStartEnd.second.Z()); 

    // ### filling the kinks for this track
    std::vector<double> kink_x, kink_y, kink_z;
    for (int iK=0; iK < (int)kinklist.size(); iK++)
    {
      // id should be the track id
      auto vtx = *kinklist[iK];
      if(vtx.ID() != tracklist[iTrack]->ID())continue;
      auto pos = vtx.position();
      std::cout<<"\n\nKINK AT "<<pos.X()<<" "<<pos.Y()<<" "<<pos.Z()<<std::endl;
      kink_x.push_back(pos.X());
      kink_y.push_back(pos.Y());
      kink_z.push_back(pos.Z());
    }
    track_kink_x.push_back(kink_x);
    track_kink_y.push_back(kink_y);
    track_kink_z.push_back(kink_z);
    
    // ### Recording the track length as calculated by the tracking module ###
    track_length.push_back(tracklist[iTrack]->Length());
    std::cout << "\ttrack length: " << tracklist[iTrack]->Length() << std::endl;
    
    // ### Grabbing the SpacePoints associated with this track ###
    std::vector< art::Ptr<recob::SpacePoint > > spts = fmsp.at(iTrack);
    ntrack_hits.push_back(fmsp.at(iTrack).size());
    std::cout << "\tntrack hits: " << spts.size() << std::endl;
    if(spts[0]->XYZ()[2] < 1){is_primary = true;}
    if(is_primary == true)
      {
      std::cout<<"\tprimary candidtate"<<std::endl;
      std::cout<<"\t\tfirst z: "<<spts[0]->XYZ()[2]<<std::endl;
      }
	if(is_primary==true){track_primary.push_back(1);}
	else                {track_primary.push_back(10);}

    // ### Looping over all the SpacePoints ###
    std::vector<double> x_spts;
    std::vector<double> y_spts;
    std::vector<double> z_spts;
    for(size_t jSpt = 0; jSpt < spts.size(); jSpt++)
      {
      x_spts.push_back(spts[jSpt]->XYZ()[0]);
      y_spts.push_back(spts[jSpt]->XYZ()[1]);
      z_spts.push_back(spts[jSpt]->XYZ()[2]);
      //if(is_primary == true)
      //  {
      //  std::cout << "\t\treco pt (x,y,z): " << spts[jSpt]->XYZ()[0] << ", "
      //                                       << spts[jSpt]->XYZ()[1] << ", "
      //                                       << spts[jSpt]->XYZ()[2] << ") " << std::endl;
      //  }
      }//<----End SpacePoint loop (j)
    track_xpos.push_back(x_spts);
    track_ypos.push_back(y_spts);
    track_zpos.push_back(z_spts);
    x_spts.clear(); y_spts.clear(); z_spts.clear();


    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							CALORIMERTY FROM THIS TRACK INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------      
    
    // ########################################################## 
    // ### Looping over Calorimetry information for the track ###
    // ########################################################## 
    if(fmcal.isValid())
      {
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(iTrack);
      for (size_t j = 0; j<calos.size(); ++j)
        {
        if(!calos[j]->PlaneID().isValid){continue;}
        // ### Grabbing this calorimetry points plane number (0 == induction, 1 == collection) ###
        int pl = calos[j]->PlaneID().Plane;
        if(pl<0||pl>1){continue;}
    
        // ### Recording the number of calorimetry points for this track in this plane ####
        if(pl==0)
          {
          //std::cout << "        number of induction calorimetry points: " << calos[j]->dEdx().size() << std::endl;
          ind_track_hits.push_back(calos[j]->dEdx().size());
          ind_track_ke.push_back(calos[j]->KineticEnergy());
          //std::vector<double> ind_wire;
          std::vector<double> ind_dedx; std::vector<double> ind_dqdx;
          std::vector<double> ind_rr;   std::vector<double> ind_pitch_hit;
          for(size_t k = 0; k<calos[j]->dEdx().size(); ++k)
            {
            if(k>=1000){continue;}
            //ind_wire.push_back(calos[j]->dEdx()[k]);
            ind_dedx.push_back(calos[j]->dEdx()[k]);
            ind_dqdx.push_back(calos[j]->dQdx()[k]);
            ind_rr.push_back(calos[j]->ResidualRange()[k]);
            ind_pitch_hit.push_back(calos[j]->TrkPitchVec()[k]);
            }//<---End calo points (k)
          //ind_track_wire.push_back(ind_wire);
          ind_track_dedx.push_back(ind_dedx); ind_track_dqdx.push_back(ind_dqdx); 
          ind_track_rr.push_back(ind_dqdx);   ind_track_pitch_hit.push_back(ind_pitch_hit);
          ind_dedx.clear(); ind_dqdx.clear(); ind_rr.clear(); ind_pitch_hit.clear(); 
          }
        if(pl==1)
          {
          std::cout << "\t\t\tnumber of collection calorimetry points: " << calos[j]->dEdx().size() << std::endl;
          col_track_hits.push_back(calos[j]->dEdx().size());
          col_track_ke.push_back(calos[j]->KineticEnergy());
          //std::vector<double> col_wire;
          //std::vector<double> col_x;std::vector<double> col_y;std::vector<double> col_z;
          std::vector<double> col_dedx; std::vector<double> col_dqdx;
          std::vector<double> col_rr;   std::vector<double> col_pitch_hit;
          std::vector<double> col_x; std::vector<double> col_y; std::vector<double> col_z;
          for(size_t k = 0; k<calos[j]->dEdx().size(); ++k)
            {
            if(k>=1000){continue;}
            //std::cout << "\t\t\t\tpoint xyz: "<<calos[j]->XYZ()[k][0]<<", "<<calos[j]->XYZ()[k][1]<<", "<<calos[j]->XYZ()[k][0]<<std::endl;
            //std::cout << "\t\t\t\tpoint xyz: " << calos[j]->xyz().X() << ", " << calos[j]->xyz().Y() << ", " << calos[j]->xyz().Z() << std::endl;
            //Double_t col_x = calos[j]->XYZ()[0];
            //auto col_x = calos[j]->XYZ()(0);
            //auto col_x = calos[j]->XYZ().x();
            //auto col_x = calos[j]->XYZ()[k][0];
            //std::cout << "\t\t\t\tcolx: " << col_x << std::endl;
            //std::cout << "              calo pt wire: " << calos[j]->wire()[k] << std::endl;
            //col_wire.push_back(calos[j]->wire()[k]);
            //if(is_primary == true)
            //  {
            //  //std::cout << "\t\t\t\tpoint xyz: "<<calos[j]->XYZ()[k][0]<<", "<<calos[j]->XYZ()[k][1]<<", "<<calos[j]->XYZ()[k][2]<<std::endl;
            //  //double min_dist_to_true = 1.;
            //  //int closest_g4pt = 999;
            //  //for(int ng4 = 0; ng4 < geant_list_size; ng4++)
            //  //  {
            //  //  if(process_primary[ng4] == 1)
            //  //    {
            //  //    for(int trpt = 0; trpt < NTrTrajPts[ng4]; trpt++)
            //  //      {
            //  //      double dist_to_true = sqrt(pow(calos[j]->XYZ()[k][0] - MidPosX[ng4][trpt],2) 
            //  //                               + pow(calos[j]->XYZ()[k][1] - MidPosY[ng4][trpt],2) 
            //  //                               + pow(calos[j]->XYZ()[k][2] - MidPosZ[ng4][trpt],2)); 
            //  //      std::cout << "dist to true: " << dist_to_true << std::endl;
            //  //      if(dist_to_true < min_dist_to_true)
            //  //        {
            //  //        std::cout << "\t\t\t\t\tdistance: " << dist_to_true << std::endl;
            //  //        min_dist_to_true = dist_to_true;
            //  //        closest_g4pt = trpt;
            //  //        }//<--- End if close
            //  //      }//<--- End spt loop
            //  //    }//<--- End if primary
            //  //  }//<--- End g4 loop
            //  //std::cout << "\t\t\t\t\t\tclosest g4 pt: (pt, d): (" << closest_g4pt << ", " << min_dist_to_true << ")" << std::endl;
            //  }//<---End if primary
            col_x.push_back(calos[j]->XYZ()[k][0]);
            col_y.push_back(calos[j]->XYZ()[k][1]);
            col_z.push_back(calos[j]->XYZ()[k][2]);
            col_dedx.push_back(calos[j]->dEdx()[k]);
            col_dqdx.push_back(calos[j]->dQdx()[k]);
            col_rr.push_back(calos[j]->ResidualRange()[k]);
            col_pitch_hit.push_back(calos[j]->TrkPitchVec()[k]);
            }//<---End calo points (k)
          //col_track_wire.push_back(col_wire);
          col_track_dedx.push_back(col_dedx); col_track_dqdx.push_back(col_dqdx); 
          col_track_rr.push_back(col_dqdx);   col_track_pitch_hit.push_back(col_pitch_hit);
          col_track_x.push_back(col_x); col_track_y.push_back(col_y); col_track_z.push_back(col_z);
          col_dedx.clear(); col_dqdx.clear(); col_rr.clear(); col_pitch_hit.clear(); 
          col_x.clear(); col_y.clear(); col_z.clear();

          }

        }//<---End looping over calo points (j)
      }//<---End checking Calo info is valid  

    }//<---End track loop (i)


    // ### trying to get all the raw hits ? ###
    nhits = hitlist.size();
    //std::cout<<"\n\nGoing to try and get Hit information for an event panel.\n";
    //std::cout<<"hitlist.size(): "<<hitlist.size()<<std::endl;;
    for(int iHit = 0; iHit < nhits; iHit++) {
      hit_time.push_back(hitlist[iHit]->PeakTime());
      hit_wire.push_back(hitlist[iHit]->Channel());
      hit_view.push_back(hitlist[iHit]->View());
      hit_amp.push_back(hitlist[iHit]->PeakAmplitude());
      hit_charge.push_back(hitlist[iHit]->Integral());
    }
    


  fTree->Fill();

}//<---End analyze()


void lariat::AnaTreeT1034UTA::beginJob()
{
   
  //std::cout<<"Check-1"<<std::endl;
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(),  fTreeName.c_str());
  fTree->Branch("run",                         &run,"run/I");
  fTree->Branch("subrun",                      &subrun,"subrun/I");
  fTree->Branch("event",                       &event,"event/I");
  fTree->Branch("evttime",                     &evttime,"evttime/D");
  fTree->Branch("efield",                      efield,"efield[3]/D");
  fTree->Branch("t0",                          &t0,"t0/I");
  fTree->Branch("ntracks_reco",                &ntracks_reco,"ntracks_reco/I");
  fTree->Branch("track_WC2TPC_match",          &track_WC2TPC_match);
  fTree->Branch("track_primary",               &track_primary);
  fTree->Branch("track_start_x",               &track_start_x);
  fTree->Branch("track_start_y",               &track_start_y);
  fTree->Branch("track_start_z",               &track_start_z);
  fTree->Branch("track_end_x",                 &track_end_x);
  fTree->Branch("track_end_y",                 &track_end_y);
  fTree->Branch("track_end_z",                 &track_end_z);
  fTree->Branch("track_length",                &track_length);
  fTree->Branch("ntrack_hits",                 &ntrack_hits);
  fTree->Branch("track_xpos",                  &track_xpos);
  fTree->Branch("track_ypos",                  &track_ypos);
  fTree->Branch("track_zpos",                  &track_zpos);
  fTree->Branch("nhit_ids",                    &nhit_ids);
  
  fTree->Branch("ind_track_hits",              &ind_track_hits);
  fTree->Branch("ind_track_ke",                &ind_track_ke);
  fTree->Branch("ind_track_wire",              &ind_track_wire);
  fTree->Branch("ind_track_dedx",              &ind_track_dedx);
  fTree->Branch("ind_track_dqdx",              &ind_track_dqdx);
  fTree->Branch("ind_track_rr",                &ind_track_rr);
  fTree->Branch("ind_track_pitch_hit",         &ind_track_pitch_hit);
  fTree->Branch("col_track_hits",              &col_track_hits);
  fTree->Branch("col_track_ke",                &col_track_ke);
  fTree->Branch("col_track_x",                 &col_track_x);
  fTree->Branch("col_track_y",                 &col_track_y);
  fTree->Branch("col_track_z",                 &col_track_z);
  fTree->Branch("col_track_wire",              &col_track_wire);
  fTree->Branch("col_track_dedx",              &col_track_dedx);
  fTree->Branch("col_track_dqdx",              &col_track_dqdx);
  fTree->Branch("col_track_rr",                &col_track_rr);
  fTree->Branch("col_track_pitch_hit",         &col_track_pitch_hit);

  fTree->Branch("nhits",                       &nhits,"nhits/I");
  fTree->Branch("hit_time",                    &hit_time);
  fTree->Branch("hit_amp",                     &hit_amp);
  fTree->Branch("hit_wire",                    &hit_wire);
  fTree->Branch("hit_view",                    &hit_view);
  fTree->Branch("hit_charge",                  &hit_charge);
  
  fTree->Branch("no_primaries",                &no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",             &geant_list_size,"geant_list_size/I");
  fTree->Branch("primary_p",                   &primary_p,"primary/D");
  fTree->Branch("PDG",                         &PDG);
  fTree->Branch("StartKE",                     &StartKE);
  fTree->Branch("LastKE",                      &LastKE);
  fTree->Branch("StartEnergy",                 &StartEnergy);
  fTree->Branch("StartPx",                     &StartPx);
  fTree->Branch("StartPy",                     &StartPy);
  fTree->Branch("StartPz",                     &StartPz);
  fTree->Branch("EndEnergy",                   &EndEnergy);
  fTree->Branch("EndPx",                       &EndPx);
  fTree->Branch("EndPy",                       &EndPy);
  fTree->Branch("EndPz",                       &EndPz);
  fTree->Branch("StartPointx",                 &StartPointx);
  fTree->Branch("StartPointy",                 &StartPointy);
  fTree->Branch("StartPointz",                 &StartPointz);
  fTree->Branch("EndPointx",                   &EndPointx);
  fTree->Branch("EndPointy",                   &EndPointy);
  fTree->Branch("EndPointz",                   &EndPointz);
  fTree->Branch("Process",                     &Process);
  fTree->Branch("NumberDaughters",             &NumberDaughters);
  //fTree->Branch("primary_simChannel_num_voxel",&primary_simChannel_num_voxel);
  //fTree->Branch("primary_simChannel_voxel_dr", &primary_simChannel_voxel_dr);
  //fTree->Branch("primary_simChannel_voxel_E",  &primary_simChannel_voxel_E);
  //fTree->Branch("primary_simChannel_voxel_e",  &primary_simChannel_voxel_e);
  //fTree->Branch("primary_simChannel_voxel_x",  &primary_simChannel_voxel_x);
  //fTree->Branch("primary_simChannel_voxel_y",  &primary_simChannel_voxel_y);
  //fTree->Branch("primary_simChannel_voxel_z",  &primary_simChannel_voxel_z);
  //fTree->Branch("primary_num_simChannel",      &primary_num_simChannel,"primary_num_simChannel/I");
  //fTree->Branch("primary_simChannel",          &primary_simChannel);
  //fTree->Branch("primary_simChannel_dr",       &primary_simChannel_dr);
  //fTree->Branch("primary_simChannel_E",        &primary_simChannel_E);
  //fTree->Branch("primary_simChannel_e",        &primary_simChannel_e);
  fTree->Branch("Mother",                      &Mother);
  fTree->Branch("TrackId",                     &TrackId);
  fTree->Branch("process_primary",             &process_primary);
  fTree->Branch("G4Process",                   &G4Process);
  fTree->Branch("G4FinalProcess",              &G4FinalProcess);  
  fTree->Branch("NTrTrajPts",                  &NTrTrajPts);
  fTree->Branch("NProtonDaughters",            &NProtonDaughters,"NProtonDaughters/I");
  fTree->Branch("NNeutronDaughters",           &NNeutronDaughters,"NNeutronDaughters/I");
  fTree->Branch("NDTrTrajPts",                 &NDTrTrajPts);
  fTree->Branch("DTrackId",                    &DTrackId);
  fTree->Branch("DPdgCode",                    &DPdgCode);
  fTree->Branch("DStartEnergy",                &DStartEnergy);
  fTree->Branch("DStartP",                     &DStartP);
  fTree->Branch("MidPosX",                     &MidPosX);
  fTree->Branch("MidPosY",                     &MidPosY);
  fTree->Branch("MidPosZ",                     &MidPosZ);
  fTree->Branch("MidPx",                       &MidPx);
  fTree->Branch("MidPy",                       &MidPy);
  fTree->Branch("MidPz",                       &MidPz);
  fTree->Branch("DMidPosX",                   &DMidPosX);
  fTree->Branch("DMidPosY",                   &DMidPosY);
  fTree->Branch("DMidPosZ",                   &DMidPosZ);
  fTree->Branch("InteractionPoint",            &InteractionPoint);
  fTree->Branch("InteractionPointType",        &InteractionPointType);

  fTree->Branch("num_tof_objects",             &num_tof_objects,"num_tof_objects/I");
  fTree->Branch("tofObject",                   tofObject,"tojObject[num_tof_objects]/D");
  fTree->Branch("num_wctracks",                &num_wctracks,"num_wctracks/I");
  fTree->Branch("wctrk_momentum",              wctrk_momentum,"wctrk_momentum[num_wctracks]/D");
  fTree->Branch("wctrk_XFace",                 wctrk_XFace,"wctrk_XFace[num_wctracks]/D");
  fTree->Branch("wctrk_YFace",                 wctrk_YFace,"wctrk_YFace[num_wctracks]/D");
  fTree->Branch("wctrk_theta",                 wctrk_theta,"wctrk_theta[num_wctracks]/D");
  fTree->Branch("wctrk_phi",                   wctrk_phi,"wctrk_phi[num_wctracks]/D");
  fTree->Branch("wctrk_missed",                wctrk_missed, "wctrk_missed[num_wctracks]/I");
  fTree->Branch("wctrk_picky",                 wctrk_picky, "wctrk_picky[num_wctracks]/I");
  fTree->Branch("wctrk_quality",               wctrk_quality, "wctrk_quality[num_wctracks]/I");
  fTree->Branch("wctrk_x_proj_3cm",            wctrk_x_proj_3cm, "wctrk_x_proj_3cm[num_wctracks]/D");
  fTree->Branch("wctrk_y_proj_3cm",            wctrk_y_proj_3cm, "wctrk_y_proj_3cm[num_wctracks]/D");
  fTree->Branch("wctrk_z_proj_3cm",            wctrk_z_proj_3cm, "wctrk_z_proj_3cm[num_wctracks]/D");

  fTree->Branch("wctrk_wc1_x", wctrk_wc1_x, "wctrk_wc1_x[num_wctracks]/D");
  fTree->Branch("wctrk_wc1_y", wctrk_wc1_y, "wctrk_wc1_y[num_wctracks]/D");
  fTree->Branch("wctrk_wc1_z", wctrk_wc1_z, "wctrk_wc1_z[num_wctracks]/D");
  fTree->Branch("wctrk_wc2_x", wctrk_wc2_x, "wctrk_wc2_x[num_wctracks]/D");
  fTree->Branch("wctrk_wc2_y", wctrk_wc2_y, "wctrk_wc2_y[num_wctracks]/D");
  fTree->Branch("wctrk_wc2_z", wctrk_wc2_z, "wctrk_wc2_z[num_wctracks]/D");
  fTree->Branch("wctrk_wc3_x", wctrk_wc3_x, "wctrk_wc3_x[num_wctracks]/D");
  fTree->Branch("wctrk_wc3_y", wctrk_wc3_y, "wctrk_wc3_y[num_wctracks]/D");
  fTree->Branch("wctrk_wc3_z", wctrk_wc3_z, "wctrk_wc3_z[num_wctracks]/D");
  fTree->Branch("wctrk_wc4_x", wctrk_wc4_x, "wctrk_wc4_x[num_wctracks]/D");
  fTree->Branch("wctrk_wc4_y", wctrk_wc4_y, "wctrk_wc4_y[num_wctracks]/D");
  fTree->Branch("wctrk_wc4_z", wctrk_wc4_z, "wctrk_wc4_z[num_wctracks]/D");

  fTree->Branch("TrackIdes_x", &TrackIdes_x);
  fTree->Branch("TrackIdes_y", &TrackIdes_y);
  fTree->Branch("TrackIdes_z", &TrackIdes_z);
  fTree->Branch("TrackIdes_e", &TrackIdes_e);

  fTree->Branch("vertex_track_ids", &vertex_track_ids);
  fTree->Branch("vertex_x", &vertex_x);
  fTree->Branch("vertex_y", &vertex_y);
  fTree->Branch("vertex_z", &vertex_z);
  fTree->Branch("track_kink_x", &track_kink_x);
  fTree->Branch("track_kink_y", &track_kink_y);
  fTree->Branch("track_kink_z", &track_kink_z);

  fTree->Branch("electron_lifetime",           &electron_lifetime, "electron_lifetime/D");

  // ### subdir for truth ###
  art::TFileDirectory truthDir = tfs->mkdir("truth");
  // ### histos! ###
  //E_vs_e   = truthDir.make<TH2D>("E_vs_e","IDE E(MeV) vs e(#)", 100, 0., 1., 600, 0., 6000.);

}

void lariat::AnaTreeT1034UTA::ResetVars()
{
  G4Process.clear();
  G4FinalProcess.clear();
  InteractionPoint.clear();
  InteractionPointType.clear();
  PDG.clear();
  Mother.clear();
  TrackId.clear();
  NumberDaughters.clear();

  TrackIdes_x.clear();
  TrackIdes_y.clear();
  TrackIdes_z.clear();
  TrackIdes_e.clear();

  track_kink_x.clear();
  track_kink_y.clear();
  track_kink_z.clear();

  vertex_track_ids.clear();
  vertex_x.clear();
  vertex_y.clear();
  vertex_z.clear();

  process_primary.clear();
  StartPointx.clear();  
  StartPointy.clear();
  StartPointz.clear();
  StartKE.clear();
  LastKE.clear();
  StartEnergy.clear();
  StartPx.clear();
  StartPy.clear();
  StartPz.clear(); 
  EndPointx.clear();
  EndPointy.clear(); 
  EndPointz.clear();  
  EndEnergy.clear();  
  EndPx.clear();  
  EndPy.clear();  
  EndPz.clear();  
  NDTrTrajPts.clear();
  DTrackId.clear();
  DPdgCode.clear();
  NTrTrajPts.clear();
  MidPosX.clear();
  MidPosY.clear();
  MidPosZ.clear();
  MidPx.clear();
  MidPy.clear();
  MidPz.clear();
  DMidPosX.clear();
  DMidPosY.clear();
  DMidPosZ.clear();
  DStartEnergy.clear();
  DStartP.clear();

  track_primary.clear();
  track_WC2TPC_match.clear();
  track_start_x.clear();
  track_start_y.clear();
  track_start_z.clear();
  track_end_x.clear();
  track_end_y.clear();
  track_end_z.clear();
  track_length.clear();
  ntrack_hits.clear();
  track_xpos.clear();
  track_ypos.clear();
  track_zpos.clear();
  nhit_ids.clear();
  ind_track_ke.clear();
  ind_track_wire.clear();
  ind_track_dedx.clear();
  ind_track_dqdx.clear();
  ind_track_rr.clear();
  ind_track_pitch_hit.clear();
  col_track_hits.clear();
  col_track_ke.clear();
  col_track_x.clear();
  col_track_y.clear();
  col_track_z.clear();
  col_track_wire.clear();
  col_track_dedx.clear();
  col_track_dqdx.clear();
  col_track_rr.clear();
  col_track_pitch_hit.clear();

  nhits = -99999;
  hit_time.clear();
  hit_wire.clear();
  hit_view.clear();
  hit_amp.clear();
  hit_charge.clear();

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  NProtonDaughters = 0;
  NNeutronDaughters = 0;
  ntracks_reco = -99999;

  no_primaries = -99999;
  geant_list_size=-999;
  primary_p=-999;

  num_tof_objects = -99999;
  for(int i = 0; i < kMaxTOF; i++){
    tofObject[i] = -99999;
  }//<---End i loop

  num_wctracks = -9999999;
  for (int i = 0; i < kMaxWCTracks; i++){
    wctrk_momentum[i] = -999999;
    wctrk_XFace[i] = -999999;
    wctrk_YFace[i] = -999999;
    wctrk_theta[i] = -999999;
    wctrk_phi[i] = -999999;
    wctrk_missed[i] = -999999;
    wctrk_picky[i] = -999999;
    wctrk_quality[i] = -999999;
    wctrk_x_proj_3cm[i] = -999999;
    wctrk_y_proj_3cm[i] = -999999;
    wctrk_z_proj_3cm[i] = -999999;
  }

  electron_lifetime = -99999;

}

DEFINE_ART_MODULE(lariat::AnaTreeT1034UTA)
