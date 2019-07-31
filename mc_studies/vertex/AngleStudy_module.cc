////////////////////////////////////////////////////////////////////////
// Class:       AngleStudy
// Module Type: analyzer
// File:        ElasticLikeVertex_module.cc
//
// Generated at Mon Jul 29 09:41:14 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"


#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

namespace piinelastic 
{
/// Container for pdg codes for charged particles (cuz I can't remember) 
struct PdgCodes_t
{
  int kElectron = 11;
  int kMuon     = 13;
  int kPion     = 211;
  int kKaon     = 321;
  int kProton   = 2212;
};

/// Vertex structure
struct Vertex_t
{
  int         point;    // int point of primary vertex
  std::string process;  // the process
  TVector3    position; // 3D point

  void Print() { std::cout << "Point " << point << "\nProcess " << process << "\n Position "; position.Print(); }
};

class AngleStudy;

class AngleStudy : public art::EDAnalyzer {
public:
  explicit AngleStudy(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AngleStudy(AngleStudy const &) = delete;
  AngleStudy(AngleStudy &&) = delete;
  AngleStudy & operator = (AngleStudy const &) = delete;
  AngleStudy & operator = (AngleStudy &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  bool inActiveRegion(const TVector3& position);
  void doAngleStudy(const sim::ParticleList& plist, const Vertex_t& vertex);
  void identifyVisibleSecondaries(const sim::ParticleList& plist, const Vertex_t& vertex);
  double toKineticEnergy(const TVector3& mom, const double& mass);
  void getFirstTpcPoint(const simb::MCParticle& mcPrimary, TVector3& position);

  /// Method to check if charged
  inline bool isCharged(int const& pdg) 
  { 
    return (pdg == fPdgcodes.kElectron || pdg == fPdgcodes.kMuon || 
            pdg == fPdgcodes.kPion     || pdg == fPdgcodes.kKaon || 
            pdg == fPdgcodes.kProton);
  }

  /// Fiducial volume definition
  float FV_X_BOUND[2] = {   1.0, 46.0 };
  float FV_Y_BOUND[2] = { -18.0, 18.0 };
  float FV_Z_BOUND[2] = {   0.0, 88.0 };

  int fVerbose = 1;
  std::vector<int> fVisSec;
  std::vector<Vertex_t> fVertices;
  PdgCodes_t fPdgcodes;

  const float SECONDARY_LENGTH_CUT = 2.0;

  TH1D* hMCSecondaryTrkLength      = nullptr;
  TH1S* hMCSecondaries             = nullptr;

  TH1D* hMCElasticAngle             = nullptr;
  TH1D* hMCElasticConeAngle         = nullptr;
  TH2D* hMCKeVsElasticAngle         = nullptr;
  TH2D* hMCKeVsElasticConeAngle     = nullptr;
  TH2D* hMCTrkLenVsElasticAngle     = nullptr;
  TH2D* hMCTrkLenVsElasticConeAngle = nullptr;

  TH1D* hMCInelasticAngle           = nullptr;

  TH1D* hMCInelasticAngleOneVisD             = nullptr;
  TH1D* hMCInelasticConeAngleOneVisD         = nullptr;
  TH2D* hMCKeVsInelasticAngleOneVisD         = nullptr;
  TH2D* hMCKeVsInelasticConeAngleOneVisD     = nullptr;
  TH2D* hMCTrkLenVsInelasticAngleOneVisD     = nullptr;
  TH2D* hMCTrkLenVsInelasticConeAngleOneVisD = nullptr;
};

//----------------------------------------------------------------------------------
AngleStudy::AngleStudy(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

//----------------------------------------------------------------------------------
void AngleStudy::analyze(art::Event const & e)
{
  if (fVerbose)
  {
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "AngleStudy processing event #" << e.id().event() << "\n";
  }
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  fVertices.clear();
  int primTrkId(0);
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // Get the true particle and its process, skip whatever is not primary 
    auto mcParticle = plist.Particle(iG4);
    if ( mcParticle->Process().find("primary") == std::string::npos ) continue;
    // Make sure the primary entered tpc
    if (mcParticle->EndPosition().Vect().Z() < FV_Z_BOUND[0]) continue;
    primTrkId = mcParticle->TrackId();

    // Store the processes for this primary, only those that occur in TPC
    simb::MCTrajectory truetraj = mcParticle->Trajectory();
    auto thisTrjProcessMap = truetraj.TrajectoryProcesses();

    // If there's nothing, check the end of the track
    if (thisTrjProcessMap.size() == 0)
    {
      if (!inActiveRegion(mcParticle->EndPosition().Vect())) continue;
      size_t thePoint = mcParticle->NumberTrajectoryPoints()-1;

      // ### If it doesn't have daughters, it's through going
      if (mcParticle->NumberDaughters())
      {
        auto theDauId = mcParticle->Daughter(0);
        for (size_t iD = 0; iD < plist.size(); iD++)
        {
          if (plist.Particle(iD)->TrackId() == theDauId) 
          {
            Vertex_t vertex;
            vertex.point    = thePoint;
            vertex.process  = plist.Particle(iD)->Process();
            vertex.position = mcParticle->Position(thePoint).Vect();
            fVertices.push_back(vertex);
          }
        }
      }//<-- End if has daughters
    }//<-- End if map is zero
    else
    {
      for (const auto& couple : thisTrjProcessMap)
      {
        // ### Get the vertex
        int thePoint = (int) couple.first;
        if (!inActiveRegion(truetraj.Position(thePoint).Vect())) continue;
        Vertex_t vertex;
        vertex.point    = thePoint;
        vertex.process  = truetraj.KeyToProcess(couple.second);
        vertex.position = truetraj.Position(thePoint).Vect();
        fVertices.push_back(vertex);
      }
    }//<-- End if map is not zero
  }
  if (!fVertices.size()) return;

  // Get first interaction in TPC
  if (fVerbose) std::cout << "\nGetting first interaction in TPC...\n";
  std::sort(fVertices.begin(), fVertices.end(), [](const auto& l, const auto& r) { return l.position.Z() < r.position.Z(); });
  Vertex_t firstVertexInTpc = fVertices[0];

  std::cout << "\nThis primary underwent the following processes...\n";
  for (const auto& v : fVertices) std::cout << v.process << std::endl;
  std::cout << "\nFirst vertex in TPC:\n";
  firstVertexInTpc.Print();

  // Begin study
  if (fVerbose) std::cout << "\nFilling track lengths...\n";
  // Make distribution of secondary track lengths
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)   
  {
    auto mcParticle = plist.Particle(iG4);

    // Only look at daughters of primary
    if (std::abs(mcParticle->Mother()) != primTrkId) continue;

    // Make sure she's charged
    if ( !isCharged(std::abs(mcParticle->PdgCode())) ) continue;
    TVector3 startPoint = mcParticle->Position().Vect();
    TVector3 endPoint   = mcParticle->EndPosition().Vect();
    hMCSecondaryTrkLength->Fill((startPoint-endPoint).Mag());
  }

  // Angle study
  if (fVerbose) std::cout << "\nStarting angle study...\n";
  doAngleStudy(plist, firstVertexInTpc);
}

//----------------------------------------------------------------------------------
void AngleStudy::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hMCSecondaryTrkLength    = tfs->make<TH1D>("hMCSecondaryTrkLength",    "Secondary Track Lengths",               100, 0, 50);
  hMCSecondaries           = tfs->make<TH1S>("hMCSecondaries",           "True number of tracks leaving vertex",  10, 0, 10);

  // should all have one secondary
  hMCElasticAngle              = tfs->make<TH1D>("hMCElasticAngle",             "Elastic Scattering Angle of Primary",                                   180, 0, 180);
  hMCElasticConeAngle          = tfs->make<TH1D>("hMCElasticConeAngle",         "Cone Angle Between Single Visible Secondary and Primary for Elastic",    90, 0, 90);
  hMCKeVsElasticAngle          = tfs->make<TH2D>("hMCKeVsElasticAngle"  ,       "Kinetic energy vs Elastic Angle",      90, 0, 90, 300, 0, 1200);
  hMCKeVsElasticConeAngle      = tfs->make<TH2D>("hMCKeVsElasticConeAngle",     "Kinetic energy vs Elastic Cone Angle", 90, 0, 90, 300, 0, 1200);
  hMCTrkLenVsElasticAngle      = tfs->make<TH2D>("hMCTrkLenVsElasticAngle",     "Track length vs Elastic Cone Angle",   90, 0, 90, 90, 0, 90);
  hMCTrkLenVsElasticConeAngle  = tfs->make<TH2D>("hMCTrkLenVsElasticConeAngle", "Track length vs Elastic Cone Angle",   90, 0, 90, 90, 0, 90);

  // general inelastic
  hMCInelasticAngle              = tfs->make<TH1D>("hMCInelasticAngle", "Angle Between Secondaries and Primary for Inelastic", 180, 0, 180);

  // inelastic one visible daughter
  hMCInelasticAngleOneVisD             = tfs->make<TH1D>("hMCInelasticAngleOneVisD",      "Angle Between Single Visible Secondary and Primary for Inelastic", 180, 0, 180);
  hMCInelasticConeAngleOneVisD         = tfs->make<TH1D>("hMCInelasticConeAngleOneVisD",  "Angle Between Single Visible Secondary and Primary for Inelastic", 180, 0, 180);
  hMCKeVsInelasticAngleOneVisD         = tfs->make<TH2D>("hMCKeVsInelasticAngleOneVisD",         "Kinetic energy vs Inelastic Angle Single Visible Secondary",       90, 0, 90, 300, 0, 1200);
  hMCKeVsInelasticConeAngleOneVisD     = tfs->make<TH2D>("hMCKeVsInelasticConeAngleOneVisD",     "Kinetic energy vs Inelastic Cone Angle Single Visible Secondary",  90, 0, 90, 300, 0, 1200);
  hMCTrkLenVsInelasticAngleOneVisD     = tfs->make<TH2D>("hMCTrkLenVsInelasticAngleOneVisD",     "Track length vs Inelastic Cone Angle Single Visible Secondary", 90, 0, 90, 90, 0, 90);
  hMCTrkLenVsInelasticConeAngleOneVisD = tfs->make<TH2D>("hMCTrkLenVsInelasticConeAngleOneVisD", "Track length vs Inelastic Cone Angle Single Visible Secondary", 90, 0, 90, 90, 0, 90);
}

//----------------------------------------------------------------------------------
void AngleStudy::endJob()
{
  // Implementation of optional member function here.
}

//----------------------------------------------------------------------------------
void AngleStudy::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

/**
 * @brief Check if in active region
 *
 */
bool AngleStudy::inActiveRegion( const TVector3& thePos  )
{
  if ( FV_X_BOUND[0] < thePos.X() && thePos.X() < FV_X_BOUND[1] &&
       FV_Y_BOUND[0] < thePos.Y() && thePos.Y() < FV_Y_BOUND[1] &&
       FV_Z_BOUND[0] < thePos.Z() && thePos.Z() < FV_Z_BOUND[1] ) return true;

  return false;
}

/**
 * @brief Area to make angle plots
 * 
 * @param plist List of mcparticles
 * @param vertex The first tpc vertex
 */
void AngleStudy::doAngleStudy(const sim::ParticleList& plist, const Vertex_t& vertex)
{
  if (vertex.process.find("none") != std::string::npos) return;

  // Identify visible secondaries
  identifyVisibleSecondaries(plist, vertex);
  if (fVerbose) std::cout << "Visible secondaries " << fVisSec.size() << std::endl;

  // Handle the case of just the primary
  bool isElastic(false);
  bool isInelastic(false);
  if (vertex.process.find("hadElastic")   != std::string::npos || vertex.process.find("CoulombScat") != std::string::npos) isElastic = true;
  if (vertex.process.find("pi-Inelastic") != std::string::npos) isInelastic = true;

  // Only looking at elastic or inelastic
  if (!isElastic && !isInelastic) return;

  int primId(0); 
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // Get the true particle and its process, skip whatever is not primary 
    auto mcParticle = plist.Particle(iG4);
    if ( mcParticle->Process().find("primary") == std::string::npos ) continue;
    // Make sure the primary entered tpc
    if (mcParticle->EndPosition().Vect().Z() < FV_Z_BOUND[0]) continue;
    primId = iG4;
  }
  auto primMcParticle = plist.Particle(primId);

  TVector3 primIncMom = primMcParticle->Momentum(vertex.point-1).Vect();
  double   primIncKe  = toKineticEnergy(1000*primIncMom, 1000*primMcParticle->Mass());
  TVector3 primFirstTpcPoint = vertex.position;
  getFirstTpcPoint(*primMcParticle, primFirstTpcPoint);
  double   primTrkLength = (vertex.position-primFirstTpcPoint).Mag();

  // If the primary doesn't end here
  if ( (vertex.point+1) < (int)(primMcParticle->NumberTrajectoryPoints()-1) )
  {
    auto trailMom = primMcParticle->Momentum(vertex.point).Vect().Unit();
    double theta  = (180/TMath::Pi())*std::acos( primIncMom.Unit().Dot(trailMom.Unit()) );
    
    if (isElastic)        
    {
      hMCElasticAngle->Fill(theta);
      hMCKeVsElasticAngle->Fill(theta, primIncKe);
      hMCTrkLenVsElasticAngle->Fill(theta, primTrkLength);
    }
    else if (isInelastic) 
    {
      // Only filling general inelastic
      hMCInelasticAngle->Fill(theta);
    }

    if (fVerbose) 
    {
      primIncMom.Unit().Print();
      trailMom.Unit().Print();
      std::cout << "Theta " << theta << "   KE " << primIncKe << "  track length " << primTrkLength << std::endl;
    }   
  }

  // Fill inelastic scattering angles for all relevant daughters
  if (isInelastic)
  {
    for (const auto& iG4 : fVisSec)
    {
      if (iG4 == primId) continue;

      // This is a charged daughter
      TVector3 dMom0 = plist.Particle(iG4)->Momentum().Vect();

      // Compute angle between daughter and incoming primary
      double theta = (180/TMath::Pi())*std::acos( primIncMom.Unit().Dot(dMom0.Unit()) );
      hMCInelasticAngle->Fill(theta);
    }//<-- End loop over visible particles
  }//<--End if inelastic

  // Do it again, except for cases in which inelastic but one visible daughter
  if (isInelastic)
  {
    // If fVisSec.size() == 1 this looks elastic
    hMCSecondaries->Fill(fVisSec.size());
    if (fVisSec.size() == 1)
    {
      int iG4 = fVisSec[0];
      TVector3 trailMom = primIncMom;

      if (iG4 == primId) trailMom  = plist.Particle(iG4)->Momentum(vertex.point).Vect();
      else               trailMom  = plist.Particle(iG4)->Momentum().Vect();

      double theta = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(trailMom.Unit()));
      hMCInelasticAngleOneVisD->Fill(theta);
      hMCKeVsInelasticAngleOneVisD->Fill(theta, primIncKe);
      hMCTrkLenVsInelasticAngleOneVisD->Fill(theta, primTrkLength);

      if (fVerbose) 
      {
        primIncMom.Unit().Print();
        trailMom.Unit().Print();
        std::cout << "Theta " << theta << "  KE " << primIncKe << std::endl;
      } 
    }//<-- End if one visible daughter
  }//<-- End if is inelastic

  // Look at the cone angle 
  if (fVisSec.size() == 1)
  {
    int iG4 = fVisSec[0];
    TVector3 trailMom = primIncMom;

    if (iG4 == primId) trailMom  = plist.Particle(iG4)->Momentum(vertex.point).Vect();
    else               trailMom  = plist.Particle(iG4)->Momentum().Vect(); 
    
    auto resultant = primIncMom + trailMom;
    double theta   = (180 / TMath::Pi()) * std::acos(primIncMom.Unit().Dot(resultant.Unit()));
      
    if (isElastic)   
    {
      hMCElasticConeAngle->Fill(theta);
      hMCKeVsElasticConeAngle->Fill(theta, primIncKe);
      hMCTrkLenVsElasticConeAngle->Fill(theta, primTrkLength);
    }
    else if (isInelastic) 
    {
      hMCInelasticConeAngleOneVisD->Fill(theta);
      hMCKeVsInelasticConeAngleOneVisD->Fill(theta, primIncKe);
      hMCTrkLenVsInelasticConeAngleOneVisD->Fill(theta, primTrkLength);
    }
  }//<-- End if one visible daughter

  // Sanity checks
  if (isElastic && fVisSec.size() != 1) 
  {
    mf::LogWarning("AngleStudy::doAngleStudy") << "Is elastic but visible secondaries does not make sense!\n";
  }
}

/**
 * @brief Method to identify visible secondaries attached to the interaction vertex
 * 
 * @param plist List of mcparticles
 * @param vertex The first tpc vertex
 */
void AngleStudy::identifyVisibleSecondaries(const sim::ParticleList& plist, const Vertex_t& vertex)
{
  fVisSec.clear();

  // If the primary doesn't end here
  int primId(0); 
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // Get the true particle and its process, skip whatever is not primary 
    auto mcParticle = plist.Particle(iG4);
    if ( mcParticle->Process().find("primary") == std::string::npos ) continue;
    // Make sure the primary entered tpc
    if (mcParticle->EndPosition().Vect().Z() < FV_Z_BOUND[0]) continue;
    primId = iG4;
  }
  if ((vertex.point+1) < (int)plist.Particle(primId)->NumberTrajectoryPoints()) fVisSec.push_back(primId);

  // Fill inelastic scattering angles for all relevant daughters
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    auto mcParticle = plist.Particle(iG4);
    // @todo Take absolute value here?
    if (mcParticle->Mother() != plist.Particle(primId)->TrackId()) continue;

    // Make sure she is charged 
    if ( !isCharged(std::abs(mcParticle->PdgCode())) ) continue;

    // This is a charged daughter
    TVector3 dPos0 = mcParticle->Position().Vect();
    TVector3 dPosf = mcParticle->EndPosition().Vect();

    // Check track length
    // @todo Do we need this cut?
    if ((dPosf-dPos0).Mag() < SECONDARY_LENGTH_CUT) continue;

    // This particle needs to be attached to the vertex
    if ((vertex.position-dPos0).Mag() < 0.01) fVisSec.push_back(iG4);
  }//<-- End loop over G4 particles
}

/**
 * @brief Get the primary's first point in the tpc
 * 
 */
void AngleStudy::getFirstTpcPoint(const simb::MCParticle& mcPrimary, TVector3& position)
{
  for (int iPt = 0; iPt < (int)mcPrimary.NumberTrajectoryPoints(); iPt++)
  {
    auto primPos = mcPrimary.Position(iPt).Vect();
    if (primPos.Z() < FV_Z_BOUND[0]) continue;
    if (!inActiveRegion(primPos)) continue;
    position = primPos;
    return;
  }
}

/// Compute kinetic energy
double AngleStudy::toKineticEnergy(const TVector3& mom, const double& mass)
{
  return std::sqrt( mom.Mag()*mom.Mag() + mass*mass ) - mass;
}


DEFINE_ART_MODULE(AngleStudy)
}
