////////////////////////////////////////////////////////////////////////
// Class:       AngularRes
// Module Type: analyzer
// File:        AngularRes_module.cc
//
// Generated at Sun Jun  2 13:17:01 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TVector3.h"
#include "TH1D.h"
#include "TTree.h"

namespace piinelastic 
{

class AngularRes;

class AngularRes : public art::EDAnalyzer {
public:
  explicit AngularRes(fhicl::ParameterSet const & p);
  AngularRes(AngularRes const &) = delete;
  AngularRes(AngularRes &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  inline void PrintVec(const TVector3& v) { std::cout << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << ")" << std::endl; }
  void HandleElastic(const sim::ParticleList& plist, const std::vector<size_t>& points);
  void HandleInelastic(const sim::ParticleList& plist, const std::vector<size_t>& points);
  bool InActiveRegion( const TVector3& thePos  );


  // tpc boundaries
  const float TPC_X_BOUND[2] = {   0.0, 47.0 };
  const float TPC_Y_BOUND[2] = { -20.0, 20.0 };
  const float TPC_Z_BOUND[2] = {   0.0, 90.0 };
  
  // fiducial volume definition
  const float FV_X_BOUND[2] = {   2.0, 45.0 };
  const float FV_Y_BOUND[2] = { -18.0, 18.0 };
  const float FV_Z_BOUND[2] = {   0.0, 88.0 };

  // histograms
  TH1D* hElasticAngle;

};



//----------------------------------------------------------------------------------------------------
// Constructor
AngularRes::AngularRes(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);  
}



//----------------------------------------------------------------------------------------------------
// Reconfigure
void AngularRes::reconfigure(fhicl::ParameterSet const & p)
{
}




//----------------------------------------------------------------------------------------------------
// Analyzer
void AngularRes::analyze(art::Event const & e)
{
  // ### Backtracker to recover truth information
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // ### get primary info
  std::map<size_t, std::string> primInterestingPoints;
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // ### Get MCParticle
    auto mcParticle = plist.Particle(iG4);
    if (mcParticle->Process().find("primary") == std::string::npos) continue;
    // ### Make sure this primary entered TPC
    if (mcParticle->EndPosition().Vect().Z() < FV_Z_BOUND[0]) continue; 

    // ### Store the processes for this primary, only those that occur in TPC
    simb::MCTrajectory truetraj = mcParticle->Trajectory();
    auto thisTrjProcessMap = truetraj.TrajectoryProcesses();
  
    // ### If there's nothing, check the end of the track
    if (thisTrjProcessMap.size() == 0)
    {
      if (!InActiveRegion(mcParticle->EndPosition().Vect())) continue; 
      size_t thePoint = mcParticle->NumberTrajectoryPoints()-1;
 
      // ### If it doesn't have daughters, it's through going
      if (mcParticle->NumberDaughters())
      {
        auto theDauId = mcParticle->Daughter(0);
        for (size_t iD = 0; iD < plist.size(); iD++) 
        {
          if (plist.Particle(iD)->TrackId() == theDauId) primInterestingPoints.emplace(thePoint, plist.Particle(iD)->Process());
        }
      }//<-- End if has daughters
    }//<-- End if map is zero
    else 
    {
      for (const auto& couple : thisTrjProcessMap)
      {
        // ### Get the vertex
        int thePoint = (int) couple.first;
        if (!InActiveRegion(truetraj.Position(thePoint).Vect())) continue; 

        primInterestingPoints.emplace( thePoint, truetraj.KeyToProcess(couple.second) );
      }
    }//<-- End if map is not zero
  }

  // ### Check for interesting points, if not quit
  if (primInterestingPoints.size() == 0) return; 

  // check for elastic and inelastic
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  std::cout << "This primary underwent the following processes...\n";
  std::vector<size_t> elasticPoints, inelasticPoints;
  for (const auto& p : primInterestingPoints) 
  {
    std::cout << p.second << std::endl;
    if (p.second.find("hadElastic") != std::string::npos || p.second.find("CoulombScat") != std::string::npos) elasticPoints.push_back(p.first);
    if (p.second.find("pi") != std::string::npos && p.second.find("Inelastic") != std::string::npos)           inelasticPoints.push_back(p.first);
  }
  std::cout << std::endl;

  // ### ignore cases in which we have neither
  if (elasticPoints.size() == 0 && inelasticPoints.size() == 0) return;
  
  // ### Handle elastic
  if (elasticPoints.size() != 0) HandleElastic(plist, elasticPoints);

  // ### Handle inelastic
  if (inelasticPoints.size() != 0) HandleInelastic(plist, inelasticPoints);

  std::cout << std::endl;
}




//----------------------------------------------------------------------------------------------------
// Begin job
void AngularRes::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  
  hElasticAngle = tfs->make<TH1D>("hElasticAngle", "Scattering Angle for Elastic", 1000, 0, 180);
}



//----------------------------------------------------------------------------------------------------
// In active region
bool AngularRes::InActiveRegion( const TVector3& thePos  )
{
  if ( FV_X_BOUND[0] < thePos.X() && thePos.X() < FV_X_BOUND[1] &&
       FV_Y_BOUND[0] < thePos.Y() && thePos.Y() < FV_Y_BOUND[1] &&
       FV_Z_BOUND[0] < thePos.Z() && thePos.Z() < FV_Z_BOUND[1] ) return true;

  return false;
}




//----------------------------------------------------------------------------------------------------
// End job
void AngularRes::endJob()
{
 
}


//----------------------------------------------------------------------------------------------------
// Handle elastic
void AngularRes::HandleElastic(const sim::ParticleList& plist, const std::vector<size_t>& points)
{
  // ### Since this is elastic, we only need the primary information here
  size_t primId(0);
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // ### Get MCParticle
    auto mcParticle = plist.Particle(iG4);
    if (mcParticle->Process().find("primary") == std::string::npos) continue; 
    primId = iG4;
  }

  // ### the primary
  auto mcPrimary = plist.Particle(primId);

  std::cout << "Elastic scatters at..." << std::endl;
  for (const auto& p : points) 
  {
    auto pos = mcPrimary->Position(p).Vect();
    PrintVec(pos);
  }
  std::cout << std::endl;

  //////////////////////////////////////////////
  // ### Loop over interesting points
  for (const auto& p : points)
  {
    // get the momentum vector before and at this point
    auto momBefore = mcPrimary->Momentum(p-1).Vect();
    auto momHere   = mcPrimary->Momentum(p).Vect();

    float theta = (180/TMath::Pi())*std::acos( momBefore.Unit().Dot(momHere.Unit()) );
    hElasticAngle->Fill(theta);
  } 
}




//----------------------------------------------------------------------------------------------------
// Handle inelastic
void AngularRes::HandleInelastic(const sim::ParticleList& plist, const std::vector<size_t>& points)
{
  // ### Get the primary id
  size_t primId(0);
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // ### Get MCParticle
    auto mcParticle = plist.Particle(iG4);
    if (mcParticle->Process().find("primary") == std::string::npos) continue; 
    primId = iG4;
  }

  // ### the primary
  auto mcPrimary = plist.Particle(primId);

  std::cout << "Inelastic interactions at..." << std::endl;
  for (const auto& p : points)
  {
    auto pos = mcPrimary->Position(p).Vect();
    PrintVec(pos);
  }
  std::cout << std::endl;

  //////////////////////////////////////////////
  // ### Loop over the interesting points
  for (const auto& p : points)
  {

  }
}

DEFINE_ART_MODULE(AngularRes)

}
