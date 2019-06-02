////////////////////////////////////////////////////////////////////////
// Class:       AngularRes
// Module Type: analyzer
// File:        AngularRes_module.cc
//
// Generated at Sun Jun  2 13:17:01 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
  void HandleElastic(const sim::ParticleList& plist);
  void HandleInelastic(const sim::ParticleList& plist);
  bool InActiveRegion( const TVector3& thePos  )

private:

  // tpc boundaries
  const float TPC_X_BOUND[2] = {   0.0, 47.0 };
  const float TPC_Y_BOUND[2] = { -20.0, 20.0 };
  const float TPC_Z_BOUND[2] = {   0.0, 90.0 };
  
  // fiducial volume definition
  const float FV_X_BOUND[2] = {   2.0, 45.0 };
  const float FV_Y_BOUND[2] = { -18.0, 18.0 };
  const float FV_Z_BOUND[2] = {   0.0, 88.0 };

};



//----------------------------------------------------------------------------------------------------
// Constructor
AngularRes::AngularRes(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure();  
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
  std::vector<std::string> primProcesses;
  std::vector<TVector3>    primInterestingPoint;
  for (size_t iG4 = 0; iG4 < plist.size(); iG4++)
  {
    // ### Get MCParticle
    auto mcParticle = plist.Particle(iG4);
    if (mcParticle->Process().find("primary") == std::string::npos) continue;
    // ### Make sure this primary entered TPC
    if (mcParticle->EndPosition().Vect().Z() < FV_BOUND_Z[0]) continue; 

    // ### Store the processes for this primary, only those that occur in TPC
    simb::MCTrajectory truetraj = mcParticle->Trajectory();
    auto thisTrjProcessMap = truetraj.TrajectoryProcesses();
  
    // ### If there's nothing, check the end of the track
    if (thisTrjProcessMap.size() == 0)
    {
      if (!InActiveRegion(mcParticle->EndPosition().Vect())) continue; 
      primInterestingPoints.push_back( mcParticle->NumberTrajectoryPoints()-1 );
 
      // ### If it doesn't have daughters, it's through going
      if (mcParticle->NumberDaughters())
      {
        auto theDauId = mcParticle->Daughter(0);
        for (size_t iD = 0; iD < plist.size(); iD++) 
        {
          if (plist.Particle(iD)->TrackId() == theDauId) primProcesses.push_back(plist.Particle(iD)->Process());
        }
      }//<-- End if has daughters
      else primProcesses.push_back("throughgoing");
    }//<-- End if map is zero
    else 
    {
      for (const auto& couple : thisTrjProcessMap)
      {
        // ### Get the vertex
        int interestingPoint = (int) couple.first;
        if (!InActiveRegion(truetraj.Position(interestingPoint).Vect())) continue; 

        primInterestingPoints.push_back( interestingPoint );
        primProcesses.push_back( truetraj.KeyToProcess(couple.second) );
      }
    }//<-- End if map is not zero
  }

  // check for elastic and inelastic
  bool isElastic   = false;
  bool isInelastic = false;
  for (const auto& p : primProcesses) 
  {
    if (p.find("hadElastic") != std::string::npos || p.find("CoulombScat") != std::string::npos) isElastic = true;
    if (p.find("pi") != std::string::npos || p.find("Inelastic") != std::string::npos) isInelastic = true;
  }

  // ### ignore cases in which we have both
  if (isElastic && is Inelastic) return;

  // ### Hopefully we have interesting points, otherwise stop
  if (!primInterestingPoints.size()) throw cet::exception("AngularRes") << "Found processes but no interesting points!\n";

  // ### Handle elastic
  if (isElastic) HandleElastic(plist);

  // ### Handle inelastic
  if (isInelastic) HandleInelastic(plist);
}




//----------------------------------------------------------------------------------------------------
// Begin job
void AngularRes::beginJob()
{
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
void AngularRes::HandleElastic(const sim::ParticleList& plist)
{
  
}


//----------------------------------------------------------------------------------------------------
// Handle inelastic
void AngularRes::HandleInelastic(const sim::ParticleList& plist)
{

}

DEFINE_ART_MODULE(AngularRes)

}
