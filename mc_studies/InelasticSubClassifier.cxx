////////////////////////////////////////////////////////////
// InelasticSubClassifier
//
// This algorithm further classifies pion inelastic events
// into inelastic, charge exchange, absorption, and pion
// production.
//
// Author: Hunter Sullivan hunter.sullivan@mavs.uta.edu
////////////////////////////////////////////////////////////

#include "LArIATFilterModule/InelasticSubClassifier.h"

#include "nusimdata/SimulationBase/MCParticle.h"

namespace piinelastic 
{

// ########################
// ### constructor
// ########################
InelasticSubClassifier::InelasticSubClassifier()
{}

// ########################
// ### destructor
// ########################
InelasticSubClassifier::~InelasticSubClassifier()
{}

// ########################
// ### classify
// ########################
const std::string InelasticSubClassifier::Classify(const sim::ParticleList& plist, 
                                                   const int& primaryTrkId)
{
  //
  // We define sub channels for inelastic processes:
  //      1) Absorption (no pions in final state)
  //      2) Inelastic (no neutral pions, 1 same charged pion in final state)
  //      3) Charge exchange (1 nuetral pion, no charged pions in final state)
  //      4) Pion production (what's left?)
  //

  bool isPrimaryChargedPion(false);
  bool foundInelastic(false);
  int nChargedPionDaughters(0), nNeutralPionDaughters(0),
      nProtonDaughters(0),      nNeutronDaughters(0),
      nMuonDaughters(0),        nKaonDaughters(0),
      nChargedParticles(0);
  TVector3 primaryFinalPosition;
  std::set<std::string> processes;
  std::vector<TVector3> vertices;

  // ### loop over the g4 particles
  for (size_t p = 0; p < plist.size(); p++)
  {
    // get the primary associated to the trackId
    auto particle = *(plist.Particle(p));
    if ( particle.TrackId() != primaryTrkId ) continue;

    // if charged pion
    if (std::abs(particle.PdgCode()) == 211)
    {
      // change flag
      isPrimaryChargedPion = true;

      // get last point
      primaryFinalPosition = particle.EndPosition().Vect();
    } 
  }

  // ### checkout if not a charged pion
  if (!isPrimaryChargedPion) return std::string("NoPrimaryChargedPion");

  // ### loop over the particles again checking daughters
  for (size_t p = 0; p < plist.size(); p++)
  {
    auto particle = *(plist.Particle(p));

    // we only care about particles that end up in the tpc
    if (!InTPC(particle)) continue;

    // get pdg code
    int pdgCode = std::abs(particle.PdgCode());
    if (pdgCode == 2212 || pdgCode == 211 ||
        pdgCode == 321  || pdgCode == 13) nChargedParticles++;

    // skip if particle is not a child of the primary
    if (particle.Mother() != primaryTrkId) continue;
 
    // we will classify the interaction based on pions, kaons, muons, protons, and neutrons
    // count the daughters
    switch(pdgCode)
    {
      // check for neutral pion
      case 111: { nNeutralPionDaughters++; break; }
      // check for charged pions
      case 211: { nChargedPionDaughters++; break; }
      // check for protons
      case 2212:{ nProtonDaughters++; break; }
      // check for neutrons
      case 2112:{ nNeutronDaughters++; break; }
      // check for kaons
      case 321: { nKaonDaughters++; break; }
      // check for muons
      case 13:  { nMuonDaughters++; break; }
      default: continue;
    }

    // add the process to our set
    processes.insert(particle.Process());

    // check for inelastic
    if (particle.Process().find("Inelastic") != std::string::npos) 
    {
      foundInelastic = true;
      // keep this vertex if we don't already have it
      if (vertices.size() == 0) vertices.push_back(particle.Position().Vect());
      else 
      {
        float min = std::numeric_limits<float>::max();
        for (const auto& v : vertices)
        {
          auto d = (v - particle.Position().Vect()).Mag();
          if (d < min) min = d;
        }
        if (min > 0.001) vertices.push_back(particle.Position().Vect());
      }
    }
  }

  // ### print out the processes in for this primary
  std::cout << "\n//////////////////////////////////////////////////"
            << "\nInelasticSubClassifier..."
            << "\nThe processes for this primary:\n";
  for (const auto& p : processes) std::cout << p << std::endl;
  std::cout << std::endl;

  std::cout << "\nNumber of protons:           " << nProtonDaughters
            << "\nNumber of neutrons:          " << nNeutronDaughters
            << "\nNumber of charged pions:     " << nChargedPionDaughters
            << "\nNumber of neutral pions:     " << nNeutralPionDaughters
            << "\nNumber of charged kaons:     " << nKaonDaughters
            << "\nNumber of muons:             " << nMuonDaughters
            << "\nNumber of charged particles: " << nChargedParticles 
            << "\nPrimary end vertex: (" 
            << primaryFinalPosition.X() << ", " 
            << primaryFinalPosition.Y() << ", "
            << primaryFinalPosition.Z() << ")"
            << "\nInteraction vertices:\n";
  for (const auto& v : vertices) std::cout << "\t(" << v.X() << ", " << v.Y() << ", " << v.Z() << ")" << std::endl;

  // ### checkout if no inelastic
  if (!foundInelastic) 
  {
    std::cout << "\nRESULT: no inelastic interaction!"
              << "\n//////////////////////////////////////////////////"
              << std::endl;
    return std::string("NoInelasticInteraction");
  }

  // ### okay... we can start to classify this event
  
  // ### we need to calculate the distance between our 
  // ### primary end point and each vertex
  std::vector<double> distances;
  for (const auto& v : vertices)
  {
    distances.push_back( (v-primaryFinalPosition).Mag() );
  }

  // ### handle absorption and inelastic
  if (nChargedPionDaughters == 0 && nNeutralPionDaughters == 0)
  {
    // if the distance = 0, then we have absorption
    for (const auto& d : distances)
    {
      if (d < 0.001)
      {
        std::cout << "\nRESULT: PionAbsorption!"
                  << "\n//////////////////////////////////////////////////"
                  << std::endl;
        return std::string("PionAbsorption");
      }
    }
    // otherwise we have inelastic
    std::cout << "\nRESULT: PionInelastic!"
              << "\n//////////////////////////////////////////////////"
              << std::endl;
    return std::string("PionInelastic");
  }
  
  // ### handle charge exchange
  if (nNeutralPionDaughters == 1 && nChargedPionDaughters == 0)
  {
    // this should do it!
    std::cout << "\nRESULT: ChargeExchange!"
              << "\n//////////////////////////////////////////////////"
              << std::endl;
    return std::string("ChargeExchange");
  }

  // ### special case
  if (nChargedPionDaughters == 1 && nNeutralPionDaughters == 0)
  {
    // this might as well be inelastic
    for (const auto& d : distances)
    {
      if (d < 0.001)
      {
        std::cout << "\nRESULT: PionInelastic!"
                  << "\n//////////////////////////////////////////////////"
                  << std::endl;
        return std::string("PionInelastic");
      }
    }
  }

  // ### everything else should be pion production
  std::cout << "\nRESULT: PionProduction!"
            << "\n//////////////////////////////////////////////////"
            << std::endl;
  return std::string("PionProduction");
}

// ########################
// ### intpc
// ########################
bool InelasticSubClassifier::InTPC(const simb::MCParticle& particle)
{
  // ### get the starting and ending position for this particle
  TVector3 pos0 = particle.Position().Vect();
  TVector3 pos1 = particle.EndPosition().Vect();

  if (pos0.Z() < 0 && pos1.Z() < 0) return false;

  // ### check that the starting point is in the tpc
  if (  0 < pos0.X() && pos0.X() < 47 &&   
      -20 < pos0.Y() && pos0.Y() < 20 &&
        0 < pos0.Z() && pos0.Z() < 90) return true;

  return false;
}

}
