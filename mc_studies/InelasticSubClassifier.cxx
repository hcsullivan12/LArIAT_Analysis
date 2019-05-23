////////////////////////////////////////////////////////////
// InelasticSubClassifier
//
// This algorithm further classifies pion inelastic events
// into inelastic, charge exchange, absorption, and pion
// production.
//
// Author: Hunter Sullivan hunter.sullivan@mavs.uta.edu
////////////////////////////////////////////////////////////

namespace piinelastic 
{

// ########################
// ### constructor
// ########################
InelasticSubClassifier::InelasticSubClassifier()
{}

// ########################
// ### classify
// ########################
std::string InelasticSubClassifier::Classify(const sim::ParticleList& plist)
{
  //
  // We define sub channels for inelastic processes:
  //      1) Absorption (no pions in final state)
  //      2) Inelastic (no neutral pions, 1 same charged pion in final state)
  //      3) Charge exchange (1 nuetral pion, no charged pions in final state)
  //      4) Pion production (what's left?)
  //

  bool isPrimaryChargedPion(false);
  int primaryPdgCode(0);
  int primaryTrkId(0);
  int nChargedPionDaughters(0);
  int nChargedPDaughters(0);

  // ### loop over the g4 particles
  for (size_t p = 0; p < plist.size(); p++)
  {
    // get the true particle and it's process, skip if not primary
    auto particle = plist.Particle();
    if ( !(particle.Process().find("primary") != std::string::npos) ) continue;

    // if charged pion
    if (std::abs(particle.PdgCode()) == 211)
    {
      // change flag
      primaryChargedPion = true;

      // get the primary pdg code
      primaryPdgCode = particle.PdgCode();

      // get trackid
      primaryTrkId = particle.TrackId();

      // get last point
      primaryFinalPosition = particle.EndPosition();
    } 
  }

  // ### checkout if not a charged pion
  if (!isPrimaryChargedPion) return std::string("NoPrimaryChargedPion");

  // ### loop over the particles again checking daughters
  for (const auto& particle : (*particleHandle))
  {
    // skip if particle is not a child of the primary
    if (particle.Mother() != primaryTrackID) continue;

    // get trackid 
    int trackId = particle.TrackId();

    // get pdg code
    int pdgCode = particle.PdgCode();

    // check for neutral pion, if so we're done
    if (pdgCode == 111) nNeutralPionDaughters++;

    // check for charged pions
    if (std::abs(pdgCode) == 211) nChargedPionDaughters++;

    // check for inelastic
    if (particle.Process().find("Inelastic") != std::string::npos) 
    {
      foundInelastic = true;
      // keep this vertex
      vertices.push_back(particle->Position());
    }

    // add the process to our set
    processes.insert(particle.Process());
  }

  // ### print out the processes in for this primary
  std::cout << "\n//////////////////////////////////////////////////"
            << "\nThe processes for this primary:\n";
  for (const auto& p : processes) std::cout << p << std::endl;
  std::cout << std::endl;

  // ### checkout if no inelastic
  if (!foundInelastic) 
  {
    std::cout << "\nRESULT: no inelastic interaction!" << std::endl;
    return std::string("NoInelasticInteraction");
  }

  // ### okay... we can start to classify this event
  
  // ### we need to calculate the distance between our 
  // ### primary end point and each vertex
  std::vector<double> distances;
  for (const auto& v : vertices)
  {
    TVector3 vertex = v.Vect();
    distances.push_back( (vertex-primaryFinalPosition).Mag() );
  }

  // ### handle absorption and inelastic
  if (nChargedPionDaughters == 0 && nNeutralPionDaughters == 0)
  {
    // if the distance = 0, then we have absorption
    for (const auto& d : distances)
    {
      if (d < 0.001)
      {
        std::cout << "\nRESULT: PionAbsorption!" << std::endl;
        return std::string("PionAbsorption");
      }
    }
    // otherwise we have inelastic
    std::cout << "\nRESULT: PionInelastic!" << std::endl;
    return std::string("PionInelastic");
  }
  
  // ### handle charge exchange
  if (nNeutralPionDaughters == 1 && nChargedPionDaughters == 0)
  {
    // this should do it!
    std::cout << "\nRESULT: ChargeExchange!" << std::endl;
    return std::string("ChargeExchange");
  }

  // ### everything else should be pion production
  std::cout << "\nRESULT: PionProduction!" << std::endl;
  return std::string("PionProduction");
}

}