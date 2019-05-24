////////////////////////////////////////////////////////////
// InelasticSubClassifier
//
// This algorithm further classifies pion inelastic events
// into inelastic, charge exchange, absorption, and pion
// production.
//
// Author: Hunter Sullivan hunter.sullivan@mavs.uta.edu
////////////////////////////////////////////////////////////

// LArSoft includes
#include "larsim/MCCheater/ParticleInventoryService.h"

// c++
#include <string>

namespace piinelastic 
{

class InelasticSubClassifier
{
  public:
    InelasticSubClassifier();
    ~InelasticSubClassifier();

    const std::string Classify(const sim::ParticleList& plist, const int& primaryTrkId);
    bool InTPC(const simb::MCParticle& particle);
};

}
