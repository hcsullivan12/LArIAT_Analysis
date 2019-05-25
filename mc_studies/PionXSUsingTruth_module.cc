////////////////////////////////////////////////////////////////////////
// Class:       PionXSUsingTruth
// Module Type: analyzer
// File:        PionXSUsingTruth_module.cc
//
// Generated at Wed May 22 11:52:52 2019 by Hunter Sullivan using artmod
// from cetpkgsupport v1_10_02.
//
// The module calculates various pion cross sections using mc truth
// information. This is supplemented with an algorithm to subclassify 
// inelastic processes, e.g. inelastic, CX, absorption, pion production.
//
////////////////////////////////////////////////////////////////////////

// framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h" 

// LArSoft includes
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// lariatsoft includes
#include "LArIATFilterModule/InelasticSubClassifier.h"

// root includes
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

namespace lariat
{

class PionXSUsingTruth : public art::EDAnalyzer {
public:
  explicit PionXSUsingTruth(fhicl::ParameterSet const & p);

  PionXSUsingTruth(PionXSUsingTruth const &) = delete;
  PionXSUsingTruth(PionXSUsingTruth &&) = delete;

  // analyzing function
  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  
private:

  bool fVerbose = false;

  // tpc boundaries
  double TPC_BOUND_X[2] = {  0, 47};
  double TPC_BOUND_Y[2] = {-20, 20};
  double TPC_BOUND_Z[2] = {  0, 90};

  // constants for xs calculation
  float RHO            = 1396; //kg/m^3
  float MOLAR_MASS     = 39.95; //g/mol
  float G_PER_KG       = 1000; 
  float AVOGADRO       = 6.022e+23; //number/mol
  float NUMBER_DENSITY = (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO;
  float SLAB_WIDTH     = 0.0047; //in m
  float M2_PER_BARN    = 1e-28;

  // histograms
  TH1D *hDE;
  TH1D *hDX;
  TH1D *hDEDX;
  TH1D *hDEUniform;
  TH1D *hDXUniform;
  TH1D *hDEDXUniform;
  TH1D *hDeltaE;
  TH1D *hSimIDEDist;
  TH1D *hUniformDistances;
  TH1D *hInteractingKEInel;
  TH1D *hInteractingKEPionAbsorp;
  TH1D *hInteractingKEChargeExch;
  TH1D *hInteractingKEPionInel;
  TH1D *hInteractingKEPionProd; 
  TH1D *hIncidentKE; 
  TH1D *hCrossSectionInel;
  TH1D *hCrossSectionPionAbsorp;
  TH1D *hCrossSectionChargeExch;
  TH1D *hCrossSectionPionInel;
  TH1D *hCrossSectionPionProd;
  TH1D *hKEAtTPCFF; 
  TH1D *hInitialKE; 
  TH1D *hInitialPz; 
  TH2D *hXZ; 
  TH2D *hYZ; 
  TH2D *hXZPre; 
  TH2D *hYZPre; 
  TH2D *hdEVsdX; 
  TH2D *hdEVsKE; 
};

//-----------------------------------------------------------------------------------------------------
// ### constructor
PionXSUsingTruth::PionXSUsingTruth(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------------------------------------
// ### reconfigure
void PionXSUsingTruth::reconfigure(fhicl::ParameterSet const & p)
{}

//-----------------------------------------------------------------------------------------------------
// ### analyze
void PionXSUsingTruth::analyze(art::Event const & evt)
{
  // ###########################
  // ### Get useful services ###
  // ###########################

  // ### Geometry
  art::ServiceHandle<geo::Geometry> geom;
  // ### Backtracker to recover truth information
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // ##################################
  // ### Loop over the g4 particles ###
  // ##################################
  for (size_t p = 0; p < plist.size(); p++)
  {
    // ### Get the true particle and it's process, skip if not primary
    auto        mcParticle = plist.Particle(p);
    if ( !(mcParticle->Process().find("primary") != std::string::npos) ) continue;
 
    // ### Get simIDE associated to the primary
    auto primSimIDE = bt->TrackIdToSimIDEs_Ps(mcParticle->TrackId(), geom->View(0));

    // ### Order them in order of increasing Z
    std::map<double, sim::IDE> orderedSimIDE;
    for (auto ide : primSimIDE ) orderedSimIDE[ide->z] = *ide;

    // ### Storing the kinetic energy and momentum most upstream
    simb::MCTrajectory trueTrj = mcParticle->Trajectory();
    auto     firstInTPCpoint = trueTrj.begin();
    double   mass            = mcParticle->Mass()*1000; // convert to MeV
    TVector3 momentum0       = firstInTPCpoint->second.Vect()*1000; // convert to MeV
    double   kinEn0          = TMath::Sqrt( momentum0.Mag()*momentum0.Mag() + mass ) - mass; 

    hInitialKE->Fill(kinEn0);
    hInitialPz->Fill(momentum0.Z());

    // ### Identify the first traj point in TPC
    for (auto itTrj = trueTrj.begin(); itTrj != std::prev(trueTrj.end()); itTrj++)
    {
      auto pos = itTrj->first;
      if (TPC_BOUND_X[0] > pos.X() || pos.X() > TPC_BOUND_X[1] ||
          TPC_BOUND_Y[0] > pos.Y() || pos.Y() > TPC_BOUND_Y[1] || 
          pos.Z() < TPC_BOUND_Z[0]) continue;
      
      firstInTPCpoint = itTrj;
      break;
    }
    //nEvtsPreTPC++;
    hXZPre->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).X());
    hYZPre->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).Y());

    // ### If we didn't enter TPC, we're done here
    if (firstInTPCpoint == trueTrj.begin()) continue;
    //nEvtsPostTPC++;
    hXZ->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).X());
    hYZ->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).Y());

    // ### Storing the interaction type
    std::string interactionLabel("");
    // ### Algorithm to identify the sub process type
    // ### Sub processes = (inelastic, ch-exch, absorp, pi-prod)
    piinelastic::InelasticSubClassifier subclassifier;
    auto lastInTPCpoint      = std::prev(trueTrj.end());
    auto interactionSubLabel = subclassifier.Classify(plist, mcParticle->TrackId());
    
  // ### Identify the last interesting point in TPC
    // The last point is a bit more complicated:
    // if there's no interaction, then it is simply the last point in the TPC
    // if there's one or more interaction points, it's the first interaction point deemed interesting (no coulomb)
    // Take the interaction Map... check if there's something there
    auto trjProcessMap =  trueTrj.TrajectoryProcesses();
    bool keepInteraction = false;
    if (trjProcessMap.size())
    {
      for (auto const& couple: trjProcessMap) 
      {
        auto tempProcess = trueTrj.KeyToProcess(couple.second); 
        std::cout << tempProcess << std::endl;
        // I'm not interested in the CoulombScat, LArVoxel, OpDetReadout (if keeping all spacepoints)
        if (tempProcess.find("CoulombScat")  != std::string::npos) continue;
        if (tempProcess.find("LArVoxel")     != std::string::npos) continue;
        if (tempProcess.find("OpDetReadout") != std::string::npos) continue;
  
        // Let's check if the interaction is in the the TPC
        auto interactionPos4D = (trueTrj.at(couple.first)).first;
        if ( TPC_BOUND_X[0] > interactionPos4D.X() || interactionPos4D.X() > TPC_BOUND_X[1] ||
             TPC_BOUND_Y[0] > interactionPos4D.Y() || interactionPos4D.Y() > TPC_BOUND_Y[1] ||
             TPC_BOUND_Z[0] > interactionPos4D.Z() || interactionPos4D.Z() > TPC_BOUND_Z[1] ) continue;
  
        // If we made it here, then this is the first interesting interaction in the TPC
        // Store the interaction label and the iterator for the final point
        interactionLabel = trueTrj.KeyToProcess(couple.second);
        lastInTPCpoint = trueTrj.begin() + couple.first; 
        keepInteraction = true;
        //nInteractInTPC++;
        break;
      }
    }

    // ### If I didn't find anything interesting in the intereaction map, let's loop back!
    if ( !keepInteraction )
    {
      std::cout << "HERERER\n";
      // Loop on the daughters 
      for (size_t d = 0; d < plist.size(); ++d) 
      {
        auto mcDaught = plist.Particle(d);
        // I'm not interested in CoulombScat
        if (mcDaught->Mother()  != 1 ) continue;
        std::cout << mcDaught->Process() << std::endl;
        if ((mcDaught->Process()).find("CoulombScat")!= std::string::npos) continue;

        // Is the daughter born inside the TPC? If yes, store the process which created it 
        simb::MCTrajectory trueDaugthTraj = mcDaught->Trajectory();	     
        TVector3 trjPoint = trueDaugthTraj.begin()->first.Vect(); 

        if ( TPC_BOUND_X[0] > trjPoint.X() || trjPoint.X() > TPC_BOUND_X[1] ||
             TPC_BOUND_Y[0] > trjPoint.Y() || trjPoint.Y() > TPC_BOUND_Y[1] ||
             TPC_BOUND_Z[0] > trjPoint.Z() || trjPoint.Z() > TPC_BOUND_Z[1] ) continue;

        interactionLabel = mcDaught->Process();
        break;
      }	  

      for (auto itTrj = std::prev(trueTrj.end()); itTrj != trueTrj.begin(); itTrj--)
      {
        auto pos = itTrj->first;
        
        if ( TPC_BOUND_X[0] > pos.X() || pos.X() > TPC_BOUND_X[1] ||
             TPC_BOUND_Y[0] > pos.Y() || pos.Y() > TPC_BOUND_Y[1] ||
             pos.Z() > TPC_BOUND_Z[1] ) continue;

        lastInTPCpoint = itTrj;
        break;
      }
    } 

    // ### Exit if the last point is the first point
    if(lastInTPCpoint == firstInTPCpoint) continue;

    // ### Exit if track is too short
    TVector3 firstPos = firstInTPCpoint->first.Vect();
    TVector3 lastPos  = lastInTPCpoint->first.Vect();
    double totalLength = (lastPos-firstPos).Mag();
    if (totalLength < (100*SLAB_WIDTH)) continue;

    std::cout << "\n//////////////////////////////////////////////////"
              << "\nPionUsingTruth...\n"
              << "\nInteraction label:     " << interactionLabel
              << "\nInteraction sublabel:  " << interactionSubLabel
              << "\nFirst position:       (" << firstPos.X() << ", " << firstPos.Y() << ", " << firstPos.Z() << ")"
              << "\nLast position:        (" << lastPos.X()  << ", " << lastPos.Y()  << ", " << lastPos.Z()  << ")"
              << "\n//////////////////////////////////////////////////"
              << std::endl;

    // ### We want to chop up the points between the first and last uniformely
    // Order them in increasing Z
    std::map<double, TVector3> orderedUniformTrjPts;
    // We want the first and last uniform point to coincide with the 
    // the first and last points we just found 
    orderedUniformTrjPts[firstPos.Z()] = firstPos;
    orderedUniformTrjPts[lastPos.Z()]  = lastPos;

    // ### Calculate the number of points between first and last position
    int nPts = (int) (totalLength/(100*SLAB_WIDTH));
    for (int iPt = 1; iPt <= nPts; iPt++)
    {
      auto newPoint = firstPos + iPt*((100*SLAB_WIDTH)/totalLength) * (lastPos - firstPos);
      orderedUniformTrjPts[newPoint.Z()] = newPoint;
    }

    // ### Calculate the initial kinetic energy
    auto   initialMom = firstInTPCpoint->second.Vect()*1000; // convert to MeV
    double initialKE  = TMath::Sqrt( initialMom.Mag()*initialMom.Mag() + mass*mass ) - mass; 
    hKEAtTPCFF->Fill(initialKE);
    double kineticEnergy = initialKE;

    // ### Start filling interacting and incident histograms
    auto old_it = orderedUniformTrjPts.begin();
    for (auto it = std::next(orderedUniformTrjPts.begin()); it != orderedUniformTrjPts.end(); it++, old_it++ )
    {
      auto oldPos        = old_it->second;
      auto currentPos    =     it->second;

      double uniformDist =  (currentPos - oldPos).Mag();
      hUniformDistances->Fill(uniformDist);

      //Calculate the energy deposited in this slice          
      auto old_iter = orderedSimIDE.begin();
      double currentDepEnergy = 0.;
      for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++)
      {
        auto currentIde = iter->second;
        if (currentIde.z < oldPos.Z()) continue;
        if (currentIde.z > currentPos.Z()) continue;
        currentDepEnergy += currentIde.energy;
      }// Determing which simIDE is within the current slice

      // avoid overfilling super tiny energy depositions
      if (currentDepEnergy/uniformDist < 0.1) continue;
      //Calculate the current kinetic energy
      kineticEnergy -= currentDepEnergy;

      hdEVsdX->Fill(currentDepEnergy,(currentPos.Z()-oldPos.Z()) );
      hdEVsKE->Fill(currentDepEnergy,kineticEnergy);
      hIncidentKE->Fill(kineticEnergy);

      hDEUniform->Fill(currentDepEnergy);
      hDXUniform->Fill(uniformDist);
      hDEDXUniform->Fill(currentDepEnergy/uniformDist);
    }//<--- Loop on OrderedPoints

    // ### fill the Inelastic and Total Interacting with the last point
    if ( interactionLabel.find("Inelastic")!= std::string::npos )
    {
      hInteractingKEInel->Fill(kineticEnergy);

           if (interactionSubLabel == "PionAbsorption") hInteractingKEPionAbsorp->Fill(kineticEnergy);
      else if (interactionSubLabel == "PionInelastic")  hInteractingKEPionInel->Fill(kineticEnergy);
      else if (interactionSubLabel == "ChargeExchange") hInteractingKEChargeExch->Fill(kineticEnergy);
      else if (interactionSubLabel == "PionProduction") hInteractingKEPionProd->Fill(kineticEnergy);
    }
      
    /*// ### fill the Elastic and Total Interacting with the last point
    if ( interactionLabel.find("Elastic")!= std::string::npos )
	  {
	    h_DeltaE ->Fill(kineticEnergy -  1000*((finTPCPoint->second).E() - mass) );
	    hInteractingKEElDep->Fill(kineticEnergy);
	    hInteractingKE->Fill(kineticEnergy);
	    auto MomentumF = finTPCPoint->second;
	    double KEF = 1000*(TMath::Sqrt(MomentumF.X()*MomentumF.X() + MomentumF.Y()*MomentumF.Y() + MomentumF.Z()*MomentumF.Z() + mass*mass ) - mass); //I want this in MeV
	    hInteractingKEEl->Fill(KEF);
	  }*/

    keepInteraction = false;
  }//<--- End loop over mc particles
}

//-----------------------------------------------------------------------------------------------------
// ### beginjob
void PionXSUsingTruth::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hDE   = tfs->make<TH1D>("h_DE","h_DE; Energy Deposited [MeV]",200, 0,100);   
  hDX   = tfs->make<TH1D>("h_DX","h_DX; Distance between points  [cm]",400, 0,20);   
  hDEDX = tfs->make<TH1D>("h_DEDEX","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
  hDEUniform        = tfs->make<TH1D>("h_DEUniform","h_DE; Energy Deposited [MeV]",200, 0,100);   
  hDXUniform        = tfs->make<TH1D>("h_DXUniform","h_DX; Distance between points  [cm]",400, 0,20);   
  hDEDXUniform      = tfs->make<TH1D>("h_DEDEXUniform","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
  hDeltaE           = tfs->make<TH1D>("h_DeltaE","h_DeltaE; dEDep - TrjDE [MeV/cm]",500, -1000,1000);   
  hSimIDEDist       = tfs->make<TH1D>("h_SimIDEDist","h_SimIDEDist; h_SimIDEDist [cm]",1000, 0,10);   
  hUniformDistances = tfs->make<TH1D>("h_UniformDistances","h_UniformDistances; Distance between uniform points  [cm]",500, 0,5);   
  hInitialPz     = tfs->make<TH1D>("hInitialPz"    , "Initial Pz [MeV/c]"    , 42, -100, 2000);
  hInitialKE     = tfs->make<TH1D>("hInitialKE"    , "Initial Kinetic Energy [MeV]"    , 42, -100, 2000);
  hKEAtTPCFF     = tfs->make<TH1D>("hKEAtTPCFF"    , "Kinetic Energy @ TPC FF [MeV]"   , 42, -100, 2000);
  hIncidentKE    = tfs->make<TH1D>("hIncidentKE"   , "Incident Kinetic Energy [MeV]"   , 42, -100, 2000); 
  hInteractingKEInel       = tfs->make<TH1D>("hInteractingKEInel", "Inelastic Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hInteractingKEPionAbsorp = tfs->make<TH1D>("hInteractingKEPionAbsorp", "Pion Absorption Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hInteractingKEPionInel   = tfs->make<TH1D>("hInteractingKEPionInel", "Pion Inelastic Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hInteractingKEChargeExch = tfs->make<TH1D>("hInteractingKEChargeExch", "Charge Exchange Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hInteractingKEPionProd   = tfs->make<TH1D>("hInteractingKEPionProd", "Pion Production Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hCrossSectionInel   = tfs->make<TH1D>("hCrossSectionInel" , "Inelastic Cross-Section [barn]"   , 42, -100, 2000);
  hCrossSectionPionAbsorp = tfs->make<TH1D>("hCrossSectionPionAbsorp" , "Pion Absorption Cross-Section [barn]"   , 42, -100, 2000);
  hCrossSectionPionInel   = tfs->make<TH1D>("hCrossSectionPionInel" , "Pion Inelastic Cross-Section [barn]"   , 42, -100, 2000);
  hCrossSectionPionProd   = tfs->make<TH1D>("hCrossSectionPionProd" , "Pion Production Cross-Section [barn]"   , 42, -100, 2000);
  hCrossSectionChargeExch = tfs->make<TH1D>("hCrossSectionChargeExch" , "Charge Exchange Cross-Section [barn]"   , 42, -100, 2000);
  hXZ      = tfs->make<TH2D>("hXZ"     , "hXZ"    , 110, -100, 10, 200, -100, 100);  
  hYZ      = tfs->make<TH2D>("hYZ"     , "hYZ"    , 110, -100, 10, 200, -100, 100); 
  hXZPre   = tfs->make<TH2D>("hXZPre"  , "hXZPre" , 110, -100, 10, 200, -100, 100); 
  hYZPre   = tfs->make<TH2D>("hYZPre"  , "hYZPre" , 110, -100, 10, 200, -100, 100); 
  hdEVsdX  = tfs->make<TH2D>("hdEVsdX"  , "hdEVsdX" , 504, -1, 50, 1100, -10, 100); 
  hdEVsKE  = tfs->make<TH2D>("hdEVsKE"  , "hdEVsKE" , 504, -1, 50, 220,  -10, 1000); 
}


//-----------------------------------------------------------------------------------------------------
// ### endjob
void PionXSUsingTruth::endJob()
{
  /*std::cout << "-------------------------------------------"                     << std::endl;
  std::cout << "True Events pre-TPC .............. "         << evtsPreTPC       << std::endl;
  std::cout << "True Events pre-TPC .............. "         << evtsInTheMiddle  << std::endl;
  std::cout << "True Events post-TPC ............. "         << evtsPostTPC      << std::endl;
  std::cout << "True Throughgoing    ............. "         << throughgoing     << std::endl;
  std::cout << "True interactingInTPC ............ "         << interactingInTPC << std::endl;
  std::cout << "-------------------------------------------"                     << std::endl;
*/
  // Calculate the Cross Section
  // ###################################################################
  // #### Looping over the exiting bins to extract the cross-section ###
  // ###################################################################

  // ### place our histograms in a single container 
  std::vector<TH1D*> ourHists;
  ourHists.push_back(hInteractingKEInel);
  ourHists.push_back(hInteractingKEPionAbsorp);
  ourHists.push_back(hInteractingKEPionInel);
  ourHists.push_back(hInteractingKEChargeExch);
  ourHists.push_back(hInteractingKEPionProd);

  // ### place our cross sections histograms in a single container
  std::vector<TH1D*> ourXSHists;
  ourXSHists.push_back(hCrossSectionInel);
  ourXSHists.push_back(hCrossSectionPionAbsorp);
  ourXSHists.push_back(hCrossSectionPionInel);
  ourXSHists.push_back(hCrossSectionChargeExch);
  ourXSHists.push_back(hCrossSectionPionProd);

  for ( int iBin = 1; iBin <= hIncidentKE->GetNbinsX(); ++iBin )
  {
    // ### If an incident bin is equal to zero then skip that bin ###
    if ( hIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
   
    // ### our cross sections
    for (size_t h = 0; h < ourHists.size(); h++)
    {
      float tempXS = (ourHists[h]->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
      ourXSHists[h]->SetBinContent(iBin, tempXS);
    }
  
    // ###########################################################
    // ### Calculating the error on the numerator of the ratio ###
    // ###########################################################

    float denomError = std::sqrt(hIncidentKE->GetBinContent(iBin));
    float denom      = hIncidentKE->GetBinContent(iBin);
    if (denom == 0) continue; 
    float term2 = denomError/denom;

    std::vector<float> ourNumError, ourNum;
    for (size_t h = 0; h < ourHists.size(); h++)
    {
      ourNumError.push_back( std::sqrt(ourHists[h]->GetBinContent(iBin)) );
      ourNum.push_back( ourHists[h]->GetBinContent(iBin) );
    }

    // ### set the errors
    for (size_t h = 0; h < ourHists.size(); h++)
    {
      // ### Putting in a protection against dividing by zero ###   
      if (ourNum[h] != 0)
      {
        float term1 = ourNumError[h]/ourNum[h];
        float xs    = ourXSHists[h]->GetBinContent(iBin);
        float totalError = xs * std::sqrt( term1*term1 + term2*term2 ); 
        ourXSHists[h]->SetBinError(iBin,totalError);
      }
    }
  }//<---End iBin Loop
}

DEFINE_ART_MODULE(PionXSUsingTruth)

}
