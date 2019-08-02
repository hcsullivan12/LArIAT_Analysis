////////////////////////////////////////////////////////////////////////
// Class:       McTruthMcTruthThinSlabPionInelasticXS
// Module Type: analyzer
// File:        McTruthThinSlabPionInelasticXS_module.cc
//
// This module is a modification of TrueXS_module for pion inelastic.
//
// Hunter Sullivan hunter.sullivan@mavs.uta.edu
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

// root includes
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

namespace lariat
{

class McTruthThinSlabPionInelasticXS : public art::EDAnalyzer {
public:
  explicit McTruthThinSlabPionInelasticXS(fhicl::ParameterSet const & p);

  McTruthThinSlabPionInelasticXS(McTruthThinSlabPionInelasticXS const &) = delete;
  McTruthThinSlabPionInelasticXS(McTruthThinSlabPionInelasticXS &&) = delete;

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
  TH1D *hDe;
  TH1D *hDx;
  TH1D *hDeDx;
  TH1D *hDeUniform;
  TH1D *hDxUniform;
  TH1D *hDeDxUniform;
  TH1D *hDeltaE;
  TH1D *hSimIdeDist;
  TH1D *hUniformDistances;
  TH1D *hInteractingKe;
  TH1D *hIncidentKe; 
  TH1D *hCrossSection;
  TH1D *hKeAtTpcFF; 
  TH1D *hInitialKe; 
  TH1D *hInitialPz; 
  TH2D *hXZ; 
  TH2D *hYZ; 
  TH2D *hXZPre; 
  TH2D *hYZPre; 
  TH2D *hDeVsDx; 
  TH2D *hDeVsKe; 
};

//-----------------------------------------------------------------------------------------------------
// constructor
McTruthThinSlabPionInelasticXS::McTruthThinSlabPionInelasticXS(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------------------------------------
// reconfigure
void McTruthThinSlabPionInelasticXS::reconfigure(fhicl::ParameterSet const & p)
{}

//-----------------------------------------------------------------------------------------------------
// analyze
void McTruthThinSlabPionInelasticXS::analyze(art::Event const & evt)
{
  // Geometry
  art::ServiceHandle<geo::Geometry> geom;
  // Backtracker to recover truth information
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // Loop over the g4 particles 
  for (size_t p = 0; p < plist.size(); p++)
  {
    // Get the true particle and it's process, skip if not primary
    auto mcParticle = plist.Particle(p);
    if ( !(mcParticle->Process().find("primary") != std::string::npos) ) continue;
 
    // Get simIDE associated to the primary
    auto primSimIDE = bt->TrackIdToSimIDEs_Ps(mcParticle->TrackId(), geom->View(0));

    // Order them in order of increasing Z
    std::map<double, sim::IDE> orderedSimIDE;
    for (auto ide : primSimIDE ) orderedSimIDE[ide->z] = *ide;

    // Storing the kinetic energy and momentum most upstream
    simb::MCTrajectory trueTrj = mcParticle->Trajectory();
    auto     firstInTPCpoint = trueTrj.begin();
    double   mass            = mcParticle->Mass()*1000; // convert to MeV
    TVector3 momentum0       = firstInTPCpoint->second.Vect()*1000; // convert to MeV
    double   kinEn0          = TMath::Sqrt( momentum0.Mag()*momentum0.Mag() + mass ) - mass; 

    hInitialKe->Fill(kinEn0);
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

    // If we didn't enter TPC, we're done here
    if (firstInTPCpoint == trueTrj.begin()) continue;
    //nEvtsPostTPC++;
    hXZ->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).X());
    hYZ->Fill((firstInTPCpoint->first).Z(), (firstInTPCpoint->first).Y());

    // Storing the interaction type
    std::string interactionLabel("");
    
    // Identify the last interesting point in TPC
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
        // I'm not interested in
        if (tempProcess.find("CoulombScat")  != std::string::npos) continue;
        if (tempProcess.find("LArVoxel")     != std::string::npos) continue;
        if (tempProcess.find("hIoni")        != std::string::npos) continue;
  
        if (fVerbose) std::cout << tempProcess << std::endl;

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
        break;
      }
    }

    // If I didn't find anything interesting in the intereaction map check the daughters
    if ( !keepInteraction )
    {
      // Loop on the daughters 
      for (size_t d = 0; d < plist.size(); ++d) 
      {
        auto mcDaught = plist.Particle(d);
        // I'm not interested in 
        if (mcDaught->Mother()  != 1 ) continue;
        if ((mcDaught->Process()).find("CoulombScat")!= std::string::npos) continue;
        if (tempProcess.find("hIoni")!= std::string::npos) continue;
        if (tempProcess.find("lastic")!= std::string::npos) continue;

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
              << "\nFirst position:       (" << firstPos.X() << ", " << firstPos.Y() << ", " << firstPos.Z() << ")"
              << "\nLast position:        (" << lastPos.X()  << ", " << lastPos.Y()  << ", " << lastPos.Z()  << ")"
              << "\n//////////////////////////////////////////////////"
              << std::endl;

    // We want to chop up the points between the first and last uniformely
    // Order them in increasing Z
    std::map<double, TVector3> orderedUniformTrjPts;
    // We want the first and last uniform point to coincide with the 
    // the first and last points we just found 
    orderedUniformTrjPts[firstPos.Z()] = firstPos;
    orderedUniformTrjPts[lastPos.Z()]  = lastPos;

    // Calculate the number of points between first and last position
    int nPts = (int) (totalLength/(100*SLAB_WIDTH));
    for (int iPt = 1; iPt <= nPts; iPt++)
    {
      auto newPoint = firstPos + iPt*((100*SLAB_WIDTH)/totalLength) * (lastPos - firstPos);
      orderedUniformTrjPts[newPoint.Z()] = newPoint;
    }

    // Calculate the initial kinetic energy
    auto   initialMom = firstInTPCpoint->second.Vect()*1000; // convert to MeV
    double initialKE  = TMath::Sqrt( initialMom.Mag()*initialMom.Mag() + mass*mass ) - mass; 
    hKeAtTpcFF->Fill(initialKE);
    double kineticEnergy = initialKE;

    // Start filling interacting and incident histograms
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

      hDeVsDx->Fill(currentDepEnergy,(currentPos.Z()-oldPos.Z()) );
      hDeVsKe->Fill(currentDepEnergy,kineticEnergy);
      hIncidentKe->Fill(kineticEnergy);

      hDeUniform->Fill(currentDepEnergy);
      hDxUniform->Fill(uniformDist);
      hDeDxUniform->Fill(currentDepEnergy/uniformDist);
    }//<--- Loop on OrderedPoints

    // fill the Inelastic and Total Interacting with the last point
    if ( interactionLabel.find("Inelastic")!= std::string::npos )
    {
      hInteractingKe->Fill(kineticEnergy);
    }

    keepInteraction = false;
  }//<--- End loop over mc particles
}

//-----------------------------------------------------------------------------------------------------
// beginjob
void McTruthThinSlabPionInelasticXS::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hDe   = tfs->make<TH1D>("h_DE","h_DE; Energy Deposited [MeV]",200, 0,100);   
  hDx   = tfs->make<TH1D>("h_DX","h_DX; Distance between points  [cm]",400, 0,20);   
  hDeDx = tfs->make<TH1D>("h_DEDEX","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
  hDeUniform        = tfs->make<TH1D>("h_DEUniform","h_DE; Energy Deposited [MeV]",200, 0,100);   
  hDxUniform        = tfs->make<TH1D>("h_DXUniform","h_DX; Distance between points  [cm]",400, 0,20);   
  hDeDxUniform      = tfs->make<TH1D>("h_DEDEXUniform","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   
  hDeltaE           = tfs->make<TH1D>("h_DeltaE","h_DeltaE; dEDep - TrjDE [MeV/cm]",500, -1000,1000);   
  hSimIdeDist       = tfs->make<TH1D>("h_SimIDEDist","h_SimIDEDist; h_SimIDEDist [cm]",1000, 0,10);   
  hUniformDistances = tfs->make<TH1D>("h_UniformDistances","h_UniformDistances; Distance between uniform points  [cm]",500, 0,5);   
  hInitialPz     = tfs->make<TH1D>("hInitialPz"    , "Initial Pz [MeV/c]"    , 42, -100, 2000);
  hInitialKe     = tfs->make<TH1D>("hInitialKe"    , "Initial Kinetic Energy [MeV]"    , 42, -100, 2000);
  hKeAtTpcFF     = tfs->make<TH1D>("hKeAtTpcFF"    , "Kinetic Energy @ TPC FF [MeV]"   , 42, -100, 2000);
  hIncidentKe    = tfs->make<TH1D>("hIncidentKe"   , "Incident Kinetic Energy [MeV]"   , 42, -100, 2000); 
  hInteractingKe       = tfs->make<TH1D>("hInteractingKe", "Inelastic Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hCrossSection   = tfs->make<TH1D>("hCrossSection" , "Inelastic Cross-Section [barn]"   , 42, -100, 2000);
  hXZ      = tfs->make<TH2D>("hXZ"     , "hXZ"    , 110, -100, 10, 200, -100, 100);  
  hYZ      = tfs->make<TH2D>("hYZ"     , "hYZ"    , 110, -100, 10, 200, -100, 100); 
  hXZPre   = tfs->make<TH2D>("hXZPre"  , "hXZPre" , 110, -100, 10, 200, -100, 100); 
  hYZPre   = tfs->make<TH2D>("hYZPre"  , "hYZPre" , 110, -100, 10, 200, -100, 100); 
  hDeVsDx  = tfs->make<TH2D>("hDeVsDx"  , "hDeVsDx" , 504, -1, 50, 1100, -10, 100); 
  hDeVsKe  = tfs->make<TH2D>("hDeVsKe"  , "hDeVsKe" , 504, -1, 50, 220,  -10, 1000); 
}

//-----------------------------------------------------------------------------------------------------
// endjob
void McTruthThinSlabPionInelasticXS::endJob()
{
  // Calculate the Cross Section

  for ( int iBin = 1; iBin <= hIncidentKe->GetNbinsX(); ++iBin )
  {
    // If an incident bin is equal to zero then skip that bin 
    if ( hIncidentKe->GetBinContent(iBin) == 0 )continue; // fix to ensure that no Infinities are propagated to pad
   
    float tempXS = (hCrossSection->GetBinContent(iBin)/hIncidentKe->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hCrossSection->SetBinContent(iBin, tempXS);

    // incident taken as poissonian
    float denomError = std::sqrt(hIncidentKe->GetBinContent(iBin));
    float denom      = hIncidentKe->GetBinContent(iBin);
    if (denom == 0) continue; 
    float term2 = denomError/denom;

    // interacting taken as binomial
    auto intCounts = ourIntHists[h]->GetBinContent(iBin);
    auto incCounts = hIncidentKe->GetBinContent(iBin);
    float var      = intCounts*( 1 - intCounts/incCounts );
    auto numError  = std::sqrt(var);
    auto num       = intCounts;

    // protection against dividing by zero 
    if (num != 0)
    {
      // errors are not independent 
      float term1 = numError/num;
      float xs    = hCrossSection[h]->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 ); 
    
      hCrossSection->SetBinError(iBin,totalError);
    }
  }//<---End iBin Loop
}

DEFINE_ART_MODULE(McTruthThinSlabPionInelasticXS)

}
