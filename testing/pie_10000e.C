// ###############################################################
// ###############################################################
// ### NOTES ON PION INELASTIC: IDEAS, QUESTIONS, AND COMMENTS ###
// 
// 1) Definition
// 	* Pion must still be present after interaction
// 	* There must be at least two charged particles
// 	  (including primary) leaving vertex 
// 	  (this should eliminate elastic scattering) 
// 	  or point like gammas (neutron capture)? 
// 	* Should not see any showers
// 	  (Should eliminate charge exchange)
// 	* Pion production and absoprtion may be the more difficult
// 	  backgrounds
//
// I am interested in applying some classification algorithm here.
// A couple choices are decision trees and nnets.
// Will first need to find some good observables for input. 
//
// ###############################################################
// ###############################################################

#define pie_10000e_cxx
#include "pie_10000e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ==================================================================================================
// ====================================     PROTOTYPES     ==========================================
// ==================================================================================================
void MakePlots();
bool InFV(const TVector3& );

// ========================================================================================================
// ====================================     USEFUL VARIABLES     ==========================================
// ========================================================================================================
// Counters
size_t nPionsSim(0);
size_t nPionsInteractFV(0);
size_t nPionsInelastic(0);

// tpc dim and fiducial volume definition
TVector3 tpcDimVec( 47.0, 40.0, 90.0 );
TVector3 fvCutsVec( 5.0, 5.0, 5.0 ); // (x,y,z)
using Vec2D = std::vector<TVector3>;
Vec2D fvVec = { TVector3 (      fvCutsVec(0),                  tpcDimVec(0)-fvCutsVec(0), 0 ),
                TVector3 ( -0.5*tpcDimVec(1)+fvCutsVec(1), 0.5*tpcDimVec(1)-fvCutsVec(1), 0 ),
                TVector3 (      fvCutsVec(2),                  tpcDimVec(2)-fvCutsVec(2), 0 ) };

// Identifiers
std::set<std::string> processes;
std::set<int>         ufos;

// ==================================================================================================
// ====================================     HISTOGRAMS     ==========================================
// ==================================================================================================
TH1F *fvCutX = new TH1F("fvCutX", "fvCutX", 10,   0, 50);
TH1F *fvCutY = new TH1F("fvCutY", "fvCutY",  8, -20, 20);
TH1F *fvCutZ = new TH1F("fvCutZ", "fvCutZ", 20,   0, 100);
bool drawFVCut(false);

TH1S *fsPiP = new TH1S("fsPiP", "fsPiP", 10, 0, 10);
TH1S *fsPiM = new TH1S("fsPiM", "fsPiM", 10, 0, 10);
TH1S *fsPi0 = new TH1S("fsPi0", "fsPi0", 10, 0, 10);
TH1S *fsPro = new TH1S("fsPro", "fsPro", 30, 0, 30);
TH1S *fsNeu = new TH1S("fsNeu", "fsNeu", 30, 0, 30);
TH1S *fsMuM = new TH1S("fsMuM", "fsMuM", 10, 0, 10);
TH1S *fsMuP = new TH1S("fsMuP", "fsMuP", 10, 0, 10);
TH1S *fsKaP = new TH1S("fsKaP", "fsKaP", 10, 0, 10);
TH1S *fsTot = new TH1S("fsTot", "fsKaP", 10, 0, 10);
bool drawFSMult(false);

TH1F *fsProAngle = new TH1F("fsProAngle", "fsProAngle", 60, 0, 180);
bool drawAngles(false);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Main loop
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    cout << "EVENT = " << event << endl;

    // ####################################################
    // ### Register processes from daughters of primary ###
    // ####################################################
    for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++) 
    { if (Mother[thisPart] == 1) processes.emplace((*G4Process)[thisPart]); }

    // #############
    // ### Flags ###
    // #############
    bool isInelastic(false);
    nPionsSim = nentries; 
    
    // #######################################
    // ### Get the vertices for this event ###
    // #######################################
    std::vector<Vertex> theVertices;
    GetVertices(theVertices);
    for (const auto& v : theVertices) cout << v << endl;

    // ######################
    // ### Appling FV Cut ###
    // ######################
    ApplyFVCut(theVertices);
    if (theVertices.size() == 0) continue;
    nPionsInteractFV++;

    // #################################
    // ### Characterize the vertices ###
    // #################################
    CharacterizeVertices(theVertices);

    // ###################################################################
    // ### Make deflection angle distributions for coloumb and elastic ###
    // ###################################################################
    //AnaColAndElastic();

    // ###########################
    // ### Apply Inelastic Cut ###
    // ###########################
    for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)
    { 
      if (isInelastic) break;
      // Tagging inelastic events:
      //    1) Mother == 1
      //    2) Process generated from = PionInelastic
      // If not, skip
      if ( Mother[thisPart] == 1 && 
           TString((*G4Process)[thisPart]).Contains("pi") &&
           TString((*G4Process)[thisPart]).Contains("Inelastic") ) 
      { isInelastic = true; nPionsInelastic++; }
    } // end inelastic check
    if (!isInelastic) continue;

    // ######################################################################
    // ### This should be an inelastic event that has passed fiducial cut ###
    // ######################################################################
    //AnaInelastic();
  }// end event loop

  std::cout << endl
            << "TRUTH LEVEL"                                                      << endl
            << "Number of pions simulated:                  " << nPionsSim        << endl
            << "Number of pions interacted in FV:           " << nPionsInteractFV << endl
            << "Number of pions with inelastic interaction: " << nPionsInelastic  << endl
            << endl;

  //for (const auto& p : processes) std::cout << p << std::endl;
  //for (const auto& u : ufos) std::cout << u << std::endl;
  MakePlots();
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Checking if in FV
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bool InFV(const TVector3& vertexVec)
{
  if ( fvVec[0](0) < vertexVec[0] && vertexVec[0] < fvVec[0](1) && 
       fvVec[1](0) < vertexVec[1] && vertexVec[1] < fvVec[1](1) &&
       fvVec[2](0) < vertexVec[2] && vertexVec[2] < fvVec[2](1) ) return true;
  return false;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%  Applying FV Cut
// %%%  We're only looking for an interaction
// %%%  in the fiducial volume
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::ApplyFVCut(std::vector<Vertex>& theVertices)
{
  // Remember, all the vertices are attached to the primary.
  // So if at least one of the vertices is in the FV, we
  // keep it
  theVertices.erase( std::remove_if(theVertices.begin(), 
                                    theVertices.end(),
                                    [](const Vertex& v) { return !InFV(v.GetVertex()); } ),
                     theVertices.end() );
  // If it's not empty, it passed
  for (const auto& v : theVertices) 
  { 
    auto vertexVec = v.GetVertex();
    fvCutX->Fill(vertexVec[0]);
    fvCutY->Fill(vertexVec[1]);
    fvCutZ->Fill(vertexVec[2]);
  }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%    Get vertices
// %%%    The point here is to construct our vertices in a more useful
// %%%    data structure. The vertices are attached to the primary
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::GetVertices(std::vector<Vertex>& theVertices)
{
  // The interaction point variable wasn't giving sensible results
  // So our vertices will be the starting points of particles
  // with Mother == primary
  for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)  
  {
    if ( Mother[thisPart] == 1 )
    {
      // Where did this particle originate?
      // If it's not in the FV, we skip
      std::string type = (*G4Process)[thisPart]; 
      TVector3 thisVertex(StartPointx[thisPart], StartPointy[thisPart], StartPointz[thisPart]);
      // Add this is we don't have any yet
      if (theVertices.size() == 0) 
      { 
        Vertex v(thisVertex, type);
        v.AddDaughter(thisPart);
        theVertices.push_back(v); 
        continue; 
      }

      // We don't want to add duplicates
      std::multimap<double, size_t> temp;
      size_t id(0);
      for (const auto& v : theVertices)
      {
        auto theVertex = v.GetVertex();
        temp.emplace( (theVertex-thisVertex).Mag(), id ); 
      }
      for (const auto& t : temp) 
      {
        // If the difference is zero, we already have this vertex
        // just add the daughter ID
        if (t.first < 0.01) theVertices[t.second].AddDaughter(thisPart);
        else
        {
          Vertex v(thisVertex, type);
          v.AddDaughter(thisPart);
          theVertices.push_back(v); 
        }
      }// <--- End loop over temp
    }// <--- End if mother is prim
  }// <-- End loop over g4 particles 
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Characterizing the vertices
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::CharacterizeVertices(std::vector<Vertex>& theVertices)
{

  
  for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)  
  {
    if ( Mother[thisPart] == 1 )
    {
      // Where did this particle originate?
      // If it's not in the FV, we skip 
      TVector3 thisVertex(StartPointx[thisPart], StartPointy[thisPart], StartPointz[thisPart]);
      if (!InFV(thisVertex)) continue;
      
      // Now this vertex should correspond to only one of the vertices
      // in our list of vertices
      std::multimap<double, size_t> temp;
      size_t id(0);
      for (const auto& v : theVertices) { temp.emplace( (thisVertex-v.GetVertex()).Mag(), id ); id++; }
      size_t count(0);
      for (const auto& v : temp) 
      { 
        if (v.first < 0.01) count++;
      }
      if (count == 0) 
      {
        cout << "------------\n" << NumberDaughters[0] << " " << theVertices.size() << endl; 
        for (const auto& v : theVertices) cout << v << endl;
        thisVertex.Print();
        cout << pdg[thisPart] << " " << (*G4Process)[thisPart] << endl;
      }
    }
  }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Looking at Inelastic
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::AnaInelastic()
{
  // Classification:
  //    1) Species of FS particles
  size_t nPiP(0), nPiM(0), nPi0(0), 
         nPro(0), nNeu(0), nGam(0),
         nMuM(0), nMuP(0), nElP(0),
         nElM(0), nKaP(0);
  for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)
  {
    // Make sure we can see this particle
    TVector3 vertexVec(StartPointx[thisPart], StartPointy[thisPart], StartPointz[thisPart]);
    if (!InFV(vertexVec)) continue;
    TVector3 momUnitVec( EndPointx[thisPart]-StartPointx[thisPart],
                         EndPointy[thisPart]-StartPointy[thisPart],
                         EndPointz[thisPart]-StartPointz[thisPart] );
    momUnitVec = momUnitVec.Unit();
    switch(pdg[thisPart])
    {
      case  211: nPiP++; break;
      case -211: nPiM++; break;
      case  111: nPi0++; break; 
      case 2212: nPro++; break; 
      case 2112: nNeu++; break;
      case   22: nGam++; break;
      case   13: nMuM++; break;
      case  -13: nMuP++; break;
      case   11: nElM++; break;
      case  -11: nElP++; break;
      case  321: nKaP++; break;
      default: ufos.emplace(pdg[thisPart]);
    }

    // angles 
    if (Mother[thisPart] == 1 && pdg[thisPart] == 2212)
    {
               
    } 
  }

  fsPiP->Fill(nPiP);
  fsPiM->Fill(nPiM);
  fsPi0->Fill(nPi0);
  fsPro->Fill(nPro);
  fsNeu->Fill(nNeu);
  fsMuM->Fill(nMuM);
  fsMuP->Fill(nMuP);
  fsKaP->Fill(nKaP);
  size_t chargedMult = nPiP+nPiM+nPro+nMuM+nMuP+nKaP;
  fsTot->Fill(chargedMult);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Looking at Coulomb and Elastic
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void pie_10000e::AnaColAndElastic()
{
  // ###############################################
  // ### Fill the deflection angle distributions ###
  // ###############################################
  std::vector<TVector3> vertexVec;
  //vertexVec.push_back
  for (const auto& v : *InteractionPoint) cout << "Int point = " <<  MidPosX[0][v] << " " << MidPosY[0][v] << " " << MidPosZ[0][v] << " " << endl;
  for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)
  {
    // Use the start position of the daughter
    if ( Mother[thisPart] == 1 && 
         //((*G4Process)[thisPart] == "hadElastic" ||
         1 )//(*G4Process)[thisPart] == "CoulombScat") ) 
    {
      cout << StartPointx[thisPart] << " " << StartPointy[thisPart] << " " << StartPointz[thisPart] << endl;
    }   
  }

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%     Making plots
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MakePlots()
{
  if (drawFVCut)
  {
    TCanvas *c1 = new TCanvas("c1", "FVCutX", 800, 800);
    fvCutX->Draw();
    TCanvas *c2 = new TCanvas("c2", "FVCutY", 800, 800);
    fvCutY->Draw();
    TCanvas *c3 = new TCanvas("c3", "FVCutZ", 800, 800);
    fvCutZ->Draw();
  }

  if (drawFSMult)
  {
    // Groups:
    //   1) Protons, neutrons
    //   2) Pions
    //   3) Everything else
    gStyle->SetOptStat(0); 
    TCanvas *c4 = new TCanvas("c4", "FSMult_P&N", 800, 800);
    c4->SetLogy();
    std::vector<TH1S*> hists;
    hists.push_back(fsPro);
    hists.push_back(fsNeu);
    std::sort(hists.begin(), hists.end(), [](const TH1S* lh, const TH1S* rh) { return lh->GetMaximum() > rh->GetMaximum(); });
    for (auto& h : hists) { h->SetLineWidth(3); h->SetFillColor(0); }
    fsPro->SetLineColor(2);
    fsNeu->SetLineColor(4);
    TLegend* leg4 = new TLegend(0.1,0.7,0.48,0.9);
    leg4->AddEntry(fsPro,"p","f");
    leg4->AddEntry(fsNeu,"n","f");
    leg4->SetBorderSize(0);
    auto iter = hists.begin();
    (*iter)->SetTitle("Final State Multiplicities"); 
    (*iter)->Draw(); iter++;
    while (iter != hists.end()) { (*iter)->Draw("same"); iter++; }
    leg4->Draw("same");

    TCanvas *c5 = new TCanvas("c5", "FSMult_Pions", 800, 800);
    c5->SetLogy();
    hists.clear();
    hists.push_back(fsPiP);
    hists.push_back(fsPiM);
    hists.push_back(fsPi0);
    std::sort(hists.begin(), hists.end(), [](const TH1S* lh, const TH1S* rh) { return lh->GetMaximum() > rh->GetMaximum(); });
    for (auto& h : hists) { h->SetLineWidth(3); h->SetFillColor(0); }
    fsPiP->SetLineColor(46);
    fsPiM->SetLineColor(36);
    fsPi0->SetLineColor(26);
    TLegend* leg5 = new TLegend(0.1,0.7,0.48,0.9);
    leg5->AddEntry(fsPiP,"#pi+","f");
    leg5->AddEntry(fsPiM,"#pi-","f");
    leg5->AddEntry(fsPi0,"#pi0","f");
    leg5->SetBorderSize(0);
    iter = hists.begin();
    (*iter)->SetTitle("Final State Multiplicities"); 
    (*iter)->Draw(); iter++;
    while (iter != hists.end()) { (*iter)->Draw("same"); iter++; }
    leg5->Draw("same");

    TCanvas *c6 = new TCanvas("c6", "FSMult_Others", 800, 800);
    c6->SetLogy();
    hists.clear();
    hists.push_back(fsMuM);
    hists.push_back(fsMuP);
    hists.push_back(fsKaP);
    std::sort(hists.begin(), hists.end(), [](const TH1S* lh, const TH1S* rh) { return lh->GetMaximum() > rh->GetMaximum(); });
    for (auto& h : hists) { h->SetLineWidth(2); h->SetFillColor(0); }
    fsMuM->SetLineColor(2);
    fsMuP->SetLineColor(4);
    fsKaP->SetLineColor(3);
    TLegend* leg6 = new TLegend(0.1,0.7,0.48,0.9);
    leg6->AddEntry(fsMuM,"#mu-","f");
    leg6->AddEntry(fsMuP,"#mu+","f");
    leg6->AddEntry(fsKaP,"K+","f");
    leg6->SetBorderSize(0);
    iter = hists.begin();
    (*iter)->SetTitle("Final State Multiplicities"); 
    (*iter)->Draw(); iter++;
    while (iter != hists.end()) { (*iter)->Draw("same"); iter++; }
    leg6->Draw("same");

    TCanvas *c7 = new TCanvas("c7", "FSMulti_All", 800, 800);
    c7->SetLogy();
    fsTot->SetTitle("Final State Multiplicity (Charged)");
    fsTot->SetLineWidth(3);
    fsTot->SetLineColor(1);
    fsTot->Draw();
  }
}

