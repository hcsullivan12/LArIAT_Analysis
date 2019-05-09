#define piInel_10000e_cxx
#include "piInel_10000e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

/////////////////////////
// USEFUL VARIABLES
size_t nPionsSim(0);
size_t nPionsInteractFV(0);
size_t nPionsInelastic(0);

// tpc dim and fiducial volume definition
TVector3 tpcDimVec( 47.0, 40.0, 90.0 );
TVector3 fvCutsVec( 5.0, 5.0, 5.0 ); // (x,y,z)
using Vec2D = std::vector<std::vector<double>>;
Vec2D fvVec = { std::vector<double>(  fvCutsVec[0],              tpcDimVec[0]-fvCutsVec[0] ),
                std::vector<double>( -tpcDimVec[1]+fvCutsVec[1], tpcDimVec[1]-fvCutsVec[1] ),
                std::vector<double>(  fvCutsVec[3],              tpcDimVec[0]-fvCutsVec[0] ) };
/////////////////////////

void piInel_10000e::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

     ////////////////////////////////////////////////////////
     // Flags
     bool passedFVCut(false);
     bool isInelastic(false);

     nPionsSim = nentries;

     // Did this pion interact within fv?
     // primary = 0, trackId = 1
     TVector3 vertexVec(EndPointx[0], EndPointy[0], EndPointz[0]);
     if ( fvVec[0][1] < vertexVec[0] && vertexVec[0] << fvVec[0][1] && 
          fvVec[1][1] < vertexVec[1] && vertexVec[1] << fvVec[1][1] &&
          fvVec[2][1] < vertexVec[2] && vertexVec[2] << fvVec[2][1] ) passedFVCut = true;

     // Tagging inelastic events:
     //    1) MotherID == 1
     //    2) Process generated from = PionInelastic
     // If not, skip
     for (size_t thisPart = 0; thisPart < geant_list_size; thisPart++)
     {
       if ( MotherId[i] == 1 && 
            TString((*G4Process)[i]).Contains("Pion") &&
            TString((*G4Process)[i]).Contains("Inelastic") ) isInelastic = true; break;
     }
     if (!isInelastic) continue;

     ////////////////////////////////////////////////////////
   }
}
