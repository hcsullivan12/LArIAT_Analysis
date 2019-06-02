#define myana_cxx
#include "myana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TH1D* hXS = new TH1D("hXS", "The XS", 23, -50, 1100);


void makeXSplot();

void myana::Loop()
{
  // make xs plots
  makeXSplot();

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;



    // would like to make some efficiency/purity plots 


    // efficiency
    if (isTrackSignal)
    {
      hSignal->Fill(trueInteractingKE);
    }
  }

  TFile g("piMinusAna.root", "RECREATE");
  hXS->Write(); 

}

void makeXSplot()
{
  TFile f("piminusanatree.root", "READ");

  TH1D* hXSInc = nullptr;
  TH1D* hXSInt = nullptr;

  f.GetObject("anatree/hRecoMCIncidentKE", hXSInc);
  f.GetObject("anatree/hRecoMCInteractingKE", hXSInt);

  if (!hXSInc || !hXSInt) {std::cout << "NOPE\n"; return;}

   // constants for xs calculation
   float RHO            = 1396; //kg/m^3
   float MOLAR_MASS     = 39.95; //g/mol
   float G_PER_KG       = 1000; 
   float AVOGADRO       = 6.022e+23; //number/mol
   float NUMBER_DENSITY = (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO;
   float SLAB_WIDTH     = 0.0047; //in m
   float M2_PER_BARN    = 1e-28; 

  for (int iBin = 0; iBin <= hXSInc->GetNbinsX(); iBin++)
  {
    if (hXSInc->GetBinContent(iBin) == 0) continue;

    // ### our cross section
    float tempXS = (hXSInt->GetBinContent(iBin)/hXSInc->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hXS->SetBinContent(iBin, tempXS);

    // ### incident taken as poissonian
    float denomError = std::sqrt(hXSInc->GetBinContent(iBin));
    float denom      = hXSInc->GetBinContent(iBin);
    if (denom == 0) continue; 
    float term2 = denomError/denom;

    auto intCounts = hXSInt->GetBinContent(iBin);
    auto incCounts = hXSInc->GetBinContent(iBin);
    float var      = intCounts*( 1 - intCounts/incCounts );
    float numError = std::sqrt(var);
    float num      = intCounts;

    if (num)
    {
      float term1 = numError/num;
      float xs    = hXS->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 ); 

      hXS->SetBinError(iBin,totalError);
    }
  } 
}
