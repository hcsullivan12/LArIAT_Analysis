#define myana_cxx
#include "myana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TH1D* hXSRaw     = new TH1D("hXSRaw", "The Raw XS", 23, -50, 1100);
TH1D* hXSEffCorrected = new TH1D("hXSEffCorrected", "The Efficiency Corrected XS", 23, -50, 1100);
TH1D* hSignalEff = new TH1D("hSignalEff", "Signal", 23, -50, 1100);
TH1D* hRecoEff   = new TH1D("hRecoEff", "Reco", 23, -50, 1100);
TH1D* hSignalPur = new TH1D("hSignalPur", "Signal", 23, -50, 1100);
TH1D* hRecoPur   = new TH1D("hRecoPur", "Reco", 23, -50, 1100);
TH1D* hEfficiencyKE = new TH1D("hEfficiencyKE", "Efficiency", 23, -50, 1100);
TH1D* hPurityKE     = new TH1D("hPurityKE", "Purity", 23, -50, 1100);
TH1D* hElaAng       = nullptr;

void makeXSplot();
void makePlots();

void myana::Loop()
{
  std::cout << "HERE\n";
 
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
      hSignalEff->Fill(trueKEFF);
      if (didDetermineSignal) hRecoEff->Fill(trueKEFF);     
    }

    // purity
    if (didDetermineSignal)
    {
      hRecoPur->Fill(trueKEFF);
      if (isTrackSignal) hSignalPur->Fill(trueKEFF);
    }
  }

  // make other plots
  makePlots();

  // make xs plot
  makeXSplot();
  
  TFile g("piMinusAna.root", "RECREATE");
  hXSRaw->Write(); 
  hEfficiencyKE->Write();
  hPurityKE->Write();
  hXSEffCorrected->Write();
  //hElaAng->Write();
}


void makePlots()
{
  // eff.
  for (int iBin = 1; iBin <= hSignalEff->GetNbinsX(); iBin++)
  {
    if (hSignalEff->GetBinContent(iBin) == 0) continue;

    auto n = hRecoEff->GetBinContent(iBin);
    auto d = hSignalEff->GetBinContent(iBin);

    hEfficiencyKE->SetBinContent(iBin, n/d);
  }

  // purity 
  for (int iBin = 1; iBin <= hSignalPur->GetNbinsX(); iBin++)
  {
    if (hRecoPur->GetBinContent(iBin) == 0) continue;

    auto n = hSignalPur->GetBinContent(iBin);
    auto d = hRecoPur->GetBinContent(iBin);

    hPurityKE->SetBinContent(iBin, n/d);
  } 


}


void makeXSplot()
{
  TFile f("piminusanatree.root", "READ");

  TH1D* hXSInc = nullptr;
  TH1D* hXSInt = nullptr;

  f.GetObject("pionxscalc/hRecoMCIncidentKE", hXSInc);
  f.GetObject("pionxscalc/hRecoMCInteractingKE", hXSInt);
  f.GetObject("pionangularres/hElasticAngle", hElaAng);

  if (!hElaAng) std::cout << "NOPE1\n";

  if (!hXSInc || !hXSInt) {std::cout << "NOPE2\n"; return;}

   // constants for xs calculation
   float RHO            = 1396; //kg/m^3
   float MOLAR_MASS     = 39.95; //g/mol
   float G_PER_KG       = 1000; 
   float AVOGADRO       = 6.022e+23; //number/mol
   float NUMBER_DENSITY = (RHO*G_PER_KG/MOLAR_MASS)*AVOGADRO;
   float SLAB_WIDTH     = 0.0047; //in m
   float M2_PER_BARN    = 1e-28; 

  for (int iBin = 1; iBin <= hXSInc->GetNbinsX(); iBin++)
  {
    if (hXSInc->GetBinContent(iBin) == 0) continue;

    // ### our cross section
    float tempXS = (hXSInt->GetBinContent(iBin)/hXSInc->GetBinContent(iBin)) * (1/NUMBER_DENSITY) * (1/SLAB_WIDTH) * (1/M2_PER_BARN);
    hXSRaw->SetBinContent(iBin, tempXS);

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
      float xs    = hXSRaw->GetBinContent(iBin);
      float totalError = xs * ( term1 + term2 ); 

      hXSRaw->SetBinError(iBin,totalError);
    }
  }

  // apply efficiency correction 
  for (int iBin = 1; iBin <= hXSEffCorrected->GetNbinsX(); iBin++)
  {
    if (hEfficiencyKE->GetBinContent(iBin) == 0) continue;
    hXSEffCorrected->SetBinContent(iBin, hXSRaw->GetBinContent(iBin)/hEfficiencyKE->GetBinContent(iBin));
    hXSEffCorrected->SetBinError(iBin, hXSRaw->GetBinError(iBin)/hEfficiencyKE->GetBinContent(iBin));
  }
}
