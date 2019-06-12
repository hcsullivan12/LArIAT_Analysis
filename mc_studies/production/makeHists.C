TFile* anaFile;

void makeXsPlots()
{
  gStyle->SetOptStat(0);

  TFile *f2 = new TFile("../hists/g4XsPredictions.root", "READ");

  TH1D *hMCXS = nullptr;
  TCanvas *temp = nullptr;
  
  anaFile->GetObject("mc/hMCXSKE", hMCXS);
  f2->GetObject("piMinInel", temp);
  TGraph *g = (TGraph*)temp->GetListOfPrimitives()->FindObject("Graph");

  if (!g || !temp || !hMCXS) {cout << "Nope\n"; return;}

  TCanvas *c1 = new TCanvas("xsplots", "c1", 800., 800.);
  
  //hMCXS->SetMinimum(0);
  //hMCXS->SetMaximum(1100);
  //hMCXS->GetXaxis()->SetRangeUser(0,1000);
  hMCXS->GetXaxis()->SetTitle("Kinetic energy (MeV)");
  hMCXS->GetYaxis()->SetTitle("Cross section (barn)");
  g->Draw("A c");
  hMCXS->Draw("same");

  hMCXS->SetLineWidth(3);
  hMCXS->SetLineColor(kRed);
  hMCXS->SetMarkerStyle(8);
  hMCXS->SetMarkerColor(kRed);

  g->SetLineWidth(3);
  g->SetLineStyle(9);
  g->SetLineColor(kRed);

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(g, "G4Prediction Inelastic XS", "l");
  legend->AddEntry(hMCXS,"MC inelastic XS","lp");
  legend->Draw("same");
}

void makeAnglePlots()
{
  gStyle->SetOptStat(0);

  TH1D* hEla     = nullptr;
  TH1D* hInelNDa = nullptr;

  anaFile->GetObject("mc/hMCElasticAngle", hEla);
  anaFile->GetObject("mc/hMCInelasticOneVisDAngle", hInelNDa);

  if (!hEla || !hInelNDa) {cout << "Nope\n"; return;}

  TCanvas *c1 = new TCanvas("angles", "c1", 800, 800);

  hEla->GetXaxis()->SetTitle("Angle (degrees)");
  hEla->SetTitle("Angle Between Primary and Secondary");
  hEla->Draw();
  hInelNDa->Draw("same");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(hEla,"Elastic (primary/primary)","l");
  legend->AddEntry(hInelNDa,"Inelastic (primary/one visible daughter)","l");
  legend->Draw("same");
  legend->SetLineColor(0);

  c1->SetLogy();

  hEla->SetLineWidth(4);
  hEla->SetLineColor(46);
  hInelNDa->SetLineWidth(4);
  hEla->SetLineColor(42);
}

void makeSmearingPlots()
{
  TH2D* hS = nullptr;

  anaFile->GetObject("reco/hSmearingMatrix", hS);
  if (!hS) {cout << "Nope\n"; return;}

  TCanvas *c1 = new TCanvas("smearing", "c1", 800, 800);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetPalette(kBlueRedYellow);
  hS->Draw("colz text");

  hS->GetXaxis()->SetTitle("Reconstructed KE (MeV)");
  hS->GetYaxis()->SetTitle("True KE (MeV)");

  // get avg smearing along diagonal
  float sum = 0;
  size_t count = 0;
  for(size_t iBinX = 0; iBinX < hS->GetXaxis()->GetNbins(); iBinX++)
  {
    for(size_t iBinY = 0; iBinY < hS->GetYaxis()->GetNbins(); iBinY++)
    {
      if (iBinX != iBinY) continue;
      auto content = hS->GetBinContent(iBinX, iBinY);
      sum = sum + content;
      count++;
    }
  }
  cout << "Average smearing along diagonal = " << sum/count << endl;
  
}

void makeEfficiencyPlots()
{
  TH1D* hMCInc = nullptr;
  TH1D* hMCInt = nullptr;
  TH1D* hRecoInc = nullptr;
  TH1D* hRecoInt = nullptr;

  anaFile->GetObject("mc/hMCIncidentKE", hMCInc);
  anaFile->GetObject("mc/hMCInteractingKE", hMCInt);
  anaFile->GetObject("reco/hRecoIncidentKE", hRecoInc);
  anaFile->GetObject("reco/hRecoInteractingKE", hRecoInt);

  if (!hMCInc || !hMCInt || !hRecoInc || !hRecoInt) {cout<<"Nope\n"; return;}
  
  size_t nBins = hMCInc->GetXaxis()->GetNbins();
  TH1D* hEffInc = new TH1D("hEffInc", "Incident Efficiency", nBins, 0, nBins*50. );
  TH1D* hEffInt = new TH1D("hEffInt", "Interacting Efficiency", nBins, 0, nBins*50. );
  
  std::vector<TH1D*> mcKE   = {hMCInt,   hMCInc};
  std::vector<TH1D*> recoKE = {hRecoInt, hRecoInc};
  std::vector<TH1D*> effKE  = {hEffInt,  hEffInc};

  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    // int
    if (mcKE[0]->GetBinContent(iBin)!=0)
    {
      auto intEff = recoKE[0]->GetBinContent(iBin)/mcKE[0]->GetBinContent(iBin);
      effKE[0]->SetBinContent(iBin, intEff);
    }
    // inc
    if (mcKE[1]->GetBinContent(iBin)!=0)
    {
      auto incEff = recoKE[1]->GetBinContent(iBin)/mcKE[1]->GetBinContent(iBin);
      effKE[1]->SetBinContent(iBin, incEff);
    }
  }

  TCanvas* c1 = new TCanvas("effInt", "effInt", 800, 800);
  effKE[0]->Draw();
  TCanvas* c2 = new TCanvas("effInc", "effInc", 800, 800);
  effKE[1]->Draw();
}

void makeHists()
{
  anaFile = new TFile("piMinusAna.root", "READ");

  makeXsPlots();
  makeAnglePlots();
  makeSmearingPlots();
  makeEfficiencyPlots();
}