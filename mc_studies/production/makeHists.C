TFile* anaFile;
TH1D* hMCInc = nullptr;
TH1D* hMCInt = nullptr;
TH1D* hRecoInc = nullptr;
TH1D* hRecoInt = nullptr;
TH1D* hEffInc  = nullptr;
TH1D* hEffInt  = nullptr;
TH1D* hG4XS    = nullptr;
TH2D* hS = nullptr;
TH1D* hRecoRawXS = nullptr;

void makeMCXsPlots()
{
  gStyle->SetOptStat(0);

  TFile *f2 = new TFile("../hists/g4XsPredictions.root", "READ");

  TH1D *hMCXS = nullptr;
  TCanvas *temp = nullptr;
  
  anaFile->GetObject("mc/hMCXSKE", hMCXS);
  f2->GetObject("piMinInel", temp);
  TGraph *g = (TGraph*)temp->GetListOfPrimitives()->FindObject("Graph");

  if (!g || !temp || !hMCXS) {cout << "Nope\n"; return;}

  // make a histogram of g4 xs
  auto nBins = hMCXS->GetXaxis()->GetNbins();
  hG4XS = new TH1D("hG4XS", "G4 Prediction", nBins, 0, 50*nBins);
  
  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    auto center = hMCXS->GetBinCenter(iBin);
    double bounds[2] = {center-25, center+25};
    double sum = 0;
    double count = 0;
    
    for (size_t iPt = 1; iPt <= g->GetN(); iPt++)
    {
      double x = 0;
      double y = 0;
      g->GetPoint(iPt, x, y);

      if (bounds[0] <= x && x < bounds[1]) {sum = sum + y; count++;}
    }
    hG4XS->SetBinContent(iBin, sum/count);
  }

  TCanvas *c1 = new TCanvas("xsplots", "c1", 800., 800.);
  
  //hMCXS->SetMinimum(0);
  //hMCXS->SetMaximum(1100);
  //hMCXS->GetXaxis()->SetRangeUser(0,1000);
  hG4XS->GetXaxis()->SetTitle("Kinetic energy (MeV)");
  hG4XS->GetYaxis()->SetTitle("Cross section (barn)");
  hG4XS->Draw("");
  //g->Draw("A c");
  hMCXS->Draw("same");

  hMCXS->SetLineWidth(3);
  hMCXS->SetLineColor(kRed);
  hMCXS->SetMarkerStyle(8);
  hMCXS->SetMarkerColor(kRed);

  g->SetLineWidth(3);
  g->SetLineStyle(9);
  g->SetLineColor(kRed);

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(hG4XS, "G4Prediction Inelastic XS", "l");
  legend->AddEntry(hMCXS,"True Inelastic XS","lp");
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
  anaFile->GetObject("reco/hSmearingMatrix", hS);
  if (!hS) {cout << "Nope\n"; return;}

  TCanvas *c1 = new TCanvas("smearing", "c1", 800, 800);
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetPalette(kBlueRedYellow);
  hS->Draw("colz text");

  hS->GetXaxis()->SetTitle("Reconstructed KE (MeV)");
  hS->GetYaxis()->SetTitle("True KE (MeV)");

  // get avg smearing along diagonal
  float sum1 = 0;
  size_t count1 = 0;
  float sum2 = 0;
  size_t count2 = 0;
  for(size_t iBinX = 0; iBinX < hS->GetXaxis()->GetNbins(); iBinX++)
  {
    for(size_t iBinY = 0; iBinY < hS->GetYaxis()->GetNbins(); iBinY++)
    {
      if (iBinX != iBinY) continue;
      auto content = hS->GetBinContent(iBinX, iBinY);
      sum1 = sum1 + content;
      count1++;

      auto center = hS->GetXaxis()->GetBinCenter(iBinX);
      if (center >= 500) {sum2 = sum2+content; count2++;}
    }
  }
  cout << "Average smearing along diagonal = " << sum1/count1 << endl;
  cout << "Average about 500 MeV = " << sum2/count2 << endl;
}

void makeEfficiencyPlots()
{
  anaFile->GetObject("mc/hMCIncidentKE", hMCInc);
  anaFile->GetObject("mc/hMCInteractingKE", hMCInt);
  anaFile->GetObject("reco/hRecoIncidentKE", hRecoInc);
  anaFile->GetObject("reco/hRecoInteractingKE", hRecoInt);
  anaFile->GetObject("reco/hRecoXSKE", hRecoRawXS);

  if (!hMCInc || !hMCInt || !hRecoInc || !hRecoInt || !hRecoRawXS) {cout<<"Nope\n"; return;}
  
  TCanvas *c3 = new TCanvas("rawxs", "c3", 800, 800);
  hRecoRawXS->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hRecoRawXS->GetYaxis()->SetTitle("Cross Section (MeV)");
  //hG4XS->Draw();
  hRecoRawXS->Draw();

  size_t nBins = hMCInc->GetXaxis()->GetNbins();
  hEffInc = new TH1D("hEffInc", "Incident Efficiency", nBins, 0, nBins*50. );
  hEffInt = new TH1D("hEffInt", "Interacting Efficiency", nBins, 0, nBins*50. );
  
  std::vector<TH1D*> mcKE   = {hMCInt,   hMCInc};
  std::vector<TH1D*> recoKE = {hRecoInt, hRecoInc};
  std::vector<TH1D*> effKE  = {hEffInt,  hEffInc};

  for (size_t iBin = 2; iBin < nBins; iBin++)
  {
    // int
    if (mcKE[0]->GetBinContent(iBin)!=0)
    {
      auto num = recoKE[0]->GetBinContent(iBin);
      auto den = mcKE[0]->GetBinContent(iBin);
      auto intEff = num/den;
      double error = intEff * ( std::sqrt(num)/num + std::sqrt(den)/den );
      effKE[0]->SetBinContent(iBin, intEff);
      effKE[0]->SetBinError(iBin, error);
    }
    // inc
    if (mcKE[1]->GetBinContent(iBin)!=0)
    {
      auto num = recoKE[1]->GetBinContent(iBin);
      auto den = mcKE[1]->GetBinContent(iBin);
      auto incEff = num/den;
      double error = incEff * ( std::sqrt(num)/num + std::sqrt(den)/den );
      effKE[1]->SetBinContent(iBin, incEff);
      effKE[1]->SetBinError(iBin, error);
    }
  }

  effKE[0]->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  effKE[0]->GetYaxis()->SetTitle("#epsilon_{int}");
  effKE[1]->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  effKE[1]->GetYaxis()->SetTitle("#epsilon_{inc}");

  TCanvas* c1 = new TCanvas("effInt", "effInt", 800, 800);
  effKE[0]->Draw();
  TCanvas* c2 = new TCanvas("effInc", "effInc", 800, 800);
  effKE[1]->Draw();
}

void makeRecoXsPlots()
{
  size_t nBins = hRecoInc->GetXaxis()->GetNbins();
  TH1D* hEffCorrXS = new TH1D("hEffCorrXS", "Efficiency Corrected XS", nBins, 0, 50*nBins);

  for (size_t iBin = 2; iBin < nBins; iBin++)
  {
    double intEff     = hEffInt->GetBinContent(iBin);
    double intErr     = hEffInt->GetBinError(iBin);
    double incEff     = hEffInc->GetBinContent(iBin);
    double incErr     = hEffInc->GetBinError(iBin);
    
    //cout << iBin << " " << intContent << " " << incContent << " " << intContent/incContent <<  endl;

    double xs = hRecoRawXS->GetBinContent(iBin);
    double xsErr = hRecoRawXS->GetBinError(iBin);
    
    double effCorrXS = (incEff/intEff) * xs;
    double error = effCorrXS * ( intErr/intEff + incErr/intEff );

    hEffCorrXS->SetBinContent( iBin, effCorrXS );
    hEffCorrXS->SetBinError( iBin, error );
  }

  TCanvas* c1 = new TCanvas("effcorrectedxs", "c1", 800, 800);
  hEffCorrXS->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hEffCorrXS->GetYaxis()->SetTitle("Cross section (MeV)");  
  hEffCorrXS->Draw();
  hG4XS->Draw("same");
}

void makeHists()
{
  anaFile = new TFile("piMinusAna.root", "READ");

  makeMCXsPlots();
  makeAnglePlots();
  makeSmearingPlots();
  makeEfficiencyPlots();
  makeRecoXsPlots();
}