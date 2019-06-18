TFile* anaFile;

// G4 prediction: histogram or graph form
TH1D* hG4XS = nullptr;
TGraph* g = nullptr;

// Pre WC to TPC matching
TH1D* hMCPreWCtoTPCInc = nullptr;
TH1D* hMCPreWCtoTPCInt = nullptr;
TH1D* hMCPreWCtoTPCXS = nullptr;

// Post WC to TPC matching
TH1D* hMCXS = nullptr;

// Post WC to TPC and after cuts
TH1D* hMCRawXS = nullptr;
TH1D* hMCInc = nullptr;
TH1D* hMCInt = nullptr;
TH1D* hRecoRawXS = nullptr;
TH1D* hRecoInc = nullptr;
TH1D* hRecoInt = nullptr;
TH1D* hEffInc  = nullptr;
TH1D* hEffInt  = nullptr;

// Smearing
TH2D* hS = nullptr;
TH1D* hRecoUnsmInc = new TH1D("hRecoUnsmInc", "hRecoUnsmInc", 22, 0, 22*50);
TH1D* hRecoUnsmInt = new TH1D("hRecoUnsmInt", "hRecoUnsmInt", 22, 0, 22*50);

// Backgrounds
TH1D* hRecoVertexBkg = nullptr;
TH1D* hRecoTypeBkg = nullptr;
TH1D* hBkgSubInt = new TH1D("hBkgSubInt", "hBkgSubInt", 22, 0, 22*50);
TH1D* hBkgSubInc = new TH1D("hBkgSubInc", "hBkgSubInc", 22, 0, 22*50);

void makeMCPreWCtoTPCXsPlots()
{
  gStyle->SetOptStat(0);

  // make a histogram of g4 xs
  auto nBins = hMCPreWCtoTPCXS->GetXaxis()->GetNbins();
  hG4XS = new TH1D("hG4XS", "G4 Prediction", nBins, 0, 50*nBins);
  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    auto center = hMCPreWCtoTPCXS->GetBinCenter(iBin);
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
  
  //hMCPreWCtoTPCXS->SetMinimum(0);
  //hMCPreWCtoTPCXS->SetMaximum(1100);
  //hMCPreWCtoTPCXS->GetXaxis()->SetRangeUser(0,1000);
  g->GetXaxis()->SetTitle("Kinetic energy [MeV]");
  g->GetYaxis()->SetTitle("#sigma [barn]");
  //hG4XS->Draw("][");
  g->Draw("A c");
  hMCPreWCtoTPCXS->Draw("same");

  hMCPreWCtoTPCXS->SetLineWidth(1);
  hMCPreWCtoTPCXS->SetLineColor(kAzure+4);
  hMCPreWCtoTPCXS->SetMarkerStyle(21);
  hMCPreWCtoTPCXS->SetMarkerColor(kAzure+4);
  hMCPreWCtoTPCXS->SetMarkerColor(kAzure+4);
  g->SetLineWidth(4);
  g->SetLineStyle(1);
  g->SetLineColor(kAzure+4);

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(g, "G4Prediction Inelastic XS", "l");
  legend->AddEntry(hMCPreWCtoTPCXS,"True Inelastic XS","lp");
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

  hEla->GetXaxis()->SetTitle("Angle [degrees]");
  hEla->SetTitle("Angle Between Primary and Secondary");
  hEla->Draw();
  hInelNDa->Draw("same");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(hEla,"Elastic (primary/primary)","l");
  legend->AddEntry(hInelNDa,"Inelastic (primary/one visible secondary)","l");
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
  for(size_t iBinX = 1; iBinX <= hS->GetXaxis()->GetNbins(); iBinX++)
  {
    for(size_t iBinY = 1; iBinY <= hS->GetYaxis()->GetNbins(); iBinY++)
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

  auto nBins = hMCPreWCtoTPCInc->GetXaxis()->GetNbins();
  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    auto recoIncCounts = hBkgSubInc->GetBinContent(iBin);
    auto recoIntCounts = hBkgSubInt->GetBinContent(iBin);
    for (size_t iTrueBin = 1; iTrueBin <= nBins; iTrueBin++)
    {
      auto weight = hS->GetBinContent(iBin, iTrueBin);
      
      auto newIntCounts  = weight * recoIntCounts;
      auto origIntCounts = hRecoUnsmInt->GetBinContent(iTrueBin);
      hRecoUnsmInt->SetBinContent(iTrueBin, origIntCounts+newIntCounts);

      auto newIncCounts  = weight * recoIncCounts;
      auto origIncCounts = hRecoUnsmInc->GetBinContent(iTrueBin);
      hRecoUnsmInc->SetBinContent(iTrueBin, origIncCounts+newIncCounts);
    }
  }

  TCanvas* c2 = new TCanvas("unsmearinginc", "unsmearinginc", 800,800);
  hRecoUnsmInc->Draw();
  TCanvas* c3 = new TCanvas("unsmearingint", "unsmearingint", 800,800);
  hRecoUnsmInt->Draw();

}

void makeEfficiencyPlots()
{
  std::vector<TH1D*> mcRawHists = {nullptr, nullptr};
  anaFile->GetObject("mc/hMCInteractingKE", mcRawHists[0]);
  anaFile->GetObject("mc/hMCIncidentKE", mcRawHists[1]);
  if (!mcRawHists[0] || !mcRawHists[1]) {cout<<"Nope\n"; return;}
  
  TCanvas *c3 = new TCanvas("rawxs", "c3", 800, 800);
  hRecoRawXS->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hRecoRawXS->GetYaxis()->SetTitle("Cross Section (MeV)");
  //hG4XS->Draw();
  //hMCRawXS->Draw();
  TH1D *htemp = new TH1D("temp", "temp", mcRawHists[0]->GetXaxis()->GetNbins(), 0, 50*mcRawHists[0]->GetXaxis()->GetNbins() );
  for (size_t iBin=1; iBin <=mcRawHists[0]->GetXaxis()->GetNbins(); iBin++)
  {
    if (mcRawHists[1]->GetBinContent(iBin) == 0) continue;
    htemp->SetBinContent(iBin, (1/0.0047) * (1/2.1029554659845297e+28) * (1/1e-28) * mcRawHists[0]->GetBinContent(iBin)/mcRawHists[1]->GetBinContent(iBin));
  }
  hRecoRawXS->SetMaximum(0.5);
  hRecoRawXS->SetMinimum(0);
  hRecoRawXS->GetXaxis()->SetRangeUser(0,1000);
  hRecoRawXS->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
  hRecoRawXS->GetYaxis()->SetTitle("#sigma [barns]");
  //htemp->Draw("][");
  hRecoRawXS->Draw("");
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(htemp, "MC raw XS", "l");
  legend->AddEntry(hRecoRawXS,"Reco raw XS","lp");
  //legend->Draw("same");
  
  size_t nBins = hMCPreWCtoTPCInc->GetXaxis()->GetNbins();
  hEffInc = new TH1D("hEffInc", "Incident Efficiency", nBins, 0, nBins*50. );
  hEffInt = new TH1D("hEffInt", "Interacting Efficiency", nBins, 0, nBins*50. );
  
  std::vector<TH1D*> mcPreWCtoTPCKE = {hMCPreWCtoTPCInt,   hMCPreWCtoTPCInc};
  std::vector<TH1D*> recoKE = {hRecoUnsmInt, hRecoUnsmInc};
  std::vector<TH1D*> effKE  = {hEffInt,  hEffInc};

  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    // int
    if (mcPreWCtoTPCKE[0]->GetBinContent(iBin)!=0)
    {
      auto num = recoKE[0]->GetBinContent(iBin);
      auto den = mcPreWCtoTPCKE[0]->GetBinContent(iBin);
      auto intEff = num/den;
      double error = intEff * ( std::sqrt(num)/num + std::sqrt(den)/den );
      effKE[0]->SetBinContent(iBin, intEff);
      effKE[0]->SetBinError(iBin, error);
    }
    // inc
    if (mcPreWCtoTPCKE[1]->GetBinContent(iBin)!=0)
    {
      auto num = recoKE[1]->GetBinContent(iBin);
      auto den = mcPreWCtoTPCKE[1]->GetBinContent(iBin);
      auto incEff = num/den;
      double error = incEff * ( std::sqrt(num)/num + std::sqrt(den)/den );
      effKE[1]->SetBinContent(iBin, incEff);
      effKE[1]->SetBinError(iBin, error);
    }
  }

  effKE[0]->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
  effKE[0]->GetYaxis()->SetTitle("#epsilon_{int}");
  effKE[1]->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
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
    double intEff = hEffInt->GetBinContent(iBin);
    double intErr = hEffInt->GetBinError(iBin);
    double incEff = hEffInc->GetBinContent(iBin);
    double incErr = hEffInc->GetBinError(iBin);
    
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
  g->Draw("same");
  //hG4XS->Draw("same");
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(g, "G4Prediction Inelastic XS", "l");
  legend->AddEntry(hEffCorrXS,"Eff. corrected Reco XS","lp");
  legend->Draw("same");
}

void makeBkgSubtraction()
{
  // add background hists
  auto nBins = hRecoVertexBkg->GetXaxis()->GetNbins();
  for (size_t iBin = 1; iBin <= nBins; iBin++)
  {
    auto term1 = hRecoVertexBkg->GetBinContent(iBin);
    auto term2 = hRecoTypeBkg->GetBinContent(iBin);

    auto rawContentInt = hRecoInt->GetBinContent(iBin);
    auto rawContentInc = hRecoInc->GetBinContent(iBin);
    hBkgSubInt->SetBinContent(iBin, rawContentInt-term1-term2);
    hBkgSubInc->SetBinContent(iBin, rawContentInc-term1-term2);
  }

  TCanvas *c1 = new TCanvas("backsubint", "backsubint", 800, 800);
  hBkgSubInt->Draw();
  TCanvas *c2 = new TCanvas("backsubinc", "backsubinc", 800, 800);
  hBkgSubInc->Draw();
}

void makeOtherPlots()
{
  TCanvas *c1 = new TCanvas("hintcomp", "hintcomp", 800, 800);
  hMCInt->Draw();
  hRecoInt->Draw("e p same");
  hMCInt->SetTitle("Interacting");
  auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry(hRecoInt, "Reco", "lp");
  legend1->AddEntry(hMCInt,"Truth","l");
  legend1->Draw("same");
  hRecoInt->SetMarkerStyle(21);
  hRecoInt->SetMarkerSize(1.5);
  hRecoInt->SetMarkerColor(kAzure+4);
  hRecoInt->SetLineWidth(2);
  hMCInt->SetLineColor(kAzure+4);
  hMCInt->SetLineWidth(4);

  TCanvas *c2 = new TCanvas("hinccomp", "hinccomp", 800, 800);
  hRecoInc->Draw("e p");
  hMCInc->Draw("same");
  hRecoInc->SetTitle("Incident");
  auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry(hRecoInc, "Reco", "lp");
  legend2->AddEntry(hMCInc,"Truth","l");
  legend2->Draw("same");
  hRecoInc->SetMarkerStyle(21);
  hRecoInc->SetMarkerSize(1.5);
  hRecoInc->SetLineWidth(2);
  hRecoInc->SetLineColor(kAzure+4);
  hRecoInc->SetMarkerColor(kAzure+4);
  hMCInc->SetLineColor(kAzure+4);
  hMCInc->SetLineWidth(4);
}

void makeHists()
{
  // get G4 prediction
  anaFile = new TFile("piMinusAna.root", "READ");
  TFile* f2 = new TFile("../hists/g4XsPredictions.root", "READ");
  TCanvas* temp = nullptr;
  anaFile->GetObject("mc/hMCPreWCtoTPCXSKE", hMCPreWCtoTPCXS);
  f2->GetObject("piMinInel", temp);
  g = (TGraph*)temp->GetListOfPrimitives()->FindObject("Graph");
  if (!g) {cout << "No G4 predictions!\n"; return;}

  anaFile->GetObject("reco/hRecoIncidentKE", hRecoInc);
  anaFile->GetObject("reco/hRecoInteractingKE", hRecoInt);
  anaFile->GetObject("reco/hRecoXSKE", hRecoRawXS);
  anaFile->GetObject("reco/hRecoSmearingMatrix", hS);
  anaFile->GetObject("mc/hMCIncidentKE", hMCInc);
  anaFile->GetObject("mc/hMCInteractingKE", hMCInt);
  anaFile->GetObject("mc/hMCXSKE", hMCRawXS);
  anaFile->GetObject("mc/hMCPreWCtoTPCIncidentKE", hMCPreWCtoTPCInc);
  anaFile->GetObject("mc/hMCPreWCtoTPCInteractingKE", hMCPreWCtoTPCInt);
  anaFile->GetObject("mc/hMCPreWCtoTPCXSKE", hMCPreWCtoTPCXS); 
  anaFile->GetObject("reco/hRecoIntVertexBkg", hRecoVertexBkg);
  anaFile->GetObject("reco/hRecoIntTypeBkg", hRecoTypeBkg);
  if (!hRecoInc)         {cout << "Nope1\n"; return;}
  if (!hRecoInt)         {cout << "Nope2\n"; return;}
  if (!hRecoRawXS)       {cout << "Nope3\n"; return;}
  if (!hS)               {cout << "Nope4\n"; return;}
  //if (!hMCRawXS)         {cout << "Nope5\n"; return;}
  if (!hMCPreWCtoTPCInc) {cout << "Nope6\n"; return;}
  if (!hMCPreWCtoTPCInt) {cout << "Nope7\n"; return;}
  if (!hMCPreWCtoTPCXS)  {cout << "Nope8\n"; return;}
  if (!hRecoVertexBkg)   {cout << "Nope9\n"; return;}
  if (!hRecoTypeBkg)     {cout << "Nope10\n"; return;}

  // temp
  //TCanvas* c1 = new TCanvas("c1", "int", 800, 800);
  //hRecoInt->SetTitle("Reconstructed Interacting");
  //hRecoInt->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
  //hRecoInt->Draw();
  //TCanvas* c2 = new TCanvas("c2", "inc", 800, 800);
  //hRecoInc->SetTitle("Reconstructed Incident");
  //hRecoInc->GetXaxis()->SetTitle("Kinetic Energy [MeV]");
  //hRecoInc->Draw();
  //TCanvas* c3 = new TCanvas("c3", "recorawxs", 800, 800);
  //hRecoRawXS->SetLineColor(kAzure+4);
  //hRecoRawXS->SetLineWidth(2);
  //hRecoRawXS->SetMarkerStyle(21);
  //hRecoRawXS->SetMarkerSize(1.5);
  //hRecoRawXS->SetMarkerColor(kAzure+4);
  //hRecoRawXS->SetMinimum(0);
  //hRecoRawXS->SetMaximum(0.5);
  //hRecoRawXS->GetXaxis()->SetRangeUser(0,1000);
  //hRecoRawXS->Draw();

  //makeMCPreWCtoTPCXsPlots();
  //makeAnglePlots();
  // bkg subtract first
  makeBkgSubtraction();
  // then unsmear
  makeSmearingPlots();
  makeEfficiencyPlots();
  //makeRecoXsPlots();
  //makeOtherPlots();
}