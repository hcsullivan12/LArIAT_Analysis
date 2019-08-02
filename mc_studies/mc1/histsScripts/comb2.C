void comb2()
{
  gStyle->SetOptStat(0);

  TFile *f1 = new TFile("temp.root", "READ");
  TFile *f2 = new TFile("../hists/g4XsPredictions.root", "READ");

  std::vector<TH1D*> hists(5, nullptr);

  //TH1D *hG4Pred = nullptr;
  TCanvas *temp = nullptr;
  
  f1->GetObject("anatree/hCrossSectionInel",       hists[0]);
  f1->GetObject("anatree/hCrossSectionPionAbsorp", hists[1]);
  f1->GetObject("anatree/hCrossSectionPionInel",   hists[2]);
  f1->GetObject("anatree/hCrossSectionPionProd",   hists[3]);
  f1->GetObject("anatree/hCrossSectionChargeExch", hists[4]);

  f2->GetObject("piMinInel", temp);

  TGraph *g = (TGraph*)temp->GetListOfPrimitives()->FindObject("Graph");

  // rescale 
  for (size_t b = 1; b <= hists[0]->GetXaxis()->GetNbins(); b++)
  {
    hists[0]->SetBinContent(b, 1000*hists[0]->GetBinContent(b)/(100.));
    hists[1]->SetBinContent(b, 1000*hists[1]->GetBinContent(b)/(100.));
    hists[2]->SetBinContent(b, 1000*hists[2]->GetBinContent(b)/(100.));
    hists[3]->SetBinContent(b, 1000*hists[3]->GetBinContent(b)/(100.));
    hists[4]->SetBinContent(b, 1000*hists[4]->GetBinContent(b)/(100.));

    hists[0]->SetBinError(b, 1000*hists[0]->GetBinError(b)/(100.));
    hists[1]->SetBinError(b, 1000*hists[1]->GetBinError(b)/(100.));
    hists[2]->SetBinError(b, 1000*hists[2]->GetBinError(b)/(100.));
    hists[3]->SetBinError(b, 1000*hists[3]->GetBinError(b)/(100.));
    hists[4]->SetBinError(b, 1000*hists[4]->GetBinError(b)/(100.));
  }
  for (int i=0;i<g->GetN();i++) g->GetY()[i] *= 1000;

  // don't know why I have to do  this
  TH1D* tempH = new TH1D("tempH", "tempH", 20, 0, 1000);
  for (size_t b = 1; b <= hists[0]->GetXaxis()->GetNbins(); b++)
  {
    tempH->SetBinContent(b, hists[0]->GetBinContent(b));
    tempH->SetBinError(b, hists[0]->GetBinError(b));
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 800., 800.);
  //c1->SetGridy();
  //c1->SetGridx();
  
  hists[0]->SetMinimum(1);
  hists[0]->SetMaximum(1100);
  hists[0]->GetXaxis()->SetRangeUser(0,1000);
  hists[0]->GetXaxis()->SetTitle("Kinetic energy [MeV]");
  hists[0]->GetYaxis()->SetTitle("#sigma [mb]");
  c1->SetLogy();
  hists[0]->Draw();
  hists[1]->Draw("same");
  hists[2]->Draw("same");
  hists[3]->Draw("same");
  hists[4]->Draw("same");
  g->Draw("same c");

  hists[0]->SetLineWidth(1);
  hists[0]->SetLineColor(kRed);
  hists[0]->SetMarkerStyle(21);
  hists[0]->SetMarkerColor(kRed);

  g->SetLineWidth(3);
  g->SetLineStyle(9);
  g->SetLineColor(kRed);

  hists[1]->SetLineWidth(1);
  hists[1]->SetLineColor(kGreen+3);
  hists[1]->SetMarkerStyle(21);
  hists[1]->SetMarkerColor(kGreen+3);

  hists[2]->SetLineWidth(1);
  hists[2]->SetLineColor(kBlue);
  hists[2]->SetMarkerStyle(21);
  hists[2]->SetMarkerColor(kBlue);

  hists[3]->SetLineWidth(1);
  hists[3]->SetLineColor(kMagenta);
  hists[3]->SetMarkerStyle(21);
  hists[3]->SetMarkerColor(kMagenta);

  hists[4]->SetLineWidth(1);
  hists[4]->SetLineColor(kAzure-8);
  hists[4]->SetMarkerStyle(21);
  hists[4]->SetMarkerColor(kAzure-8);


  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(g, "G4Prediction Inelastic XS", "l");
  legend->AddEntry(hists[0],"Total inelastic XS","lp");
  legend->AddEntry(hists[1],"Pion absorption","lp");
  legend->AddEntry(hists[2],"Inelastic", "lp");
  legend->AddEntry(hists[3],"Pion production","lp");
  legend->AddEntry(hists[4],"Charge exchange","lp");
  legend->Draw("same");
}
