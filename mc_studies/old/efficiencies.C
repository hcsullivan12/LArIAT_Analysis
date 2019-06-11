void efficiencies()
{
  TFile* f1 = new TFile("piMinusAna.root", "READ");
  TFile* f2 = new TFile("../hists/piMinusMCAna.root", "READ");

  TH1D* hMCInt = nullptr;
  TH1D* hMCInc = nullptr;
  TH1D* hRecoInt = nullptr;
  TH1D* hRecoInc = nullptr;

  f1->GetObject("hRecoMCIncidentKE", hRecoInc);
  f1->GetObject("hRecoMCInteractingKE", hRecoInt);
  f2->GetObject("hXsG4IncidentKinEn", hMCInc);
  f2->GetObject("hXsG4InteractingKinEn", hMCInt);

  size_t nBins = 20;

  TH1D* hIntEff = new TH1D("hIntEff", "Interacting Eff", nBins, 0, 1000);
  TH1D* hIncEff = new TH1D("hIncEff", "Incident Eff", nBins, 0, 1000);

  if (!hRecoInt || !hRecoInc) {cout << "Nope1\n"; return;}
  if (!hMCInt || !hMCInc) {cout << "NOPE2\n"; return;}  

  std::vector<std::pair<double, double>> h0, h1, h2, h3;
  std::vector<TH1D*> hists;
  hists.push_back(hRecoInc);
  hists.push_back(hRecoInt);
  hists.push_back(hMCInc);
  hists.push_back(hMCInt);

  for (size_t iBin = 1; iBin <= hists[0]->GetXaxis()->GetNbins(); iBin++)
  {
    if (hists[0]->GetBinCenter(iBin) < 0) continue;
    if (hists[0]->GetBinCenter(iBin) > 1000) continue;
    h0.emplace_back(hists[0]->GetBinContent(iBin), hists[0]->GetBinError(iBin));
  }
  for (size_t iBin = 1; iBin <= hists[1]->GetXaxis()->GetNbins(); iBin++)
  {
    if (hists[1]->GetBinCenter(iBin) < 0) continue;
    if (hists[1]->GetBinCenter(iBin) > 1000) continue;
    h1.emplace_back(hists[1]->GetBinContent(iBin), hists[1]->GetBinError(iBin));
  }
  for (size_t iBin = 1; iBin <= hists[2]->GetXaxis()->GetNbins(); iBin++)
  {
    if (hists[2]->GetBinCenter(iBin) < 0) continue;
    if (hists[2]->GetBinCenter(iBin) > 1000) continue;
    h2.emplace_back(hists[2]->GetBinContent(iBin), hists[2]->GetBinError(iBin));
  }
  for (size_t iBin = 1; iBin <= hists[3]->GetXaxis()->GetNbins(); iBin++)
  {
    if (hists[3]->GetBinCenter(iBin) < 0) continue;
    if (hists[3]->GetBinCenter(iBin) > 1000) continue;
    h3.emplace_back(hists[3]->GetBinContent(iBin), hists[3]->GetBinError(iBin));
  }

  //cout << h0.size() << " " << h1.size() << " " << h2.size() << " " << h3.size() << endl;
  for (size_t i = 0; i<h0.size(); i++)
  {
    hIncEff->SetBinContent(i+1,h0[i].first/h2[i].first);
  }
  hIncEff->Draw();
  for (size_t i = 0; i<h1.size(); i++)
  {
    hIntEff->SetBinContent(i+1,h1[i].first/h3[i].first);
  }
  hIntEff->Draw();

  for (size_t i = 0; i < hists.size(); i++)
  {
    TCanvas *c = new TCanvas(std::to_string(i).c_str(), "temp", 800, 800);
    hists[i]->Draw();
  }

}