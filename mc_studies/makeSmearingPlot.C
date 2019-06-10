void makeSmearingPlot()
{
  TFile* f = new TFile("piMinusAna.root", "READ");

  TH2D* hS = nullptr;

  f->GetObject("hSmearingMatrix", hS);
  if (!hS) {cout << "Nope\n"; return;}

  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetPalette(kBlueRedYellow);
  hS->Draw("colz text");
}