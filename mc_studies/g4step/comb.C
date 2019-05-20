void comb()
{
  TFile *f1 = new TFile("piMinusAna.root", "READ");
  TFile *f2 = new TFile("../../hists/g4XsPredictions.root", "READ");

  TH1D *hInelXs = nullptr;
  TH1D *hG4Pred = nullptr;
  TCanvas *temp = nullptr;
  
  f1->GetObject("hXsG4", hInelXs);
  f2->GetObject("piMinInel", temp);

  TGraph *g = (TGraph*)temp->GetListOfPrimitives()->FindObject("Graph");

  if (!g) cout << "Error1\n";
  if (!hInelXs) cout << "Error2\n";

  TCanvas *c1 = new TCanvas("c1", "c1", 800., 800);
  c1->SetGridy();
  c1->SetGridx();
  g->SetMinimum(0);
  g->SetMaximum(1.2);
  g->GetXaxis()->SetTitle("Kinetic energy (MeV)");
  g->GetYaxis()->SetTitle("Cross section (barn)");
  g->SetLineColor(kBlue-3);
  g->SetLineWidth(3);
  g->Draw("AC");
  hInelXs->SetMarkerSize(1.2);
  hInelXs->SetMarkerStyle(8);
  hInelXs->SetMarkerColor(kBlue-3);
  hInelXs->SetLineColor(kBlue-3);
  hInelXs->Draw("p same e1");
}
