void combAngles()
{
  TFile* f = new TFile("piMinusAna.root", "READ");
  gStyle->SetOptStat(0);
  TH1D* hEla     = nullptr;
  TH1D* hInelNDa = nullptr;

  f->GetObject("hElasticAngle", hEla);
  f->GetObject("hInelasticOneVisDAngle", hInelNDa);

  if (!hEla || !hInelNDa) {cout << "Nope\n"; return;}

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

  hEla->GetXaxis()->SetTitle("Angle (degrees)");
  hEla->SetTitle("Angle Between Primary and Secondary");
  hEla->Draw();
  hInelNDa->Draw("same");

   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(hEla,"Elastic (primary/primary)","l");
   legend->AddEntry(hInelNDa,"Inelastic (primary/one visible daughter)","l");
   legend->Draw("same");
   legend->SetLineColor(0);
}