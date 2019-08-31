import ROOT

f1 = ROOT.TFile.Open('XS_ANA_5deg_nearestInelastic.root', 'READ')
f2 = ROOT.TFile.Open('XS_ANA_10deg_nearestInelastic.root', 'READ')

h5 = f1.Get('hTruVsRecoLength')
h10 = f2.Get('hTruVsRecoLength')

ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
c5 = ROOT.TCanvas('c5', 'c5', 800, 800)
c5.SetLogz()
fit5 = ROOT.TF1('fit5', '[0]+[1]*x',40,60)
h5.Fit(fit5, 'QR')
h5.Draw('colz')
fit5.Draw('same')
c10 = ROOT.TCanvas('c10', 'c10', 800, 800)
c10.SetLogz()
fit10 = ROOT.TF1('fit10', '[0]+[1]*x',40,60)
h10.Fit(fit10, 'QR')
h10.Draw('colz')
fit10.Draw('same')

wait = input(' ')