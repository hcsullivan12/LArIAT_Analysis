import ROOT

f1 = ROOT.TFile.Open('plots/g4PredAndFirstAttempt.root', 'READ')
f2 = ROOT.TFile.Open('plots/g4PredAndCurrentBest.root',  'READ')

g4p = f1.Get('Graph')
hold = f1.Get('hCrossSectionKe')
hnew = f2.Get('hCrossSectionKe')

c = ROOT.TCanvas('c','c',800,800)

hold.Draw()
g4p.Draw('same l')
hnew.Draw('same')
hold.SetMarkerStyle(1)
hnew.SetMarkerStyle(1)

leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(g4p,  'G4Prediction Inelastic XS', 'l')
leg.AddEntry(hnew, 'New MC Reco (only pions)', 'pl')
leg.AddEntry(hold, 'MC Reco (only pions)', 'pl')
leg.Draw('same')

wait = input(' ')