import ROOT

g4f = ROOT.TFile.Open('root_files/g4XsPredictionPiMinusInel.root', 'READ')
mcf = ROOT.TFile.Open('myana_mc_full.root', 'READ')

hG4Xs = g4f.Get('Graph') 
hMcRecoXs = mcf.Get('hCrossSectionKe')


hMcRecoXs.Draw('p')
hMcRecoXs.SetMinimum(0)
hMcRecoXs.SetMaximum(1.2)
hMcRecoXs.SetMarkerStyle(8)
hMcRecoXs.SetMarkerColor(4)
hMcRecoXs.GetXaxis().SetTitle('KE [MeV]')
hMcRecoXs.GetYaxis().SetTitle('Cross Section [barns]')
hG4Xs.Draw('same l')
hG4Xs.SetLineColor(2)
hG4Xs.SetLineWidth(2)

leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(hMcRecoXs, 'MC Reco (only pions)', 'pl')
leg.AddEntry(hG4Xs, 'MC Thin Target', 'l')
leg.Draw('same')


wait = input(' ')