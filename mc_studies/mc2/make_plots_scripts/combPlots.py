import ROOT

f1 = ROOT.TFile.Open('XS_ANA.root', 'READ')

hInel = f1.Get('hInteractingKe')
hEl   = f1.Get('hInteractingKeElasticBkg')

c = ROOT.TCanvas('c','c',1000,1000)

temp = hEl.Divide(hInel)
#hInel.Draw()
#hEl.Draw('sames')
hEl.Draw()

leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(hInel, 'Inelastic', 'pl')
leg.AddEntry(hEl, 'Elastic', 'pl')
leg.Draw('same')

wait = input(' ')