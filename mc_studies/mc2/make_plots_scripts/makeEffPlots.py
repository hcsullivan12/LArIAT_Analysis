import ROOT

f1 = ROOT.TFile.Open('XS_ANA.root', 'READ')
f2 = ROOT.TFile.Open('XS_ANA_TRUE_WELL_RECO.root', 'READ')

hInc = f1.Get('hTrueIncidentKe')
hInt = f1.Get('hTrueInteractingKe')

hIncWlR = f2.Get('hTrueIncidentKe')
hIntWlR = f2.Get('hTrueInteractingKe')

hIncEff = hIncWlR #ROOT.TH1D('hIncEff', 'hIncEff', 24, 0, 1200)
hIntEff = hIntWlR #ROOT.TH1D('hIntEff', 'hIntEff', 24, 0, 1200)

hIncEff.Divide(hInc)
hIntEff.Divide(hInt)

hIncEff.SetLineColor(1)
hIntEff.SetLineColor(2)
hIncEff.SetLineWidth(4)
hIntEff.SetLineWidth(4)

leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(hIncEff, 'Incident', 'l')
leg.AddEntry(hIntEff, 'Interacting', 'l')

c = ROOT.TCanvas('c','c',800,800)
hIncEff.Draw()
hIntEff.Draw('same')
leg.Draw('same')



wait = input(' ')