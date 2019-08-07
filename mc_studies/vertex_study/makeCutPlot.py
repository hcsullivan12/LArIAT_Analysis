import ROOT

f = ROOT.TFile.Open('results.root', 'READ')
elas = f.Get('mc/hMCIdeVertexElastic')
inel = f.Get('mc/hMCIdeVertexInelastic')

ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
cEl = ROOT.TCanvas('cEl', 'cEl', 800,800)
elas.Draw('colz')
cIn = ROOT.TCanvas('cIn', 'cIn', 800,800)
inel.Draw('colz')

wait = input(' ')
