import ROOT
import matplotlib.pyplot as plt
from matplotlib  import cm

f = ROOT.TFile.Open('anaTree.root')

hElasticHists = [0 for x in range(0, 6)]
hInelasticHists = [0 for x in range(0, 6)]

hElasticAngle            = f.Get('angleana/hMCElasticAngle')
hElasticConeAngle        = f.Get('angleana/hMCElasticConeAngle')
hKeVsElasticAngle        = f.Get('angleana/hMCKeVsElasticAngle')
hKeVsElasticConeAngle    = f.Get('angleana/hMCKeVsElasticConeAngle')
TrkLenVsElasticAngle     = f.Get('angleana/hMCTrkLenVsElasticAngle')
TrkLenVsElasticConeAngle = f.Get('angleana/hMCTrkLenVsElasticConeAngle')

hInelasticAngle            = f.Get('angleana/hMCInelasticAngle')

hInelasticAngleOneVisD             = f.Get('angleana/hMCInelasticAngleOneVisD')
hInelasticConeAngleOneVisD         = f.Get('angleana/hMCInelasticConeAngleOneVisD')
hKeVsInelasticAngleOneVisD         = f.Get('angleana/hMCKeVsInelasticAngleOneVisD')
hKeVsInelasticConeAngleOneVisD     = f.Get('angleana/hMCKeVsInelasticConeAngleOneVisD')
hTrkLenVsInelasticAngleOneVisD     = f.Get('angleana/hMCTrkLenVsInelasticAngleOneVisD')
hTrkLenVsInelasticConeAngleOneVisD = f.Get('angleana/hMCTrkLenVsInelasticConeAngleOneVisD')

elEntries = hKeVsElasticConeAngle.GetEntries()
inEntries = hKeVsInelasticConeAngle.GetEntries()

hKeVsInelasticConeAngle.SetMarkerColor(ROOT.kOrange+8)
hKeVsInelasticConeAngle.SetMarkerSize(0.7)
hKeVsInelasticConeAngle.SetMarkerStyle(8)
hKeVsInelasticConeAngle.GetXaxis().SetTitle('Angle [degrees]')
hKeVsInelasticConeAngle.GetYaxis().SetTitle('KE[MeV]')
hKeVsElasticConeAngle.SetMarkerColor(ROOT.kBlue)
hKeVsElasticConeAngle.SetMarkerSize(0.7)
hKeVsElasticConeAngle.SetMarkerStyle(8)
hKeVsElasticConeAngle.GetXaxis().SetTitle('Angle [degrees]')
hKeVsElasticConeAngle.GetYaxis().SetTitle('KE[MeV]')


ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)

cEl = ROOT.TCanvas('cEl', 'Elastic', 800, 800)
hKeVsElasticConeAngle.Draw('colz')
fit = ROOT.TF1('fitToElasticEdge', '[0] + [1]*TMath::Exp(-[2]*x)', 0, 25)
fit.SetParameters(80, 2500, 0.25)
fit.SetLineColor(ROOT.kYellow)
fit.SetLineWidth(7)
fit.Draw('same')

cIn = ROOT.TCanvas('cIn', 'Inelastic', 800, 800)
hKeVsInelasticConeAngle.Draw('colz')
fit.Draw('same')

#cBoth = ROOT.TCanvas('cBoth', 'Both', 800, 800)
#hKeVsElasticConeAngle.Draw()
#hKeVsInelasticConeAngle.Draw('same')


#ex1 = ROOT.TExec('ex1', 'gStyle->SetPalette(kBlueRedYellow);')
#ex2 = ROOT.TExec('ex2', 'gStyle->SetPalette(53);')
#hKeVsElasticConeAngle.Draw('colz')
#ex1.Draw()
#hKeVsElasticConeAngle.Draw('col same')
#ex2.Draw()
#hKeVsInelasticConeAngle.Draw('col same')

wait = input(' ')