import ROOT
import matplotlib.pyplot as plt
from matplotlib  import cm

f = ROOT.TFile.Open('anaTree.root')

##############################################################
el_hists = [0 for x in range(0, 6)]
in_hists = [0 for x in range(0, 6)]

el_hists[0] = f.Get('angleana/hMCElasticAngle')
el_hists[1] = f.Get('angleana/hMCElasticConeAngle')
el_hists[2] = f.Get('angleana/hMCKeVsElasticAngle')
el_hists[3] = f.Get('angleana/hMCKeVsElasticConeAngle')
el_hists[4] = f.Get('angleana/hMCTrkLenVsElasticAngle')
el_hists[5] = f.Get('angleana/hMCTrkLenVsElasticConeAngle')

hInelasticAngle = f.Get('angleana/hMCInelasticAngle')

in_hists[0] = f.Get('angleana/hMCInelasticAngleOneVisD')
in_hists[1] = f.Get('angleana/hMCInelasticConeAngleOneVisD')
in_hists[2] = f.Get('angleana/hMCKeVsInelasticAngleOneVisD')
in_hists[3] = f.Get('angleana/hMCKeVsInelasticConeAngleOneVisD')
in_hists[4] = f.Get('angleana/hMCTrkLenVsInelasticAngleOneVisD')
in_hists[5] = f.Get('angleana/hMCTrkLenVsInelasticConeAngleOneVisD')
##############################################################

# draw elastic and inelastic angles
#cv1 = ROOT.TCanvas('cv1', 'Elastic', 800, 800)
#cv1.SetLogy()
#el_hists[0].Draw()
el_hists[0].GetXaxis().SetTitle('Angle [degrees]')
el_hists[0].GetYaxis().SetTitle('Entries/1 degree')
el_hists[0].SetLineWidth(4)
el_hists[0].SetLineColor(ROOT.kBlue+1)
#in_hists[0].Draw('sames')
in_hists[0].SetLineWidth(4)
in_hists[0].SetLineColor(ROOT.kOrange+7)

in_hists[2].SetMarkerColor(ROOT.kOrange+8)
in_hists[2].SetMarkerSize(0.7)
in_hists[2].SetMarkerStyle(8)
in_hists[2].GetXaxis().SetTitle('Angle [degrees]')
in_hists[2].GetYaxis().SetTitle('KE[MeV]')
el_hists[2].SetMarkerColor(ROOT.kBlue)
el_hists[2].SetMarkerSize(0.7)
el_hists[2].SetMarkerStyle(8)
el_hists[2].GetXaxis().SetTitle('Angle [degrees]')
el_hists[2].GetYaxis().SetTitle('KE[MeV]')


ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)

cEl = ROOT.TCanvas('cEl', 'Elastic', 800, 800)
el_hists[2].Draw('colz')
fit = ROOT.TF1('fitToElasticEdge', '[0] + [1]*TMath::Exp(-[2]*x)', 0, 20)
fit2 = ROOT.TF1('fit2', '[0]+[1]*x', 19, 40)
#fit.SetParameters(80, 2500, 0.25)
fit.SetParameters(200, 4000, 0.25)
fit2.SetParameters(350, -1*6)
#fit = ROOT.TF1('fitToElasticEdge', '[0] + [1]*TMath::Power([2],[3]*x)', 0, 30)
#fit.SetParameters(0,1000,2,-0.01)
fit.SetLineColor(ROOT.kYellow)
fit.SetLineWidth(7)
fit2.SetLineColor(ROOT.kYellow)
fit2.SetLineWidth(7)
fit.Draw('same')
fit2.Draw('same')

cIn = ROOT.TCanvas('cIn', 'Inelastic', 800, 800)
in_hists[2].Draw('colz')
#fit.Draw('same')
fit2.Draw('same')

nIn, nEl, nEl_in_In = 0., 0., 0.
nInEntries = in_hists[2].GetEntries()
nElEntries = el_hists[2].GetEntries()
for xBin in range(0, el_hists[2].GetNbinsX()):
    for yBin in range(0, el_hists[2].GetNbinsY()):
        x = el_hists[2].GetXaxis().GetBinCenter(xBin)
        y = el_hists[2].GetYaxis().GetBinCenter(yBin)
        
        above = False
        if x > 30: 
            above = True

        if above: 
            nIn += in_hists[2].GetBinContent(xBin, yBin)
            nEl_in_In += el_hists[2].GetBinContent(xBin, yBin)
        else: nEl += el_hists[2].GetBinContent(xBin, yBin)

print 'Inelastic efficiency:', nIn/nInEntries
print 'Inelastic purity:    ', nIn/(nIn+nEl_in_In)
wait = input(' ')