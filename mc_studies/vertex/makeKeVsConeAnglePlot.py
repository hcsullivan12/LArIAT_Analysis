import ROOT
import matplotlib.pyplot as plt
from matplotlib  import cm

f = ROOT.TFile.Open('anaTree.root')

hKeVsElasticConeAngle    = f.Get('angleana/hMCKeVsElasticConeAngle')
hKeVsInelasticConeAngle  = f.Get('angleana/hMCKeVsInelasticConeAngle')

el_x, el_y, el_z = [], [], []
in_x, in_y, in_z = [], [], []

nXbins = hKeVsElasticConeAngle.GetNbinsX()
nYbins = hKeVsInelasticConeAngle.GetNbinsY()
for xBin in range(1, nXbins+1):
    for yBin in range(1, nYbins+1):
        x = hKeVsElasticConeAngle.GetXaxis().GetBinCenter(xBin)
        y = hKeVsElasticConeAngle.GetYaxis().GetBinCenter(yBin)

        elContent = hKeVsElasticConeAngle.GetBinContent(xBin, yBin)
        inContent = hKeVsInelasticConeAngle.GetBinContent(xBin, yBin)

        if elContent != 0:
            el_x.append(x)
            el_y.append(y)
            el_z.append(elContent)
        if inContent != 0:
            in_x.append(x)
            in_y.append(y)
            in_z.append(inContent)

elEntries = hKeVsElasticConeAngle.GetEntries()
inEntries = hKeVsInelasticConeAngle.GetEntries()
#hKeVsElasticConeAngle.Scale(1/elEntries)
#hKeVsInelasticConeAngle.Scale(1/inEntries)


#for xBin in range(1, nXbins+1):
#    for yBin in range(1, nYbins+1):
#        hKeVsElasticConeAngle.SetBinContent(xBin, yBin, hKeVsElasticConeAngle.GetBinContent(xBin, yBin)/elEntries)
#        hKeVsInelasticConeAngle.SetBinContent(xBin, yBin, hKeVsInelasticConeAngle.GetBinContent(xBin, yBin)/inEntries)

#plt.scatter(in_x, in_y, c=in_z, s=20, marker='o', cmap=cm.hot)
#plt.scatter(el_x, el_y, c=el_z, s=20, marker='o', cmap=cm.jet)
#plt.hist2d(in_x, in_y, in_z, cmap=plt.cm.Reds)
#plt.show()

hKeVsInelasticConeAngle.SetMarkerColor(ROOT.kOrange+8)
hKeVsInelasticConeAngle.SetMarkerSize(0.7)
hKeVsInelasticConeAngle.SetMarkerStyle(8)
hKeVsElasticConeAngle.SetMarkerColor(ROOT.kBlue)
hKeVsElasticConeAngle.SetMarkerSize(0.7)
hKeVsElasticConeAngle.SetMarkerStyle(8)



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