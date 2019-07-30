import ROOT
import matplotlib.pyplot as plt

f = ROOT.TFile.Open('anaTree.root')

hKeVsElasticConeAngle    = f.Get('angleana/hMCKeVsElasticConeAngle')
hKeVsInelasticConeAngle  = f.Get('angleana/hMCKeVsInelasticConeAngle')

el_x, el_y = [], []
in_x, in_y = [], []

nXbins = hKeVsElasticConeAngle.GetNbinsX()
nYbins = hKeVsInelasticConeAngle.GetNbinsY()
for xBin in range(1, nXbins):
    for yBin in range(1, nYbins):
        x = hKeVsElasticConeAngle.GetXaxis().GetBinCenter(xBin)
        y = hKeVsElasticConeAngle.GetYaxis().GetBinCenter(yBin)

        elContent = hKeVsElasticConeAngle.GetBinContent(xBin, yBin)
        inContent = hKeVsInelasticConeAngle.GetBinContent(xBin, yBin)

        if elContent != 0:
            el_x.append(x)
            el_y.append(y)
        if inContent != 0:
            in_x.append(x)
            in_y.append(y)

plt.plot(el_x, el_y, c='b', marker='.')
plt.plot(in_x, in_y, c='r', marker='s')
#plt.show()

hKeVsElasticConeAngle.Draw()
hKeVsElasticConeAngle.SetMarkerColor(ROOT.kBlue)
hKeVsElasticConeAngle.SetMarkerSize(0.7)
hKeVsElasticConeAngle.SetMarkerStyle(8)
hKeVsInelasticConeAngle.Draw('same')
hKeVsInelasticConeAngle.SetMarkerColor(ROOT.kOrange+8)
hKeVsInelasticConeAngle.SetMarkerSize(0.7)
hKeVsInelasticConeAngle.SetMarkerStyle(8)

wait = input(' ')