import ROOT
import matplotlib.pyplot as plt 

f = ROOT.TFile.Open("anaTree.root", "READ")
hElaHist  = f.Get("angleana/hMCElasticAngle")
hInelHist = f.Get("angleana/hMCInelasticAngleOneVisD")

# renormalize the histograms
nBins = hElaHist.GetNbinsX()
nElaEntries  = hElaHist.GetEntries()
nInelEntries = hInelHist.GetEntries()
nTotEntries  = nElaEntries + nInelEntries
hElaHist.Scale(1/nTotEntries)
hInelHist.Scale(1/nTotEntries)

cv1 = ROOT.TCanvas('cv1', 'Angles combined', 800, 800)
cv1.SetLogy()
hElaHist.Draw()
hInelHist.Draw('sames')
hElaHist.GetXaxis().SetTitle('Angle [degrees]')
hElaHist.GetYaxis().SetTitle('Fraction of events')

# make plot of integral versus cut line
cuts = [c for c in range(0, 30, 1) ]

elaIntegrals = []
inelIntegrals = []
for c in cuts:
    iBin = hElaHist.GetXaxis().FindBin(c)
    jBin = hInelHist.GetXaxis().FindBin(c)
    assert(iBin == jBin)

    elaInt = 0
    inelInt = 0
    for b in range(iBin, nBins+1):
        elaInt  += hElaHist.GetBinContent(b)
        inelInt += hInelHist.GetBinContent(b)

    elaIntegrals.append(elaInt*nTotEntries)
    inelIntegrals.append(inelInt*nTotEntries)

# percentage of elastic like events passing cuts
plt.figure(0)
yEl_1   = [x/nElaEntries for x in elaIntegrals]
yInel_1 = [x/nInelEntries for x in inelIntegrals]
plt.plot(cuts, yEl_1,   c='mediumblue', linewidth=3, label='Elastic eff.')
plt.plot(cuts, yInel_1, c='darkorange', linewidth=3, label='Inelastic eff.')

plt.xlabel('Cut [degrees]', fontsize=25)
#plt.ylabel('Passing cut [%]', fontsize=25)
plt.subplot(111).legend(prop={'size': 20})
plt.ylim([0, 1.1])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(linestyle='dotted', linewidth=2)

# efficiency/purity
plt.figure(1)
yInel_purity = [x/(x+y) for x,y in zip(inelIntegrals, elaIntegrals)]
yInel_eff    = [x/nInelEntries for x in inelIntegrals]
plt.plot(cuts, yInel_eff,    c='darkorange', linewidth=3, linestyle='-', label='Efficiency')
plt.plot(cuts, yInel_purity, c='darkorange', linewidth=3, linestyle='--', label='Purity')
plt.xlabel('Cut [degrees]', fontsize=25)
plt.ylabel('', fontsize=25)
plt.ylim([0, 1.1])
plt.subplot(111).legend(prop={'size': 40})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(linestyle='dotted', linewidth=2)

plt.show()
wait = input('')