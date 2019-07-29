import ROOT
import matplotlib.pyplot as plt 

f = ROOT.TFile.Open("oldWithSkip.root", "READ")
hElaHist  = f.Get("mc/hMCElasticAngle")
hInelHist = f.Get("mc/hMCInelasticOneVisDAngle")

# renormalize the histograms
nBins = hElaHist.GetNbinsX()
nElaEntries  = hElaHist.GetEntries()
nInelEntries = hInelHist.GetEntries()
nTotEntries  = nElaEntries + nInelEntries
for b in range(1, nBins+1):
    hElaHist.SetBinContent(b, hElaHist.GetBinContent(b)/nTotEntries)
    hInelHist.SetBinContent(b, hInelHist.GetBinContent(b)/nTotEntries)

hElaHist.Draw()
stat1 = hElaHist.GetListOfFunctions().FindObject('stats')
hInelHist.Draw('sames')
stat2 = hInelHist.GetListOfFunctions().FindObject('stats')
hElaHist.GetXaxis().SetTitle('Angle [degrees]')
hElaHist.GetYaxis().SetTitle('Fraction of events')
#ROOT.gPad.Update()
#stat2.Draw('sames')

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
plt.plot(cuts, yEl_1, 'r' , label='Elastic')
plt.plot(cuts, yInel_1, 'b', label='Inelastic')
plt.xlabel('Cut [degrees]', fontsize=25)
plt.ylabel('Passing cut [%]', fontsize=25)
plt.subplot(111).legend(prop={'size': 40})
plt.ylim([0, 1.1])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(linestyle='dotted', linewidth=2)

# efficiency/purity
plt.figure(1)
yInel_purity = [x/(x+y) for x,y in zip(inelIntegrals, elaIntegrals)]
yInel_eff    = [x/nInelEntries for x in inelIntegrals]
plt.plot(cuts, yInel_eff, 'b', label='Efficiency')
plt.plot(cuts, yInel_purity, 'r', label='Purity')
plt.xlabel('Cut [degrees]', fontsize=25)
plt.ylabel('', fontsize=25)
plt.ylim([0, 1.1])
plt.subplot(111).legend(prop={'size': 40})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(linestyle='dotted', linewidth=2)

plt.show()
wait = input('')