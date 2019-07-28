import ROOT
import matplotlib.pyplot as plt 

f = ROOT.TFile.Open("piMinusAna.root", "READ")
hElaHist  = f.Get("mc/hMCElasticAngle")
hInelHist = f.Get("mc/hMCInelasticOneVisDAngle")

# make plot of integral versus cut line
cuts = [c for c in range(0, 90, 1) ]

elaCounts  = hElaHist.GetEntries()
inelCounts = hInelHist.GetEntries()
nBins = hElaHist.GetNbinsX()

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

    elaIntegrals.append(elaInt/(elaCounts+inelCounts))
    inelIntegrals.append(inelInt/(inelCounts+elaCounts))

plt.plot(cuts, elaIntegrals, 'r' , label='Elastic')
plt.plot(cuts, inelIntegrals, 'b', label='Inelastic')
plt.xlabel('Cut [degrees]', fontsize=25)
plt.ylabel('Passing cut [%]', fontsize=25)
plt.title('Events with one visible secondary', fontsize=25)
ax = plt.subplot(111)
ax.legend(prop={'size': 20})
plt.show()

# renormalize the histograms
for b in range(1, nBins+1):
    hElaHist.SetBinContent(b, hElaHist.GetBinContent(b)/(elaCounts+inelCounts))
    hInelHist.SetBinContent(b, hInelHist.GetBinContent(b)/(elaCounts+inelCounts))

hElaHist.Draw()
hInelHist.Draw('same')

wait = input('')