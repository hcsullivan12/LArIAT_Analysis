import ROOT
import matplotlib.pyplot as plt 

f = ROOT.TFile.Open("piMinusAna.root", "READ")
hElaHist  = f.Get("mc/hMCElasticAngle")
hInelHist = f.Get("mc/hMCInelasticOneVisDAngle")

# make plot of integral versus cut line
cuts = [c for c in range(1, 90, 1) ]

elaCounts = 0
inelCounts = 0
nBins = hElaHist.GetNbinsX()
for b in range(1, nBins+1):
    elaCounts  += hElaHist.GetBinContent(b)
    inelCounts += hInelHist.GetBinContent(b)

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