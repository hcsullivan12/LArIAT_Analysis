import ROOT
import matplotlib.pyplot as plt 
import matplotlib.colors as colors

f = ROOT.TFile.Open("piMinusAna.root", "READ")
hElaHist  = f.Get("mc/hMCIdeVertexElastic")
hInelHist = f.Get("mc/hMCIdeVertexInelastic")

# normalize xbins
nBins = hElaHist.GetNbinsX()
integrals = [1 for x in range(0, nBins)]
#for xBin in range(0, nBins):
#    for yBin in range(0, hElaHist.GetNbinsY()):
#        integrals[xBin] += hElaHist.GetBinContent(xBin+1, yBin+1) + hInelHist.GetBinContent(xBin+1, yBin+1)

for xBin in range(0, nBins):
    if integrals[xBin] == 0: continue
    for yBin in range(0, hElaHist.GetNbinsY()):
        hElaHist.SetBinContent(xBin+1,  yBin+1, hElaHist.GetBinContent (xBin+1, yBin+1)/integrals[xBin])
        hInelHist.SetBinContent(xBin+1, yBin+1, hInelHist.GetBinContent(xBin+1, yBin+1)/integrals[xBin])

ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
cEla = ROOT.TCanvas("cEla", "Elastic", 1000, 1000)
cEla.SetLogz()
hInelHist.Draw('colz')
cInel = ROOT.TCanvas("cInela", "Elastic", 1000, 1000)
cInel.SetLogz()
hElaHist.Draw('colz')

rcuts = [0.5 * x for x in range(0, 10)]
ecuts = [x for x in range(5, 20, 5)]
COLORS = plt.cm.get_cmap('hsv', len(ecuts))
yEl   = {}
yInel = {}
for e in ecuts:
    yEl[e]   = []
    yInel[e] = []

# integrate for each of the cuts
for r in rcuts:
    xBin = hElaHist.GetXaxis().FindBin(r)

    for e in ecuts:
        elBin   = hElaHist.GetYaxis().FindBin(e)
        inelBin = hInelHist.GetYaxis().FindBin(e)
        s = 0
        for iBin in range(elBin, hElaHist.GetNbinsY()+1):
            s += hElaHist.GetBinContent(xBin, iBin)
        yEl[e].append(s)
        s = 0
        for iBin in range(inelBin, hInelHist.GetNbinsY()+1):
            s += hInelHist.GetBinContent(xBin, iBin)
        yInel[e].append(s)

# plot the results for each cut
counter = [x for x in range(1, 5)]
for e, c in zip(ecuts, counter):
    print yEl[e]
    plt.plot(rcuts, yEl[e],   color=COLORS(c), label='Elastic_'+str(e)+'cut')
    plt.plot(rcuts, yInel[e], color=COLORS(c), label='Inelastic_'+str(e)+'cut')


plt.xlabel('Bubble radius [cm]', fontsize=25)
plt.ylabel('Events', fontsize=25)
plt.title('Elastic like events', fontsize=25)
ax = plt.subplot(111)
ax.legend(prop={'size': 20})
plt.show()
wait = input(' ')
