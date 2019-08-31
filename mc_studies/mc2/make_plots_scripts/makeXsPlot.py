
import ROOT
import argparse

# define histograms
hCrossSectionKe = ROOT.TH1D('hCrossSectionKe', 'PiMinus MC Inelastic Cross Section', 24, 0, 1200)
g4f = ROOT.TFile.Open('../mc_data/root_files/g4XsPredictionPiMinusInel.root', 'READ')
hG4Xs = g4f.Get('Graph') 
# define constants
SLAB_WIDTH     = 0.0047 
CONVERSION     = 2.1043084 # NUMBER_DENSITY * 1e-28

########################################################
def doXsCalculation(file):
    f = ROOT.TFile.Open(file, 'READ')
    hIntKe = f.Get('hWellRecoInteractingKe')
    hIncKe = f.Get('hWellRecoIncidentKe')

    nBins = hIncKe.GetNbinsX()
    for iBin in range(1, nBins+1):
        if hIncKe.GetBinContent(iBin) == 0:
            continue
        if hIntKe.GetBinContent(iBin) == 0:
            continue

        int_ke = hIntKe.GetBinContent(iBin)
        inc_ke = hIncKe.GetBinContent(iBin)
        temp_xs = (int_ke/inc_ke) * (1/CONVERSION) * (1/SLAB_WIDTH)
        hCrossSectionKe.SetBinContent(iBin, temp_xs)

        var = int_ke * (1 - int_ke/inc_ke)
        num_err = var**0.5
        den_err = (inc_ke)**0.5
        
        term1 = num_err/int_ke
        term2 = den_err/inc_ke

        total_err = temp_xs * ( term1 + term2 )
        hCrossSectionKe.SetBinError(iBin, total_err)

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    parser.add_argument("-t", "--type", help="The input type (data/mc)", required=True)
    args = parser.parse_args()

    doXsCalculation(str(args.source))

    cXs = ROOT.TCanvas('cXs', 'Cross section', 1000, 1000)
    hCrossSectionKe.SetLineWidth(2)
    hCrossSectionKe.SetMarkerStyle(8)
    hCrossSectionKe.SetMarkerSize(1.5)
    if args.type == 'data':
        hCrossSectionKe.SetMarkerColor(1)
        hCrossSectionKe.SetLineColor(1)
    else:
        hCrossSectionKe.SetFillColor(ROOT.kAzure-3)
        hCrossSectionKe.SetLineColor(ROOT.kViolet+4)
        hCrossSectionKe.SetMarkerColor(ROOT.kViolet+4)

    if args.type == 'mc':
        hCrossSectionKe.Draw('e')
    else:
        hCrossSectionKe.Draw()

    hCrossSectionKe.GetXaxis().SetRangeUser(0, 1000)
    #cXs.SetGrid()
    #ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetGridStyle(2)
    #ROOT.gStyle.SetGridWidth(2)

    hCrossSectionKe.SetMinimum(0)
    hCrossSectionKe.SetMaximum(1.2)
    hCrossSectionKe.GetXaxis().SetTitle('KE [MeV]')
    hCrossSectionKe.GetYaxis().SetTitle('Cross Section [barns]')
    hG4Xs.Draw('same l')
    hG4Xs.SetLineColor(ROOT.kOrange+7)
    hG4Xs.SetLineWidth(4)

    leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
    name = 'Data Reco'
    if args.type == 'mc': name = 'MC Reco (only pions)'
    leg.AddEntry(hCrossSectionKe, name, 'pl')
    leg.AddEntry(hG4Xs, 'G4Prediction Inelastic XS', 'l')
    leg.Draw('same')

    ##########
    output = ROOT.TFile.Open('histograms.root', 'RECREATE')
    hG4Xs.Write()
    hCrossSectionKe.Write()

    wait = input(' ')