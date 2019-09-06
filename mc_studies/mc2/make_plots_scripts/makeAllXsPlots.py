
import ROOT
import argparse

# define histograms
hTrueXsKe = ROOT.TH1D('hTrueXsKe', 'PiMinus MC Inelastic Cross Section', 24, 0, 1200)
hTrueWellRecoXsKe = ROOT.TH1D('hTrueWellRecoXsKe', 'PiMinus MC Inelastic Cross Section', 24, 0, 1200)
hWellRecoXsKe = ROOT.TH1D('hWellRecoXsKe', 'PiMinus MC Inelastic Cross Section', 24, 0, 1200)
g4f = ROOT.TFile.Open('../mc_data/root_files/g4XsPredictionPiMinusInel.root', 'READ')
hG4Xs = g4f.Get('Graph') 
# define constants
SLAB_WIDTH     = 0.0047 
CONVERSION     = 2.1043084 # NUMBER_DENSITY * 1e-28

########################################################
def doXsCalculation(hIncKe, hIntKe, type):
    nBins = hIncKe.GetNbinsX()
    for iBin in range(1, nBins+1):
        if hIncKe.GetBinContent(iBin) == 0:
            continue
        if hIntKe.GetBinContent(iBin) == 0:
            continue

        int_ke = hIntKe.GetBinContent(iBin)
        inc_ke = hIncKe.GetBinContent(iBin)
        temp_xs = (int_ke/inc_ke) * (1/CONVERSION) * (1/SLAB_WIDTH)
        if type=='true':           hTrueXsKe.SetBinContent(iBin, temp_xs)
        if type=='true_well_reco': hTrueWellRecoXsKe.SetBinContent(iBin, temp_xs)
        if type=='well_reco':      hWellRecoXsKe.SetBinContent(iBin, temp_xs)

        var = int_ke * (1 - int_ke/inc_ke)
        num_err = var**0.5
        den_err = (inc_ke)**0.5
        
        term1 = num_err/int_ke
        term2 = den_err/inc_ke

        total_err = temp_xs * ( term1 + term2 )
        if type=='true':hTrueXsKe.SetBinError(iBin, total_err)
        if type=='true_well_reco':hTrueWellRecoXsKe.SetBinError(iBin, total_err)
        if type=='well_reco':hWellRecoXsKe.SetBinError(iBin, total_err)

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    args = parser.parse_args()

    f = ROOT.TFile.Open(args.source, 'READ')
    hIntKe = f.Get('hTrueInteractingKe')
    hIncKe = f.Get('hTrueIncidentKe')
    doXsCalculation(hIncKe, hIntKe, 'true')
    hIntKe = f.Get('hWellRecoTrueInteractingKe')
    hIncKe = f.Get('hWellRecoTrueIncidentKe')
    doXsCalculation(hIncKe, hIntKe, 'true_well_reco')
    hIntKe = f.Get('hWellRecoInteractingKe')
    hIncKe = f.Get('hWellRecoIncidentKe')
    doXsCalculation(hIncKe, hIntKe, 'well_reco')

    cXs = ROOT.TCanvas('cXs', 'Cross section', 1000, 1000)
    hTrueXsKe.SetLineWidth(3)
    hTrueXsKe.SetLineColor(4)
    hTrueXsKe.GetXaxis().SetRangeUser(0,1000)
    hTrueXsKe.SetMinimum(0)
    hTrueXsKe.SetMaximum(2.0)
    hTrueXsKe.GetXaxis().SetTitle('KE [MeV]')
    hTrueXsKe.GetYaxis().SetTitle('Cross Section [barn]')
    hTrueWellRecoXsKe.SetLineWidth(3)
    hTrueWellRecoXsKe.SetLineColor(2)
    hTrueWellRecoXsKe.GetXaxis().SetRangeUser(0,1000)
    hTrueWellRecoXsKe.SetMinimum(0)
    hTrueWellRecoXsKe.SetMaximum(2.0)
    hTrueWellRecoXsKe.GetXaxis().SetTitle('KE [MeV]')
    hTrueWellRecoXsKe.GetYaxis().SetTitle('Cross Section [barn]')
    hWellRecoXsKe.SetLineWidth(3)
    hWellRecoXsKe.SetLineColor(6)
    hWellRecoXsKe.GetXaxis().SetRangeUser(0,1000)
    hWellRecoXsKe.SetMinimum(0)
    hWellRecoXsKe.SetMaximum(2.0)
    hWellRecoXsKe.GetXaxis().SetTitle('KE [MeV]')
    hWellRecoXsKe.GetYaxis().SetTitle('Cross Section [barn]')

    #cXs.SetGrid()
    #ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetGridStyle(2)
    #ROOT.gStyle.SetGridWidth(2)

    hTrueXsKe.Draw('l')
    hTrueWellRecoXsKe.Draw('l same')
    hWellRecoXsKe.Draw('l same')

    leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
    leg.AddEntry(hTrueXsKe, 'True XS')
    leg.AddEntry(hTrueWellRecoXsKe, 'True XS Well Reco')
    leg.AddEntry(hWellRecoXsKe, 'Reco XS Well Reco')
    leg.Draw('same')

    wait = input(' ')