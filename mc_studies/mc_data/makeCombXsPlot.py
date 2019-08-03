import ROOT

df = ROOT.TFile.Open('myana_data.root', 'READ')
mf = ROOT.TFile.Open('myana_mc.root', 'READ')

# cross section plots
hXsData = df.Get('hCrossSectionKe')
hXsMc   = mf.Get('hCrossSectionKe')
# reco plots
hDeDxData  = df.Get('hDeDx')
hLastZData = df.Get('hLastZ')
hPitchData = df.Get('hPitch')
hDeDxMc  = mf.Get('hDeDx')
hLastZMc = mf.Get('hLastZ')
hPitchMc = mf.Get('hPitch')

#########################################
def makeDataPretty(h, t, u):
    h.SetMarkerStyle(8)
    h.SetMarkerColor(1)
    h.SetLineColor(1)
    h.GetXaxis().SetTitle(t+' ['+u+']')
    b = h.GetXaxis().GetBinWidth(1)
    h.GetYaxis().SetTitle('Entries/'+str(b)+' '+u)
    return h

#########################################
def makeMcPretty(h, t, u):
    h.SetLineColor(4)
    h.SetLineWidth(3)
    return h

#########################################
def getOrderedHists(hD, hM):
    if hD.GetMaximum() > hM.GetMaximum():
        return zip([hD, hM], ['pe','hist same'])
    else:
        return zip([hM, hD], ['hist', 'same pe'])

cXs = ROOT.TCanvas('cXs', 'Cross Section', 1000, 1000)
hXsData = makeDataPretty(hXsData, 'Cross Section', 'MeV')
hXsData.GetXaxis().SetTitle('KE [MeV]')
hXsData.GetYaxis().SetTitle('Cross Section/50 MeV [barns]')
hXsData.GetXaxis().SetRangeUser(0, 1000)
hXsMc = makeMcPretty(hXsMc, '','')
hXsData.Draw('p')
hXsMc.Draw('same')

data_int = [hDeDxData.Integral(), hLastZData.Integral(), hPitchData.Integral()]
mc_int   = [hDeDxMc.Integral(), hLastZMc.Integral(), hPitchMc.Integral()]

hDeDxMc.Scale (data_int[0]/mc_int[0])
hLastZMc.Scale(data_int[1]/mc_int[1])
hPitchMc.Scale(data_int[2]/mc_int[2])

hDeDxData  = makeDataPretty(hDeDxData,  'dEdX',   'MeV/cm')
hLastZData = makeDataPretty(hLastZData, 'Last Z', 'cm')
hPitchData = makeDataPretty(hPitchData, 'Pitch',  'cm')

hDeDxMc  = makeMcPretty(hDeDxMc,  'dEdX',   'MeV/cm')
hLastZMc = makeMcPretty(hLastZMc, 'Last Z', 'cm')
hPitchMc = makeMcPretty(hPitchMc, 'Pitch',  'cm')

cDeDex = ROOT.TCanvas('cDeDex', 'dEdx', 1000, 1000)
for h,o in getOrderedHists(hDeDxData, hDeDxMc):
    h.Draw(o)
cLastZ = ROOT.TCanvas('cLastZ', 'Last z', 1000, 1000)
for h,o in getOrderedHists(hLastZData, hLastZMc):
    h.Draw(o)
cPitch = ROOT.TCanvas('cPitch', 'Pitch', 1000, 1000)
for h,o in getOrderedHists(hPitchData, hPitchMc):
    h.Draw(o)


wait = input(' ')