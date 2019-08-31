import ROOT

path = '/home/hunter/Desktop/piInelasticStudies/hists'
f_true_inc = ROOT.TFile.Open(path+'/hMCIncidentKE.root', 'READ')
f_true_int = ROOT.TFile.Open(path+'/hMCInteractingKE.root', 'READ')
f_reco = ROOT.TFile.Open('XS_ANA.root', 'READ')

h_true_inc = f_true_inc.Get('Canvas_1').GetPrimitive('hXsG4IncidentKinEn') or 0
h_true_int = f_true_int.Get('Canvas_1').GetPrimitive('hXsG4InteractingKinEn') or 0
h_reco_inc = f_reco.Get('hTrueIncidentKe') or 0
h_reco_int = f_reco.Get('hTrueInteractingKe') or 0

h_true_inc.Scale(1/h_true_inc.Integral(0, 1000))
h_true_int.Scale(1/h_true_int.Integral(0, 1000))
h_reco_inc.Scale(1/h_reco_inc.Integral(0, 1000))
h_reco_int.Scale(1/h_reco_int.Integral(0, 1000))

c_inc = ROOT.TCanvas('c_inc', 'c_inc', 800, 800)
h_reco_inc.SetMarkerStyle(8)
h_reco_inc.SetMarkerColor(ROOT.kAzure-7)
h_reco_inc.SetLineWidth(4)
h_reco_inc.SetMarkerSize(1.5)
h_reco_inc.SetLineColor(ROOT.kAzure-7)
h_reco_inc.GetXaxis().SetRangeUser(0,1000)
h_reco_inc.GetXaxis().SetTitle('KE [MeV]')
h_reco_inc.GetYaxis().SetTitle('Normalized events/50 MeV')
h_reco_inc.Draw('ep')
h_true_inc.SetLineWidth(4)
h_true_inc.SetLineColor(ROOT.kAzure-7)
h_true_inc.Draw('same hist ][')
leg_inc = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg_inc.AddEntry(h_true_inc, 'Thin Slab True EDep', 'l')
leg_inc.AddEntry(h_reco_inc, 'MC Reco', 'pl')
leg_inc.SetLineColor(0)
leg_inc.Draw('same')

c_int = ROOT.TCanvas('c_int', 'c_int', 800, 800)
h_reco_int.SetMarkerStyle(8)
h_reco_int.SetMarkerSize(1.5)
h_reco_int.SetMarkerColor(ROOT.kAzure-7)
h_reco_int.SetLineWidth(4)
h_reco_int.SetLineColor(ROOT.kAzure-7)
h_reco_int.GetXaxis().SetRangeUser(0,1000)
h_reco_int.GetXaxis().SetTitle('KE [MeV]')
h_reco_int.GetYaxis().SetTitle('Normalized events/50 MeV')
h_reco_int.Draw('ep')
h_true_int.SetLineWidth(4)
h_true_int.SetLineColor(ROOT.kAzure-7)
h_true_int.Draw('same hist ][')
leg_int = ROOT.TLegend(0.1,0.7,0.48,0.9)
leg_int.AddEntry(h_true_int, 'Thin Slab True EDep', 'l')
leg_int.AddEntry(h_reco_int, 'MC Reco', 'pl')
leg_int.SetLineColor(0)
leg_int.Draw('same')

wait = input(' ')

