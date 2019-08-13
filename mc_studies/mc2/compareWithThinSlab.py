import ROOT

path = '/home/hunter/Desktop/piInelasticStudies/hists'
f_true_inc = ROOT.TFile.Open(path+'/hMCIncidentKE.root', 'READ')
f_true_int = ROOT.TFile.Open(path+'/hMCInteractingKE.root', 'READ')
f_reco = ROOT.TFile.Open('XS_ANA.root', 'READ')

h_true_inc = f_true_inc.Get('Canvas_1').GetPrimitive('hXsG4IncidentKinEn') or 0
h_true_int = f_true_int.Get('Canvas_1').GetPrimitive('hXsG4InteractingKinEn') or 0
h_reco_inc = f_reco.Get('hIncidentKe') or 0
h_reco_int = f_reco.Get('hInteractingKe') or 0

h_true_inc.Scale(h_reco_inc.Integral()/h_true_inc.Integral())
h_true_int.Scale(h_reco_int.Integral()/h_true_int.Integral())

c_inc = ROOT.TCanvas('c_inc', 'c_inc', 800, 800)
h_true_inc.SetLineWidth(3)
h_reco_inc.Draw()
h_true_inc.Draw('same')

c_int = ROOT.TCanvas('c_int', 'c_int', 800, 800)
h_true_int.SetLineWidth(3)
h_reco_int.Draw()
h_true_int.Draw('same')

wait = input(' ')

