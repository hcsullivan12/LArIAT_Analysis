
import ROOT
import argparse
import math

# define histograms
hTruthInteractingKe  = ROOT.TH1D('hInteractingKe',  'Interacting',   24, 0, 1200)
hTruthIncidentKe     = ROOT.TH1D('hIncidentKe',     'Incident',      24, 0, 1200)
hTruthCrossSectionKe = ROOT.TH1D('hCrossSectionKe', 'Cross Section', 24, 0, 1200)

hLastZ          = ROOT.TH1D('hLastZ', 'Last z in tpc', 50, 0, 100)
hDeDx           = ROOT.TH1D('hDeDx', 'dEdX', 200, 0, 50)
hPitch          = ROOT.TH1D('hPitch', 'Track pitch', 100, 0, 5)

# define counters
n_events = 0
n_beam_match = 0
elastic_events = []
inelastic_events = []

# define cuts
VERTEX_CUT = 3.0
ANGLE_CUT  = 10
DOWNSTREAM_Z_CUT = 85

# define constants
SLAB_WIDTH     = 0.0047 
CONVERSION     = 2.1043084 # NUMBER_DENSITY * 1e-28
PARTICLE_MASS  = 139.57
ENERGY_LOSS    = 40

########################################################
def makePlots(dtype):
    name = 'myana_'+str(dtype)+'.root'
    f = ROOT.TFile.Open(name, 'RECREATE')
    hInteractingKe.Write()
    hIncidentKe.Write()
    hCrossSectionKe.Write()
    hLastZ.Write()
    hDeDx.Write()
    hPitch.Write()

    f.Close()

########################################################
def doXsCalculation(dtype):
    nBins = hCrossSectionKe.GetNbinsX()
    for iBin in range(1, nBins+1):
        if hIncidentKe.GetBinContent(iBin) == 0:
            continue
        if hInteractingKe.GetBinContent(iBin) == 0:
            continue

        int_ke = hInteractingKe.GetBinContent(iBin)
        inc_ke = hIncidentKe.GetBinContent(iBin)
        temp_xs = (int_ke/inc_ke) * (1/CONVERSION) * (1/SLAB_WIDTH)
        hCrossSectionKe.SetBinContent(iBin, temp_xs)

        var = int_ke * (1 - int_ke/inc_ke)
        num_err = var**0.5
        den_err = (inc_ke)**0.5
        
        term1 = num_err/int_ke
        term2 = den_err/inc_ke

        total_err = temp_xs * ( term1 + term2 )
        if dtype == 'data':
            hCrossSectionKe.SetBinError(iBin, total_err)


########################################################
def inActiveVolume(pos):
    b = False
    if FV_X_BOUND[0] < pos.X() and pos.X() < FV_X_BOUND[1]:
        if FV_Y_BOUND[0] < pos.Y() and pos.Y() < FV_Y_BOUND[1]:
            if FV_Z_BOUND[0] < pos.Z() and pos.Z() < FV_Z_BOUND[1]:
                b = True
    return b

########################################################
def getFirstIntInTpc(ent):
    first_p = 0
    first_proc = 'none'
    for p, proc in zip(ent.InteractionPoint, ent.InteractionPointType):
        pos = ROOT.TVector3(ent.MidPosX[0][p], ent.MidPosY[0][p], ent.MidPosZ[0][p])
        if not inActiveVolume(pos):
            continue;
        first_p = p
        first_proc = proc
    return p,proc

########################################################
def getFirstAndLastPtInTpc(ent):
    enteredTpc = False
    sP = -1
    for c in range(0, len(MidPosX[0])):
        point = ROOT.TVector3(MidPosX[0][c],MidPosY[0][c],MidPosZ[0][c])
        if not inActiveVolume(point):
            continue
        # she entered!
        enteredTpc = True
        sP = c
        break
    
    # get last primary point
    # we have it easy; the last tpc point is either
    # 1) last visible pt of primary (no interaction)
    # 2) last physical pt of primary (non elastic interaction)
    eP = len(MidPosX[0])
    for c in range(len(MidPosX[0]), sP, -1):
        point = ROOT.TVector3(MidPosX[0][c],MidPosY[0][c],MidPosZ[0][c])
        if not inActiveVolume(point):
            continue
        # she entered!
        eP = c
        break
    
    assert(sP != eP)
    return enteredTpc, sP, eP

########################################################
def checkForMidpoints(ent, sp, ep):
    mp = []
    for ip, iproc in zip(InterestingPoint, InterstingPointType):
        if ip <= sp or ip >= ep:
            continue
        point = ROOT.TVector3(MidPosX[0][ip],MidPosY[0][ip],MidPosZ[0][ip])
        if not inActiveVolume(point):
            continue
        mp.append([ip, iproc])
        

########################################################
def getOrderedPoints(ent, sp, ep):
    # start point should be first in tpc,
    # end point should be last in tpc point
    # or last interaction point
    mp = checkForMidpoints(ent, sp, ep)
    
    start_point = ROOT.TVector3(MidPosX[0][sp],MidPosY[0][sp],MidPosZ[0][sp])
    end_point   = ROOT.TVector3(MidPosX[0][ep],MidPosY[0][ep],MidPosZ[0][ep])
    
    # order the points
    orderedPoints = {}
    orderedPoints.append(start_point)

    current_pt = sp
    current_point = start_point
    for intpt in mp:
        next_vertex = ROOT.TVector3(MidPosX[0][intpt],MidPosY[0][intpt],MidPosZ[0][intpt])
        length = getDistance(current_point, next_vertex)
        increment = 100*SLAB_WIDTH
        n_points = int(length/increment)

        for pt in range(1, n_points+1):
            next_point = current_point + pt * (increment/length) * (next_vertex - current_point)
            orderedPoints.append(next_point)
        # set the current pt as the vertex pt
        current_pt = intpt
        current_point = next_vertex
    
    orderedPoints.append(end_point)
    return orderedPoints

########################################################
def getIdes(ent, sp, ep):
    start_point = ROOT.TVector3(MidPosX[0][sp],MidPosY[0][sp],MidPosZ[0][sp])
    end_point   = ROOT.TVector3(MidPosX[0][ep],MidPosY[0][ep],MidPosZ[0][ep])

    for tid,e in zip(track_spt_idarr, track_spt_earr):
        

########################################################
def doAna(file, stop):
    global n_events
    global n_beam_match

    file = ROOT.TFile.Open(str(args.source), 'READ')
    anatree = file.Get('anatree/anatree')
    n_entries = anatree.GetEntries()

    # loop over the events
    for ent in anatree:
        n_events += 1
        if n_events%500 == 0:
            print n_events, '/', n_entries
        if stop > 0 and n_events == stop:
            break

        # get first and last tpc point
        entered, sp, ep = getFirstAndLastPtInTpc(ent)
        if not entered:
            return
        
        # get ordered points along primary track
        orderedPoints = getOrderedPoints(ent, sp, ep)

        # get ides
        ides = getIdes(ent, sp, ep)

        # intial ke
        in_mom = ROOT.TVector3(MidPx[0][sp], MidPy[0][sp], MidPz[0][sp])
        kin_en = (in_mom**2 + PARTICLE_MASS**2)**0.5
        hTruthKinEnFF.Fill(kin_en)

        # fill interacting and incident histograms
        old_point = 
        for pt range(1, len(orderedPoints)+1):
            this_point = orderedPoints[pt]
            hUniformDist.Fill((this_point-old_point).Mag())

            for dsf:
                


        # fill interacting and incident histograms
        assert(ent.num_wctracks == 1)
        wc_mom = ent.wctrk_momentum[0]
        kin_en = (wc_mom**2 + PARTICLE_MASS**2)**0.5 - PARTICLE_MASS
        kin_en -= ENERGY_LOSS
        for slb in range(0, len(trk_dedx)):
            hIncidentKe.Fill(kin_en)
            #print kin_en, trk_dedx[slb], trk_pitch[slb]
            kin_en -= trk_dedx[slb] * trk_pitch[slb]
        
        kin_en += trk_dedx[-1] * trk_pitch[-1]
        if is_inelastic:
            hInteractingKe.Fill(kin_en)

    # finish up
    doXsCalculation(dtype)
    makePlots(dtype)

########################################################
def getDistance(p1, p2):
    return (p1-p2).Mag()

########################################################
def getUnitDirection(s, e):
    return (e-s).Unit()

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for truth pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    parser.add_argument("-n", "--stop", help="The number of events to do (default=all)", default=-1, type=int)
    args = parser.parse_args()

    doAna(str(args.source), args.stop)