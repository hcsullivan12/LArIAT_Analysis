
import ROOT
import argparse
import math

# define histograms
hInteractingKe  = ROOT.TH1D('hInteractingKe',  'Interacting',   24, 0, 1200)
hIncidentKe     = ROOT.TH1D('hIncidentKe',     'Incident',      24, 0, 1200)
hCrossSectionKe = ROOT.TH1D('hCrossSectionKe', 'Cross Section', 24, 0, 1200)
hTruthInteractingKe  = ROOT.TH1D('hInteractingKe',  'Interacting',   24, 0, 1200)
hTruthIncidentKe     = ROOT.TH1D('hIncidentKe',     'Incident',      24, 0, 1200)
hTruthCrossSectionKe = ROOT.TH1D('hCrossSectionKe', 'Cross Section', 24, 0, 1200)

hLastZ          = ROOT.TH1D('hLastZ', 'Last z in tpc', 50, 0, 100)
hDeDx           = ROOT.TH1D('hDeDx', 'dEdX', 200, 0, 50)
hPitch          = ROOT.TH1D('hPitch', 'Track pitch', 100, 0, 5)

'''
Maybe show plot of largest background or fault in selection
'''

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
def doFilter(id_matched_trk, ent):
    return

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
def updateForEvtDisplay(is_inelastic, ent, dtype):
    global elastic_events
    global inelastic_events

    # elastic
    if not is_inelastic and len(elastic_events)<10:
        elastic_events.append(ent.event)
    # inelastic
    if is_inelastic and len(inelastic_events)<10:
        inelastic_events.append(ent.event)
    # inelastic missed 
    #if dtype is 'mc':
    #    p, proc = getFirstIntInTpc(ent)
    #    if proc == 'pi-Inelastic':
    #    

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
def doAna(file, dtype, stop):
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

        # loop over tracks
        n_trks_matched = 0
        id_matched_trk = -1
        for iTrk in range(0, ent.ntracks_reco):
            # count matched tracks
            if ent.track_WC2TPC_match[iTrk]:
                id_matched_trk = iTrk
                n_trks_matched += 1
        if n_trks_matched != 1: 
            continue
        n_beam_match += 1

        # apply upstream filter
        doFilter(id_matched_trk, ent)
            
        # check if inelastic
        is_inelastic = isInelastic(id_matched_trk, ent)

        # fill example event displays
        updateForEvtDisplay(is_inelastic, ent, dtype)
        
        # we should have one trk matched
        # get the containers for calo information
        trk_x = ent.col_track_x[id_matched_trk]
        trk_y = ent.col_track_y[id_matched_trk]
        trk_z = ent.col_track_z[id_matched_trk]
        trk_dedx  = ent.col_track_dedx[id_matched_trk]
        trk_pitch = ent.col_track_pitch_hit[id_matched_trk]
        for z,d,p in zip(trk_z, trk_dedx, trk_pitch):
            hLastZ.Fill(z)
            hDeDx.Fill(d)
            hPitch.Fill(p)

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
# check if track is inverted
def checkInversion(s, e):
    if s.Z() < e.Z():
        return s,e
    else: 
        return e,s

########################################################
def getDistance(p1, p2):
    return (p1-p2).Mag()

########################################################
def getUnitDirection(s, e):
    return (e-s).Unit()

########################################################
def getAngle(trk1_id, trk2_id, ent):
    trk1_start_point = ROOT.TVector3( ent.track_start_x[trk1_id],
                                      ent.track_start_y[trk1_id],
                                      ent.track_start_z[trk1_id] )
    trk1_end_point   = ROOT.TVector3( ent.track_end_x[trk1_id],
                                      ent.track_end_y[trk1_id],
                                      ent.track_end_z[trk1_id] )
    trk1_start_point, trk1_end_point = checkInversion(trk1_start_point, trk1_end_point)
    trk1_dir = getUnitDirection(trk1_start_point, trk1_end_point)

    trk2_start_point = ROOT.TVector3( ent.track_start_x[trk2_id],
                                      ent.track_start_y[trk2_id],
                                      ent.track_start_z[trk2_id] )
    trk2_end_point   = ROOT.TVector3( ent.track_end_x[trk2_id],
                                      ent.track_end_y[trk2_id],
                                      ent.track_end_z[trk2_id] )
    trk2_start_point, trk2_end_point = checkInversion(trk2_start_point, trk2_end_point)
    trk2_dir = getUnitDirection(trk2_start_point, trk2_end_point)
    return math.acos(trk1_dir.Dot(trk2_dir))
    

########################################################
# test if inelastic
def isInelastic(id_matched_trk, ent):
    # three cases to consider
    # 1) If more than one track attached to vertex, yes
    # 2) If one, check angle
    # 3) If none, where is the end point? 

    # get the start and enpoint
    mtrk_n_points = len(ent.col_track_x[id_matched_trk])
    mtrk_start_point = ROOT.TVector3( ent.col_track_x[id_matched_trk][0],
                                      ent.col_track_y[id_matched_trk][0],
                                      ent.col_track_z[id_matched_trk][0] )
    mtrk_end_point   = ROOT.TVector3( ent.col_track_x[id_matched_trk][mtrk_n_points-1],
                                      ent.col_track_y[id_matched_trk][mtrk_n_points-1],
                                      ent.col_track_z[id_matched_trk][mtrk_n_points-1] )
    
    # check if inverted
    mtrk_start_point, mtrk_end_point = checkInversion(mtrk_start_point, mtrk_end_point)

    # loop over tracks to see if any are attached to vertex
    sec_ids = []
    for iTrk in range(0, ent.ntracks_reco):
        # skip the matched trk
        if iTrk == id_matched_trk:
            continue
        # make sure this track has points
        n_points = len(ent.col_track_x[iTrk])
        if n_points == 0:
            continue;
        
        start_point = ROOT.TVector3( ent.col_track_x[iTrk][0],
                                     ent.col_track_y[iTrk][0],
                                     ent.col_track_z[iTrk][0] )
        end_point   = ROOT.TVector3( ent.col_track_x[iTrk][n_points-1],
                                     ent.col_track_y[iTrk][n_points-1],
                                     ent.col_track_z[iTrk][n_points-1] )
        dist1 = getDistance(mtrk_end_point, start_point)
        dist2 = getDistance(mtrk_end_point, end_point)
        dist_min = min(dist1, dist2)

        # checl inversion
        if dist2 < dist1:
            temp = start_point
            start_point = end_point
            end_point = temp
        
        # we're done with this track if it's not close enough
        if dist_min > VERTEX_CUT:
            continue

        # we got one!
        sec_ids.append(iTrk)

    # we should have all secondaries connected to vertex
    # Case 1)
    if len(sec_ids) > 1:
        return True

    if len(sec_ids) == 1:
        assert(id_matched_trk != sec_ids[0])
        angle = getAngle(id_matched_trk, sec_ids[0], ent)
        if angle > ANGLE_CUT:
            return True
        else: 
            return False

    # '''
    # @todo Think about what to do if sec = 0 but the track endpoint 
    #       is in the middle of the tpc. It's not obvious to me right 
    #       what's better to do: assume no interaction, assume inelastic
    #       'or some other third thing' - SS
    #       For now, assume inelastic.
    # '''
    assert(len(sec_ids) is 0)
    if mtrk_end_point.Z() < DOWNSTREAM_Z_CUT:
        return True
    else:
        return False
    return False

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    parser.add_argument("-t", "--type", help="The input type (data/mc)", required=True)
    parser.add_argument("-n", "--stop", help="The number of events to do (default=all)", default=-1, type=int)
    args = parser.parse_args()
    assert(args.type == 'mc' or args.type == 'data')

    doAna(str(args.source), str(args.type), args.stop)