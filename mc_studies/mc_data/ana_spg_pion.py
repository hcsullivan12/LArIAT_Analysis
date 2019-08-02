
import ROOT
import argparse
import math

# define histograms
hInteractingKe = ROOT.TH1D('hInteractingKe',  'Interacting',   24, 0, 1200)
hIncidentKe    = ROOT.TH1D('hIncidentKe',     'Incident',      24, 0, 1200)
hInteractingKe = ROOT.TH1D('hCrossSectionKe', 'Cross Section', 24, 0, 1200)

# define counters

# define cuts
VERTEX_CUT = 3.0

# define constants
SLAB_WIDTH     = 0.0047 
CONVERSION     = 2.1043084 # NUMBER_DENSITY * 1e-28

########################################################
def doXsCalculation():
    nBins = hCrossSection.GetNbinsX()
    for iBin in range(1, nBins+1):
        if hIncidentKe.GetBinContent(iBin) == 0:
            continue
        if hInteractingKe.GetBinContent(iBin) == 0:
            continue

        int_ke = hInteractingKe.GetBinContent(iBin)
        inc_ke = hIncidentKe.GetBinContent(iBin)
        temp_xs = (int_ke/inc_ke) * (1/CONVERSION) * (1/SLAB_WIDTH)
        hCrossSection.SetBinContent(iBin, temp_xs)

        var = int_ke * (1 - int_ke/inc_ke)
        num_err = var**0.5
        den_err = (inc_ke)**0.5
        
        term1 = num_err/int_ke
        term2 = den_err/inc_ke

        total_err = temp_xs * ( term1 + term2 )
        hCrossSection.SetBinError(iBin, total_err)

########################################################
# do ana
def doAna(file):
    file = ROOT.TFile.Open(str(args.source), 'READ')
    anatree = file.Get('anatree/anatree')
    n_entries = anatree.GetEntries()
    n_events = 0

    # loop over the events
    for ent in anatree:
        n_events += 1
        print n_events, '/', n_entries

        # loop over tracks
        n_trks_matched = 0
        id_matched_trk = -1
        for iTrk in range(0, ent.ntracks_reco):
            # count matched tracks
            if ent.track_WC2TPC_match[iTrk]:
                id_matched_trk = iTrk
                n_trks_matched += 1
        if n_trks_matched != 0: 
            continue
            
        # check if inelastic
        is_inelastic = isInelastic(id_matched_trk, ent)
        
        # we should have one trk matched
        # get the containers for calo information
        trk_x = ent.col_track_x[id_matched_trk]
        trk_y = ent.col_track_y[id_matched_trk]
        trk_z = ent.col_track_z[id_matched_trk]
        trk_dedx  = ent.col_track_dedx[id_matched_trk]
        trk_pitch = ent.col_track_pitch_hit[ide_matched_trk]

        # fill interacting and incident histograms
        assert(num_wctracks == 1)
        wc_mom = ent.wctrk_momentum[0]
        kin_en = (wc_mom**2 + PARTICLE_MASS**2)**0.5 - PARTICLE_MASS
        kin_en -= ENERGY_LOSS
        for slb in range(0, len(trk_dedx)):
            hIncidentKe.Fill(kin_en)
            kin_en -= trk_dedx[slb] * trk_pitch[slb]
        
        kin_en += trk_dedx[-1] * trk_pitch[-1]
        if is_inelastic:
            hInteractingKe.Fill(kin_en)

    doXsCalculation()

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
    trk2_id = trk2_id[0]
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

    trk1_dir.Print()
    print trk1_dir.Mag()
    trk2_dir.Print()
    print trk2_dir.Mag()
    print trk1_dir.Dot(trk2_dir)
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
        n_points = len(ent.col_track_x[iTrk])
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
        angle = getAngle(id_matched_trk, sec_ids, ent)
        if angle > ANGLE_CUT:
            return True

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
    return False

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    args = parser.parse_args()
    doAna(str(args.source))