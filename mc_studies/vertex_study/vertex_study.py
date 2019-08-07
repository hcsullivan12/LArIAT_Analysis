import ROOT
from sets import Set
import argparse
import math

g_do_plotting = False
if g_do_plotting:
    import matplotlib.pyplot as plt
    from matplotlib  import cm

# histograms
hIdeEnVsRadElastic   = ROOT.TH2D('hIdeEnVsRadElastic',   'Energy vs Radius cut elastic', 10, 0, 10, 100, 0, 100)
hIdeEnVsRadInelastic = ROOT.TH2D('hIdeEnVsRadInelastic', 'Energy vs Radius cut inelastic', 10, 0, 10, 100, 0, 100)

# counters and varibles
g_n_events = 0
g_n_entered_tpc = 0
g_n_elastic_like = 0
g_processes = Set([])
g_pdg_codes = {11:'electron', 13:'muon', 211:'pion', 321:'kaon', 2212:'proton'}

EVENT_FOR_EVT_DSP = 97
g_plot_x, g_plot_y, g_plot_z, g_plot_e = [], [], [], []
g_plot_sub_x, g_plot_sub_y, g_plot_sub_z, g_plot_sub_e = [], [], [], []

# cuts
SECONDARY_LENGTH_CUT = 2.0
CONE_DISTANCE_CUT = 2.4
CONE_ANGLE_CUT = 120
CYLINDER_DISTANCE_CUT = 3

########################################################
def inActiveVolume(pos):
    b = False
    if 2. < pos.X() and pos.X() < 45:
        if -18 < pos.Y() and pos.Y() < 18:
            if 0 < pos.Z() and pos.Z() < 88:
                b = True
    return b

########################################################
def doFilter(ent):
    in_tpc = False
    # make sure primary entered tpc 
    for x,y,z in zip(ent.MidPosX[0], ent.MidPosY[0], ent.MidPosZ[0]):
        v = ROOT.TVector3(x, y, z)
        if not inActiveVolume(v):
            continue
        in_tpc = True
        break
    
    if not in_tpc:
        return False
    # entered!
    global g_n_entered_tpc
    g_n_entered_tpc += 1
    return True

########################################################
def getPrimaryPoint(ent, position):
    for ipt in range(0, len(ent.MidPosX[0])):
            v = ROOT.TVector3(ent.MidPosX[0][ipt], ent.MidPosY[0][ipt], ent.MidPosZ[0][ipt])
            if (v-position).Mag() < 0.01:
                return ipt

########################################################
def getInteractionsInTpc(ent):
    vertices = []
    # check if g4 gave us anything
    if len(ent.InteractionPoint) != 0:
        for ipt, iproc in zip(ent.InteractionPoint, ent.InteractionPointType):
            v = ROOT.TVector3(ent.MidPosX[0][ipt], ent.MidPosY[0][ipt], ent.MidPosZ[0][ipt])
            if not inActiveVolume(v):
                continue
            # in tpc, add it
            _vtx = {}
            _vtx['point']    = ipt
            _vtx['process']  = iproc
            _vtx['position'] = v
            vertices.append(_vtx)
    else:
        # check the daughters
        for ig4 in range(0, ent.geant_list_size):
            if ent.Mother[ig4] != 1:
                continue
            # make sure it began in tpc
            v = ROOT.TVector3(ent.StartPointx[ig4], ent.StartPointy[ig4], ent.StartPointz[ig4])
            if not inActiveVolume(v):
                continue
            # check the processes
            # we don't care about elastic or inelastic
            proc = ent.G4Process[ig4]
            if proc == 'hadElastic' or proc == 'pi-Inelastic' or proc == 'hIoni':
                continue

            global g_processes
            g_processes.add(proc)
            # in tpc, add it
            _vtx = {}
            _vtx['point']    = getPrimaryPoint(ent, v)
            _vtx['process']  = proc
            _vtx['position'] = v
            vertices.append(_vtx)
            # we're done!
            break 
    return vertices

########################################################
def dump():
    print('\n')
    print '--------------------------------------'
    print 'Events', g_n_events
    print 'Entered tpc', g_n_entered_tpc
    print 'Unknown process:'
    print(g_processes)
    print 'Elastic like', g_n_elastic_like

########################################################
def isCharged(pdg):
    pdg = abs(pdg)
    if pdg in g_pdg_codes:
        return True
    return False

########################################################
def getVisSecondaries(ent, vertex):
    vis_sec = []
    # get prim g4 id
    prim_ig4 = -1
    prim_trk_id = -1
    for ig4 in range(0, ent.geant_list_size):
        if ent.process_primary[ig4] == 0:
            continue
        prim_ig4 = ig4
        prim_trk_id = ent.TrackId[ig4]
        break

    # if prim doesn't end here
    if (vertex['point']+1) < (ent.NTrTrajPts[0]-1):
        vis_sec.append(prim_ig4)

    # daughters
    for ig4 in range(0, ent.geant_list_size):
        if ent.Mother[ig4] != prim_trk_id:
            continue
        
        # make sure she is charged
        if not isCharged(ent.PDG[ig4]):
            continue

        # this is a charged daughter
        p0 = ROOT.TVector3(ent.StartPointx[ig4], ent.StartPointy[ig4], ent.StartPointz[ig4])
        pf = ROOT.TVector3(ent.EndPointx[ig4], ent.EndPointy[ig4], ent.EndPointz[ig4])
        length = (pf-p0).Mag()
        if length < SECONDARY_LENGTH_CUT:
            continue

        # make sure she's attached
        if (vertex['position']-p0).Mag() < 0.01:
            vis_sec.append(ig4)
    
    return vis_sec

########################################################
def getIdes(ent):
    _ides = []
    for iT in range(0, len(ent.TrackIdes_x)):
        tid = ent.TrackId[iT]
        for x,y,z,e in zip(ent.TrackIdes_x[iT], ent.TrackIdes_y[iT], ent.TrackIdes_z[iT], ent.TrackIdes_e[iT]):
            ide = {}
            ide['tid'] = tid
            ide['x'] = x
            ide['y'] = y
            ide['z'] = z
            ide['e'] = e
            _ides.append(ide)
    return _ides

########################################################
def checkVicinity(ent, ide, vis_sec, vertex):
    new_tide = ide

    # get direction vectors
    position = ROOT.TVector3(ide['x'], ide['y'], ide['z'])
    pt       = vertex['point']
    vtx      = vertex['position']
    prim_dir = ROOT.TVector3(ent.MidPosX[0][pt-1], ent.MidPosY[0][pt-1], ent.MidPosZ[0][pt-1])
    prim_dir = (prim_dir - vtx).Unit()
    sec_dir  = ROOT.TVector3(ent.EndPointx[vis_sec[0]], ent.EndPointy[vis_sec[0]], ent.EndPointz[vis_sec[0]])
    sec_dir  = (sec_dir - vtx).Unit()

    # check if this is within some small radius of the vertex
    if (position-vtx).Mag() < 0.1:
        new_tide['e'] = None

    # first check primary
    u = position - vtx
    u_dot_prim  = u.Dot(prim_dir)
    u_cos_prim  = u_dot_prim/u.Mag()
    u_sin_prim  = (1 - u_cos_prim*u_cos_prim)**0.5
    u_dist_prim = u.Mag() * u_sin_prim

    if u_dot_prim < CONE_DISTANCE_CUT:
        cone_dist = u_dot_prim * math.tan( 0.5*CONE_ANGLE_CUT*math.pi/180. )
        if u_dist_prim < cone_dist:
            new_tide['e'] = None
    # check if in cylinder
    elif u_dist_prim < CYLINDER_DISTANCE_CUT:
        new_tide['e'] = None

    # Now check for secondary
    u_dot_sec   = u.Dot(sec_dir)
    u_cos_sec   = u_dot_sec/u.Mag()
    u_sin_sec   = (1 - u_cos_sec*u_cos_sec)**0.5
    u_dist_sec  = u.Mag() * u_sin_sec
  
    # check if in cone first
    if u_dot_sec < CONE_DISTANCE_CUT:
        cone_dist = u_dot_sec * math.tan( 0.5*CONE_ANGLE_CUT*math.pi/180. )
        if u_dist_sec < cone_dist:
            new_tide['e'] = None
  
    # check if in cylinder
    else:
        #if ent.event == EVENT_FOR_EVT_DSP and new_tide['z'] > 60:
        #    print u_dist_sec, CYLINDER_DISTANCE_CUT
        if u_dist_sec < CYLINDER_DISTANCE_CUT:
            new_tide['e'] = None

    return new_tide


########################################################
def subtractIdes(ent, _ides, vis_sec, vertex):
    sub_ides = []
    for ide in _ides:
        side = checkVicinity(ent, ide, vis_sec, vertex)
        if side['e'] is not None:
            sub_ides.append(side) 
    return sub_ides

########################################################
def doVertexStudy(ent, vis_sec, vertex):
    # make copy of ides
    _ides = getIdes(ent)

    # plotting example before
    if ent.event == EVENT_FOR_EVT_DSP and g_do_plotting:
        global g_plot_x
        global g_plot_y
        global g_plot_z
        global g_plot_e
        for ide in _ides:
            g_plot_x.append(ide['x'])
            g_plot_y.append(ide['y'])
            g_plot_z.append(ide['z'])
            g_plot_e.append(ide['e'])

    # subtract ides
    _ides = subtractIdes(ent, _ides, vis_sec, vertex)
    if ent.event == EVENT_FOR_EVT_DSP and g_do_plotting:
        global g_plot_sub_x
        global g_plot_sub_y
        global g_plot_sub_z
        global g_plot_sub_e
        for ide in _ides:
            g_plot_sub_x.append(ide['x'])
            g_plot_sub_y.append(ide['y'])
            g_plot_sub_z.append(ide['z'])
            g_plot_sub_e.append(ide['e'])

    # fill plots
    for xbin in range(1, hIdeEnVsRadElastic.GetNbinsX()+1):
        for ide in _ides:
            pos = ROOT.TVector3(ide['x'], ide['y'], ide['z'])
            r_cut = hIdeEnVsRadElastic.GetXaxis().GetBinCenter(xbin)
            diff = (pos-vertex['position']).Mag()
            if diff < r_cut:
                ybin = hIdeEnVsRadElastic.GetYaxis().FindBin(1000*ide['e'])
                if vertex['process'] == 'hadElastic':
                    hIdeEnVsRadElastic.Fill(xbin, ybin)
                elif vertex['process'] == 'pi-Inelastic':
                    hIdeEnVsRadInelastic.Fill(xbin, ybin)
                else:
                    print 'WOAH!'

########################################################
def doAna(file, stop, debug):
    file = ROOT.TFile.Open(str(args.source), 'READ')
    anatree = file.Get('anatree/anatree')
    n_entries = anatree.GetEntries()

    # loop over the events
    for ent in anatree:
        global g_n_events
        g_n_events += 1
        if debug:
            print '\nEvent', ent.event
        if g_n_events%50 == 0:
            print g_n_events, '/', n_entries
        if stop > 0 and g_n_events == stop:
            break

        # apply upstream filter
        if not doFilter(ent):
            if debug:
                print 'Did not enter TPC!'
            continue
            
        # get interactions in tpc
        vertices = getInteractionsInTpc(ent)
        # leave if we have no interaction
        if len(vertices) == 0:
            continue
        
        if debug:
            for v in vertices:
                print v['process'], v['position'].X(), v['position'].Y(), v['position'].Z()
        first_vtx_in_tpc = vertices[0]

        # get visible secondaries
        vis_sec = getVisSecondaries(ent, first_vtx_in_tpc)
        if debug:
            print 'Visible secondaries', len(vis_sec)

        # check if elastic like
        # we only care about elastic and inelastic
        if len(vis_sec) != 1:
            continue
        if first_vtx_in_tpc['process'] != 'hadElastic' and first_vtx_in_tpc['process'] != 'pi-Inelastic':
            continue
        if debug:
            print 'Elastic like!'

        # elastic like!
        global g_n_elastic_like
        g_n_elastic_like += 1

        # vertex study!
        doVertexStudy(ent, vis_sec, first_vtx_in_tpc)

    # finish up
    dump()



########################################################
def writeHists():
    cIdesElastic = ROOT.TCanvas('cIdesElastic', 'cIdesElastic', 800, 800)
    hIdeEnVsRadElastic.Draw('colz')
    cIdesInelastic = ROOT.TCanvas('cIdesInelastic', 'cIdesInelastic', 800, 800)
    hIdeEnVsRadInelastic.Draw('colz')
    wait = input(' ')

########################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ana code for pi inelastic")
    parser.add_argument("-s", "--source", help="The input root file", required=True)
    parser.add_argument("-n", "--stop", help="The number of events to do (default=all)", default=-1, type=int)
    parser.add_argument("-d", "--debug", help="Debug mode (default=0)", default=0, type=int)
    args = parser.parse_args()

    doAna(str(args.source), args.stop, args.debug)
    writeHists()

    # plot 
    if g_do_plotting:
        fig0 = plt.figure(0)
        ax0 = fig0.add_subplot(111)
        s0 = plt.scatter(g_plot_z, g_plot_x, c=g_plot_e, cmap=cm.gnuplot, vmin=0, vmax=0.5)
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)
        s1 = plt.scatter(g_plot_sub_z, g_plot_sub_x, c=g_plot_sub_e, cmap=cm.gnuplot, vmin=0, vmax=0.5)
        fig0.colorbar(s0)
        fig1.colorbar(s1)
        ax0.set_xlim([0,100])
        ax0.set_ylim([0,47])
        ax1.set_xlim([0,100])
        ax1.set_ylim([0,47])
        plt.show()