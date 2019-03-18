import numpy as np

from hashwrap.hash_utils import get_sta_coords


def read_sp_ratios(ampfile=None,
                   stfile=None,
                   corfile=None,
                   fpfit_events=None,
                   ratmin=3.0,
                   delmax=120.):
    """
    """

    fname = 'read_sp_ratios'

    sta_coords = get_sta_coords(stfile)

    events = []

    s_to_p_ratios = {}

    sta_cors = read_sta_cor_file(corfile)
    ratios = read_amp_ratio_file(ampfile)

    for evid, event_dict in ratios.items():

        sta_ratios = {}

        found = False
        for event in fpfit_events:
            if event['event']['icusp'] == evid:
                found = True
                break

        if not found:
            print(">>> %s: No matching evid found in phase file for evid:%s --> Skip this event!" % (fname, evid))
            continue

        qlat = event['event']['qlat']
        qlon = event['event']['qlon']
        qdep = event['event']['qdep']

        event = {}
        event['icusp'] = int(evid)
        event['qdep'] = float(qdep)

        event['sname'] = []
        event['qazi'] = []
        event['dist'] = []
        event['sp_ratio'] = []

        aspect = np.cos(qlat*np.pi/180.)

        print("%s: Process icusp:%d" % (fname, evid))
        iload=1
        for key, rd in event_dict.items():

            sname = "%-4s" % key.split(".")[0].strip()

            if sname not in sta_coords:
                print("%s: sname=[%s] --> NOT FOUND in sta_coords" % (fname, sname))
                for k, v in sta_coords.items():
                    print("k=[%s] --> v=[%s]" % (k,v))
                exit()
                continue

            flat,flon,felv = sta_coords[sname]

            dx = (flon - qlon) * 111.2 * aspect
            dy = (flat - qlat) * 111.2
            dist = np.sqrt(dx**2 + dy**2)
            qazi = 90. - np.arctan2(dy,dx) * 180./np.pi
            if (qazi < 0.):
                qazi = qazi + 360.
            if (dist > delmax):
                print("dist=%.1f > delmax=%.1f --> Skip" % (dist, delmax))
                continue

            qpamp = rd['qpamp']
            qsamp = rd['qsamp']
            qns1 = rd['qns1']
            qns2 = rd['qns2']
            #print("key:%s: qns1:%f qns2:%f qpamp:%f qsamp:%f" % (key, qns1, qns2, qpamp, qsamp))
            scn = key.split(".")
            # For some reason the net codes in sta_corrections file are all set to 'XX':
            key = "%s.%s.%s" % (scn[0], scn[1], 'XX')
            qcor = 0
            if key in sta_cors:
                #print("Found --> %f" % (sta_cors[key]))
                qcor = sta_cors[key]

            if qcor == -999.:
                print("%s: key=%s has qcorr = %f !!! ==> Skip" % (fname, key, qcor))
                continue

            if qpamp == 0.:
                print("%s: key=%s has qpamp = %f !!! ==> Skip" % (fname, key, qpamp))
                continue

            if qns1 == 0. or qns2 == 0:
                print("%s: key=%s has qns1 (%f) or qns2(%f) = 0 !!! ==> Skip" % (fname, key, qns1, qns2))
                continue

            s2n1 = np.abs(qpamp)/qns1
            s2n2 = qsamp/qns2

            if s2n1 < ratmin or s2n2 < ratmin:
                print("%s: key:%s s2n1=%f or s2n2=%f < ratmin:%f --> Skip" % (fname, key, s2n1, s2n2, ratmin))
                continue

            spin = qsamp/np.abs(qpamp)
            sp_ratio = np.log10(spin) - qcor

            event['sname'].append(sname)
            event['qazi'].append(qazi)
            event['dist'].append(dist)
            event['sp_ratio'].append(sp_ratio)


            sta_ratios[key] = sp_ratio
            print("%s: Add sta_ratios[%s]=%.2f iload:%d" % (fname, key, sp_ratio, iload))
            #print("spin=%f qcor:%f --> sp_ratio:%f" % (spin, qcor, sp_ratio))
            iload += 1

        events.append(event)

        #s_to_p_ratios[evid] = sta_ratios


    return events
    #return s_to_p_ratios



def read_sta_cor_file(corfile):
    sta_cors = {}
    with open(corfile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            sname, scomp, snet, correction = line.split()
            key = "%s.%s.%s" % (sname, scomp, snet)

            sta_cors[key] = float(correction)

    return sta_cors


"""
         1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890
                           x qns1     x     qns2 x    qpamp x   qsamp
2148509     12
GRH  EHZ CI   36.53   28.30      0.715     35.705     16.989    124.245
NHL  ELZ CI    6.16   45.41      0.265      1.731      1.234     31.955
SYL  ELZ CI   51.13   50.76      0.111      0.159     -0.757     61.141
BRCY EHZ NO   10.89   23.49      0.075      0.102      8.449     78.074
"""
def read_amp_ratio_file(ampfile):
    ratios = {}
    with open(ampfile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            evid, nratios = line.split()
            d = {}

            nratios = int(nratios)
            #print("evid=%s nratios=%s" % (evid, nratios))
            for i in range(nratios):
                line = f.readline()
                sname, scomp, snet, x1, x2, qns1, qns2, qpamp, qsamp = line.split()
                qns1 = float(qns1)
                qns2 = float(qns2)
                qpamp = float(qpamp)
                qsamp = float(qsamp)
                #print("Got sname:%s scomp:%s qns1:%f qns2:%f" % (sname, scomp, qns1, qns2))

                key = "%s.%s.%s" % (sname, scomp, snet)
                d[key] = {}
                d[key]['qns1'] = qns1
                d[key]['qns2'] = qns2
                d[key]['qpamp'] = qpamp
                d[key]['qsamp'] = qsamp

            ratios[int(evid)] = d

    return ratios


