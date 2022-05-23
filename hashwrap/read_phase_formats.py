
import numpy as np

from hashwrap.libhashpy import check_pol
from hashwrap.hash_utils import get_sta_coords

degrad = 180./np.pi

"""
    Collection of functions to read the various FPFIT-like
    formats used in HASH examples:
        example 1 - north1.phase = FPFIT + the/az uncertainties
        example 2 - north2.phase = SCEC-like 
        example 3 - uses same phase as example 2
        example 4 - north4.phase = SCEC new 5-char stations
"""


"""
Example 1: north1.phase, Similar to FPFIT input file
           but with takeoff & azim uncertainties added:

  Event line:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-10      5i2     origin time, year, month, day, hour, minute
  11-14     f4.2    origin time, seconds
  15-16     i2      latitude, degrees
  17        a1      latitude, 'S'=south
  18-21     f4.2    latitude, minutes
  22-24     i3      longitude, degrees
  25        a1      longitude, 'E'=east
  26-29     f4.2    longitude, minutes
  30-34     f5.2    depth, km
  35-36     f2.1    magnitude
  81-88     2f4.2   horizontal and vertical uncertainty, km
  123-138   a16     event ID

  Polarity lines:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-4       a4      station name
  7         a1      polarity:U,u,+,D,d,or-
  8         i1      quality: 0=good quality, 1=lower quality, etc 
  59-62     f4.1    source-station distance (km)

  ** MTH:  The following specs from HASH_manual_v1.2.pdf DO NOT MATCH north1.phase:
  66-68     i3      takeoff angle 
  79-81     i3      azimuth
  83-85     i3      takeoff angle uncertainty ** NOT in standard FPFIT files
  87-89     i3      azimuth uncertainty ** NOT in standard FPFIT files

  ** MTH:  Use these instead:
  63-65     i3      takeoff angle 
  75-78     i3      azimuth
  81-83     i3      takeoff angle uncertainty ** NOT in standard FPFIT files
  86-88     i3      azimuth uncertainty ** NOT in standard FPFIT files
"""

'''
  94 1211104155034 1455118 3706 181323 31  0  0   1                  23              7  10                                           3143312 230
           1         2         3         4         5         6         7         8         9        100
 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
  IR2 IPD0                   0                               258121           51  10   1     X   VHZ
  SWM IPU0                   0                               528103            3  10   1     X   VHZ
  PYR IPU0                   0                               379110          342  10   1     X   VHZ
  SSN EPU1                   0                              1378 86          216  10   1     X   VHZ
'''


def read_fpfit_file(fpfile=None, plfile=None, stfile=None, delmax=120.,
                    fpfit_phase_format=1):

    """
    Reads an FPFIT input phase file to get polarities
      fpfit_phase_format = 1: FPFIT polarity lines + takeoff/azim uncertainties
                              e.g., Hardebeck & Shearer Example 1
                         = 2: Another FPFIT format wo/takeoff & azimuth. Use velocity
                              files to estimate range of takeoff and azimuth
                              e.g., Hardebeck & Shearer Example 2

    :param fpfile: path to FPFIT file
    :type fpfile: str
    :param plfile: path to Polarity Reversal file
    :type plfile: str
    :param delmax: max dist (km) to use polarity
    :type delmax: float
    :returns: events: list of event_dicts (one per event) with polarity info
    :rtype: list
    """

    fname = "read_fpfit_file"
    print("%s: fpfit_format:%d" % (fname, fpfit_phase_format))

    if fpfit_phase_format == 2:
        if stfile is not None:
            return read_fpfit_file_2(fpfile=fpfile, plfile=plfile, stfile=stfile,
                                     delmax=delmax)
        else:
            print("%s: Error: fpfit_phase_format=2 but no stfile given!!" % fname)
            return None
    elif fpfit_phase_format == 4:
        if stfile is not None:
            return read_fpfit_file_4(fpfile=fpfile, plfile=plfile, stfile=stfile,
                                     delmax=delmax)
        else:
            print("%s: Error: fpfit_phase_format=4 but no stfile given!!" % fname)
            return None

    events = []

    with open(fpfile) as f:
        while True:
            line = f.readline()
            if not line:
                break

            event_dict = {}

            try:
                iyr    = int(line[0:2])
                imon   = int(line[2:4])
                idy    = int(line[4:6])
                ihr    = int(line[6:8])
                imn    = int(line[8:10])
                qsec   = float(line[10:14])/100.
                ilatd  = int(line[14:16])
                cns    = line[16]
                qlatm  = float(line[17:21])/100.
                ilond  = int(line[21:24])
                cew    = line[24]
                qlonm  = float(line[25:29])/100.
                qdep   = float(line[29:34])/100.
                qmag   = float(line[34:36])/10.

                seh    = float(line[80:84])/100.
                sez    = float(line[84:88])/100.

                icusp  = int(line[122:138])
            except:
                print("Caught exception: processing event --> SKIP")

            if iyr < 50:    # fix Y2K problem
                iyr = iyr+2000
            else:
                iyr = iyr+1900

            qlat = ilatd + (qlatm / 60.0)
            if ('S' in cns):
                qlat *= -1
            qlon = -(ilond + (qlonm / 60.0))
            if ('E' in cew):
                qlon *= -1


            print("MTH: Process event: %4d-%02d-%02d %02d:%02d <%.3f, %.3f> h:%f [cusp id:%d]" % \
                (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp))

            event_dict['event_info'] = '%4d-%02d-%02d %02d:%02d <%.3f, %.3f> h:%f [cusp id:%d]' % \
                                        (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp)
            event_dict['event'] = {} 
            event_dict['event']['qlat'] = qlat
            event_dict['event']['qlon'] = qlon
            event_dict['event']['qdep'] = qdep
            event_dict['event']['sez'] = sez
            event_dict['event']['seh'] = seh
            event_dict['event']['qmag'] = qmag
            event_dict['event']['icusp'] = icusp


            sname_list = []
            p_pol_list = []
            p_qual_list = []
            qdist_list = []
            qthe_list = []
            qazi_list = []
            sazi_list = []
            sthe_list = []

            while True:
                ph = f.readline()
                sname     = ph[0:4]
                if not sname.strip():
                    #print("MTH: Encountered empty sta_code --> break")
                    break

                pickpol   = ph[6]
                p_qual    = int(ph[7])
                qdist     = float(ph[58:62])/10.
                ith       = int(ph[62:65])
                iaz       = int(ph[75:78])
                #print("MTH: sname:%s pickpol:%s p_qual:%s qdist:%f" % (sname, pickpol,p_qual,qdist))

                if qdist > delmax:
                    print("sname:%s dist:%.1f > delmax(%.1f) --> Skip" % (sname, qdist, delmax))
                    continue
                if p_qual > 1:
                    print("sname:%s P qual [%d] > 1 --> Skip" % (sname, p_qual))
                    continue

                try:
                    isthe     = int(ph[80:83])
                    isazi     = int(ph[85:88])
                except ValueError as e:
                    isthe = 10
                    isazi = 1

                #print("sta:%s pickpol:%s p_qual:%s qdist:%.1f ith:%d iaz:%d isthe:%d isazi:%d" % \
                     #(sname, pickpol, p_qual, qdist, ith, iaz, isthe, isazi))

                if pickpol in 'Uu+':
                    p_pol = 1
                elif pickpol in 'Dd-':
                    p_pol = -1
                else:
                    print("Unkown pick polarity=[%s]" % pickpol)

                if plfile:
                    spol = check_pol(plfile,sname,iyr,imon,idy,ihr)
                    if spol < 0:
                        #print("MTH: flip polarity for sname=%s" % (sname))
                        pass
                    p_pol *= spol

                qazi=float(iaz)
                if qazi < 0:
                    qazi = qazi+360.
        # MTH: Looks like the thetas in this input file are wrt Down and we want wrt Up:
                qthe=180.-float(ith)

                sname_list.append(sname)
                p_pol_list.append(p_pol)
                p_qual_list.append(p_qual)
                qdist_list.append(qdist)
                qthe_list.append(qthe)
                qazi_list.append(qazi)
                sazi_list.append(float(isazi))
                sthe_list.append(float(isthe))


            event_dict['sname'] = sname_list
            event_dict['p_pol'] = p_pol_list
            event_dict['p_qual'] = p_qual_list
            event_dict['qdist'] = qdist_list
            event_dict['qthe'] = qthe_list
            event_dict['qazi'] = qazi_list
            event_dict['sazi'] = sazi_list
            event_dict['sthe'] = sthe_list

            events.append(event_dict)


    return events


"""
Example 2: north2.phase, Similar to old SCEDC phase format:

  Event line:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-4       i4      origin time, year
  5-12      4i2     origin time, month, day, hour, minute
  13-17     f5.2    origin time, seconds
  18-19     i2      latitude, degrees
  20        a1      latitude, 'S'=south
  21-25     f5.2    latitude, minutes
  26-28     i3      longitude, degrees
  29        a1      longitude, 'E'=east
  30-34     f5.2    longitude, minutes
  35-39     f5.2    depth, km
  89-93     f5.2    horizontal location uncertainty, km
  95-99     f5.2    vertical location uncertainty, km
  140-143   f4.2    magnitude
  150-165   a16     event ID

  Polarity lines:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-4       a4      station name
  6-7       a2      station network
  10-12     a3      station component 
  14        a1      Ponset,IorE 
  16        a1      Ppolarity:U,u,+,D,d,or-

"""

def read_fpfit_file_2(fpfile=None, plfile=None, stfile=None, delmax=120.):
    """
    Reads a FPFIT format 2 input phase file to get polarities

    Note that this phase format does NOT have station azim,theta - these must be calculated
       from velocity models for given source/station location

    :param fpfile: path to FPFIT file
    :type fpfile: str
    :param plfile: path to Polarity Reversal file
    :type plfile: str
    :param sta_coords: dict of station coords: sta_coords['sta_name'] = <st_lat,st_lon,st_dep>
    :type sta_coords: dict
    :param delmax: max dist (km) to use polarity
    :type delmax: float
    :returns: events: list of event_dicts (one per event) with polarity info
    :rtype: list
    """

    events = []

    sta_coords = get_sta_coords(stfile)

    with open(fpfile) as f:
        while True:
            line = f.readline()
            if not line:
                print("Reached End of File")
                break

            event_dict = {}

            try:
                iyr    = int(line[0:4])
                imon   = int(line[4:6])
                idy    = int(line[6:8])
                ihr    = int(line[8:10])
                imn    = int(line[10:12])
                qsec   = float(line[12:17])
                ilatd  = int(line[17:19])
                cns    = line[19]
                qlatm  = float(line[20:25])
                ilond  = int(line[25:28])
                cew    = line[28]
                qlonm  = float(line[29:34])
                qdep   = float(line[34:39]) # then 49x
                seh    = float(line[88:93]) # then 1x
                sez    = float(line[94:99]) # then 40x
                qmag   = float(line[139:143]) # then 6x
                icusp  = int(line[149:165])
            except:
                print("Caught exception: processing event --> SKIP")

            qlat = ilatd + (qlatm / 60.0)
            if ('S' in cns):
                qlat *= -1
            qlon = -(ilond + (qlonm / 60.0))
            if ('E' in cew):
                qlon *= -1

            print("MTH: Process event: %4d-%02d-%02d %02d:%02d <%.1f, %.1f> h:%f [cusp id:%d]" % \
                (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp))

            event_dict['event_info'] = '%4d-%02d-%02d %02d:%02d <%.3f, %.3f> h:%f [cusp id:%d]' % \
                                        (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp)

            event_dict['event'] = {} 
            event_dict['event']['qlat'] = qlat
            event_dict['event']['qlon'] = qlon
            event_dict['event']['qdep'] = qdep
            event_dict['event']['sez'] = sez
            event_dict['event']['seh'] = seh
            event_dict['event']['qmag'] = qmag
            event_dict['event']['icusp'] = icusp


            sname_list = []
            p_pol_list = []
            p_qual_list = []
            qdist_list = []
            qthe_list = []
            qazi_list = []
            sazi_list = []
            sthe_list = []

            aspect = np.cos(qlat*np.pi/180.)

            npol_lines = 0
            while True:
                ph = f.readline()
                npol_lines += 1
                sname     = ph[0:4]
                if not sname.strip():
                    #print("MTH: Encountered empty sta_code --> break")
                    break
                snet      = ph[5:7]
                scomp     = ph[9:12]
                pickonset = ph[13]
                pickpol   = ph[15]

                if pickpol in 'Uu+':
                    p_pol = 1
                elif pickpol in 'Dd-':
                    p_pol = -1
                else:
                    print("Unkown pick polarity=[%s]" % pickpol)

                if pickonset in 'Ii':
                    p_qual = 0
                else:
                    p_qual = 1
                    #print("pickonset:%s ==> Skip" % pickonset)
                    continue

                #print("sta:%s net:%s comp:%s pickonset:%s pickpol:%s [p_pol:%d]" % \
                    #(sname, snet, scomp, pickonset, pickpol, p_pol))

                if sname in sta_coords:
                    flat,flon,felv = sta_coords[sname]
                else:
                    print("MTH: sname=%s not in sta_coords" % (sname))
                    continue

                dx = (flon - qlon) * 111.2 * aspect
                dy = (flat - qlat) * 111.2
                dist = np.sqrt(dx**2 + dy**2)
                qazi = 90. - np.arctan2(dy,dx) * 180./np.pi
                if (qazi < 0.):
                    qazi = qazi + 360.
                if (dist > delmax):
                    #print("dist=%.1f > delmax=%.1f --> Skip" % (dist, delmax))
                    continue

                if plfile:
                    spol = check_pol(plfile,sname,iyr,imon,idy,ihr)
                    if spol < 0:
                        #print("MTH: flip polarity for sname=%s" % (sname))
                        pass
                    p_pol *= spol


                sname_list.append(sname)
                p_pol_list.append(p_pol)
                p_qual_list.append(p_qual)
                qdist_list.append(dist)
                qazi_list.append(qazi)


            event_dict['sname'] = sname_list
            event_dict['p_pol'] = p_pol_list
            event_dict['p_qual'] = p_qual_list
            event_dict['qdist'] = qdist_list
            event_dict['qthe'] = qthe_list
            event_dict['qazi'] = qazi_list
            event_dict['sazi'] = sazi_list
            event_dict['sthe'] = sthe_list

            events.append(event_dict)
            print("read_phase icusp:%d npha:%d nkeep:%d" % (icusp, npol_lines, len(p_pol_list)))

    return events

"""
Example 4: north4.phase, Similar to *NEW* SCEDC phase format:

  Event line:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-4       i4      origin time, year
  5-12      4i2     origin time, month, day, hour, minute
  13-16     f4.2    origin time, seconds
  17-18     i2      latitude, degrees
  19        a1      latitude, 'S'=south
  20-23     f4.2    latitude, minutes
  24-26     i3      longitude, degrees
  27        a1      longitude, 'E'=east
  28-31     f4.2    longitude, minutes
  32-36     f5.2    depth, km
  131-146   a16     event ID
  148-150   f3.2    magnitude

  Polarity lines:
  -------------------------------------------------------------
  columns   format  value
  -------------------------------------------------------------
  1-5       a5      station name
  6-7       a2      station network
  10-12     a3      station component 
  14        a1      Ponset,IorE 
  16        a1      Ppolarity:U,u,+,D,d,or-
"""

def read_fpfit_file_4(fpfile=None, plfile=None, stfile=None, delmax=120.):
    """
    Reads a FPFIT format 4 (=Similar to *NEW* SCEDC format) input phase file to get polarities

    Note that this phase format does NOT have station azim,theta - these must be calculated
       from velocity models for given source/station location

    :param fpfile: path to FPFIT file
    :type fpfile: str
    :param plfile: path to Polarity Reversal file
    :type plfile: str
    :param sta_coords: dict of station coords: sta_coords['sta_name'] = <st_lat,st_lon,st_dep>
    :type sta_coords: dict
    :param delmax: max dist (km) to use polarity
    :type delmax: float
    :returns: events: list of event_dicts (one per event) with polarity info
    :rtype: list
    """

    events = []

    sta_coords = get_sta_coords(stfile)

    print("MTH: read from fpfile:%s" % fpfile)
    exit()
    with open(fpfile) as f:
        while True:
            line = f.readline()
            if not line:
                print("Reached End of File")
                break

            event_dict = {}

            try:
                iyr    = int(line[0:4])
                imon   = int(line[4:6])
                idy    = int(line[6:8])
                ihr    = int(line[8:10])
                imn    = int(line[10:12])
                qsec   = float(line[12:16])
                ilatd  = int(line[16:18])
                cns    = line[18]
                qlatm  = float(line[19:23])
                ilond  = int(line[23:26])
                cew    = line[26]
                qlonm  = float(line[27:31])
                qdep   = float(line[31:36]) 
                icusp  = int(line[130:146])
                qmag   = float(line[147:150])
            except:
                print("Caught exception: processing event --> SKIP")

            qlat = ilatd + (qlatm / 60.0)
            if ('S' in cns):
                qlat *= -1
            qlon = -(ilond + (qlonm / 60.0))
            if ('E' in cew):
                qlon *= -1

            print("MTH: Process event: %4d-%02d-%02d %02d:%02d <%.1f, %.1f> h:%f [cusp id:%d]" % \
                (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp))

            event_dict['event_info'] = '%4d-%02d-%02d %02d:%02d <%.3f, %.3f> h:%f [cusp id:%d]' % \
                                        (iyr, imon, idy, ihr, imn, qlat, qlon, qdep, icusp)

            event_dict['event'] = {} 
            event_dict['event']['qlat'] = qlat
            event_dict['event']['qlon'] = qlon
            event_dict['event']['qdep'] = qdep
            event_dict['event']['sez'] = sez
            event_dict['event']['seh'] = seh
            event_dict['event']['qmag'] = qmag
            event_dict['event']['icusp'] = icusp

            sname_list = []
            p_pol_list = []
            p_qual_list = []
            qdist_list = []
            qthe_list = []
            qazi_list = []
            sazi_list = []
            sthe_list = []

            aspect = np.cos(qlat*np.pi/180.)

            while True:
                ph = f.readline()
                sname     = ph[0:5]
                if not sname.strip():
                    #print("MTH: Encountered empty sta_code --> break")
                    break
                snet      = ph[5:7]
                scomp     = ph[9:12]
                pickonset = ph[13]
                pickpol   = ph[15]

                if pickpol in 'Uu+':
                    p_pol = 1
                elif pickpol in 'Dd-':
                    p_pol = -1
                else:
                    print("Unkown pick polarity=[%s]" % pickpol)

                if pickonset in 'Ii':
                    p_qual = 0
                else:
                    p_qual = 1

                #print("sta:%s net:%s comp:%s pickonset:%s pickpol:%s [p_pol:%d]" % \
                    #(sname, snet, scomp, pickonset, pickpol, p_pol))

                if sname in sta_coords:
                    flat,flon,felv = sta_coords[sname]
                else:
                    print("MTH: sname=%s not in sta_coords" % (sname))
                    continue

                dx = (flon - qlon) * 111.2 * aspect
                dy = (flat - qlat) * 111.2
                dist = np.sqrt(dx**2 + dy**2)
                qazi = 90. - np.arctan2(dy,dx) * 180./np.pi
                if (qazi < 0.):
                    qazi = qazi + 360.
                if (dist > delmax):
                    print("dist=%.1f > delmax=%.1f --> Skip" % (dist, delmax))
                    #continue

                if plfile:
                    spol = check_pol(plfile,sname,iyr,imon,idy,ihr)
                    if spol < 0:
                        #print("MTH: flip polarity for sname=%s" % (sname))
                        pass
                    p_pol *= spol


                sname_list.append(sname)
                p_pol_list.append(p_pol)
                p_qual_list.append(p_qual)
                qdist_list.append(dist)
                qazi_list.append(qazi)


            event_dict['sname'] = sname_list
            event_dict['p_pol'] = p_pol_list
            event_dict['p_qual'] = p_qual_list
            event_dict['qdist'] = qdist_list
            event_dict['qthe'] = qthe_list
            event_dict['qazi'] = qazi_list
            event_dict['sazi'] = sazi_list
            event_dict['sthe'] = sthe_list

            events.append(event_dict)


    return events


