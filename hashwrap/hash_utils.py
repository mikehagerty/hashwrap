# -*- coding: utf-8 -*-
"""
  hash_utils.py

 -by Mark Williams (2013), NSL
 
 utilities for running HASH using the python version of hash_driver2
 (hash_driver2.py)
"""

table = """oklahoma=> select a.sta, a.seedchan, a.iphase, a.qual, a.fm, s.arid, s.seaz, s.ema, m.oridin from arrival a, assocaro s, mec m where s.orid=m.oridin and s.arid=a.arid and s.orid=908308 and a.iphase='P' order by a.sta;
 sta  | seedchan | iphase | qual | fm |  arid  |       seaz       | ema | oridin
------+----------+--------+------+----+--------+------------------+-----+--------
 LR01 | HHZ      | P      | x    | c. | 968643 |  340.52360234072 | 116 | 908308
 LR02 | HHZ      | P      | x    | c. | 968653 | 43.4639287121458 | 109 | 908308
 LR03 | HHZ      | P      | i    | d. | 968833 | 247.030183457485 | 123 | 908308
 LR04 | HHZ      | P      | i    | c. | 968828 | 61.3944048473663 | 136 | 908308
 LR05 | HHZ      | P      | x    | d. | 968683 |   96.39565050853 | 113 | 908308
 LR06 | HHZ      | P      | i    | d. | 968838 |  263.81024246661 | 104 | 908308
 LR07 | HHZ      | P      | i    | c. | 968818 | 122.031669147261 |  98 | 908308
 LR08 | HHZ      | P      | x    | c. | 968713 | 225.321512616947 | 102 | 908308
 LR09 | HHZ      | P      | x    | c. | 968723 |  163.92592544886 | 115 | 908308
 LR10 | HHZ      | P      | x    | c. | 968733 | 132.366526048533 | 107 | 908308
 LR11 | HHZ      | P      | x    | c. | 968743 | 217.587589990806 |  96 | 908308
 LR12 | HHZ      | P      | x    | c. | 968753 | 169.440171201542 | 102 | 908308
(12 rows)
"""

def main():

    azimuths = []
    takeoffs = []
    polarities = []

    #print(table)
    lines = table.split('\n')
    for line in lines[3:-2]:
        print(line)
        fields = line.split('|')
        sta = fields[0]
        cha = fields[1]
        fm = fields[4]
        azi = float(fields[6])
        ema = float(fields[7])
        #ema = 0.
        polarity = 1 if 'c' in fm else -1
        print(sta, cha, fm, azi, ema)
        azimuths.append(azi)
        takeoffs.append(ema)
        polarities.append(polarity)

    import numpy as np
    polarities = np.array(polarities)
    takeoffs = np.array(takeoffs)
    azimuths = np.array(azimuths)
    sdr = [120, 60, -50]
    test_stereo(azimuths,takeoffs,polarities,sdr)

    return



def parameter(**kwargs):
	'''returns variables inside a fortran 'parameter' call'''
	# FUTURE: could just make them all globals and import from namespace
	return [kwargs[key] for key in kwargs]

def fortran_include(fname):
	'''functions similar to a fortran include'''
	vars = []
	f_inc = open(fname)
	for line in f_inc:
		if 'c' in line[0]:
			pass
		elif 'parameter' in line:
			vars.extend(eval(line))
		else:
			pass
	f_inc.close()
	return vars
	
def get_sta_coords(stfile):
	'''Stub which reads in the HASH example station file'''
	f = open(stfile)
	sta_coords = {}
	for s in f:
		sta = s[0:4]
		lat = float(s[42:50])
		lon = float(s[51:61])
		elv = float(s[62:67])/1000.
		if sta in sta_coords:
			pass
		else:
			sta_coords[sta] = [lat,lon,elv]
	f.close()
	return sta_coords

def test_stereo(azimuths,takeoffs,polarities,sdr=[]):
    '''
	Plots points with given azimuths, takeoff angles, and 
	polarities on a stereonet. Will also plot both planes
	of a double-couple given a strike/dip/rake
    '''
    import matplotlib.pyplot as plt
    import mplstereonet
    from obspy.imaging.beachball import aux_plane

    print("azimuths:")
    print(azimuths)
    print("takeoffs:")
    print(takeoffs)
    print("polarities:")
    print(polarities)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='stereonet')
    print(type(ax))
    up = polarities > 0
    print("up:")
    print(up)
    dn = polarities < 0
    print("dn:")
    print(dn)
    #h_rk = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'ro')
    # MTH: this seems to put the observations in the right location
    #  We're plotting a lower-hemisphere focal mech, and we have to convert
    #  the up-going rays to the right az/dip quadrants:
    #h_rk = ax.rake(azimuths[up]-90.+180.,90.-takeoffs[up],90, 'ro')
    #h_rk = ax.rake(azimuths[dn]-90.+180.,90.-takeoffs[dn],90, 'b+')
    #ax.rake(strike-90., 90.-dip, rake, 'ro', markersize=14)
    # aqms event:
    h_rk = ax.rake(azimuths[up]-90 +180.,takeoffs[up]-90,90, 'ro')
    h_rk = ax.rake(azimuths[dn]-90 +180.,takeoffs[dn]-90,90, 'b+')

    #h_rk = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'b+')
    print("HERE we ARE")
    if sdr:
        s2,d2,r2 = aux_plane(*sdr)
        h_rk = ax.plane(sdr[0],sdr[1],'g')
        h_rk = ax.rake(sdr[0],sdr[1],-sdr[2], 'go')
        h_rk = ax.plane(s2,d2, 'g')
        plt.show()



if __name__ == '__main__':
    main()

