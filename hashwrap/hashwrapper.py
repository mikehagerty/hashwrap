# -*- coding: utf-8 -*-
"""
Python wrapper around HASH focalmech code,
for use with AQMS db

:copyright:
    Mike Hagerty (m.hagerty@isti.com), 2022
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)

Originally forked from Mark Williams' HASHPy 2013
"""


from hashwrap.hash_utils import fortran_include, get_sta_coords, test_stereo

from hashwrap.libhashpy import (mk_table_add, angtable, ran_norm, get_tts, get_gap, check_pol,
                              focalmc, mech_prob, get_misf, focalamp_mc, get_misf_amp)

from hashwrap.read_phase_formats import read_fpfit_file
from hashwrap.read_sp_ratios import read_sp_ratios

import numpy as np
import sys, os
import toml

import logging
logger = logging.getLogger(__name__)

from pathlib import Path
base_dir = Path(__file__).parent


def calc_focal_mechanisms(events, hash_settings, phase_format='FPFIT',
                          events_sp=None, logger_in=None):
    """

    :param events: list of event_dictionaries
    :type events: list
    :param settings: hash settings
    :type settings: toml dict
    :param phase_format: {'SCEDC','FPFIT'} - indicates format of input polarity/phase file
    :type phase_format: str

    """

    global sname
    global p_pol, p_qual
    global p_azi_mc, p_the_mc
    global index, qdep2
    global sp_ratio

    global f1norm, f2norm
    global strike2, dip2, rake2
    global str_avg, dip_avg, rak_avg
    global var_est, var_avg
    global mfrac, stdr, mavg, prob, qual

    fname = 'calc_focal_mechanisms'

    global logger
    if logger_in is not None:
        logger = logger_in


    src_dir = os.path.join(base_dir, 'src')

    npick0, nmc0, nmax0 = fortran_include(os.path.join(src_dir,'param.inc'))
    dang0, ncoor        = fortran_include(os.path.join(src_dir,'rot.inc'))


    npolmin = hash_settings['npolmin']
    max_agap = hash_settings['max_agap']
    max_pgap = hash_settings['max_pgap']
    dang = hash_settings['dang']
    nmc = hash_settings['nmc']
    maxout = hash_settings['maxout']
    badfrac = hash_settings['badfrac']
    delmax = hash_settings['delmax']
    cangle = hash_settings['cangle']
    prob_max = hash_settings['prob_max']

    # Optional settings to use S/P ratios:
    qbadfrac = hash_settings.get('qbadfrac', None)
    ratmin = hash_settings.get('ratmin', None)



    if 'velocity_models' in hash_settings:
        velocity_models = hash_settings['velocity_models']
        for i in range(len(velocity_models)):
            ntab = mk_table_add(i+1,velocity_models[i])


    dang2=max(dang0,dang) # don't do finer than dang0

# initialize arrays

    #initialize_arrays(npick0, nmc0, nmax0)

    outputs = []

# List of events, each is an event_dict with lists of polarity, azi, the, etc

    #print("Initialize arrays: npick0=%d nmc0=%d nmax0=%d" % (npick0, nmc0, nmax0))
    for event in events:

        initialize_arrays(npick0, nmc0, nmax0)

        icusp = event['event']['icusp']
        #print(event)

        """
NHL  ELZ CI    5.77   45.31      0.620      2.068      1.937     16.295
SMF  ELZ CI  147.40   62.82      0.472      0.473      1.545      0.834
BRCY EHZ NO   10.09   23.22      2.328      2.901     26.992     54.247
CPCP EHZ NO  165.73   11.72      0.712      3.756     -2.314     13.117
LA00 EHZ NO  135.25   53.15      0.163      1.749     -1.248     15.408
LA01 EHZ NO  126.94   51.88      0.871      8.648      2.875     97.636
SMIP EHZ NO  296.77   16.57      0.396      0.447     -6.291     19.795
SMIP ELZ NO  296.77   16.57      0.086      0.123     -0.865      2.477
        """
        event_sp = None
        if events_sp is not None:
            for evt in events_sp:
                if evt['icusp'] == icusp:
                    event_sp = evt
                    #print(event_sp)
                    break

        sez = event['event']['sez']

        # Load velocity files used to calculate takeoff angles:
        if phase_format == 'SCEDC':
            qdep = event['event']['qdep']
            qdep2[0] = qdep
            index[0] = 1
            for nm in range(1,nmc):
                val = ran_norm()
                qdep2[nm] = qdep + sez * val # randomly perturbed source depth
                index[nm] = (nm % ntab) + 1  # index used to choose velocity model


      # Loop over station polarities for this event

        for k,sta in enumerate(event['sname']):

            sname[k] = sta

            # In example 1 FPFIT, pol, qual, qazi, qthe and qdist are directly read in on phase line
            # In example 2 FPFIT, pol, qual are read in, qazi and qdist are calculated from sta_coords

            p_pol[k]  = event['p_pol'][k]
            p_qual[k] = event['p_qual'][k]
            qazi = event['qazi'][k]
            dist = event['qdist'][k]

            if phase_format == 'SCEDC':
            # Calculate nmc trial takeoff angles: nmc = n_vel_models * n_perturbed_source_depths
            # For example 2 we try a range of takeoff angles (one for each vel model + source dep)
            #               with azi fixed
                for nm in range(nmc):
                    p_azi_mc[k,nm] = qazi
                    p_the_mc[k,nm], iflag = get_tts(index[nm],dist,qdep2[nm])
                    if iflag != 0:
                        logger.warning("%s: get_tts returned iflag=%d (should be 0 !)" % (fname, iflag))

                #print("sname[%2d]:%4s dist:%6.1f   az:%6.1f   qthe:%4.1f   pol:[%d]" % \
                    #(k, sta, dist, qazi, p_the_mc[k,0], p_pol[k]))


            elif phase_format == 'FPFIT':
                qthe = event['qthe'][k]
                sthe = event['sthe'][k]    # Calculate nmc trial takeoff angles
                sazi = event['sazi'][k]    # using input error estimates in theta and azi

                p_azi_mc[k,0] = qazi
                p_the_mc[k,0] = qthe

                #print("sname[%2d]:%4s dist:%6.1f   az:%6.1f   qthe:%4.1f   pol:[%d]" % \
                    #(k, sta, dist, qazi, p_the_mc[k,0], p_pol[k]))

            # np.random.normal returns size=n randomly selected vals from normal distribution 
            # with mean=qazi(or qthe) and std=sazi(or sthe)
                azims  = np.random.normal(loc=qazi, scale=sazi, size=nmc-1)
                thetas = np.random.normal(loc=qthe, scale=sthe, size=nmc-1)

                for nm in range(1,nmc):
                    #val = ran_norm()
                    #p_azi_mc[k,nm] = qazi + sazi*val
                    #val = ran_norm()
                    #p_the_mc[k,nm] = qthe + sthe*val
                    p_azi_mc[k,nm] = azims[nm-1]
                    p_the_mc[k,nm] = thetas[nm-1]

        nppl = len(event['sname'])

        nspr = 0

        if event_sp:
            #print("Loop over sp_ratio: %d amps" % len(event_sp['sp_ratio']))
            for k in range(len(event_sp['sp_ratio'])):

                j = k + nppl
                sname[j]    = event_sp['sname'][k]
                sp_ratio[j] = event_sp['sp_ratio'][k]
                qazi        = event_sp['qazi'][k]
                dist        = event_sp['dist'][k]

                nspr += 1

                if phase_format == 'SCEDC':
                # Calculate nmc trial takeoff angles: nmc = n_vel_models * n_perturbed_source_depths
                    for nm in range(nmc):
                        p_azi_mc[j,nm] = qazi
                        p_the_mc[j,nm], iflag = get_tts(index[nm],dist,qdep2[nm])
                        if iflag != 0.:
                            logger.warning("%s: get_tts returned iflag=%d for SP ratio observation (should be 0 !)" \
                                  % (fname, iflag))

                elif phase_format == 'FPFIT':
                    #qthe = event['qthe'][k]
                    #sthe = event['sthe'][k]    # Calculate nmc trial takeoff angles
                    #sazi = event['sazi'][k]    # using input error estimates in theta and azi

                    qthe = event_sp['qthe'][k]
                    sthe = event_sp['sthe'][k]    # Calculate nmc trial takeoff angles
                    sazi = event_sp['sazi'][k]    # using input error estimates in theta and azi
                    if k < len(event['qthe']):
                        if event['qthe'][k] != event_sp['qthe'][k]:
                            logger.warning("k=%d event['sname']=%s qthe=%.2f -vs- event_sp['sname']=%s qthe=%.2f" %
                                  (k, event['sname'][k], event['qthe'][k], event_sp['sname'][k], event_sp['qthe'][k]))

                    p_azi_mc[j,0] = qazi
                    p_the_mc[j,0] = qthe
                    azims  = np.random.normal(loc=qazi, scale=sazi, size=nmc-1)
                    thetas = np.random.normal(loc=qthe, scale=sthe, size=nmc-1)

                    for nm in range(1,nmc):
                        p_azi_mc[j,nm] = azims[nm-1]
                        p_the_mc[j,nm] = thetas[nm-1]


        npol = nppl + nspr

        logger.info("%s: icusp=[%d] npol:%d = nppl:%d + nspr:%d" % \
                    (fname, event['event']['icusp'], npol, nppl, nspr))

        logger.info(" i: sname    azi      the  pol   s/p")
        for i in range(npol):
            if i == nppl:
                logger.info("=====================================================")
            logger.info("%2d: %s %8.2f %8.2f %3d %6.3f" % (i, sname[i].decode("utf-8"), p_azi_mc[i,0], p_the_mc[i,0], p_pol[i], sp_ratio[i]))


        # stop if there aren't enough polarities
        if (npol < npolmin):
            logger.warning("%s: npol=%d < npolmin (%d) NOT Enough Polarities --> Skip" % (fname, npol, npolmin))
            #qual[0] = 'F'
            break
            #continue

        use_amplitudes = False
        if nspr > 0:
            use_amplitudes = True

        # determine maximum azimuthal and takeoff gap in polarity observations and stop if either gap is too big
        magap,mpgap = get_gap(p_azi_mc[:npol,0],p_the_mc[:npol,0],npol)

        if ((magap > max_agap) or (mpgap > max_pgap)):
            logger.warning("%s: Azimuthal/takeoff gap too large !!!", fname)
            logger.warning("%s: magap:%f max_agap:%f mpgap:%f max_pgap:%f" % (fname, magap, max_agap, mpgap, max_pgap))
            #qual[0] = 'E'
            break

        # determine maximum acceptable number misfit polarities
        #nmismax = max(int(npol * badfrac),2)        # nint
        #nextra  = max(int(npol * badfrac * 0.5),2)  # nint
        nmismax = max(int(nppl * badfrac),2)        # nint
        nextra  = max(int(nppl * badfrac * 0.5),2)  # nint

        # find the set of acceptable focal mechanisms for all trials

        logger.info("%s: nppl:%d badfrac:%f --> ntotal=nmismax:%d nextra:%d" % (fname, nppl, badfrac, nmismax, nextra))

        if use_amplitudes:
            # determine maximum acceptable number misfit polarities
            #nmismax = max(int(nppl * self.badfrac),2)        # nint
            #nextra  = max(int(nppl * self.badfrac * 0.5),2)  # nint
            qmismax = max(int(nspr * qbadfrac),2)        # nint
            qextra  = max(int(nspr * qbadfrac * 0.5),2)  # nint
            logger.info("%s: nspr:%d qbadfrac:%f --> qmismax:%d" % (fname, nspr, qbadfrac, qmismax))
            # find the set of acceptable focal mechanisms for all trials
            #for i in range(npol):
                #print("%2d: %s %5.2f %5.2f %3d %6.3f" % (i, sname[i], p_azi_mc[i,0], p_the_mc[i,0], p_pol[i], sp_ratio[i]))

            #exit()
            nf2,strike2,dip2,rake2,f1norm,f2norm = focalamp_mc(p_azi_mc, p_the_mc, sp_ratio[:npol], p_pol[:npol],\
                                                        nmc, dang2, nmax0, nextra, nmismax, qextra, qmismax, npol)
        else:
            # determine maximum acceptable number misfit polarities
            #nmismax = max(int(npol * badfrac),2)        # nint
            #nextra  = max(int(npol * badfrac * 0.5),2)  # nint
            # find the set of acceptable focal mechanisms for all trials

            # MTH: nf2 = nfaults returned from FOCALMC --> will typically = nmax0 (e.g., 500) sampled
            #            randomly from the nrot acceptable mechanisms FOCALMC finds

            nf2,strike2,dip2,rake2,f1norm,f2norm = focalmc(p_azi_mc, p_the_mc, p_pol[:npol], p_qual[:npol], \
                                                           nmc, dang2, nmax0, nextra, nmismax, npol)


        # nout2 will typically = nf2=nmax0=500
        nout2 = min(nmax0,nf2)  # number mechs returned from sub
        # nout1 is not currently used
        nout1 = min(maxout,nf2) # number mechs to return

        logger.debug("%s: nmax0:%d nf2:%d --> sub returned nout2:%d  maxout:%d --> return nout1:%d mechs" %
                    (fname, nmax0,nf2, nout2, maxout, nout1))

        # find the probable mechanism from the set of acceptable solutions
#      subroutine MECH_PROB(nf,norm1in,norm2in,cangle,prob_max,
#                           nsltn,str_avg,dip_avg,rak_avg,prob,rms_diff)

        nmult,str_avg,dip_avg,rak_avg,prob,var_est = mech_prob(f1norm[:,:nout2], f2norm[:,:nout2], cangle,\
                                                               prob_max, nout2)

        logger.info("%s: mech_prob(prob_max=%f nout2=%d) returned: nmult:%d" % (fname, prob_max, nout2, nmult))


        # MTH: decide what to do
        if nmult > 1:
            logger.warning("%s: cuspid=%d --> nmult=%d : Only returning First preferred solution" % (fname, icusp, nmult))
            #exit()

        for imult in range(nmult):
        #for imult in range(1):
            var_avg[imult] = (var_est[0,imult] + var_est[1,imult]) / 2.
            #print("MTH: imult=%d nmult=%d Here is the pref soln:" % (imult,nmult))
            #print(('cid = {0} {1}  mech = {2} {3} {4}'.format(icusp,imult,str_avg[imult],dip_avg[imult],rak_avg[imult])))
            # find misfit for prefered solution

            mavge = "N/A"
            if use_amplitudes:
                #print("npol=%d type(p_azi_mc)=%s p_azi_mc.size=%d" % (npol,type(p_azi_mc), p_azi_mc.size))
                #for x in sp_ratio[:npol]:
                    #print(x)
                mfrac[imult], mavg[imult], stdr[imult] =  get_misf_amp(p_azi_mc[:npol,0], p_the_mc[:npol,0],
                                                                       sp_ratio[:npol], p_pol[:npol],
                                                                       str_avg[imult], dip_avg[imult], rak_avg[imult])
                """
                mfrac[imult], mavg[imult], stdr[imult] =  get_misf_amp(p_azi_mc[:npol,0], p_the_mc[:npol,0],
                                                                       sp_ratio[:npol], p_pol[:npol],
                                                                       109, 67, -86)
            # MTH: Doesn't seem to matter if you pass in npol or not
                                                                       #str_avg[imult], dip_avg[imult], rak_avg[imult],
                                                                       #npol)
                exit()
                """
                mavge = "%.2f" % mavg[imult]

            else:
                mfrac[imult], stdr[imult] =  get_misf(p_azi_mc[:npol,0], p_the_mc[:npol,0],
                                                      p_pol[:npol], p_qual[:npol],
                                                      str_avg[imult], dip_avg[imult], rak_avg[imult],
                                                      npol)

            # solution quality rating  ** YOU MAY WISH TO DEVELOP YOUR OWN QUALITY RATING SYSTEM **
            if ((prob[imult] > 0.8) and (var_avg[imult] < 25) and (mfrac[imult] <= 0.15) and (stdr[imult] >= 0.5)):
                qual[imult]='A'
            elif ((prob[imult] > 0.6) and (var_avg[imult] <= 35) and (mfrac[imult] <= 0.2) and (stdr[imult] >= 0.4)):
                qual[imult]='B'
            elif ((prob[imult] > 0.5) and (var_avg[imult] <= 45) and (mfrac[imult] <= 0.3) and (stdr[imult] >= 0.3)):
                qual[imult]='C'
            else:
                qual[imult]='D'

            logger.info("%s: imult:%d str:%d dip:%d rake:%d rms_diff:%.2f prob:%.2f mfrac:%.2f mavg:%s stdr:%.2f Q:%s" %
                        (fname, imult, round(str_avg[imult]), round(dip_avg[imult]), round(rak_avg[imult]),
                                      var_avg[imult], prob[imult], mfrac[imult], mavge, stdr[imult], qual[imult].decode('UTF-8')))

            #print("Solution Quality:[%s]" % qual[imult])

# Things that obspy event.FocalMechanism wants:
      # Turn strike,dip,rake into (2) obspy nodal planes
        strike = str_avg[0]
        dip    = dip_avg[0]
        rake   = rak_avg[0]

    # Extras:
        rms_angle = var_avg[0]

        out_dict = {}
        out_dict['icusp'] = icusp
        out_dict['event_info'] = event['event_info']
        out_dict['strike'] = strike
        out_dict['dip'] = dip
        out_dict['rake'] = rake
        out_dict['azim_gap'] = magap
        out_dict['theta_gap'] = mpgap
        out_dict['stdr'] = stdr[0]
        out_dict['misfit'] = mfrac[0]
        out_dict['use_amplitudes'] = use_amplitudes
        if use_amplitudes:
            out_dict['S_P_log10_misfit'] = mavg[0]
            out_dict['number_SP_amp_ratios'] = nspr
            out_dict['mavg'] = mavg[0]

        out_dict['number_P_polarities'] = nppl
        out_dict['station_polarity_count'] = npol  # Used by obspy npol = nppl + nspr

        out_dict['quality'] = qual[0]
        out_dict['rms_angle'] = rms_angle

        outputs.append(out_dict)
        #print("append icusp:%s to outputs" % icusp)

    return outputs



def initialize_arrays(npick0, nmc0, nmax0):

# input arrays
    global sname, scomp, snet
    global pickpol, pickonset, p_pol, p_qual, sp_ratio
    global spol, p_azi_mc, p_the_mc
    global index, qdep2

# output arrays
    global f1norm, f2norm
    global strike2, dip2, rake2
    global str_avg, dip_avg, rak_avg
    global var_est, var_avg
    global mfrac, stdr, mavg, prob, qual

    sname     = np.empty(npick0, 'a4', 'F')
    scomp     = np.empty(npick0, 'a3', 'F')
    snet      = np.empty(npick0, 'a2', 'F')
    pickpol   = np.empty(npick0, 'a1', 'F')
    pickonset = np.empty(npick0, 'a1', 'F')
    p_pol     = np.zeros(npick0, int, 'F')
    #p_pol     = np.empty(npick0, int, 'F')
    p_qual    = np.empty(npick0, int, 'F')
    spol      = np.empty(npick0, int, 'F')
    p_azi_mc  = np.empty((npick0,nmc0), float, 'F')
    p_the_mc  = np.empty((npick0,nmc0), float, 'F')
    index     = np.empty(nmc0, int, 'F')
    qdep2     = np.empty(nmc0, float, 'F')

    sp_ratio  = np.zeros(npick0, float, 'F')
    #sp_ratio  = np.empty(npick0, float, 'F')

    f1norm  = np.empty((3,nmax0), float, 'F')
    f2norm  = np.empty((3,nmax0), float, 'F')
    strike2 = np.empty(nmax0, float, 'F')
    dip2    = np.empty(nmax0, float, 'F')
    rake2   = np.empty(nmax0, float, 'F')
    str_avg = np.empty(5, float, 'F')
    dip_avg = np.empty(5, float, 'F')
    rak_avg = np.empty(5, float, 'F')
    var_est = np.empty((2,5), float, 'F')
    var_avg = np.empty(5, float, 'F')
    mfrac   = np.empty(5, float, 'F')
    stdr    = np.empty(5, float, 'F')
    mavg    = np.empty(5, float, 'F')
    prob    = np.empty(5, float, 'F')
    qual    = np.empty(5, 'a', 'F')

    return



"""
HashDriver1.f - reads in north1.phase = FPFIT format with takeoff/azim uncertainties added
                In example 1, we assume that we have already computed the azimuth and takeoff angle to each station,
                and have estimated the uncertainty of both of these angles. The input format (file: north1.phase) is
                a modified version of the standard FPFIT input format, expanded to include the azimuth and takeoff angle
                uncertainties. In this case much of the data can be read directly into the input arrays. The azimuth
                and takeoff angle for each trial are computed by randomly selecting values from normal distributions
                with the given average value and standard deviation.

HashDriver2.f - reads in north2.phase = SCEDC format (?) with no takeoff/azim info
                In example 2, we assume that we don’t have a good idea of the uncertainty in the azimuth and takeoff angles.
                A 1D ray-tracing routine is used to compute the takeoff angles for plausible velocity models 
                (the azimuth to a station is the same for any 1D model.) Each trial uses a different combination of
                (gradient) velocity model and perturbed source depth to compute a set of takeoff angles.
                Because the azimuth and takeoff angle to each station is being computed, a file with station names and
                locations is needed (file: scsn.stations)

HashDriver3.f - reads in north2.phase = SCEDC format (?) with no takeoff/azim info & reads in S/P_amp observations from north3.amp
                In example 3, S/P amplitude ratios are also used. An additional input file with the P and S amplitudes,
                as well as the noise levels prior to the arrival of the P and S waves, is included (file: north3.amp.)
                We apply a station correction of a constant shift in log10(S/P) (Hardebeck & Shearer, 2003), and
                therefore need to input a file with these corrections (file: north3.statcor)

                Two additional control parameters are also needed. We limit the S/P observations to those waveforms with
                signal to noise ratio of at least ratmin. In order to estimate the allowed log10(S/P) misfit
                (qextra and qtotal, above), an estimate of the log10(S/P) uncertainty is given by qbadfrac. For this example,
                it is very important that each event have a unique ID number, as the ID is used to match the P-polarity and
                S/P amplitudes from different files.

"""

from obspy.imaging.beachball import aux_plane

class Examples():

    def __init__(self, **kwargs):
        events = None
        events_sp = None
        outputs = None
        n = None

    def __repr__(self):
        return '{0} n={1}'.format(self.__class__.__name__, self.n)

    def printout(self, maxrows=None, filename=None):
        return printout(self.outputs, maxrows, filename)

    def run_example(self, n, show_focal_plots=False):

        sample_dir = os.path.join(base_dir, 'sample')

        if n not in [1,2,3]:
            logger.warning("hashwrapper.Examples.run_example(n=%s) --> n must be in [1,2,3]" % n)
            return None

        self.n = n

        example_file = 'example%d.toml' % n
        toml_file = os.path.join(sample_dir, example_file)
        # Chek if toml_file exists

        settings = toml.load(toml_file)
        plfile = os.path.join(sample_dir, settings['input_files']['plfile'])
        fpfile = os.path.join(sample_dir, settings['input_files']['fpfile'])
        delmax = settings['hash_parameters']['delmax']

    # Example 1
        if n == 1:
            events = read_fpfit_file(fpfile=fpfile, plfile=plfile, delmax=delmax)
            outputs = calc_focal_mechanisms(events, settings['hash_parameters'], phase_format='FPFIT')

    # Example 2, 3
        elif n in [2, 3]:
            stfile = os.path.join(sample_dir, settings['input_files']['stfile'])
            events = read_fpfit_file(fpfile=fpfile, plfile=plfile, stfile=stfile, delmax=delmax,
                                     fpfit_phase_format=2)

            new_list = []
            for fname in settings['hash_parameters']['velocity_models']:
                new_list.append(os.path.join(sample_dir, fname))

            settings['hash_parameters']['velocity_models'] = new_list

            events_sp = None
            if n == 3:
                ampfile = os.path.join(sample_dir, settings['input_files']['ampfile'])
                corfile = os.path.join(sample_dir, settings['input_files']['corfile'])
                ratmin = settings['hash_parameters']['ratmin']
                events_sp = read_sp_ratios(ampfile=ampfile,
                                           corfile=corfile,
                                           fpfit_events=events,
                                           stfile=stfile,
                                           ratmin=ratmin,
                                           delmax=delmax)

            outputs = calc_focal_mechanisms(events, settings['hash_parameters'],
                                            phase_format='SCEDC',
                                            events_sp=events_sp)

        self.events = events
        self.outputs = outputs
        #return outputs


def printout(outputs, maxrows=None, filename=None):

    usage = "Usage: printout([maxrows : default=10, use 'all' to get all rows], [filename : default=stdout]"

    if maxrows is None:
        if filename is None:
            maxrows = 10   # first 10 rows to stdout
        else:              # all to file
            maxrows = len(outputs)
    elif maxrows=='all':
        maxrows = len(outputs)

    try:
        n = int(maxrows)
    except ValueError as e:
        print(usage)
        return

    if filename is None:
        filename = "/dev/tty"
    with open(filename, 'a') as f:
        f.write("%6s %66s %3s %3s %3s [%3s %3s %3s]  %s %s %5s %5s %6s %5s %5s %5s    Q  useAmp\n" % \
                ("icusp", "event_info", "str", "dip", "rke",
                 "str2", "dip2", "rke2",
                 "agap", "pgap", "stdr", "msf", "rms_ang", "npol", "nspr", "mavg"))

        for out_dict in outputs[:n]:
            icusp = out_dict['icusp']
            #out_dict['event_info'] = event['event_info']
            strike = out_dict['strike']
            dip = out_dict['dip']
            rake = out_dict['rake']
            s2,d2,r2 = aux_plane(strike,dip,rake)
            magap = out_dict['azim_gap']
            mpgap = out_dict['theta_gap']

            stdr = int(100. * out_dict['stdr'])
            misfit = out_dict['misfit']
            use_amplitudes = out_dict['use_amplitudes']
            useAmp = ''
            nspr = 0
            mavg = "N/A"
            if use_amplitudes:
                S_P_log10_misfit = out_dict['S_P_log10_misfit']
                nspr = out_dict['number_SP_amp_ratios']
                useAmp = 'T'
                mavg = "%5.2f" % out_dict['mavg']

            #nppl = out_dict['number_P_polarities']
            #npol = out_dict['station_polarity_count']
            npol = out_dict['number_P_polarities']

            quality = out_dict['quality'].decode('UTF-8')
            rms_angle = out_dict['rms_angle']
            #f.write("%d %3d %3d %3d %.1f %.1f %5.3f %5.3f %6.2f [Q: %s] %1s\n" % \
            f.write("%d %s %3d %3d %3d  [%3d %3d %3d]    %.1f %.1f %5d %5.3f %6.2f %5d %5d %5s [Q: %s] %1s\n" % \
                (icusp, out_dict['event_info'], int(strike), int(dip), int(rake),
                int(s2), int(d2), int(r2),
                magap, mpgap, stdr, misfit, rms_angle, npol, nspr, mavg, quality, useAmp) )

    return


from obspy.core.event.source import FocalMechanism, NodalPlanes, NodalPlane
from obspy.imaging.beachball import aux_plane
from obspy.core.event.base import Comment

import sys, os
def main():

    #pathname = os.path.dirname(sys.argv[0])
    #print('path =', pathname)
    #print('full path =', os.path.abspath(pathname))
    n = 1
    n = 2
    n = 3
    exs = Examples()
    exs.run_example(n)
    outfile = 'mth_outputfile.%d' % n
    #exs.printout()
    exs.printout(filename=outfile)
    exit()

   # if outputs is not None:
        #for out in outputs:
            #print(out['strike'], out['dip'], out['rake'])
            #s,d,r = aux_plane(out['strike'], out['dip'], out['rake'])
            #print(s,d,r)
            #exit()
            #print(out)


    for i,out in enumerate(outputs):
        p1 = NodalPlane(strike=out['strike'], dip=out['dip'], rake=out['rake'])
        s,d,r = aux_plane(out['strike'], out['dip'], out['rake'])
        p2 = NodalPlane(strike=s, dip=d, rake=r)

        fc = FocalMechanism(nodal_planes = NodalPlanes(nodal_plane_1=p1, nodal_plane_2=p2),
                            azimuthal_gap = out['azim_gap'],
                            station_polarity_count = out['station_polarity_count'],
                            station_distribution_ratio = out['stdr'],
                            misfit = out['misfit'],
                            evaluation_mode = 'automatic',
                            evaluation_status = 'preliminary',
                            comments = [Comment(text="HASH v1.2 Quality=[%s]" % out['quality'])]
                           )

if __name__ == '__main__':
    main()

