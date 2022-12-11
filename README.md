
# hashwrap

Wrap of HASHv1.2 Fortran library (see Hardebeck & Shearer, 2002)

## Getting Started

The HASH f77 subroutines use fixed dimension arrays for efficiency.
These array dimensions must be known at compile time, and hence when this
module is installed.  The dimensions must also be shared with the
python side in order to ensure the numpy arrays are correctly sized
to match their fortran counterparts.

The relevant parameters are:
 - npick0: maximum number of observations (polarities + amp ratios) for an
   event
 - nmc0:   maximum number of trials (ray takeoff angles and azims are
   perturbed for each trial)
 - nmax0:  maximum number of acceptable mechanisms to output for an event
 - dang0:  minimum angular spacing (degrees) for mechanism grid search
 - ncoor:  maximum number of gridpoints

If you do nothing, these values will automatically get set to the
default values indicated in hashwrap/fortran_arrays.toml:

parameter(npick0=500, nmc0=500, nmax0=500, dang0=5.0, ncoor=31032)

If you clone the repo (see below), then you can edit fortran_arrays.toml
directly to change the array sizes before you install.

However, if you install directly without cloning, then you won't
ever see the fortran_arrays.toml file.
You can still modify them by setting env vars *before* running pip install, e.g.,

    >export npick0=500
    >export nmc0=500
    ...
    >export ncoor=31032

The env vars will override the fortran_arrays.toml values.

### Installing

The easiest way is to install directly from either pypi or github:

If you wish to change the default max array sizes (see above),
you should set these vals in shell environment vars first.

#### Install toml

This needs to be installed before you install hashwrap

    >pip install toml

#### Install directly from pypi

    >pip install hashwrap

#### Or, Install from github

    >pip install git+https://github.com/mikehagerty/hashwrap.git

Alternatively you can clone the repo and install it:

If you wish to change the default max array sizes (see above),
you can either set these vals in shell environment vars first,
or change the vals directly in hashwrap/fortran_arrays.toml.

    >git clone https://github.com/mikehagerty/hashwrap.git
    >cd hashwrap
    >pip install .


## Using hashwrap

In general, hashwrap exists to support other modules (e.g., the
focalmech app in the aqms-pdl repo) and there is no need to
interact directly with it.

Nevertheless, you may want to check out the shared Fortran lib
and/or run the equivalent of the hash_driver{1,2,3}.f tests
that come with HASH v1.2.

### Checking everything is kosher

While you won't generally want/need to interact directly with the
compiled Fortran shared object library, you can see what it contains:

    > python
    >>> from hashwrap import libhashpy
    >>> print(libhashpy.__doc__)
```
This module 'libhashpy' is auto-generated with f2py (version:1.23.1).
Functions:
    nf,strike,dip,rake,faults,slips = focalamp_mc(p_azi_mc,p_the_mc,sp_amp,p_pol,nmc,dang,maxout,nextra,ntotal,qextra,qtotal,npsta=shape(sp_amp, 0))
    mfrac,mavg,stdr = get_misf_amp(p_azi_mc,p_the_mc,sp_ratio,p_pol,str_avg,dip_avg,rak_avg,npol=shape(p_azi_mc, 0))
    nf,strike,dip,rake,faults,slips = focalmc(p_azi_mc,p_the_mc,p_pol,p_qual,nmc,dang,maxout,nextra,ntotal,npsta=shape(p_pol, 0))
    nsltn,str_avg,dip_avg,rak_avg,prob,rms_diff = mech_prob(norm1in,norm2in,cangle,prob_max,nf=shape(norm1in, 1))
    norm1_avg,norm2_avg = mech_avg(norm1,norm2,nf=shape(norm1, 1))
    rota = mech_rot(norm1,norm2,slip1,slip2)
    v3 = cross(v1,v2)
    to_car(the,phi,r,x,y,z)
    fpcoor(strike,dip,rake,fnorm,slip,idir)
    fran = ran_norm()
    mfrac,stdr = get_misf(p_azi_mc,p_the_mc,p_pol,p_qual,str_avg,dip_avg,rak_avg,npol=shape(p_azi_mc, 0))
    magap,mpgap = get_gap(p_azi_mc,p_the_mc,npol=shape(p_azi_mc, 0))
    sort(ra,n=shape(ra, 0))
    ntab = mk_table(ntab)
    tt,iflag = get_tts(ip,del,qdep)
    layertrace(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
    flat,flon,felev = getstat_tri(stlfile,snam,scom,snet)
    stpol = check_pol(polfile,snam,evyr,evmon,evdy,evhr)
    qcor = get_cor(stlfile,snam,scom,snet)
    ntab = mk_table_add(ind,vmodel)
COMMON blocks:
  /angtable/ table(101,14,10),delttab(101),deptab(14),ndel,ndep
```


### Run the hash_driver{1,2,3} examples

    >>> from hashwrap import hashwrapper
    >>> ex = hashwrapper.Examples()        # Init the Examples class
    >>> ex.run_example(1)                  # Equivalent to hash_driver1.f

You can printout the results in a format very similar to what
hash_driver1.f prints out, for direct comparison.
By default Examples.printout() only prints the first 10 events:

    >>> ex.printout()

```
 icusp                                                         event_info str dip rke [str2 dip2 rke2]  agap pgap  stdr   msf rms_ang  npol  nspr  mavg    Q  useAmp
3143312 1994-01-21 11:04 <34.242, -118.618> h:18.130000 [cusp id:3143312] 134  48 143  [250  63  47]    84.0 17.0    66 0.090  24.51    30     0   N/A [Q: A]
3145744 1994-01-25 10:05 <34.241, -118.621> h:18.540000 [cusp id:3145744] 280  47  54  [147  53 122]    44.0 15.0    63 0.172  31.77    33     0   N/A [Q: B]
3146815 1994-01-28 07:44 <34.239, -118.621> h:18.960000 [cusp id:3146815] 137  45 132  [264  58  55]    32.0 11.0    50 0.130  15.44    73     0   N/A [Q: A]
3146907 1994-01-28 15:47 <34.238, -118.611> h:17.790000 [cusp id:3146907] 297  39  99  [105  51  82]    67.0 20.0    61 0.069  33.02    23     0   N/A [Q: B]
3147167 1994-01-29 07:52 <34.238, -118.625> h:19.200000 [cusp id:3147167] 141  54 107  [292  38  67]    36.0 16.0    57 0.097  21.14    55     0   N/A [Q: A]
3148047 1994-01-31 07:13 <34.242, -118.615> h:17.920000 [cusp id:3148047] 139  49 111  [288  45  66]    35.0 15.0    61 0.070  28.07    39     0   N/A [Q: B]
3149674 1994-02-05 07:11 <34.240, -118.617> h:18.200000 [cusp id:3149674] 131  47 112  [280  47  67]    47.0 15.0    59 0.146  25.57    50     0   N/A [Q: B]
3150936 1994-02-10 04:19 <34.242, -118.620> h:18.310000 [cusp id:3150936] 143  59 130  [264  48  42]    43.0 15.0    49 0.053  18.96    57     0   N/A [Q: B]
3150947 1994-02-10 04:46 <34.241, -118.620> h:18.450000 [cusp id:3150947] 144  54 130  [268  51  47]    35.0 15.0    54 0.084  19.93    50     0   N/A [Q: A]
3151649 1994-02-13 09:04 <34.238, -118.613> h:17.740000 [cusp id:3151649] 129  47 110  [279  46  68]    42.0 16.0    54 0.101  24.09    33     0   N/A [Q: A]
```

By default, when sending the hash_driver output to a file,
all events are included.
To run hash_driver3.f and send all rows to a file, do:

    >>> ex.run_example(3)
    >>> ex.printout('example3.out')


## HASH references

  - Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first- motion focal mechanisms, Bulletin of the Seismological Society of America, 92, 2264-2276, 2002.
  - Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the Seismological Society of America, 93, 2434-2444, 2003.


## Acknowledgments

* This code began as a fork of Mark Williams' HASHpy repo.  It has
  diverged quite a bit but still retains several of his ideas and
utility functions.

