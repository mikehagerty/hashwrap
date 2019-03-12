#
# setup.py file for compiling HASH and installing hashpy
#
import sys
import os
from numpy.distutils.core import setup, Extension

# MTH: 2019
# This project was forked from Marc William's hashpy at:
# https//github.com/markcwill/hashpy
# I rewrote most of the python mods but kept Marc's
#  wrapping of the Fortran itself, hence I'm keeping his
#  name for the extension (libhashpy)

#--- libhashpy Fortran extension --------------------------------------------#
#
# Build extension from FORTRAN source of HASH subroutines
# (based on the fucntion numpy.f2py.f2py2e.run_compile)
#srcdir = os.path.join('hashpy', 'src')
srcdir = os.path.join('.', 'src')
srcf = ['fmamp_subs.f', 'fmech_subs.f', 'uncert_subs.f', 'util_subs.f',
        'pol_subs.f', 'vel_subs.f', 'station_subs.f', 'vel_subs2.f']
src_list = [os.path.join(srcdir, src) for src in srcf]
ext_args = {'sources': src_list}

requirements = [
    'numpy',
    'matplotlib',
]

### SETUP ####################################################################
s_args = {'name'         : 'hashwrap',
          'version'      : '0.0.1',
          'description'  : 'Routines for running HASH algorithms for earthquake focal mechanism',
          'author'       : 'Mike Hagerty',
          'author_email' : 'm.hagerty@isti.com',
          'url'          : 'https//github.com/mhagerty/hashwrap',
          'packages'     : ['hashwrap'],
          'install_requires' : requirements,
          #'package_data' : {'hashpy': ['src/*.inc','src/Makefile','data/*',
                                       #'scripts/*', 'src/*.f']},
          'add_data_dir' : 'data',
          'license': 'MIT',
          'include_package_data': True,
          'ext_modules'  : [Extension('hashwrap.libhashpy', **ext_args)],
}
##############################################################################


setup(**s_args)
