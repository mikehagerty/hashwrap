#
# setup.py file for compiling HASH and installing hashpy
#
import sys
import os
# By importing setuptools *before* numpy.distutils.core
#  numpy will fork and use the (more modern) setuptools
import setuptools
from numpy.distutils.core import setup, Extension

# MTH: 2019
# This project was forked from Marc William's hashpy at:
# https//github.com/markcwill/hashpy
# I rewrote most of the python mods but kept Marc's
#  wrapping of the Fortran itself, hence I'm keeping his
#  name for the extension (libhashpy)

# MTH: This is a hack to set the fortran arrays at install time
#      The fortran and python array dimensions must match
#      The fortran array dimensions are read from inc files in src/
#      This hack will read in hashwrap/fortran_arrays.toml
#      and use it to create the hashwrap/src/*inc files
#      Later, the python code will read in these same src/*inc files

#      So only set the fortran array dimensions in hashwrap/fortran_arrays.toml
#        and then re-run this setup!

def dict_to_string(inc_dict):
    string = "      parameter("
    for k,v in inc_dict.items():
        string += "%s=%s," % (k,v)
    string = string[:-1] + ")"
    return string

def write_inc_files(inc_dict, filename):
    foo = inc_dict[ list(inc_dict.keys())[0]]
    if isinstance(foo, dict):
        string=""
        for k,v in inc_dict.items():
            string += dict_to_string(v) + "\n"
    else:
        string = dict_to_string(inc_dict)

    with open(filename, 'w') as f:
        f.write(string)

    return

import toml
settings = toml.load('hashwrap/fortran_arrays.toml')['fortran_array_dimensions']
write_inc_files(settings['param_inc'], 'hashwrap/src/param.inc')
write_inc_files(settings['rot_inc'], 'hashwrap/src/rot.inc')
write_inc_files(settings['vel_inc'], 'hashwrap/src/vel.inc')


#--- libhashpy Fortran extension --------------------------------------------#
#
# Build extension from FORTRAN source of HASH subroutines
# (based on the function numpy.f2py.f2py2e.run_compile)
srcdir = os.path.join('hashwrap', 'src')
srcf = ['fmamp_subs.f', 'fmech_subs.f', 'uncert_subs.f', 'util_subs.f',
        'pol_subs.f', 'vel_subs.f', 'station_subs.f', 'vel_subs2.f']
src_list = [os.path.join(srcdir, src) for src in srcf]
ext_args = {'sources': src_list}

requirements = [
    'numpy',
    'matplotlib',
    'toml',
]

setup(name='hashwrap',
      packages=['hashwrap'],
      #packages=setuptools.find_packages(),
      version='0.0.2',
      description='Routines for running HASH algorithms for earthquake focal mechanism',
      author='Mike Hagerty',
      author_email='m.hagerty@isti.com',
      url='https://github.com/mikehagerty/hashwrap',
      # src code will NOT be copied to site-packages/hashwrap without this:
      package_data={'hashwrap': ['*.toml', 'src/*', 'docs/*', 'sample/*', 'sample/data/*']},
      include_package_data=True,
      install_requires=requirements,
      ext_modules=[Extension('hashwrap.libhashpy', **ext_args)],
      classifiers=(
          "Programming Language :: Python :: 3.7",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
      ),
    )
