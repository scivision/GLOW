#!/usr/bin/env python
from numpy.distutils.core import Extension, setup
from glob import glob
from os.path import join

src = ('basic.f90 cglow.f90 glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90'
       'etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90'
        'qback.f90 rcolum.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f'
        'nrlmsise00.f').split(' ')

ext = Extension(name='glow', sources=src,
                f2py_options=['--quiet'],)

iridata = glob(join('data/iri90', '*.asc'))
fortdata = glob(join('data', '*.dat'))

setup(ext_modules=[ext],
      data_files=[('glow/data', fortdata),
                  ('glow/data/iri90', iridata)]
      )
