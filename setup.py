#!/usr/bin/env python
req = ['nose','numpy','matplotlib','seaborn']
pipreq=['sciencedates']
# %%
import pip
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception:
    pip.main(['install'] + req)
pip.main(['install'] + pipreq)
# %%
import setuptools # enables develop
from numpy.distutils.core import Extension, setup

src = 'cglow.f90 glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 rout.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f'.split(' ')

ext = Extension(name='glow', sources=src, 
                 f2py_options=['--quiet'])


setup(name='glowiono',
      packages=['glowiono'],
      version='1.0',
      author=['Stan Solomon','Liam Kilcommons','Michael Hirsch, Ph.D.'],
      url = 'https://github.com/NCAR/GLOW',
      description='GLobal airglOW model',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 5 - Production/Stable',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3',
      ],
      ext_modules=[ ext ],
      data_files=[],
    )
