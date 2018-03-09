#!/usr/bin/env python
install_requires = ['numpy',
                    'sciencedates']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import find_packages # enables develop
from numpy.distutils.core import Extension, setup

src = 'cglow.f90 glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 rout.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f'.split(' ')

ext = Extension(name='glow', sources=src, 
                 f2py_options=['--quiet'])


setup(name='glowiono',
      packages=find_packages(),
      version='0.1.0',
      author=['Stan Solomon','Liam Kilcommons','Michael Hirsch, Ph.D.'],
      url = 'https://github.com/scivision/glow',
      description='GLobal airglOW model',
      long_description=open('README.rst').read(),
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 4 - Beta',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python :: 3.6',
      ],
      ext_modules=[ ext ],
      data_files=[],
      python_requires=">=3.6",
      tests_require=tests_require,
      extras_require={'tests':tests_require},
    )
