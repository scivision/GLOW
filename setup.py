#!/usr/bin/env python
install_requires = ['numpy','xarray',
                    'sciencedates']
tests_require=['pytest','nose','coveralls']
# %%
from setuptools import find_packages # enables develop
from numpy.distutils.core import Extension, setup
from glob import glob
from os.path import join

src = 'basic.f90 cglow.f90 glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f'.split(' ')

ext = Extension(name='glow', sources=src,
                 f2py_options=['--quiet'],)

iridata = glob(join('data/iri90','*.asc'))
fortdata = glob(join('data','*.dat'))

setup(name='glowiono',
      packages=find_packages(),
      version='0.2.0',
      author=['Stan Solomon','Liam Kilcommons','Michael Hirsch, Ph.D.'],
      url = 'https://github.com/scivision/glow',
      description='GLobal airglOW model 0.981',
      long_description=open('README.rst').read(),
      classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
       'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3.6',
      'Topic :: Scientific/Engineering :: Atmospheric Science',

      ],
      ext_modules=[ ext ],
      data_files=[('glow/data',fortdata),
		           ('glow/data/iri90',iridata)],
      install_requires=install_requires,
      python_requires=">=3.6",
      tests_require=tests_require,
      extras_require={'tests':tests_require,
		      'plot':['matplotlib']},
    )

