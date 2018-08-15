#!/usr/bin/env python
"""
default parameter values like those of Stan's fortran examples--yield rather similar output

Note, this simulation uses a specific input differential number flux spectrum
"""
from matplotlib.pyplot import show
import glowiono
from glowiono.plots import plotaurora
from argparse import ArgumentParser


def main():
    p = ArgumentParser(description="Stan Solomon's GLOW auroral model")
    p.add_argument('-t', '--simtime', help='yyyy-mm-ddTHH:MM:SSZ time of sim', default='2013-04-14T15:54Z')
    p.add_argument('-c', '--latlon', help='geodetic latitude/longitude (deg)', type=float, nargs=2, default=(65., -148.))
#    p.add_argument('-n','--nbins',help='Number energy bins in incident diff num flux',type=int,default=190) #hard-coded in cglow.h
    p.add_argument('-q', '--flux', help='overall incident flux [erg ...]', type=float, default=1.)
    p.add_argument('--e0', help='characteristic energy [eV]', type=float, default=1e3)
    p.add_argument('-m', '--makeplot', help='show to show plots, png to save pngs of plots', nargs='+', default=['show'])
    p.add_argument('-zlim', help='plot limits of altitude axis [km]', nargs=2, type=float)
    p = p.parse_args()

    params = {'time': p.simtime,
              'glat': p.latlon[0],
              'glon': p.latlon[1],
              'flux': p.flux,
              'E0': p.e0,
              'makeplot': p.makeplot,
              'zlim': p.zlim,
              'plotformat': 'png',
              }

    sim = E0aurora(params)
    print(sim)

    show()


def E0aurora(params: dict):

    sim = glowiono.runglowaurora(params)

    plotaurora(params, sim)

    return sim


if __name__ == '__main__':
    main()
