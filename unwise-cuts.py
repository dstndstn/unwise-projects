from __future__ import print_function
import os
import numpy as np

from astrometry.util.fits import *

def brandt():
    T = fits_table('sdss-dr10d-tiles.fits')
    print(len(T), 'tiles')
    print('unique Decs', np.unique(T.dec))
    T.cut((T.dec > -2.) * (T.dec < 2.))
    print(len(T), 'tiles in Dec range')
    #print('unique RAs', np.unique(T.ra))
    T.cut(np.logical_or(T.ra < 62., T.ra > (360-52)))
    print(len(T), 'tiles in RA & Dec range')

    WW = []
    for t in T:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % t.coadd_id
        W = fits_table(fn)
        print(len(W), 'from', fn)
        W.cut((W.dec > -1.25) * (W.dec < 1.25))
        W.cut(np.logical_or(W.dec < 60., W.dec > (350-50)))
        print(len(W), 'in range')
        WW.append(W)
    W = merge_tables(WW)
    print(len(W), 'total')
    W.writeto('brandt.fits')



if __name__ == '__main__':
    brandt()

#sdss-dr10d-phot/phot-
