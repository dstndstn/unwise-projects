from __future__ import print_function
import os
import numpy as np

from astrometry.util.fits import *
from astrometry.sdss.common import cas_flags, photo_flags1_map, photo_flags2_map
from astrometry.libkd.spherematch import match_radec
#tree_build_radec, tree_free, tree_search_radec, trees_match

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


def nunez():
    tag = 'nunez'

    '''
    Read all SDSS photoObj primary objects, pulling out select
    columns and applying cuts.
    '''

    W = fits_table('window_flist-cut.fits', columns=['run','camcol','field'])
    print('Read', len(W), 'fields')

    istart = 0
    k = 0
    ilast = 0

    # This is for resuming a crashed run:
    # look in the log file for the last nunez-###.fits file written:
    #
    #   312014 of 544402 : Run,camcol,field 4849 3 849
    #   ...
    #   Wrote 143434 to out/nunez-311.fits
    #
    istart = 312014
    k = 312
    ilast = istart

    TT = []
    for i,(run,camcol,field) in enumerate(zip(W.run, W.camcol, W.field)):
        if i < istart:
            continue
        print()
        print(i+1, 'of', len(W), ': Run,camcol,field', run,camcol,field)
        fn = os.path.join(os.environ['BOSS_PHOTOOBJ'], '301', '%i'%run,
                          '%i'%camcol, 'photoObj-%06i-%i-%04i.fits' % (run,camcol,field))
    
        T = fits_table(fn, columns=['objc_type', 'clean', 'flags', 'flags2',
                                    'psfmagerr', 'resolve_status',])
        if T is None:
            continue
        T.index = np.arange(len(T))
        T.flags_r = T.flags[:,2]
        T.flags_i = T.flags[:,3]
        T.flags2_r = T.flags2[:,2]
        T.flags2_i = T.flags2[:,3]
        T.cut((T.resolve_status & 256) > 0)
        T.cut(T.objc_type == 3)
        T.cut(T.clean == 1)
    
        def casToPhoto(casval):
            flag1val = 0
            flag2val = 0
            casvals = 0
            for k,v in cas_flags.items():
                if v & casval:
                    print('CAS flag', k)
                    casvals |= v
                    if k in photo_flags1_map:
                        flag1val |= photo_flags1_map[k]
                    elif k in photo_flags2_map:
                        flag2val |= photo_flags2_map[k]
                    else:
                        print('Flag not found:', k)
                        assert(False)
            assert(casvals == casval)
            print('Flag values 0x%x 0x%x' % (flag1val,flag2val))
            return flag1val,flag2val
    
        f1,f2 = casToPhoto(0x10000000)
        assert(f2 == 0)
        T.cut((T.flags_r & f1) > 0)
        T.cut((T.flags_i & f1) > 0)
    
        f1,f2 = casToPhoto(0x800a0)
        assert(f2 == 0)
        T.cut((T.flags_r & f1) == 0)
        T.cut((T.flags_i & f1) == 0)
    
        f1,f2 = casToPhoto(0x400000000000)
        assert(f1 == 0)
        T.cut(np.logical_or((T.flags2_r & f2) == 0, T.psfmagerr[:,2] <= 0.2))
        T.cut(np.logical_or((T.flags2_i & f2) == 0, T.psfmagerr[:,3] <= 0.2))
    
        print(len(T), 'pass cut')
        if len(T) == 0:
            continue
    
        # Re-read SDSS photoObjs and grab more columns.
        inds = T.index
        T = fits_table(fn, rows=inds, columns=['ra','dec','raerr','decerr',
                                            'cmodelmag', 'cmodelmagerr',
                                            'psfmag','psfmagerr', 'flags', 'flags2', 'objid'])
    
        unwdir = '/project/projectdirs/cosmo/data/unwise/unwise-phot/sdss-collab/sdss-dr13-pobj'
        fn = os.path.join(unwdir, '%i'%run, '%i'%camcol, 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field))
        U = fits_table(fn, rows=inds, columns=['treated_as_pointsource', 'pointsource',
                                               'w1_mag', 'w1_mag_err',
                                               'w2_mag', 'w2_mag_err',
                                               'w3_mag', 'w3_mag_err',
                                               'w4_mag', 'w4_mag_err',
                                               'has_wise_phot', 'objid'])
    
        print('has_wise_phot:', np.unique(U.has_wise_phot))
        # *almost* all the time, has_wise_phot=True; exception is
        # sdss-dr13-pobj/4135/4/photoWiseForced-004135-4-0169.fits
        # but I haven't investigated.
        T.cut(U.has_wise_phot)
        U.cut(U.has_wise_phot)
    
        assert(np.all(T.objid == U.objid))
    
        U.delete_column('objid')
        T.add_columns_from(U)
    
        TT.append(T.copy())
        print('Total of', sum([len(x) for x in TT]), 'sources')
        if i - ilast >= 1000:
            ilast = i
            T = merge_tables(TT)
            print('Total of', len(T), 'pass cut')
            N = 10000
            if len(T) > N:
                fn = 'out/%s-%03i.fits' % (tag,k)
                k += 1
                T.writeto(fn)
                print('Wrote', len(T), 'to', fn)
                TT = []
            else:
                T.writeto('out/%s-x.fits' % tag)
                TT = [T]
    
    T = merge_tables(TT)
    fn = 'out/%s-%03i.fits' % (tag,k)
    k += 1
    T.writeto(fn)
    print('Wrote', len(T), 'to', fn)


def manga():
    T = fits_table('sdss-dr13-atlas.fits')
    print(len(T), 'tiles')

    # https://data.sdss.org/sas/mangawork/manga/target/v1_2_12/MaNGA_targets_extNSA.fits
    M = fits_table('MaNGA_targets_extNSA.fits')
    print(len(M), 'MaNGA')

    I,J,d = match_radec(M.ra, M.dec, T.ra, T.dec, 1.2)
    keep = np.zeros(len(M), bool)
    keep[I] = True
    assert(np.all(keep))
    keep = np.zeros(len(T), bool)
    keep[J] = True
    T.cut(keep)
    print(len(T), 'tiles near MaNGA targets')

    # Attempted to find nearest...
    #mw = [None for i in range(len(M))]
    #md = np.zeros(len(M))
    #md[:] = 100.
    WW = []

    for t in T:
        fn = 'sdss-dr13-phot/phot-%s.fits' % t.coadd_id
        W = fits_table(fn)
        print(len(W), 'from', fn)
        if len(W) == 0:
            continue
        #I,J,d = match_radec(M.ra, M.dec, W.ra, W.dec, 1./3600., nearest=True)
        I,J,d = match_radec(M.ra, M.dec, W.ra, W.dec, 60./3600.)
        print(len(I), 'matched')
        if len(I) == 0:
            continue
        W.cut(J)
        WW.append(W)

        # nearest = np.flatnonzero(d < md[I])
        # print(len(nearest), 'nearest')
        # if len(nearest) == 0:
        #     continue
        # I = I[nearest]
        # J = J[nearest]
        # D = d[nearest]
        # for i,j,d in zip(I,J,D):
        #     md[i] = d
        #     mw[i] = W[j]
        # print('Now', np.sum([m is not None for m in mw]), 'of', len(M), 'matched')

    #W = merge_tables(mw)
    W = merge_tables(WW)
    print(len(W), 'total matches')
    W.writeto('manga-unwise.fits')

if __name__ == '__main__':
    #brandt()
    nunez()
    #manga()



