from __future__ import print_function
import os
import numpy as np

from astrometry.util.fits import *
from astrometry.sdss.common import cas_flags, photo_flags1_map, photo_flags2_map
from astrometry.libkd.spherematch import match_radec, tree_build_radec, tree_free, tree_search_radec, trees_match

def gupta():
    T = fits_table('sdss-dr10d-tiles.fits')
    O = fits_table('gupta.fits')
    #O.rename('dr8objid', 'objid')
    Omap = dict([(objid,i) for i,objid in enumerate(O.objid)])
    PP = []
    J = []
    for tile in T.coadd_id:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % tile
        if not os.path.exists(fn):
            print('File not found:', fn)
            continue
        P = fits_table(fn, columns=['objid'])
        I = []
        for i,objid in enumerate(P.objid):
            objid = int(objid)
            try:
                j = Omap[objid]
            except KeyError:
                #print('Objid not found:', objid)
                continue
            I.append(i)
            J.append(j)
        print(len(I), 'in', tile)
        if len(I) == 0:
            continue
        P = fits_table(fn, rows=np.array(I))
        PP.append(P)
    P = merge_tables(PP)
    P.matched = np.ones(len(P), bool)
    del PP
    N = len(P)
    print('Total of', N, 'read')
    J = np.array(J)
    Jinv = np.empty(len(O), int)
    # no-match
    Jinv[:] = N
    Jinv[J] = np.arange(N)

    jj = np.flatnonzero(Jinv < N)
    print('First match:', jj[0])

    empty = fits_table()
    empty.matched = np.zeros(1, bool)
    Pplus = merge_tables([P, empty], columns='fillzero')
    print('Pplus:', len(Pplus))
    Pplus = Pplus[Jinv]
    #O.rename('objid', 'dr8objid')
    Pplus.add_columns_from(O)
    Pplus.writeto('/global/cscratch1/sd/dstn/gupta.fits')


def fengshuai():
    T = fits_table('sdss-dr10d-tiles.fits')
    O = fits_table('maingalDR7_fengshuai0210.fit')
    O.rename('dr8objid', 'objid')

    Omap = dict([(objid,i) for i,objid in enumerate(O.objid)])

    PP = []
    J = []
    for tile in T.coadd_id:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % tile
        if not os.path.exists(fn):
            print('File not found:', fn)
            continue
        P = fits_table(fn, columns=['objid'])
        I = []
        for i,objid in enumerate(P.objid):
            objid = int(objid)
            try:
                j = Omap[objid]
            except KeyError:
                #print('Objid not found:', objid)
                continue
            I.append(i)
            J.append(j)
        print(len(I), 'in', tile)
        if len(I) == 0:
            continue
        P = fits_table(fn, rows=np.array(I))
        PP.append(P)
    P = merge_tables(PP)
    P.matched = np.ones(len(P), bool)
    del PP
    N = len(P)
    print('Total of', N, 'read')
    J = np.array(J)
    Jinv = np.empty(len(O), int)
    # no-match
    Jinv[:] = N
    Jinv[J] = np.arange(N)

    jj = np.flatnonzero(Jinv < N)
    print('First match:', jj[0])

    empty = fits_table()
    empty.matched = np.zeros(1, bool)
    Pplus = merge_tables([P, empty], columns='fillzero')
    print('Pplus:', len(Pplus))
    Pplus = Pplus[Jinv]
    O.rename('objid', 'dr8objid')
    Pplus.add_columns_from(O)
    Pplus.writeto('/global/cscratch1/sd/dstn/fengshuai.fits')
    

def leslie():
    T = fits_table('sdss-dr10d-tiles.fits')
    #O = fits_table('dr12_objIDs_leslie.fits')
    O = fits_table('leslie_SDSS13_objids.fits')
    #O.rename('objid_1', 'objid')
    
    Omap = dict([(objid,i) for i,objid in enumerate(O.objid)])

    PP = []
    J = []
    for tile in T.coadd_id:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % tile
        if not os.path.exists(fn):
            print('File not found:', fn)
            continue
        P = fits_table(fn, columns=['objid'])
        I = []
        for i,objid in enumerate(P.objid):
            objid = int(objid)
            try:
                j = Omap[objid]
            except KeyError:
                #print('Objid not found:', objid)
                continue
            I.append(i)
            J.append(j)
        print(len(I), 'in', tile)
        if len(I) == 0:
            continue
        P = fits_table(fn, rows=np.array(I))
        PP.append(P)
    P = merge_tables(PP)
    P.matched = np.ones(len(P), bool)
    del PP
    N = len(P)
    print('Total of', N, 'read')
    J = np.array(J)
    Jinv = np.empty(len(O), int)
    # no-match
    Jinv[:] = N
    Jinv[J] = np.arange(N)

    jj = np.flatnonzero(Jinv < N)
    print('First match:', jj[0])

    empty = fits_table()
    empty.matched = np.zeros(1, bool)
    Pplus = merge_tables([P, empty], columns='fillzero')
    print('Pplus:', len(Pplus))
    Pplus = Pplus[Jinv]
    O.rename('objid', 'objid_orig')
    Pplus.add_columns_from(O)
    #Pplus.writeto('/scratch1/scratchdirs/dstn/leslie.fits')
    Pplus.writeto('/global/cscratch1/sd/dstn/leslie2.fits')


def ogrady():
    T = fits_table('sdss-dr10d-tiles.fits')

    O = fits_table('ogrady_dr9_objects.fits')
    #O.rename('radeg', 'ra')
    #O.rename('dedeg', 'dec')
    #okd = tree_build_radec(O.ra, O.dec)

    #print('Objid', O.objid[0], 'RA,Dec', O.ra[0], O.dec[0])
    #I,J,d = match_radec(np.array([O.ra[0]]), np.array([O.dec[0]]), T.ra, T.dec, 1.,
    #                    nearest=True)
    #tile = T.coadd_id[J[0]]
    #print('Tile', tile)

    Omap = dict([(objid,i) for i,objid in enumerate(O.objid)])

    PP = []
    J = []
    for tile in T.coadd_id: #[:100]:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % tile
        if not os.path.exists(fn):
            print('File not found:', fn)
            continue
        P = fits_table(fn, columns=['objid'])
        I = []
        for i,objid in enumerate(P.objid):
            objid = int(objid)
            try:
                j = Omap[objid]
            except KeyError:
                continue
            I.append(i)
            J.append(j)
        print(len(I), 'in', tile)
        if len(I) == 0:
            continue
        P = fits_table(fn, rows=np.array(I))
        PP.append(P)
    P = merge_tables(PP)
    P.matched = np.ones(len(P), bool)
    del PP
    N = len(P)
    print('Total of', N, 'read')
    J = np.array(J)
    Jinv = np.empty(len(O), int)
    # no-match
    Jinv[:] = N
    Jinv[J] = np.arange(N)

    jj = np.flatnonzero(Jinv < N)
    print('First match:', jj[0])

    empty = fits_table()
    empty.matched = np.zeros(1, bool)
    Pplus = merge_tables([P, empty], columns='fillzero')
    print('Pplus:', len(Pplus))
    Pplus = Pplus[Jinv]
    O.rename('objid', 'objid_anna')
    Pplus.add_columns_from(O)
    Pplus.writeto('/scratch1/scratchdirs/dstn/ogrady.fits')


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
    #istart = 312014
    #k = 312
    #ilast = istart

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
                                               'modelmag', 'modelmagerr',
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



def all_matches_near(catfn, radius, atlasfn, photdir, outfn, projname, **kwargs):
    '''
    radius in arcsec
    '''
    T = fits_table(atlasfn)
    print(len(T), 'tiles')

    M = fits_table(catfn)
    print(len(M), 'targets')

    M.set('%s_index' % projname, np.arange(len(M)))

    I,J,d = match_radec(M.ra, M.dec, T.ra, T.dec, 1.2)
    #keep = np.zeros(len(M), bool)
    #keep[I] = True
    #assert(np.all(keep))
    keep = np.zeros(len(T), bool)
    keep[J] = True
    T.cut(keep)
    print(len(T), 'tiles near targets')

    WW = []
    for t in T:
        fn = os.path.join(photdir, 'phot-%s.fits' % t.coadd_id)
        if not os.path.exists(fn):
            print('Does not exist:', fn)
            continue
        W = fits_table(fn)
        print(len(W), 'from', fn)
        if len(W) == 0:
            continue

        I,J,d = match_radec(M.ra, M.dec, W.ra, W.dec, radius/3600., **kwargs)
        print(len(I), 'matched')
        if len(I) == 0:
            continue
        W.cut(J)
        MI = M[I]
        MI.rename('ra',  'ra_%s'  % projname)
        MI.rename('dec', 'dec_%s' % projname)
        W.add_columns_from(MI)
        WW.append(W)

    W = merge_tables(WW)
    print(len(W), 'total matches')
    W = W[np.argsort(W.get('%s_index' % projname))]

    print(len(np.unique(W.objid)), 'unique objids')

    W.writeto(outfn)


def manga():
    # https://data.sdss.org/sas/mangawork/manga/target/v1_2_12/MaNGA_targets_extNSA.fits
    all_matches_near('MaNGA_targets_extNSA.fits', 60., 'sdss-dr13-atlas.fits',
                    'sdss-dr13-phot', 'manga-unwise.fits', 'manga')


def sien():
    all_matches_near('ozdes-unwise.fits', 3., 'sdss-dr10d-tiles.fits',
                    'sdss-dr10d-phot', 'ozdes-unwise-matches.fits', 'ozdes')

def hoyle():
    all_matches_near('target_list_SDSS_DR7.fits', 1., 'sdss-dr10d-tiles.fits',
                    'sdss-dr10d-phot', 'hoyle-matches.fits', 'hoyle', nearest=True)

def amaral():
    all_matches_near('farnes.fits', 3., 'sdss-dr10d-tiles.fits',
                     'sdss-dr10d-phot', 'farnes-nearest-matches.fits', 'farnes-nearest',
                     nearest=True)
    all_matches_near('farnes.fits', 3., 'sdss-dr10d-tiles.fits',
                     'sdss-dr10d-phot', 'farnes-all-matches.fits', 'farnes-all')

    all_matches_near('taylor.fits', 3., 'sdss-dr10d-tiles.fits',
                     'sdss-dr10d-phot', 'taylor-nearest-matches.fits', 'taylor-nearest',
                     nearest=True)
    all_matches_near('taylor.fits', 3., 'sdss-dr10d-tiles.fits',
                     'sdss-dr10d-phot', 'taylor-all-matches.fits', 'taylor-all')

def salvato():
    '''
    > fitscopy RADEC_forDustin.fits+1"[col ra=xray_ra; dec=xray_dec]" RADEC_forDustin_x.fits
    '''
    all_matches_near('RADEC_forDustin_x.fits', 30., 'sdss-dr10d-tiles.fits',
                     'sdss-dr10d-phot', 'salvato-matches.fits', 'salvato')

def fan():
    T = fits_table('sdss-dr10d-tiles.fits')

    TT = []
    for tile in T.coadd_id:
        fn = 'sdss-dr10d-phot/phot-%s.fits' % tile
        if not os.path.exists(fn):
            print('File not found:', fn)
            continue
        P = fits_table(fn)
        I = np.flatnonzero(np.logical_or(
            P.w3_nanomaggies * np.sqrt(P.w3_nanomaggies_ivar) > 5,
            P.w4_nanomaggies * np.sqrt(P.w4_nanomaggies_ivar) > 5))
        print(len(P), 'from', fn, 'and', len(I), 'with W3 or W4 S/N>5')
        if len(I) == 0:
            continue
        TT.append(P[I])
    T = merge_tables(TT)
    T.writeto('/global/cscratch1/sd/dstn/fan.fits')


def mccleary():
    all_matches_near('/global/cscratch1/sd/dstn/tmp/redmagic_dr8_public_v6.3_faint.fits',
                     1., 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot',
                     'mccleary-matches.fits', 'mccleary')


if __name__ == '__main__':
    mccleary()
    #fan()
    #hoyle()
    #amaral()
    #ogrady()
    #leslie()
    #fengshuai()
    #gupta()
    #salvato()
    sys.exit()

    # Match to the AllWISE catalog for unmatched targets.
    from wise.allwisecat import allwise_catalog_dec_range
    from astrometry.libkd.spherematch import *
    T1 = fits_table('RADEC_forDustin_x.fits')
    T2 = fits_table('salvato-matches.fits')
    unmatched = np.ones(len(T1), bool)
    unmatched[np.unique(T2.salvato_index)] = False
    T = T1[unmatched]
    print(len(T), 'of', len(T1), 'targets did not have SDSS matches')

    print('Dec range', T.dec.min(), T.dec.max())

    dd = np.array(allwise_catalog_dec_range)
    I = np.flatnonzero((dd[:,0] < T.dec.max()) * (dd[:,1] > T.dec.min()))
    print('AllWISE slices atalogs in range:', I)
    MM = []
    for i in I:
        rdfn = '/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i-radec.fits' % i
        print('Reading', rdfn)
        RD = fits_table(rdfn)
        II,JJ,d = match_radec(T.ra, T.dec, RD.ra, RD.dec, 30./3600.)
        print(len(II), 'matches')
        if len(II) == 0:
            continue
        catfn = '/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i.fits' % i
        print('Reading', catfn)
        M = fits_table(catfn, rows=JJ)
        M.salvato_index = II
        M.ra_salvato = T.ra[II]
        M.dec_salvato = T.dec[II]
        MM.append(M)
    M = merge_tables(MM)
    M.writeto('salvato-allwise.fits')


    import sys
    sys.exit(0)

    T = fits_table('target_list_SDSS_DR7.fits')
    print(len(T), 'targets')
    W = fits_table('hoyle-matches.fits')
    I,J,d = match_radec(T.ra, T.dec, W.ra, W.dec, 1.0, nearest=True)
    print(len(I), 'matches')
    print(len(np.unique(I)), 'unique targets')
    print(len(np.unique(T.ra[I])), 'unique RAs')
    print(len(np.unique(J)), 'unique WISE')
    print(len(np.unique(W.objid[J])), 'unique WISE objids')
    K = np.argsort(I)
    I = I[K]
    J = J[K]
    W.cut(J)
    #assert(np.all(W.hoyle_index == np.arange(len(T))))

    for k in ['index', 'hoyle_index', 'ra_hoyle', 'dec_hoyle']:
        W.delete_column(k)
    W.writeto('hoyle.fits')

    #brandt()
    #manga()
    #nunez()
    #sien()


