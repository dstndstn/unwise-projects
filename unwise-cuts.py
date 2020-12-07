from __future__ import print_function
import os
import numpy as np
from glob import glob

from astrometry.util.fits import *
from astrometry.sdss.common import cas_flags, photo_flags1_map, photo_flags2_map
from astrometry.libkd.spherematch import match_radec, tree_build_radec, tree_free, tree_search_radec, trees_match

def li():
    #(W1, W2, W3 and W4) for SDSS (galaxy) sources with 12 < Ra < 40 deg,  -2 < dec < 2 deg and redshift 0 < z < 0.2.
    if False:
        T = fits_table('sdss-dr13-tiles.fits')
        photdir = 'sdss-dr13-phot'
        margin = 0.8
        T.cut((T.ra > 12-margin) * (T.ra < (40+margin)) * (T.dec > (-2-margin)) * (T.dec < (2+margin)))
        print(len(T), 'tiles')
        WW = []
        for tile in T.coadd_id:
            fn = os.path.join(photdir, 'phot-%s.fits' % tile)
            if not os.path.exists(fn):
                print('Does not exist:', fn)
                continue
            W = fits_table(fn)
            print(len(W), 'from', fn)
            if len(W) == 0:
                continue
            W.cut((W.ra > 12) * (W.ra < 40) * (W.dec > -2) * (W.dec < 2))
            print('Cut to', len(W), 'in box')
            WW.append(W)
        W = merge_tables(WW)
        W.writeto('li.fits')

    S = fits_table('/global/cfs/cdirs/cosmo/data/sdss/dr13/sdss/spectro/redux/specObj-dr13.fits')
    print(len(S), 'spectra')
    S.about()
    T = fits_table('li.fits.gz')
    print(len(T), 'WISE')
    T.about()
    I,J,d = match_radec(S.plug_ra, S.plug_dec, T.ra, T.dec, 1., nearest=True)
    print(len(I), 'match')
    S.cut(I)
    T.cut(J)
    S.add_columns_from(T)
    S.writeto('li-match.fits')
    

def kovlakas():
    M = fits_table('SDSS12_objids.fits')
    M.rename('radeg_sdss', 'ra')
    M.rename('dedeg_sdss', 'dec')
    # all_matches_near(None, 5., 'sdss-dr10d-tiles.fits',
    #                  'sdss-dr10d-phot', 'kovlakas-unwise-dr10.fits', 'kovlakas', M=M, closest=True)
    all_matches_near(None, 5., 'sdss-dr13-tiles.fits',
                     'sdss-dr13-phot', 'kovlakas-unwise-dr13.fits', 'kovlakas', M=M, closest=True)

def pizarro():
    T = fits_table('cv-sources.fits')
    T.rename('ra_deg_best', 'ra')
    T.rename('dec_deg_best', 'dec')

    # all_matches_near(None, 5.0, 'sdss-dr13-atlas.fits', 'sdss-dr13-phot',
    #                  'cv-sdss-matches.fits', 'pizarro', M=T)
                     
    
    cat_all_matches_near(T, 5., 'cv-catdr2-matches.fits', 'pizarro',
                         catdir='/global/cfs/cdirs/cosmo/work/wise/unwise_catalog/dr2-cat/band-merged')

def finn4():
    all_matches_near('vf_north_v1_main.fits', 10.,
                     'sdss-dr13-tiles.fits',
                     'sdss-dr13-phot',
                     'vf_north_v1_main-matched.fits', 'vf')

def finn3():
    if False:
        all_matches_near('a100.SDSSObjID.191001.fits', 1.,
                         'sdss-dr10d-tiles.fits', 'sdss-dr10d-phot',
                         'a100.SDSSObjID.191001.match-dr10.fits', 'finn')

    T = fits_table('a100.SDSSObjID.191001.fits')
    M = fits_table('a100.SDSSObjID.191001.match.fits')
    M.objid = np.array([int(o.strip()) for o in M.objid])
    print('Targets:', len(T))
    print('Matches:', len(M))

    print('Unique objids in Targets:', len(set(T.objid)))
    
    I,J,d = match_radec(T.ra, T.dec, M.ra, M.dec, 1./3600.)
    
    K = np.flatnonzero(M.objid[J] == T.objid[I])
    Tmatched = T[I[K]]
    Mmatched = M[J[K]]
    print(len(Tmatched), 'matched with objid')

    print('Targets:')
    Tmatched.about()
    print('Matches:')
    Mmatched.about()

    Mmatched.writeto('a100.SDSSObjID.191001-matched-A.fits')

    matched = np.zeros(len(T))
    matched[I[K]] = True
    TU = T[matched == False]
    print(len(TU), 'un-matched')

    I,J,d = match_radec(TU.ra, TU.dec, M.ra, M.dec, 2./3600, nearest=True)
    print('Matched', len(I))

    MU = M[J]
    MU.writeto('a100.SDSSObjID.191001-matched-B.fits')
    print(len(set(MU.objid).intersection(Mmatched.objid)), 'objids in common A-B')

    matched = np.zeros(len(TU))
    matched[I] = True
    TUU = TU[matched == False]
    print(len(TUU), 'un-matched')
    TUU.writeto('a100.SDSSObjID.191001-unmatched.fits')

    

def finn2():
    all_matches_near('vf_main.fits', 5.,
                     'sdss-dr13-atlas.fits', 'sdss-dr13-phot',
                     'vf_main-match.fits', 'vf')
    
def finn():
    if False:
        all_matches_near('a100.SDSSObjID.191001.fits', 1.,
                         'sdss-dr13-atlas.fits', 'sdss-dr13-phot',
                         'a100.SDSSObjID.191001.match.fits', 'finn')

    T = fits_table('a100.SDSSObjID.191001.fits')
    M = fits_table('a100.SDSSObjID.191001.match.fits')
    print('Targets:', len(T))
    print('Matches:', len(M))
    
    M.objid = np.array([int(o.strip()) for o in M.objid])
    M.agc = T.agc[M.finn_index]
    M.finn_objid = T.objid[M.finn_index]

    M.writeto('a100.SDSSObjID.191001.match2.fits')
    o1 = set(T.objid)
    o2 = set(M.objid)
    print('Unique targets:', len(o1))
    print('Unique matches:', len(o2))
    print('Union:', len(o1.union(o2)))
    print('Intersection:', len(o1.intersection(o2)))

    name = 'finn'
    from collections import Counter
    dups = Counter(zip(M.objid, M.get('%s_index' % name)))
    print('Duplicate matches:', dups.most_common(5))
    keep = np.ones(len(M), bool)
    for (objid,ind),n in dups.most_common():
        if n <= 1:
            break
        I = np.flatnonzero((M.objid == objid) * (M.get('%s_index' % name) == ind))
        #print(len(I), 'for', objid, ind)
        dists = np.hypot(M.x[I] - 1024.5, M.y[I] - 1024.5)
        imin = np.argmin(dists)
        ii = I[imin]
        keep[I] = False
        keep[ii] = True
    M.cut(keep)
    print('Cut to matches:', len(M))
    o2 = set(M.objid)
    print('Unique targets:', len(o1))
    print('Unique matches:', len(o2))
    print('Union:', len(o1.union(o2)))
    print('Intersection:', len(o1.intersection(o2)))
    K = np.flatnonzero(M.objid == M.finn_objid)
    print(len(K), 'have matching objids')

    #c = Counter(M.objid)
    c = Counter(M.finn_objid)
    hist_nmatches = Counter(c.values())
    print('Number of matches per target:', hist_nmatches.most_common())
    
    M.writeto('a100.SDSSObjID.191001.match3.fits')

    
    

def secrest():
    # text2fits -H "id ra dec" -f jdd unWISE_request.ipac unWISE_request.fits
    T = fits_table('unWISE_request.fits')

    cat_all_matches_near(T, 5., 'unWISE_request_match.fits', 'secrest')
    



def cat_all_matches_near(M, radius, outfn, projname,
                         catdir = '/global/project/projectdirs/cosmo/work/wise/unwise_catalog/dr1/band-merged',
                         closest=False, **kwargs):
    '''
    radius in arcsec
    '''
    atlasfn = 'allsky-atlas.fits'
    T = fits_table(atlasfn)
    print(len(T), 'tiles')

    #M = fits_table(catfn)
    print(len(M), 'targets')

    M.set('%s_index' % projname, np.arange(len(M)))

    I,J,d = match_radec(M.ra, M.dec, T.ra, T.dec, 1.2)
    keep = np.zeros(len(T), bool)
    keep[J] = True
    T.cut(keep)
    print(len(T), 'tiles near targets')

    WW = []
    for t in T:
        fn = os.path.join(catdir, '%s.cat.fits' % t.coadd_id)
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

    print(len(np.unique(W.unwise_objid)), 'unique objids')

    if closest:
        from collections import Counter
        dups = Counter(zip(W.unwise_objid, W.get('%s_index' % projname)))
        print('Duplicate matches:', dups.most_common(5))
        keep = np.ones(len(W), bool)
        for (objid,ind),n in dups.most_common():
            if n <= 1:
                break
            I = np.flatnonzero((W.unwise_objid == unwise_objid) * (W.get('%s_index' % projname) == ind))
            print(len(I), 'for', objid, ind)
            dists = np.hypot(W.x[I] - 1024.5, W.y[I] - 1024.5)
            imin = np.argmin(dists)
            ii = I[imin]
            keep[I] = False
            keep[ii] = True
        W.cut(keep)

    W.writeto(outfn)

    

def hamadouche():
    if False:
        T = fits_table('hamadouche.fits')
        T.rename('ras', 'ra')
        T.rename('decs', 'dec')
        T.rename('objid', 'hamadouche_objid')
        T.writeto('hamadouche2.fits')
        all_matches_near('hamadouche2.fits', 1.0, 'sdss-dr13-atlas.fits',
                         'sdss-dr13-phot', 'table1forunWISE-matched.fits', 'hamadouche')#, closest=True)
    M = fits_table('table1forunWISE-matched.fits')
    M.objid = np.array([int(o) for o in M.objid])
    print(len(M), 'matches,', len(np.unique(M.objid)), 'objids')
    M.cut(M.hamadouche_objid == M.objid)
    print(len(M), 'with matching objids')
    o,I = np.unique(M.hamadouche_objid, return_index=True)
    M.cut(I)
    print(len(M), 'unique matching objids')
    M.writeto('table1forunWISE-matched2.fits')




def frasermckelvie():
    all_matches_near('MWAG_coords.fits', 1.0, 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'fm.fits', 'frasermckelvie', closest=True)

def frasermckelvie2():
    all_matches_near('mwa.fits', 1.0, 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'fm2.fits', 'frasermckelvie', closest=True)

def mccleary2():
    # A random subsample of 200k galaxies over the sky with SDSS & WISE mags/fluxes
    np.random.seed(42)

    # # 8k tiles
    # A = fits_table('sdss-dr13-atlas.fits')
    # # Choose tiles randomly
    # I = np.random.permutation(len(A))
    # A.cut(I)
    # for tile in A.coadd_id:

    W = fits_table('window_flist-cut.fits', columns=['run','camcol','field'])
    print('Read', len(W), 'fields')
    # Choose SDSS fields randomly
    I = np.random.permutation(len(W))
    W.cut(I)

    os.environ['BOSS_PHOTOOBJ'] = '/global/project/projectdirs/cosmo/data/sdss/dr13/eboss/photoObj/'

    ntotal = 0
    TT = []
    for i,(run,camcol,field) in enumerate(zip(W.run, W.camcol, W.field)):
        print()
        print(i+1, 'of', len(W), ': Run,camcol,field', run,camcol,field)
        fn = os.path.join(os.environ['BOSS_PHOTOOBJ'], '301', '%i'%run,
                          '%i'%camcol, 'photoObj-%06i-%i-%04i.fits' % (run,camcol,field))
    
        T = fits_table(fn, columns=['objc_type', 'resolve_status',])
        if T is None:
            continue
        T.index = np.arange(len(T))
        T.cut((T.resolve_status & 256) > 0)
        T.cut(T.objc_type == 3)
        print('Reading', len(T), 'from', fn)
        if len(T) == 0:
            continue
    
        # Re-read SDSS photoObjs and grab more columns.
        inds = T.index
        T = fits_table(fn, rows=inds, columns=['run','camcol','field',
                                               'ra','dec','raerr','decerr',
                                               'cmodelmag', 'cmodelmagerr',
                                               'modelmag', 'modelmagerr',
                                               'psfmag','psfmagerr',
                                               'flags', 'flags2', 'objid',
                                               'modelflux', 'modelflux_ivar',
                                               'cmodelflux', 'cmodelflux_ivar'])
    
        unwdir = '/global/project/projectdirs/cosmo/data/unwise/allwise/unwise-phot/sdss-dr13-pobj'
        fn = os.path.join(unwdir, '%i'%run, '%i'%camcol, 'photoWiseForced-%06i-%i-%04i.fits' % (run, camcol, field))
        print('Reading from', fn)
        U = fits_table(fn, rows=inds, columns=['treated_as_pointsource', 'pointsource',
                                               'w1_mag', 'w1_mag_err',
                                               'w2_mag', 'w2_mag_err',
                                               'w3_mag', 'w3_mag_err',
                                               'w4_mag', 'w4_mag_err',
                                               'w1_nanomaggies', 'w1_nanomaggies_ivar',
                                               'w2_nanomaggies', 'w2_nanomaggies_ivar',
                                               'w3_nanomaggies', 'w3_nanomaggies_ivar',
                                               'w4_nanomaggies', 'w4_nanomaggies_ivar',
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
    
        TT.append(T)
        ntotal += len(T)
        print('Total:', ntotal)
        if ntotal > 200000:
            break
    
    T = merge_tables(TT)
    T[:200000].writeto('mccleary2.fits')



def seuss():
    #all_matches_near('psbs.fits', 1.0, 'sdss-dr13-atlas.fits',
    #                 'sdss-dr13-phot', 'psbs-unwise.fits', 'seuss')
    name = 'seuss'
    T = fits_table('psbs-unwise.fits')
    from collections import Counter
    dups = Counter(zip(T.objid, T.get('%s_index' % name)))
    print('Duplicate matches:', dups.most_common(5))
    keep = np.ones(len(T), bool)
    for (objid,ind),n in dups.most_common():
        if n <= 1:
            break
        I = np.flatnonzero((T.objid == objid) * (T.get('%s_index' % name) == ind))
        print(len(I), 'for', objid, ind)
        dists = np.hypot(T.x[I] - 1024.5, T.y[I] - 1024.5)
        imin = np.argmin(dists)
        ii = I[imin]
        keep[I] = False
        keep[ii] = True

    T.cut(keep)
    T.writeto('psbs-unwise-2.fits')

def lodieu():
    #for fn in ['mytable1.fits', 'mytable2.fits', 'mytable3.fits']:
    #T = fits_table(fn)
    #all_matches_near('mytable1.fits', 1.0, 'sdss-dr13-atlas.fits',
    #                 'sdss-dr13-phot', 'mytable1-unwise.fits', 'lodieu')
    all_matches_near('mytable2.fits', 1.0, 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'mytable2-unwise.fits', 'lodieu')
    all_matches_near('mytable3.fits', 1.0, 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'mytable3-unwise.fits', 'lodieu')


def cerulo():
    T = fits_table('SDSS_DR13_galaxy_table_for_unWISE_query_pcerulo_v2.fit.gz')
    #T = fits_table('SDSS_DR13_galaxy_table_for_unWISE_query_pcerulo.fit.gz')
    print(len(T), 'objids')

    if True:
        jj = []
        SS = []
        fns = glob('sdss-dr13-phot/phot-*.fits')
        fns.sort()
        for fn in fns:
            print(fn)
            S = fits_table(fn)
            print(len(S), 'sources')
            #S.objid = np.array([int(o) for o in S.objid])
            #sobjid = set(S.objid)
            #u = sobjid.intersection(objids)
            #print(len(u), 'in intersection')
            inds = dict([(int(o), i) for i,o in enumerate(S.objid)])
            I = []
            for j,o in enumerate(T.objid):
                try:
                    i = inds[o]
                    I.append(i)
                    jj.append(j)
                except KeyError:
                    continue
            print('Matched', len(I))
            if len(I) == 0:
                continue
            S.cut(np.array(I))
            SS.append(S)
        S = merge_tables(SS)
        S.writeto('cerulo-matched.fits')

    S = fits_table('cerulo-matched.fits')
    #sind = dict([(int(o), i) for i,o in enumerate(S.objid)])
    S.dist = np.hypot(S.x - 1024, S.y - 1024)

    print('x', S.x.min(), S.x.max(), 'y', S.y.min(), S.y.max())
    T.sindex = np.empty(len(T), int)
    T.sindex[:] = -1
    T.sdist = np.empty(len(T), np.float32)
    T.sdist[:] = 1e6
    tind = dict([(o, i) for i,o in enumerate(T.objid)])
    for i,o in enumerate(S.objid):
        ti = tind[int(o)]
        if S.dist[i] < T.sdist[ti]:
            T.sdist[ti] = S.dist[i]
            T.sindex[ti] = i

    print(np.sum(T.sindex > -1), 'of', len(T), 'match')
    #assert(np.all(T.sindex > -1))
    ## HACK
    T.cut(T.sindex > -1)

    for c in S.columns():
        if c in ['objid', 'dist']:
            continue
        T.set(c, S.get(c)[T.sindex])
    T.delete_column('sindex')
    T.delete_column('sdist')
    T.writeto('cerulo.fits')


def wang():
    #text2fits objID_GSWLC_sdss.dat wang.fits -f k
    T = fits_table('wang.fits')
    print(len(T), 'objids')
    #objids = set(T.objid)

    if False:
        jj = []
        SS = []
        fns = glob('sdss-dr13-phot/phot-*.fits')
        fns.sort()
        for fn in fns:
            print(fn)
            S = fits_table(fn)
            print(len(S), 'sources')
            #S.objid = np.array([int(o) for o in S.objid])
            #sobjid = set(S.objid)
            #u = sobjid.intersection(objids)
            #print(len(u), 'in intersection')
            inds = dict([(int(o), i) for i,o in enumerate(S.objid)])
            I = []
            for j,o in enumerate(T.objid):
                try:
                    i = inds[o]
                    I.append(i)
                    jj.append(j)
                except KeyError:
                    continue
            print('Matched', len(I))
            if len(I) == 0:
                continue
            S.cut(np.array(I))
            SS.append(S)
        S = merge_tables(SS)
        S.writeto('wang-matched.fits')

    S = fits_table('wang-matched.fits')
    #sind = dict([(int(o), i) for i,o in enumerate(S.objid)])
    S.dist = np.hypot(S.x - 1024, S.y - 1024)

    print('x', S.x.min(), S.x.max(), 'y', S.y.min(), S.y.max())
    T.sindex = np.empty(len(T), int)
    T.sindex[:] = -1
    T.sdist = np.empty(len(T), np.float32)
    T.sdist[:] = 1e6
    tind = dict([(o, i) for i,o in enumerate(T.objid)])
    for i,o in enumerate(S.objid):
        ti = tind[int(o)]
        if S.dist[i] < T.sdist[ti]:
            T.sdist[ti] = S.dist[i]
            T.sindex[ti] = i

    print(np.sum(T.sindex > -1), 'of', len(T), 'match')
    #assert(np.all(T.sindex > -1))
    ## HACK
    T.cut(T.sindex > -1)

    for c in S.columns():
        if c in ['objid', 'dist']:
            continue
        T.set(c, S.get(c)[T.sindex])
    T.delete_column('sindex')
    T.delete_column('sdist')
    T.writeto('wang-m2.fits')


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



def all_matches_near(catfn, radius, atlasfn, photdir, outfn, projname,
                     closest=False, M=None, **kwargs):
    '''
    radius in arcsec
    '''
    T = fits_table(atlasfn)
    print(len(T), 'tiles')

    if M is None:
        M = fits_table(catfn)
    print(len(M), 'targets')

    if closest:
        kwargs.update(nearest=True)
    
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

    if closest:
        from collections import Counter
        dups = Counter(zip(W.objid, W.get('%s_index' % projname)))
        print('Duplicate matches:', dups.most_common(5))
        keep = np.ones(len(W), bool)
        for (objid,ind),n in dups.most_common():
            if n <= 1:
                break
            I = np.flatnonzero((W.objid == objid) * (W.get('%s_index' % projname) == ind))
            print(len(I), 'for', objid, ind)
            dists = np.hypot(W.x[I] - 1024.5, W.y[I] - 1024.5)
            imin = np.argmin(dists)
            ii = I[imin]
            keep[I] = False
            keep[ii] = True
        W.cut(keep)

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


def lai():
    all_matches_near('lai.fits', 4., 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'lai-matches.fits', 'lai')

def razink():
    all_matches_near('razink.fits', 3., 'sdss-dr13-atlas.fits',
                     'sdss-dr13-phot', 'razink-matches.fits', 'razink', nearest=True)


if __name__ == '__main__':
    import sys
    li()
    #finn4()
    #finn3()
    #finn2()
    #kovlakas()
    #pizarro()
    #finn()
    #secrest()
    #hamadouche()
    #frasermckelvie2()
    #mccleary2()
    #seuss()
    #cerulo()
    #razink()
    #lai()
    #lodieu()
    #wang()
    #mccleary()
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


