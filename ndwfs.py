from wise.forcedphot import unwise_tiles_touching_wcs
from tractor import *
from astrometry.util.fits import *
from astrometry.util.util import Tan
from astrometry.util.multiproc import multiproc

def main():
    mp = multiproc(4)

    T = fits_table('NDWFS.fits')
    #T = T[:100]
    T.rename('alpha_j2000', 'ra')
    T.rename('delta_j2000', 'dec')
    
    ra0, ra1 = T.ra.min(), T.ra.max()
    dec0,dec1 = T.dec.min(), T.dec.max()
    roiradec = [ra0, ra1, dec0, dec1]
    
    print('RA,Dec bounds', roiradec)
    
    ra_center = (ra0 + ra1)/2.
    dec_center = (dec0 + dec1) / 2.
    dra = ra1 - ra0
    ddec = dec1 - dec0
    cd = 1./3600.
    w = int(dra * np.cos(np.deg2rad(dec_center)) / cd)
    h = int(ddec / cd)
    print('W,H', w, h)
    
    fakewcs = Tan(ra_center, dec_center, w/2.+0.5, h/2.+0.5,
                  -cd, 0., 0., cd, float(w), float(h))
    
    tiles = unwise_tiles_touching_wcs(fakewcs)
    print('Cut to', len(tiles), 'unWISE tiles')
    
    wcat = []
    for r,d in zip(T.ra, T.dec):
        wcat.append(PointSource(RaDecPos(r, d),
                                NanoMaggies(w=1.)))
    
    # PSF broadening in post-reactivation data, by band.
    # Newer version from Aaron's email to decam-chatter, 2018-06-14.
    broadening = { 1: 1.0405, 2: 1.0346, 3: None, 4: None }
    
    unwise_dir = os.environ['UNWISE_COADDS_DIR']
    wise_ceres = True
    args = []
    for band in [1,2,3,4]:
        args.append((wcat, tiles, band, roiradec, unwise_dir, wise_ceres,
                     broadening[band]))
    
    
    phots = mp.map(_unwise_phot, args)
    
    WISE = phots[0]
    
    for i,p in enumerate(phots[1:len(args)]):
        if p is None:
            (wcat,tiles,band) = args[i+1][:3]
            print('"None" result from WISE forced phot:', tiles, band)
            continue
        WISE.add_columns_from(p)
    WISE.rename('tile', 'wise_coadd_id')
    
    # Look up mask bits
    ra  = T.ra
    dec = T.dec
    WISE.wise_mask = np.zeros((len(T), 2), np.uint8)
    for tile in tiles.coadd_id:
        # unwise_dir can be a colon-separated list of paths
        found = False
        for d in unwise_dir.split(':'):
            fn = os.path.join(d, tile[:3], tile,
                              'unwise-%s-msk.fits.gz' % tile)
            print('Looking for unWISE mask file', fn)
            if os.path.exists(fn):
                found = True
                break
        if not found:
            print('unWISE mask file for tile', tile, 'does not exist')
            continue
        # read header to pull out WCS (because it's compressed)
        M,hdr = fitsio.read(fn, header=True)
        wcs = Tan(*[float(hdr[k]) for k in
                    ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                     'CD1_1', 'CD1_2', 'CD2_1','CD2_2','NAXIS2','NAXIS1']])
        #print('Read WCS header', wcs)
        ok,xx,yy = wcs.radec2pixelxy(ra, dec)
        hh,ww = wcs.get_height(), wcs.get_width()
        #print('unWISE image size', hh,ww)
        xx = np.round(xx - 1).astype(int)
        yy = np.round(yy - 1).astype(int)
        I = np.flatnonzero(ok * (xx >= 0)*(xx < ww) * (yy >= 0)*(yy < hh))
        print(len(I), 'sources are within tile', tile)
        if len(I) == 0:
            continue
        # Reference the mask image M at yy,xx indices
        Mi = M[yy[I], xx[I]]
        # unpack mask bits
        # The WISE mask files have:
        #  bit 0: W1 bright star, south-going scan
        #  bit 1: W1 bright star, north-going scan
        #  bit 2: W2 bright star, south-going scan
        #  bit 3: W2 bright star, north-going scan
        WISE.wise_mask[I, 0] = ( Mi       & 3)
        WISE.wise_mask[I, 1] = ((Mi >> 2) & 3)
    
    print(WISE)
    WISE.about()

    WISE.writeto('wise.fits')

    T.add_columns_from(WISE)
    T.writeto('ndwfs-wise.fits')



def _unwise_phot(X):
    from wise.forcedphot import unwise_forcedphot
    (wcat, tiles, band, roiradec, unwise_dir, wise_ceres, broadening) = X
    try:
        W = unwise_forcedphot(wcat, tiles, roiradecbox=roiradec, bands=[band],
            unwise_dir=unwise_dir, use_ceres=wise_ceres,
            psf_broadening=broadening)
    except:
        import traceback
        print('unwise_forcedphot failed:')
        traceback.print_exc()

        if wise_ceres:
            print('Trying without Ceres...')
            W = unwise_forcedphot(wcat, tiles, roiradecbox=roiradec,
                                  bands=[band], unwise_dir=unwise_dir,
                                  use_ceres=False, psf_broadening=broadening)
        else:
            W = None
    return W


if __name__ == '__main__':
    main()
