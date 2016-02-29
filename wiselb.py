from __future__ import print_function
import tempfile
import sys

###
#sys.path.insert(0, 'django-1.7')
###

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import fitsio
from astrometry.util.util import *
from astrometry.util.fits import *
from astrometry.util.resample import *
from astrometry.util.starutil_numpy import *

def _read_tan_wcs(sourcefn, ext, hdr=None, W=None, H=None):
    from astrometry.util.util import Tan
    wcs = Tan(sourcefn, ext)
    return wcs

def halfsize(sourcefn, halffn, read_wcs=None):
    I,hdr = fitsio.read(sourcefn, header=True)
    H,W = I.shape
    # make even size; smooth down
    if H % 2 == 1:
        I = I[:-1,:]
    if W % 2 == 1:
        I = I[:,:-1]

    # ??
    #im = gaussian_filter(I, 1.)
    im = I

    # bin (excluding NaN)
    q1 = im[::2,::2]
    q2 = im[1::2,::2]
    q3 = im[1::2,1::2]
    q4 = im[::2,1::2]

    f1 = np.isfinite(q1)
    f2 = np.isfinite(q2)
    f3 = np.isfinite(q3)
    f4 = np.isfinite(q4)

    I2 = (np.where(f1, q1, 0) +
          np.where(f2, q2, 0) +
          np.where(f3, q3, 0) +
          np.where(f4, q4, 0)) / np.maximum(1, f1+f2+f3+f4)
    #I2 = (im[::2,::2] + im[1::2,::2] + im[1::2,1::2] + im[::2,1::2])/4.
    I2 = I2.astype(np.float32)
    # shrink WCS too
    if read_wcs is None:
        read_wcs = _read_tan_wcs
    wcs = read_wcs(sourcefn, 0, hdr=hdr, W=W, H=H)
    # include the even size clip; this may be a no-op
    H,W = im.shape
    wcs = wcs.get_subimage(0, 0, W, H)
    subwcs = wcs.scale(0.5)
    hdr = fitsio.FITSHDR()
    subwcs.add_to_header(hdr)
    dirnm = os.path.dirname(halffn)
    f,tmpfn = tempfile.mkstemp(suffix='.fits.tmp', dir=dirnm)
    os.close(f)
    # To avoid overwriting the (empty) temp file (and fitsio
    # printing "Removing existing file")
    os.unlink(tmpfn)
    fitsio.write(tmpfn, I2, header=hdr, clobber=True)
    os.rename(tmpfn, halffn)
    print('Wrote', halffn)


ver = 1
patdata = dict(ver=ver)

# basescale = 5
# tilesize = 256
# tiles = 2**basescale
# side = tiles * tilesize

if False:
    H,W = 4096,8192
    scalelevel = 5
    wcs = anwcs_create_allsky_hammer_aitoff2(0., 0., W, H)
    #wcs.write_to('hammertime.wcs')
    basepat = 'w%i-lb-%i-%i.fits' % (band, scalelevel, H)
    jpegfn = 'wlb.jpg'
    
else:
    # ~2 arcmin per pixel
    # |l| < 60
    # |b| < 30
    width = 120.
    W = int(width * 60.) / 2
    H = W/2
    zoom = 360. / width
    wcs = anwcs_create_hammer_aitoff(0., 0., zoom, W, H, FALSE)

    scale_allwise = False

    # AllWISE coadds
    #basepat = 'allwise-w%i-lbzoom-%i-%i.fits'
    #jpegfn = 'allwise-wlbzoom.png'

    # unWISE
    basepat = 'unwise-w%i-lbzoom-%i-%i.fits'
    jpegfn = 'unwise-wlbzoom.png'

    # unWISE 2
    #basepat = 'unwise2-w%i-lbzoom-%i-%i.fits'
    #jpegfn  = 'unwise2-wlbzoom.png'

    #basepat = 'w%i-lbzoom-%i-%i.fits'
    #jpegfn = 'wlbzoom.jpg'
    scalelevel = 5
    
from decals import settings
from map.views import _unwise_to_rgb

#w1bfn = 'w1lb-%i.fits' % basescale
#w2bfn = 'w2lb-%i.fits' % basescale

#T = fits_table('allsky-atlas.fits')
T = fits_table('wisex-atlas.fits')

imgs = []

#for band in [3,4, 1,2]:
for band in [1,2]:
    outfn = basepat % (band, scalelevel, H)
    if os.path.exists(outfn):
        print('Exists:', outfn)
        img = fitsio.read(outfn)

        #outfn = outfn.replace('.fits', '-u.fits')
        # fn = 'hammertime.wcs'
        # wcs.writeto(fn)
        # hdr = fitsio.read_header(fn)
        # hdr['CTYPE1'] = 'GLON-AIT'
        # hdr['CTYPE2'] = 'GLAT-AIT'
        # fitsio.write(outfn.replace('.fits','-wcs.fits'), img, header=hdr, clobber=True)

        imgs.append(img)
        continue

    img = np.zeros((H,W), np.float32)
    uimg = np.zeros((H,W), np.float32)
    nimg = np.zeros((H,W), np.uint32)

    for i,brick in enumerate(T.coadd_id):
        # unWISE viewer-scaled
        # fn = os.path.join('data/scaled/unwise/%iw%i' %
        #                   (scalelevel, band), brick[:3],
        #                   'unwise-%s-w%i.fits' % (brick, band))

        # unWISE
        fn = os.path.join('unwise-coadds', brick[:3], brick,
                          'unwise-%s-w%i-img-u.fits' % (brick, band))
        qfn = os.path.join('unwise-coadds-quarter',
                           'unwise-%s-w%i.fits' % (brick, band))
        hfn = os.path.join('unwise-coadds-half',
                           'unwise-%s-w%i.fits' % (brick, band))

        # unWISE-R
        # fn = os.path.join('unwise2-coadds', brick[:3], brick,
        #                   'unwise-%s-w%i-img-u.fits' % (brick, band))
        # qfn = os.path.join('unwise2-coadds-quarter',
        #                    'unwise-%s-w%i.fits' % (brick, band))
        # hfn = os.path.join('unwise2-coadds-half',
        #                    'unwise-%s-w%i.fits' % (brick, band))


        # AllWISE
        # qfn = os.path.join('wise-coadds-quarter',
        #                    '%s_ac51-w%i-int-3.fits' % (brick, band))
        # hfn = os.path.join('wise-coadds-half',
        #                    '%s_ac51-w%i-int-3.fits' % (brick, band))
        # fn = os.path.join('wise-coadds', brick[:2], brick[:4], '%s_ac51' % brick,
        #                   '%s_ac51-w%i-int-3.fits' % (brick, band))

        if not os.path.exists(qfn):
            if not os.path.exists(hfn):
                print('Reading', fn)
                halfsize(fn, hfn)
                #print('Wrote', hfn)
            halfsize(hfn, qfn)
            #print('Wrote', qfn)
        fn = qfn

        #fn = os.path.join('wise-coadds', brick[:2], brick[:4], '%s_ac51' % brick,
        #                  '%s_ac51-w%i-int-3.fits' % (brick, band))
                 
        print('Reading', fn)
        I = fitsio.read(fn)
        bwcs = Tan(fn, 0)
        bh,bw = I.shape
        print('Image shape', bh,bw)

        assert(np.all(np.isfinite(I)))
    
        xx,yy = np.meshgrid(np.arange(bw), np.arange(bh))
        rr,dd = bwcs.pixelxy2radec(xx, yy)
        #print('RA,Dec range', rr.min(), rr.max(), dd.min(), dd.max())
        ll,bb = radectolb(rr.ravel(), dd.ravel())
        #print('L,B range', ll.min(), ll.max(), bb.min(), bb.max())
        ll = ll.reshape(rr.shape)
        bb = bb.reshape(rr.shape)
    
        ok,ox,oy = wcs.radec2pixelxy(ll, bb)
        #print('Unique ok:', np.unique(ok))
        ox = np.round(ox - 1).astype(int)
        oy = np.round(oy - 1).astype(int)
        K = (ox >= 0) * (ox < W) * (oy >= 0) * (oy < H) * ok
        n1 = np.sum(K)
        K *= np.isfinite(I)
        n2 = np.sum(K)
        print('Non-finite input pixels:', n1-n2)
        if np.sum(K) == 0:
            # no overlap
            print('No overlap')
            continue

        # img [oy[K], ox[K]] += I[K]
        # uimg[oy[K], ox[K]] += (I[K] * (nimg[oy[K], ox[K]] == 0))
        # nimg[oy[K], ox[K]] += 1

        np.add.at( img, (oy[K], ox[K]), I[K])
        np.add.at(nimg, (oy[K], ox[K]), 1)

        # for y,x,im in zip(oy[K],ox[K],I[K]):
        #     img [y,x] += im
        #     uimg[y,x] += im * (nimg[y,x] == 0)
        #     nimg[y,x] += 1

        #if i % 10 == 0:
        #    fitsio.write('interim.fits', img / np.maximum(nimg, 1), clobber=True)

            
    img /= np.maximum(nimg, 1)

    fn = 'hammertime.wcs'
    wcs.writeto(fn)
    hdr = fitsio.read_header(fn)
    hdr['CTYPE1'] = 'GLON-AIT'
    hdr['CTYPE2'] = 'GLAT-AIT'

    fitsio.write(outfn, img, header=hdr, clobber=True)
    #fitsio.write(outfn.replace('.fits', '-u.fits'), uimg, header=hdr, clobber=True)
    fitsio.write(outfn.replace('.fits', '-n.fits'), nimg, header=hdr, clobber=True)
    imgs.append(img)


sys.exit(0)

w1,w2 = imgs[:2]

if scale_allwise:
    # AllWISE zeropoints: W1 = 20.5, W2 = 19.5
    w1 /= 10.**((20.5 - 22.5)/2.5)
    w2 /= 10.**((19.5 - 22.5)/2.5)

    w1 /= 4
    w2 /= 4

    # Just for kicks?
    w1 /= 2
    w2 /= 2

fns = ['w1-lbzoom-5-1800-u-wcs.fits',
       'w2-lbzoom-5-1800-u-wcs.fits',]
u1,u2 = [fitsio.read(fn) for fn in fns]

plt.clf()
plt.loglog(np.maximum(1000, w1.ravel()), np.maximum(1000, u1.ravel()), 'k.', alpha=0.01)
plt.savefig('w1.png')
plt.clf()
plt.loglog(np.maximum(1000, w2.ravel()), np.maximum(1000, u2.ravel()), 'k.', alpha=0.01)
plt.savefig('w2.png')


S,Q = 3000,25
rgb = _unwise_to_rgb([w1, w2], S=[S,S], Q=Q)
plt.imsave(jpegfn, rgb, origin='lower')


#w1 -= np.median(w1[600:1200, 1200:2400])
#w2 -= np.median(w2[600:1200, 1200:2400])

w1 -= np.percentile(w1[600:1200, 1200:2400], 25)
w2 -= np.percentile(w2[600:1200, 1200:2400], 25)


S,Q = 3000,25
rgb = _unwise_to_rgb([w1, w2], S=[S,S], Q=Q)
plt.imsave(jpegfn.replace('.png','-2.png'), rgb, origin='lower')
