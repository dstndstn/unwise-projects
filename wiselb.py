
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
    basepat = 'w%i-lbzoom-%i-%i.fits'
    jpegfn = 'wlbzoom.jpg'
    scalelevel = 5
    
from decals import settings
from map.views import _unwise_to_rgb

#w1bfn = 'w1lb-%i.fits' % basescale
#w2bfn = 'w2lb-%i.fits' % basescale

T = fits_table('allsky-atlas.fits')

imgs = []

#for band in [1,2,3,4]:
for band in [1,2]:
    outfn = basepat % (band, scalelevel, H)
    if os.path.exists(outfn):
        outfn = outfn.replace('.fits', '-u.fits')
        img = fitsio.read(outfn)

        fn = 'hammertime.wcs'
        wcs.writeto(fn)
        hdr = fitsio.read_header(fn)
        hdr['CTYPE1'] = 'GLON-AIT'
        hdr['CTYPE2'] = 'GLAT-AIT'
        fitsio.write(outfn.replace('.fits','-wcs.fits'), img, header=hdr, clobber=True)

        imgs.append(img)
        continue

    img = np.zeros((H,W), np.float32)
    uimg = np.zeros((H,W), np.float32)
    nimg = np.zeros((H,W), np.uint8)

    for i,brick in enumerate(T.coadd_id):
        fn = os.path.join('data/scaled/unwise/%iw%i' %
                          (scalelevel, band), brick[:3],
                          'unwise-%s-w%i.fits' % (brick, band))
        print 'Reading', fn
        I = fitsio.read(fn)
        bwcs = Tan(fn, 0)
        bh,bw = I.shape
    
        xx,yy = np.meshgrid(np.arange(bw), np.arange(bh))
        rr,dd = bwcs.pixelxy2radec(xx, yy)
        ll,bb = radectolb(rr.ravel(), dd.ravel())
        ll = ll.reshape(rr.shape)
        bb = bb.reshape(rr.shape)
    
        ok,ox,oy = wcs.radec2pixelxy(ll, bb)
        ox = np.round(ox - 1).astype(int)
        oy = np.round(oy - 1).astype(int)
        K = (ox >= 0) * (ox < W) * (oy >= 0) * (oy < H) * (ok == 0)
        if np.sum(K) == 0:
            # no overlap
            continue
    
        img [oy[K], ox[K]] += I[K]
        uimg[oy[K], ox[K]] += (I[K] * (nimg[oy[K], ox[K]] == 0))
        nimg[oy[K], ox[K]] += 1
            
    img /= np.maximum(nimg, 1)
    fitsio.write(outfn, img, clobber=True)
    fitsio.write(outfn.replace('.fits', '-u.fits'), uimg, clobber=True)
    fitsio.write(outfn.replace('.fits', '-n.fits'), nimg, clobber=True)
    imgs.append(img)


w1,w2 = imgs
S,Q = 3000,25
rgb = _unwise_to_rgb([w1, w2], S=[S]*len(imgs), Q=Q)
plt.imsave(jpegfn, rgb, origin='lower')

