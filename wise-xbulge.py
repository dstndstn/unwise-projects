from __future__ import print_function

import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='computer modern roman')
matplotlib.rc('font', **{'sans-serif': 'computer modern sans serif'})
import pylab as plt
import fitsio
from astrometry.util.util import *
from astrometry.util.plotutils import *
from astrometry.util.starutil_numpy import *
from astrometry.util.fits import *
from astrometry.libkd.spherematch import *

if False:
    T = fits_table('allsky-atlas.fits')
    print(len(T), 'tiles')
    
    # Aim for ~2 arcmin per pixel
    # |l| < 50
    # |b| < 30
    
    
    width = 60
    W = int(width * 60.) / 2
    H = W/2
    zoom = 360. / width
    wcs = anwcs_create_hammer_aitoff(0., 0., zoom, W, H, FALSE)
    
    xx,yy = np.meshgrid(np.linspace(1, W, 100), np.linspace(1, H, 100))
    xx = xx.ravel()
    yy = yy.ravel()
    ok,ll,bb = wcs.pixelxy2radec(xx, yy)
    ra,dec = lbtoradec(ll, bb)
    
    plt.clf()
    plt.plot(ra, dec, 'b.')
    plt.savefig('1.png')
        
    I,J,d = match_radec(T.ra, T.dec, ra, dec, 2., nearest=True)
    print(len(I), 'matches')
    
    plt.clf()
    plt.plot(T.ra[I], T.dec[I], 'b.')
    plt.savefig('2.png')
    
    if False:
        for t in T[I]:
            brick = t.coadd_id
            for band in [1,2]:
                dirnm = '%s/%s/%s_ac51' % (brick[:2], brick[:4], brick)
                #print('mkdir -p %s' % dirnm)
                print('wget --continue -r -nH --cut-dirs 5 http://irsa.ipac.caltech.edu/ibe/data/wise/merge/merge_p3am_cdd/%s/%s_ac51-w%i-int-3.fits'
                      % (dirnm, brick, band))


fns = ['w1-lbzoom-5-1800-u-wcs.fits',
       'w2-lbzoom-5-1800-u-wcs.fits',]
imgs = [fitsio.read(fn) for fn in fns]
wcs = anwcs(fns[0])
print('WCS header', wcs)


#from decals import settings
from map.views import _unwise_to_rgb

S,Q = 3000,25
rgb = _unwise_to_rgb(imgs, S=[S]*len(imgs), Q=Q)

H,W = imgs[0].shape
print('Image size', W, 'x', H)
ok,l1,b1 = wcs.pixelxy2radec(1, (H+1)/2.)
ok,l2,b2 = wcs.pixelxy2radec(W, (H+1)/2.)
ok,l3,b3 = wcs.pixelxy2radec((W+1)/2., 1)
ok,l4,b4 = wcs.pixelxy2radec((W+1)/2., H)

print('L,B', (l1,b1), (l2,b2), (l3,b3), (l4,b4))

llo,lhi = l2,l1+360
blo,bhi = b3,b4

ps = PlotSequence('xbulge', suffix='pdf')

plt.figure(figsize=(10,5))
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)

plt.clf()
plt.imshow(rgb, origin='lower', interpolation='nearest',
           extent=[lhi,llo,blo,bhi], aspect=1.)#float(H)/W)
xt = 300 + np.arange(5)*30
plt.xticks(xt, ['%i' % (x % 360) for x in xt])
plt.yticks([-30,-15,0,15,30])
plt.xlim(lhi,llo)
plt.xlabel('Galactic longitude $\ell$ (deg)')
plt.ylabel('Galactic latitude $b$ (deg)')
ps.savefig()

w1,w2 = imgs
medy1 = np.median(w1, axis=1)
medy2 = np.median(w2, axis=1)

rgb = _unwise_to_rgb([w1 - medy1[:,np.newaxis],
                      w2 - medy2[:,np.newaxis]], S=[S]*len(imgs), Q=Q)

plt.clf()
plt.imshow(rgb, origin='lower')
plt.savefig('2.png')

wm1 = w1 - medy1[:,np.newaxis]
wm2 = w2 - medy2[:,np.newaxis]

medx1 = np.median(wm1, axis=0)
medx2 = np.median(wm2, axis=0)

#medx1 = np.median(w1, axis=0)
#medx2 = np.median(w2, axis=0)

plt.clf()
plt.plot(np.median(w1 - medy1[:,np.newaxis], axis=0))
plt.savefig('3.png')

rgb = _unwise_to_rgb([wm1 - medx1[np.newaxis,:],
                      wm2 - medx2[np.newaxis,:]], S=[S]*len(imgs), Q=Q)

plt.clf()
plt.imshow(rgb, origin='lower')
plt.savefig('4.png')
