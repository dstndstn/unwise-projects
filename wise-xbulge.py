from __future__ import print_function
import sys
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

#from colormaps import viridis

from scipy.ndimage import median_filter, gaussian_filter

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

plt.figure(1, figsize=(10,5))
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)

plt.figure(2, figsize=(5,5))
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)

plt.figure(1)
plt.clf()
plt.imshow(rgb, origin='lower', interpolation='nearest',
           extent=[lhi,llo,blo,bhi], aspect=1.)
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
plt.imshow(rgb, origin='lower', interpolation='nearest',
           extent=[lhi,llo,blo,bhi], aspect=1.)
#plt.axis([390, 330, -15, 15])
#xt = 330 + np.arange(5)*15
#plt.yticks([-15,-10,-5,0,5,10,15])
plt.axis([400, 320, -20, 20])
xt = 320 + np.arange(5)*20
plt.yticks(np.arange(-20, 20+0.1, 10))
plt.xticks(xt, ['%i' % (x % 360) for x in xt])
plt.xlabel('Galactic longitude $\ell$ (deg)')
plt.ylabel('Galactic latitude $b$ (deg)')
ps.savefig()


w1 -= medy1[:,np.newaxis]
w2 -= medy2[:,np.newaxis]

#lhi,llo,blo,bhi = 400, 320, -20, 20
lhi,llo,blo,bhi = 15, 345, -15, 15

ok,x1,y1 = wcs.radec2pixelxy(llo, blo)
ok,x2,y2 = wcs.radec2pixelxy(llo, bhi)
ok,x3,y3 = wcs.radec2pixelxy(lhi, blo)
ok,x4,y4 = wcs.radec2pixelxy(lhi, bhi)

xlo = int(np.floor(min(x1,x2,x3,x4)))
xhi = int(np.ceil (max(x1,x2,x3,x4)))
ylo = int(np.floor(min(y1,y2,y3,y4)))
yhi = int(np.ceil (max(y1,y2,y3,y4)))

w1 = w1[ylo:yhi, xlo:xhi]
w2 = w2[ylo:yhi, xlo:xhi]


def lbticks(wcs, xlo, ylo):
    ax = plt.axis()
    tv = [-15, -10, -5, 0, 5, 10, 15]
    xx = [wcs.radec2pixelxy(v, 0)[-2] - xlo for v in tv]
    yy = [wcs.radec2pixelxy(0, v)[-1] - ylo for v in tv]
    plt.xticks(xx, tv)
    plt.yticks(yy, tv)
    plt.xlabel('Galactic longitude $\ell$ (deg)')
    plt.ylabel('Galactic latitude $b$ (deg)')
    plt.axis(ax)


# fn = 'resid1.fits'
# if os.path.exists(fn):
#     resid1 = fitsio.read(fn)
#     rm = np.percentile(np.abs(resid1), 95)
#     plt.figure(2)
#     for fsize in [5, 9, 17, 25, 40]:
#         fw1 = median_filter(resid1, size=fsize)
#         #fw1 = gaussian_filter(fw1, 1.)
#         plt.clf()
#         plt.imshow(fw1, origin='lower', interpolation='nearest',
#                    vmin=-rm, vmax=rm, cmap='gray')
#         plt.title('%s-pixel Median-filtered W1 Residual' % fsize)
#         lbticks(wcs, xlo,ylo)
#         ps.savefig()
# 
#         #contours = np.percentile(resid1, [25, 50, 70, 85, 90, 95])
#         #for c,cc in zip(contours, ['r', (1.,0.6,0.), 'y', 'g', 'b', 'm']):
# 
#         #contours = np.percentile(resid1, [50, 60, 70, 80, 85])
#         #contours = np.percentile(resid1, [55, 60, 65, 70, 80])
#         contours = np.percentile(resid1, [70, 80])
#         #for c,cc in zip(contours, ['r', (1.,0.6,0.), 'y', 'g', 'b', 'm']):
#         #plt.contour(fw1, [c], colors=[cc], linestyles='solid')
#         plt.contour(fw1, contours, colors='r', linestyles='solid')
#         ps.savefig()
#     sys.exit(0)
        
                


# print('WCS:', wcs)
# subwcs = wcs.get_subimage(xlo, ylo, xhi-xlo, yhi-ylo)
# print('Sub wcs:', subwcs)
# print('subwcs ll,ur:', subwcs.pixelxy2radec(1,1),
#       subwcs.pixelxy2radec(xhi-xlo, yhi-ylo))
print('Subimage ll,ur:', wcs.pixelxy2radec(1+xlo,1+ylo), wcs.pixelxy2radec(xhi,yhi))

H,W = w1.shape
print('Image size', W,'x',H)
cx,cy = W/2, H/2

# Look in 5-degree wedges
xx,yy = np.meshgrid(np.arange(W), np.arange(H))
rr = (xx - cx)**2 + (yy - cy)**2
theta = np.rad2deg(np.arctan2(yy-cy, xx-cx))

print('range of theta:', theta.min(), theta.max())

tstep = 10
tt = np.arange(theta.min(), theta.max()+tstep/2., tstep)
th = []
meanw1,medw1 = [],[]
meanw2,medw2 = [],[]

rrmax = min((W/2)**2, (H/2)**2)
print('rrmax', rrmax)

# plt.clf()
# plt.imshow(rr, interpolation='nearest', origin='lower')
# plt.colorbar()
# ps.savefig()
# 
# plt.clf()
# plt.imshow(rr < rrmax, interpolation='nearest', origin='lower')
# plt.colorbar()
# ps.savefig()


for it,(tlo,thi) in enumerate(zip(tt, tt+tstep)):
    iy,ix = np.nonzero((rr < rrmax) *
                       (theta >= tlo) * (theta < thi))
    th.append((tlo+thi)/2.)
    pix = w1[iy,ix]
    meanw1.append(np.mean(pix))
    medw1.append(np.median(pix))
    pix = w2[iy,ix]
    meanw2.append(np.mean(pix))
    medw2.append(np.median(pix))

    if False and it < 10:
        plt.clf()
        w1m = np.zeros_like(w1)
        w1m[iy,ix] = w1[iy,ix]
        w2m = np.zeros_like(w1)
        w2m[iy,ix] = w2[iy,ix]
        rgb = _unwise_to_rgb([w1m, w2m], S=[S,S], Q=Q)
        plt.imshow(rgb, origin='lower', interpolation='nearest')
        plt.title('Wedge %.1f to %.1f deg' % (tlo,thi))
        ps.savefig()

rmax = np.sqrt(rrmax)
plt.clf()
rgb = _unwise_to_rgb([w1 * (rr < rrmax), w2 * (rr < rrmax)], S=[S,S], Q=Q)
plt.imshow(rgb, origin='lower', interpolation='nearest')
ax = plt.axis()
plt.plot(np.vstack([cx+np.zeros_like(tt), cx + np.cos(np.deg2rad(tt))*rmax]),
         np.vstack([cy+np.zeros_like(tt), cy + np.sin(np.deg2rad(tt))*rmax]),
         'b-')
plt.axis(ax)
plt.title('Wedges')
ps.savefig()

plt.clf()
plt.plot(th, meanw1, 'b-', label='W1 mean')
plt.plot(th, medw1, 'b--', label='W1 median')
plt.plot(th, meanw2, 'r-', label='W2 mean')
plt.plot(th, medw2, 'r--', label='W2 median')
plt.xlabel('Wedge angle (deg)')
plt.ylabel('WISE flux')
plt.legend()
plt.title('WISE flux in wedges')
ps.savefig()

mn1,mx1 = np.percentile(w1, [25,98])
# plt.clf()
# plt.imshow(w1 * (rr < rrmax), origin='lower', interpolation='nearest',
#            vmin=mn1, vmax=mx1)
# # ax = plt.axis()
# # plt.plot(np.vstack([cx+np.zeros_like(tt), cx + np.cos(np.deg2rad(tt))*rmax]),
# #          np.vstack([cy+np.zeros_like(tt), cy + np.sin(np.deg2rad(tt))*rmax]),
# #          'b-')
# # plt.axis(ax)
# ps.savefig()
# 
# plt.colorbar()
# ps.savefig()

from tractor import *
from tractor.galaxy import *

tim = Image(data=np.zeros_like(w1), inverr=np.ones_like(w1),
            psf=NCircularGaussianPSF([1.],[1.]))
gal = ExpGalaxy(PixPos(cx, cy), Flux(w1.sum()),
                GalaxyShape(200, 0.75, 90.))
tractor = Tractor([tim],[gal])
# mod = tractor.getModelImage(0)
# 
# plt.clf()
# plt.imshow(mod, interpolation='nearest', origin='lower')
# plt.colorbar()
# ps.savefig()
# 
# plt.clf()
# plt.imshow(mod, vmin=mn1, vmax=mx1, interpolation='nearest', origin='lower')
# plt.colorbar()
# ps.savefig()

w1masked = w1.copy()
slc = slice(cy - 25, cy + 25)
w1masked[slc, :] = 0.
ie = np.ones_like(w1)
ie[slc,:] = 0.

tim.data = w1masked
tim.inverr = ie


plt.figure(2)
plt.clf()
plt.imshow(w1masked, origin='lower', interpolation='nearest',
           vmin=mn1, vmax=mx1, cmap='gray')
#plt.colorbar()
plt.title('W1 data')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.figure(1)



mn2,mx2 = np.percentile(w1, [50,99])

plt.figure(2)
for fsize in [9, 25]: #[5,9,13,17,25]:
    fw1 = median_filter(w1, size=fsize)
    #fw1 = gaussian_filter(fw1, 1.)
    plt.clf()
    #pix = fw1 * (rr < rrmax)
    pix = fw1
    plt.imshow(pix, origin='lower', interpolation='nearest',
               vmin=mn2, vmax=mx2, cmap='gray')
    #plt.contour(pix, contours, color='k')
    plt.title('%s-pixel Median filtered W1' % fsize)

    lbticks(wcs, xlo,ylo)

    ps.savefig()

contours = np.percentile(w1, [50, 70, 85, 90, 95])
#plt.contour(pix, contours, colors='k')
plt.contour(pix, contours, colors='r')
ps.savefig()




tractor.freezeParam('images')
for step in range(50):
    dlnp,x,alpha = tractor.optimize()
    print('dlnp', dlnp)
    print('x', x)
    print('alpha', alpha)
    print('Galaxy', gal)
    if dlnp == 0:
        break

mod1 = tractor.getModelImage(0)
resid1 = w1 - mod1

# plt.clf()
# plt.imshow(mod1, interpolation='nearest', origin='lower')
# #plt.colorbar()
# ps.savefig()

plt.figure(2)

plt.clf()
plt.imshow(mod1 * ie, vmin=mn1, vmax=mx1, interpolation='nearest', origin='lower', cmap='gray')
plt.title('W1 model')
lbticks(wcs, xlo,ylo)
ps.savefig()

plt.figure(1)

#rkwa = dict(vmin=-1000, vmax=1000, origin='lower', interpolation='nearest')

rm = np.percentile(np.abs(resid1), 95)

plt.figure(2)
plt.clf()
#plt.imshow(resid1, **rkwa)
#plt.imshow(resid1, vmin=mn1, vmax=mx1, interpolation='nearest', origin='lower')
plt.imshow(resid1, vmin=-rm, vmax=rm, interpolation='nearest',
           origin='lower', cmap='gray')
plt.title('W1 Residual')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.figure(1)


fitsio.write('resid1.fits', resid1, clobber=True)

plt.figure(2)
#for fsize in [5, 9, 17, 25]:
for fsize in [40]:
    fw1 = median_filter(resid1, size=fsize)
    plt.clf()
    plt.imshow(fw1, origin='lower', interpolation='nearest',
               vmin=-rm, vmax=rm, cmap='gray')
    plt.title('%s-pixel Median-filtered W1 Residual' % fsize)

    lbticks(wcs, xlo,ylo)

    ps.savefig()

    # contours = np.percentile(fw1, [25, 50, 70, 85, 90, 95])
    # #plt.contour(pix, contours, colors='k')
    # plt.contour(fw1, contours, colors='r')

    #contours = np.percentile(resid1, [70, 80])
    contours = np.percentile(resid1, [60, 70, 80])
    plt.contour(fw1, contours, colors='r', linestyles='solid')

    ps.savefig()
plt.figure(1)

sys.exit(0)
    

# plt.clf()
# plt.imshow(resid2, **rkwa)
# plt.title('W2 Residual')
# ps.savefig()

w2masked = w2.copy()
slc = slice(cy - 25, cy + 25)
w2masked[slc, :] = 0.
tim.data = w2masked

for step in range(50):
    dlnp,x,alpha = tractor.optimize()
    print('dlnp', dlnp)
    print('x', x)
    print('alpha', alpha)
    print('Galaxy', gal)
    if dlnp == 0:
        break

mod2 = tractor.getModelImage(0)
resid2 = w2 - mod2

RR = np.linspace(0, rmax, 6)**2
residslices = []
for rlo,rhi in zip(RR, RR[1:]):
    meanw1,medw1 = [],[]
    meanw2,medw2 = [],[]
    for it,(tlo,thi) in enumerate(zip(tt, tt+tstep)):
        iy,ix = np.nonzero((rr >= rlo) * (rr < rhi) *
                           (theta >= tlo) * (theta < thi))
        pix = resid1[iy,ix]
        meanw1.append(np.mean(pix))
        medw1.append(np.median(pix))
        pix = resid2[iy,ix]
        meanw2.append(np.mean(pix))
        medw2.append(np.median(pix))
    residslices.append((meanw1, medw1, meanw2, medw2))
plt.clf()
for i,(meanw1,medw1, meanw2,medw2) in enumerate(residslices):
    cc = 'bcgmr'
    plt.plot(th, np.array(medw1) + 1000*i, '-', label='W1 R%i' % i, color=cc[i])
    #plt.plot(th, np.array(medw2) + 1000*i, '--', color=cc[i])
plt.xlabel('Wedge angle (deg)')
plt.ylabel('WISE flux')
#plt.title('R slices, resid, medians')
plt.title('R slices, W1 resid, medians')
plt.legend()
ps.savefig()




w1mean = np.zeros_like(w1)
w1med  = np.zeros_like(w1)
w2mean = np.zeros_like(w2)
w2med  = np.zeros_like(w2)

RR = np.linspace(0, rmax, 6)**2
rslices = []
for rlo,rhi in zip(RR, RR[1:]):
    meanw1,medw1 = [],[]
    meanw2,medw2 = [],[]
    
    for it,(tlo,thi) in enumerate(zip(tt, tt+tstep)):
        iy,ix = np.nonzero((rr >= rlo) * (rr < rhi) *
                           (theta >= tlo) * (theta < thi))
        pix = w1[iy,ix]
        meanw1.append(np.mean(pix))
        medw1.append(np.median(pix))
        pix = w2[iy,ix]
        meanw2.append(np.mean(pix))
        medw2.append(np.median(pix))

        w1mean[iy,ix] = meanw1[-1]
        w1med [iy,ix] =  medw1[-1]
        w2mean[iy,ix] = meanw2[-1]
        w2med [iy,ix] =  medw2[-1]
        
    rslices.append((meanw1, medw1, meanw2, medw2))

# plt.clf()
# plt.imshow(w1mean, vmin=mn1, vmax=mx1, origin='lower', interpolation='nearest')
# plt.title('wedge mean w1')
# ps.savefig()
# plt.clf()
# plt.imshow(w1med, vmin=mn1, vmax=mx1, origin='lower', interpolation='nearest')
# plt.title('wedge median w1')
# ps.savefig()

plt.clf()
plt.imshow(w1mean, origin='lower', interpolation='nearest')
plt.title('wedge mean W1')
ps.savefig()
plt.clf()
plt.imshow(w1med, origin='lower', interpolation='nearest')
plt.title('wedge median W1')
ps.savefig()

plt.clf()
plt.imshow(w2mean, origin='lower', interpolation='nearest')
plt.title('wedge mean W2')
ps.savefig()
plt.clf()
plt.imshow(w2med, origin='lower', interpolation='nearest')
plt.title('wedge median W2')
ps.savefig()


plt.clf()
rgb = _unwise_to_rgb([w1mean, w2mean], S=[S,S], Q=Q)
plt.imshow(rgb, origin='lower', interpolation='nearest')
plt.title('wedge mean')
ps.savefig()
plt.clf()
rgb = _unwise_to_rgb([w1med, w2med], S=[S,S], Q=Q)
plt.imshow(rgb, origin='lower', interpolation='nearest')
plt.title('wedge median')
ps.savefig()

plt.figure(1)


for i,(meanw1,medw1, meanw2,medw2) in enumerate(rslices):
    plt.clf()
    #plt.plot(th, meanw1, 'b-')
    #plt.plot(th, medw1, 'b--')
    #plt.plot(th, meanw2, 'r-')
    #plt.plot(th, medw2, 'r--')
    plt.plot(th, medw1, 'b-')
    plt.plot(th, medw2, 'r-')
    plt.xlabel('Wedge angle (deg)')
    plt.ylabel('WISE flux')
    plt.title('R slice %i, medians' % (i))
    ps.savefig()


plt.clf()
for i,(meanw1,medw1, meanw2,medw2) in enumerate(rslices):
    plt.plot(th, medw1, 'b-')
    plt.plot(th, medw2, 'r-')
plt.xlabel('Wedge angle (deg)')
plt.ylabel('WISE flux')
plt.title('R slice %i, medians' % (i))
ps.savefig()




import sys
sys.exit(0)

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
