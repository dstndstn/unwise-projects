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



def resid_rgb(resid1, resid2):
    S,Q = 3000,25
    alpha = 5.
    
    w1,w2 = resid1,resid2
    S1,S2 = S,S
    b = w1 / S1
    r = w2 / S2
    g = (r + b) / 2.

    I = (r+g+b)/3.
    fI = np.arcsinh(alpha * Q * I) / np.sqrt(Q)

    R = fI * r / I
    G = fI * g / I
    B = fI * b / I

    RGB = np.dstack([R,G,B])
    RGB = (np.clip((RGB + 1.) / 2., 0., 1.) * 255.99).astype(np.uint8)
    return RGB


from tractor import Tractor
class AsymmetricTractor(Tractor):
    def getLogLikelihood(self):
        chisq = 0.
        for i,chi in enumerate(self.getChiImages()):
            # chi = (img - mod) / error
            #chisq += (chi.astype(float) ** 2).sum()

            # positive residuals: unexplained flux penalized less
            I = (chi > 0)
            chisq += ((0.5 * chi[I]) ** 2).sum()
            I = (chi < 0)
            chisq += ((1.0 * chi[I]) ** 2).sum()
        return -0.5 * chisq

if False:
    from astrometry.util.file import *
    from tractor import *
    from tractor.galaxy import *

    disable_galaxy_cache()
    ps = PlotSequence('xb')

    X = unpickle_from_file('sample.pickle')

    tr = X['tractor']
    #at = AsymmetricTractor(tr.images, tr.catalog)
    #at.freezeParam('images')
    #tractor = at

    tractor = tr
    
    tractor.printThawedParams()

    import emcee

    w1 = tractor.images[0].data
    w2 = tractor.images[1].data
    wcs = anwcs('wcs.fits')
    xlo,ylo = 1362, 450
    

    # Create emcee sampler
    nw = 30
    ndim = len(tractor.getParams())
    print('N dim:', ndim)
    
    sampler = emcee.EnsembleSampler(nw, ndim, tractor)

    p0 = tractor.getParams()
    #std = at.getStepSizes()
    std = np.array([1.0, 1.0, 1e7, 1e7, 0.01, 0.01, 0.01])
    pp = emcee.utils.sample_ball(p0, std, size=nw)

    print('Fitting params: (%i):' % len(p0))
    tractor.printThawedParams()

    print('Step sizes:', std)
    
    nsteps = 100
    
    allpp = np.zeros((nsteps, nw, ndim), np.float32)
    alllnp = np.zeros((nsteps, nw), np.float32)

    #pickle_to_file(dict(allpp=allpp, alllnp=alllnp, tractor=tractor),
    #               'sample2.pickle')

    rstate = None
    lnp = None
    for step in range(nsteps):
        print('Taking step', step)
        pp,lnp,rstate = sampler.run_mcmc(pp, 1, lnprob0=lnp, rstate0=rstate)
        print('Max lnprob:', np.max(lnp))
        i = np.argmax(lnp.ravel())
        pbest = pp[i,:]
        print('Best params:', pbest)

        tractor.setParams(pbest)
        mod1 = tractor.getModelImage(0)
        resid1 = w1 - mod1
        mod2 = tractor.getModelImage(1)
        resid2 = w2 - mod2
        
        rgb = resid_rgb(resid1, resid2)
        plt.clf()
        dimshow(rgb)
        plt.title('Residuals')
        lbticks(wcs, xlo,ylo)
        plt.savefig('resid.png')

        chi1 = resid1 * tractor.images[0].getInvError()
        chi2 = resid2 * tractor.images[1].getInvError()
        
        plt.clf()
        plt.hist(chi1[chi1 != 0].ravel(), range=(-10, 10), bins=100,
                 histtype='step', color='b')
        plt.hist(chi2[chi2 != 0].ravel(), range=(-10, 10), bins=100,
                 histtype='step', color='r')
        plt.savefig('residhist.png')
        
        allpp[step,:,:] = pp
        alllnp[step,:] = lnp

        plt.figure(1)
        plt.clf()
        plt.plot(alllnp[:step+1,:], 'b', alpha=0.5)
        #mx = alllnp.max()
        #plt.ylim(mx-20, mx+5)
        plt.title('log posterior')
        plt.savefig('lnp.png')

        if step % 10 == 9:
            import triangle
            X = allpp[:step+1, :,:].reshape(((step+1) * nw, ndim))
            print('std in X:', np.std(X, axis=0))
            plt.figure(2)
            plt.clf()
            triangle.corner(X, labels=tractor.getParamNames(), plot_contours=False)
            plt.savefig('corner.png')
            plt.clf()
            plt.figure(1)

    pickle_to_file(dict(allpp=allpp, alllnp=alllnp, tractor=tractor),
                   'sample3.pickle')

    
    
    sys.exit(0)
    
        
if False:
    import emcee
    from tractor import *
    from tractor.galaxy import *

    disable_galaxy_cache()
    
    ps = PlotSequence('xb')
    
    # ExpGalaxy at pixel (423.64, 444.09) with Fluxes: w1=9.06254e+08, w2=1.02335e+09 and Galaxy Shape: re=400.31, ab=0.38, phi=89.8
    wcs = anwcs('wcs.fits')
    xlo,ylo = 1362, 450
    mask   = fitsio.read('mask.fits').astype(np.float32)
    w1 = fitsio.read('data-w1.fits')
    w2 = fitsio.read('data-w2.fits')

    ie1 = np.zeros_like(w1)
    ie2 = np.zeros_like(w2)

    for img,ie in [(w1,ie1), (w2,ie2)]:
        # # Estimate per-pixel noise via Blanton's 5-pixel MAD
        slice1 = (slice(0,-5,10),slice(0,-5,10))
        slice2 = (slice(5,None,10),slice(5,None,10))
        mad = np.median(np.abs(img[slice1] - img[slice2]).ravel())
        sig1 = 1.4826 * mad / np.sqrt(2.)
        print('sig1 estimate:', sig1)
        ie += 1. / sig1
        ie[mask == 0] = 0
    
    tim1 = Image(data=w1, inverr=ie1,
                 psf=NCircularGaussianPSF([1.],[1.]),
                 photocal=LinearPhotoCal(1., 'w1'))
    tim2 = Image(data=w2, inverr=ie2,
                 psf=NCircularGaussianPSF([1.],[1.]),
                 photocal=LinearPhotoCal(1., 'w2'))

    # ExpGalaxy at pixel (423.64, 444.09) with Fluxes: w1=9.06254e+08, w2=1.02335e+09 and Galaxy Shape: re=400.31, ab=0.38, phi=89.8
    cx,cy = 423.64, 444.09
    flux1 = 9.06e8
    flux2 = 1.02e9

    #shape = GalaxyShape(400.31, 0.38, 89.8)
    shape = EllipseESoft.fromRAbPhi(400.31, 0.38, 89.8)
    
    gal = ExpGalaxy(PixPos(cx, cy), Fluxes(w1=flux1, w2=flux2), shape)

    tractor = Tractor([tim1, tim2],[gal])
    tractor.freezeParam('images')

    # Create emcee sampler
    nw = 30
    ndim = len(tractor.getParams())
    print('N dim:', ndim)
    
    sampler = emcee.EnsembleSampler(nw, ndim, tractor)

    p0 = tractor.getParams()
    #std = tractor.getStepSizes()
    std = np.array([1.0, 1.0, 1e7, 1e7, 0.01, 0.01, 0.01])
    pp = emcee.utils.sample_ball(p0, std, size=nw)

    print('Fitting params: (%i):' % len(p0))
    tractor.printThawedParams()

    print('Step sizes:', std)
    
    nsteps = 100
    
    allpp = np.zeros((nsteps, nw, ndim), np.float32)
    alllnp = np.zeros((nsteps, nw), np.float32)

    from astrometry.util.file import *
    pickle_to_file(dict(allpp=allpp, alllnp=alllnp, tractor=tractor),
                   'sample.pickle')

    rstate = None
    lnp = None
    for step in range(nsteps):
        print('Taking step', step)
        pp,lnp,rstate = sampler.run_mcmc(pp, 1, lnprob0=lnp, rstate0=rstate)
        print('Max lnprob:', np.max(lnp))
        i = np.argmax(lnp.ravel())
        pbest = pp[i,:]
        print('Best params:', pbest)

        tractor.setParams(pbest)
        mod1 = tractor.getModelImage(0)
        resid1 = w1 - mod1
        mod2 = tractor.getModelImage(1)
        resid2 = w2 - mod2
        
        rgb = resid_rgb(resid1, resid2)
        plt.clf()
        dimshow(rgb)
        plt.title('Residuals')
        lbticks(wcs, xlo,ylo)
        plt.savefig('resid.png')
        
        allpp[step,:,:] = pp
        alllnp[step,:] = lnp

        plt.figure(1)
        plt.clf()
        plt.plot(alllnp[:step+1,:], 'b', alpha=0.5)
        #mx = alllnp.max()
        #plt.ylim(mx-20, mx+5)
        plt.title('log posterior')
        plt.savefig('lnp.png')

        if step % 10 == 9:
            import triangle
            X = allpp[:step+1, :,:].reshape(((step+1) * nw, ndim))
            print('std in X:', np.std(X, axis=0))
            plt.figure(2)
            plt.clf()
            triangle.corner(X, labels=tractor.getParamNames(), plot_contours=False)
            plt.savefig('corner.png')
            plt.clf()
            plt.figure(1)

    from astrometry.util.file import *
    pickle_to_file(dict(allpp=allpp, alllnp=alllnp, tractor=tractor),
                   'sample.pickle')
            
    sys.exit(0)

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

    if True:
        for t in T[I]:
            brick = t.coadd_id
            for band in [1,2,3,4]:
                dirnm = '%s/%s/%s_ac51' % (brick[:2], brick[:4], brick)
                #print('mkdir -p %s' % dirnm)
                print('wget --continue -r -nH --cut-dirs 4 http://irsa.ipac.caltech.edu/ibe/data/wise/merge/merge_p3am_cdd/%s/%s_ac51-w%i-int-3.fits'
                      % (dirnm, brick, band))

    T.cut(I)
    T.writeto('wisex-atlas.fits')

    sys.exit(0)


#fns = ['w1-lbzoom-5-1800-u-wcs.fits',
#       'w2-lbzoom-5-1800-u-wcs.fits',]

#fns = ['unwise2-w1-lbzoom-5-1800.fits',
#       'unwise2-w2-lbzoom-5-1800.fits',]

fns = ['unwise-w1-lbzoom-5-1800.fits',
       'unwise-w2-lbzoom-5-1800.fits',]

#fns = ['unwise2-w1-lbzoom-5-1800.fits',
#       'unwise2-w2-lbzoom-5-1800.fits',]

imgs = [fitsio.read(fn) for fn in fns]
wcs = anwcs(fns[0])
print('WCS header', wcs)

# HACK -- colormap scalings
for img in imgs:
    img /= 10.

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
plt.subplots_adjust(left=0.11, right=0.96, bottom=0.1, top=0.95)

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
w1orig = w1.copy()
w2orig = w2.copy()

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

print('xlo,ylo', xlo, ylo)

w1 = w1[ylo:yhi, xlo:xhi]
w2 = w2[ylo:yhi, xlo:xhi]



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

plt.figure(2)

from tractor import *
from tractor.galaxy import *

w1masked = w1orig.copy()[ylo:yhi, xlo:xhi]
w2masked = w2orig.copy()[ylo:yhi, xlo:xhi]

print('masked size W1:', w1masked.shape)
print('masked size W2:', w2masked.shape)

w1mag = -2.5*(np.log10(w1masked) - 9.)
w2mag = -2.5*(np.log10(w2masked) - 9.)

print(np.sum(np.isfinite(w1mag)), 'pixels finite')
print(np.sum(np.isfinite(w2mag)), 'pixels finite')

plt.clf()
plt.hist((w1mag - w2mag).ravel(), range=(-1,1), bins=100)
plt.xlabel('W1 - W2')
ps.savefig()

# plt.clf()
# loghist((w1mag - w2mag).ravel(), w1mag.ravel(), nbins=100,
#         range=((-1,1),np.percentile(w1mag.ravel(), [10,90])))
# plt.xlabel('W1 - W2')
# plt.ylabel('W1')
# ps.savefig()

# plt.clf()
# plt.imshow(np.abs((w1mag - w2mag) - 0) > 0.5)
# plt.title('Color mask')
# ps.savefig()

from astrometry.util.util import median_smooth
cc = w1mag - w2mag
goodcolor = np.isfinite(cc)

#cf = median_filter(w1mag - w2mag, size=5)
#cf = np.zeros_like(cc)
#median_smooth(cc, np.logical_not(goodcolor), 2, cf)
#mlo,mhi = np.percentile(cf, [5,95])
#mlo,mhi = np.percentile(cf[goodcolor * np.isfinite(cf)], [5,95])
mlo,mhi = np.percentile(cc[goodcolor], [5,95])

print('W1 - W2 color masks:', mlo,mhi)

print('Bad colors:', np.sum(np.logical_not(goodcolor)))
print('Colors cut:', np.sum(goodcolor * np.logical_or(cc < mlo, cc > mhi)))

#mask = np.isfinite(cf) * (cf > mlo) * (cf < mhi)
mask = goodcolor * (cc > mlo) * (cc < mhi)
#mask = (cf > mlo) * (cf < mhi)
#mask = (np.abs(cf) < 0.5)

# plt.clf()
# plt.imshow(np.abs(cf) > 0.5)
# plt.title('Color mask (2)')
# ps.savefig()

plt.clf()
rgb = _unwise_to_rgb([w1masked, w2masked], S=[S,S], Q=Q)
dimshow(rgb)
plt.title('Data')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-data.pdf')

plt.clf()
rgb = _unwise_to_rgb([w1masked * mask, w2masked * mask], S=[S,S], Q=Q)
dimshow(rgb)
lbticks(wcs, xlo,ylo)
plt.title('Data (masked)')
ps.savefig()
plt.savefig('xbulge-fit-masked.pdf')

# slc = slice(cy - 25, cy + 25)
# w1masked[slc, :] = 0.
# ie1 = np.ones_like(w1masked)
# ie1[slc,:] = 0.
# 
# slc = slice(cy - 25, cy + 25)
# w2masked[slc, :] = 0.
# ie2 = np.ones_like(w2masked)
# ie2[slc,:] = 0.

ie = mask.astype(np.float32)
ie1 = ie2 = ie

fitsio.write('w1masked.fits', w1masked * mask, clobber=True)
fitsio.write('w2masked.fits', w2masked * mask, clobber=True)
#plt.imsave('rgb-masked.jpg', rgb)


tim1 = Image(data=w1masked, inverr=ie1,
             psf=NCircularGaussianPSF([1.],[1.]),
             photocal=LinearPhotoCal(1., 'w1'))
tim2 = Image(data=w2masked, inverr=ie2,
             psf=NCircularGaussianPSF([1.],[1.]),
             photocal=LinearPhotoCal(1., 'w2'))
gal = ExpGalaxy(PixPos(cx, cy), Fluxes(w1=w1.sum(), w2=w2.sum()),
                GalaxyShape(200, 0.75, 90.))
tractor = Tractor([tim1, tim2],[gal])

fitsio.write('data-w1.fits', w1masked, clobber=True)
fitsio.write('data-w2.fits', w2masked, clobber=True)
wcs.writeto('wcs.fits')



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

# plt.clf()
# plt.imshow(w1masked, origin='lower', interpolation='nearest',
#            vmin=mn1, vmax=mx1, cmap='gray')
# #plt.colorbar()
# plt.title('W1 data')
# lbticks(wcs, xlo,ylo)
# ps.savefig()

# rgb = _unwise_to_rgb([w1masked, w2masked], S=[S,S], Q=Q)
# plt.clf()
# dimshow(rgb)
# plt.title('Masked data')
# lbticks(wcs, xlo,ylo)
# ps.savefig()
#plt.figure(1)


# mn2,mx2 = np.percentile(w1, [50,99])
# 
# #plt.figure(2)
# for fsize in [9, 25]: #[5,9,13,17,25]:
#     fw1 = median_filter(w1, size=fsize)
#     #fw1 = gaussian_filter(fw1, 1.)
#     plt.clf()
#     #pix = fw1 * (rr < rrmax)
#     pix = fw1
#     plt.imshow(pix, origin='lower', interpolation='nearest',
#                vmin=mn2, vmax=mx2, cmap='gray')
#     #plt.contour(pix, contours, color='k')
#     plt.title('%s-pixel Median filtered W1' % fsize)
# 
#     lbticks(wcs, xlo,ylo)
# 
#     ps.savefig()
# 
# contours = np.percentile(w1, [50, 70, 85, 90, 95])
# #plt.contour(pix, contours, colors='k')
# plt.contour(pix, contours, colors='r')
# ps.savefig()


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
resid1 = w1masked - mod1
mod2 = tractor.getModelImage(1)
resid2 = w2masked - mod2


#plt.figure(2)

rgb = _unwise_to_rgb([mod1, mod2], S=[S,S], Q=Q)
plt.clf()
dimshow(rgb)
plt.title('Model')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-model.pdf')

rgb = resid_rgb(resid1, resid2)
plt.clf()
dimshow(rgb)
plt.title('Residuals')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-resid.pdf')

rgb = resid_rgb(resid1*mask, resid2*mask)
plt.clf()
dimshow(rgb)
plt.title('Residuals (masked)')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-residmasked.pdf')

print('Mask:', mask.dtype)
fitsio.write('resid1.fits', resid1, clobber=True)
fitsio.write('resid2.fits', resid2, clobber=True)
fitsio.write('mask.fits', mask.astype(np.uint8), clobber=True)

fr1 = median_filter(resid1*mask, size=50)
fr2 = median_filter(resid2*mask, size=50)

rgb = resid_rgb(fr1, fr2)
plt.clf()
dimshow(rgb)
plt.title('Residuals (smoothed)')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-smooth.pdf')

from astrometry.util.util import median_smooth
fr1 = np.zeros_like(resid1)
fr2 = np.zeros_like(resid2)
median_smooth(resid1, np.logical_not(mask), 25, fr1)
median_smooth(resid2, np.logical_not(mask), 25, fr2)
#median_smooth(resid1, np.logical_not(mask), 15, fr1)
#median_smooth(resid2, np.logical_not(mask), 15, fr2)

rgb = resid_rgb(fr1, fr2)
plt.clf()
dimshow(rgb)
plt.title('Residuals (smoothed)')
lbticks(wcs, xlo,ylo)
ps.savefig()
plt.savefig('xbulge-fit-smooth2.pdf')




# plt.clf()
# plt.imshow(mod1, interpolation='nearest', origin='lower')
# #plt.colorbar()
# ps.savefig()

# plt.clf()
# plt.imshow(mod1 * ie1, vmin=mn1, vmax=mx1, interpolation='nearest', origin='lower', cmap='gray')
# plt.title('W1 model')
# lbticks(wcs, xlo,ylo)
# ps.savefig()
# 
# plt.figure(1)

#rkwa = dict(vmin=-1000, vmax=1000, origin='lower', interpolation='nearest')

resid1 *= mask

rm = np.percentile(np.abs(resid1), 95)

#plt.figure(2)
plt.clf()
#plt.imshow(resid1, **rkwa)
#plt.imshow(resid1, vmin=mn1, vmax=mx1, interpolation='nearest', origin='lower')
plt.imshow(resid1, vmin=-rm, vmax=rm, interpolation='nearest',
           origin='lower', cmap='gray')
plt.title('W1 Residual')
lbticks(wcs, xlo,ylo)
ps.savefig()
#plt.figure(1)


fitsio.write('resid1.fits', resid1, clobber=True)

#plt.figure(2)
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

#plt.figure(1)
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
