import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits

## Open fits file
hdul = fits.open('TotalDat.fits')

"""
Unpack the data - I am sure there is a more elegant way
But I am a dinosaur
"""

hdr = hdul[1].data

### Read database ID
DBID = hdr['DBID']

### Read redshift
z    = hdr['Z']

### Read bolometric luminosities and errors
lbol = hdr['log_LBOL']
lbol_err = hdr['log_LBOL_ERR']

### Read in g, r, i data

### g = 472.0
### r = 641.5
### i = 783.5

### Read in rest frame wavelengths
lr_g = hdr['LAMBDA_REST_g']
lr_r = hdr['LAMBDA_REST_r']
lr_i = hdr['LAMBDA_REST_i']

### Read in the observation dates
mjd_g = hdr['MJD_g']
mjd_r = hdr['MJD_r']
mjd_i = hdr['MJD_i']

### Read in the photometric data and errors
mag_g = hdr['MAG_g']
mag_r = hdr['MAG_r']
mag_i = hdr['MAG_i']

mag_err_g = hdr['MAG_ERR_g']
mag_err_r = hdr['MAG_ERR_r']
mag_err_i = hdr['MAG_ERR_i']

## Read in the DRW timescales

tau_o_g = hdr['log_TAU_OBS_g']
tau_o_r = hdr['log_TAU_OBS_r']
tau_o_i = hdr['log_TAU_OBS_i']

tau_o_gl= hdr['log_TAU_OBS_g_ERR_L']
tau_o_rl= hdr['log_TAU_OBS_r_ERR_L']
tau_o_il= hdr['log_TAU_OBS_i_ERR_L']

tau_o_gu= hdr['log_TAU_OBS_g_ERR_U']
tau_o_ru= hdr['log_TAU_OBS_r_ERR_U']
tau_o_iu= hdr['log_TAU_OBS_i_ERR_U']

tau_g_err = np.vstack((tau_o_gl,tau_o_gu))
tau_r_err = np.vstack((tau_o_rl,tau_o_ru))
tau_i_err = np.vstack((tau_o_il,tau_o_iu))

tau_r_g = hdr['log_TAU_REST_g']
tau_r_r = hdr['log_TAU_REST_r']
tau_r_i = hdr['log_TAU_REST_i']

tau_r_gl= hdr['log_TAU_REST_g_ERR_L']
tau_r_rl= hdr['log_TAU_REST_r_ERR_L']
tau_r_il= hdr['log_TAU_REST_i_ERR_L']

tau_r_gu= hdr['log_TAU_REST_g_ERR_U']
tau_r_ru= hdr['log_TAU_REST_r_ERR_U']
tau_r_iu= hdr['log_TAU_REST_i_ERR_U']

## Close fits file
hdul.close()

print('Data read and file closed')

### Make masks for the nans
ind_o_g = ~np.isnan( mjd_g )
ind_o_r = ~np.isnan( mjd_r )
ind_o_i = ~np.isnan( mjd_i )

"""
Let's make some test plots
"""

###
### Diagnostic plots
###

### g = 472.0
### r = 641.5
### i = 783.5

### Plot 9
lr_g_o = 4720/(1+z)

plt.figure()
plt.plot(z,(lr_g - lr_g_o)*100/lr_g_o,'o')

plt.xlabel(r'Redshift')
plt.ylabel(r'$(\lambda_r - \lambda_o) * 100 / \lambda_o$')
plt.title(r'g-band : $\lambda_o = 4720/(1+z)$')
plt.savefig('z_g_rest.png',dpi=300)
plt.close('all')

### Plot 10
lr_r_o = 6415/(1+z)

plt.figure()
plt.plot(z,(lr_r - lr_r_o)*100/lr_r_o,'o')

plt.xlabel(r'Redshift')
plt.ylabel(r'$(\lambda_r - \lambda_o) * 100 / \lambda_o$')
plt.title(r'r-band : $\lambda_o = 6415/(1+z)$')
plt.savefig('z_r_rest.png',dpi=300)
plt.close('all')

### Plot 11
lr_i_o = 7835/(1+z)

plt.figure()
plt.plot(z,(lr_i - lr_i_o)*100/lr_i_o,'o')

plt.xlabel(r'Redshift')
plt.ylabel(r'$(\lambda_r - \lambda_o) * 100 / \lambda_o$')
plt.title(r'i-band : $\lambda_o = 7835/(1+z)$')
plt.savefig('z_i_rest.png',dpi=300)
plt.close('all')


### Plot 1
plt.figure()
plt.errorbar(lr_r_o,lbol,yerr=lbol_err,fmt='r.')
plt.errorbar(lr_g_o,lbol,yerr=lbol_err,fmt='g.')
plt.errorbar(lr_i_o,lbol,yerr=lbol_err,fmt='m.')

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('RestFrame.png',dpi=300)
plt.close('all')

### Plot 2
plt.figure()
plt.scatter(lr_r_o,lbol,c=z,vmin=0,vmax=4.2,cmap='rainbow')
plt.scatter(lr_g_o,lbol,c=z,vmin=0,vmax=4.2,cmap='rainbow')
plt.scatter(lr_i_o,lbol,c=z,vmin=0,vmax=4.2,cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('RestFrame_z.png',dpi=300)
plt.close('all')

### Plot 3
plt.figure()
plt.scatter(lr_r_o,lbol,c=tau_o_r,vmin=tau_o_r.min(),vmax=tau_o_r.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_r.png',dpi=300)
plt.close('all')

### Plot 4
plt.figure()
plt.scatter(lr_g_o,lbol,c=tau_o_g,vmin=tau_o_g.min(),vmax=tau_o_g.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_g.png',dpi=300)
plt.close('all')

### Plot 5
plt.figure()
plt.scatter(lr_i_o,lbol,c=tau_o_i,vmin=tau_o_i.min(),vmax=tau_o_i.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_i.png',dpi=300)
plt.close('all')

### Plot 6
plt.figure()
plt.scatter(lr_r_o,lbol,c=tau_r_r,vmin=tau_r_r.min(),vmax=tau_r_r.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_r_rest.png',dpi=300)
plt.close('all')

### Plot 7
plt.figure()
plt.scatter(lr_g_o,lbol,c=tau_r_g,vmin=tau_r_g.min(),vmax=tau_r_g.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_g_rest.png',dpi=300)
plt.close('all')

### Plot 8
plt.figure()
plt.scatter(lr_i_o,lbol,c=tau_r_i,vmin=tau_r_i.min(),vmax=tau_r_i.max(),cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Tau_i_rest.png',dpi=300)
plt.close('all')

### Plot 2
plt.figure()
plt.scatter(lr_r_o,lbol,c=tau_o_r,vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_g_o,lbol,c=tau_o_g,vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_i_o,lbol,c=tau_o_i,vmin=1.5,vmax=4.5,cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('All_tau_o.png',dpi=300)
plt.close('all')

### This is Plot 1 in the paper
plt.figure()
plt.scatter(lr_r_o,lbol,c=tau_o_r,vmin=2.2,vmax=4.4,cmap='rainbow',zorder=10,alpha=0.8,edgecolors='k')
plt.scatter(lr_g_o,lbol,c=tau_o_g,vmin=2.2,vmax=4.4,cmap='rainbow',zorder=10,alpha=0.8,edgecolors='k')
plt.scatter(lr_i_o,lbol,c=tau_o_i,vmin=2.2,vmax=4.4,cmap='rainbow',zorder=10,alpha=0.8,edgecolors='k')

cbar = plt.colorbar(pad=0.01,ticks=[2.5,3.0,3.5,4.0])
cbar.ax.set_yticklabels([2.5,3.0,3.5,4.0],rotation=270,va='center')
cbar.set_label(r'$log_{10} (  \tau_{DRW} / days )$', rotation=270, labelpad=17,fontsize=12)



"""
OK - Let's do some analysis
"""

lbol_l = 46.7
lbol_u = 47.2

lam_l =  900.0
lam_u = 1900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g1 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r1 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i1 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 46.2
lbol_u = 46.7

lam_l =  900.0
lam_u = 1900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g2 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r2 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i2 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 46.2
lbol_u = 46.7

lam_l = 1900.0
lam_u = 2900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g3 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r3 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i3 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 46.7
lbol_u = 47.2

lam_l = 1900.0
lam_u = 2900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g4 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r4 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i4 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 46.2
lbol_u = 46.7

lam_l = 2900.0
lam_u = 3900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g5 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r5 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i5 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.7
lbol_u = 46.2

lam_l =  900.0
lam_u = 1900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g6 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r6 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i6 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.7
lbol_u = 46.2

lam_l = 1900.0
lam_u = 2900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g7 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r7 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i7 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.7
lbol_u = 46.2

lam_l = 2900.0
lam_u = 3900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g8 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r8 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i8 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.7
lbol_u = 46.2

lam_l = 3900.0
lam_u = 4900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_g9 = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_r9 = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_i9 = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.2
lbol_u = 45.7

lam_l = 1900.0
lam_u = 2900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_ga = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_ra = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_ia = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.2
lbol_u = 45.7

lam_l = 2900.0
lam_u = 3900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_gb = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_rb = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_ib = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%

lbol_l = 45.2
lbol_u = 45.7

lam_l = 3900.0
lam_u = 4900.0

xx = np.array([lam_l,lam_u,lam_u,lam_l])
yy = np.array([lbol_l,lbol_l,lbol_u,lbol_u])

plt.fill(xx,yy,alpha=0.3,facecolor='lightsalmon', edgecolor='orangered', linewidth=3, zorder=5)

ind_gc = (lr_g_o > lam_l) & (lr_g_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_rc = (lr_r_o > lam_l) & (lr_r_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)
ind_ic = (lr_i_o > lam_l) & (lr_i_o < lam_u) & (lbol>lbol_l) & (lbol<lbol_u)

#%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)',fontsize=15)
plt.ylabel(r'$log_{10} ( L_{Bol} / L_\odot )$',fontsize=15)
plt.tight_layout()
plt.savefig('All_tau_r.png',dpi=300)
plt.close('all')

"""
### Plot 2
plt.figure()
plt.scatter(lr_r_o[ind_r],lbol[ind_r],c=tau_r_r[ind_r],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_g_o[ind_g],lbol[ind_g],c=tau_r_g[ind_g],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_i_o[ind_i],lbol[ind_i],c=tau_r_i[ind_i],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Selection_rr.png',dpi=300)
plt.close('all')

### Plot 2
plt.figure()
plt.scatter(lr_r_o[ind_r],lbol[ind_r],c=tau_o_r[ind_r],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_g_o[ind_g],lbol[ind_g],c=tau_o_g[ind_g],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.scatter(lr_i_o[ind_i],lbol[ind_i],c=tau_o_i[ind_i],vmin=1.5,vmax=4.5,cmap='rainbow')
plt.colorbar()

plt.xlabel(r'Rest Wavelength (${\rm \AA}$)')
plt.ylabel(r'$L_{Bol}$')
plt.savefig('Selection_ro.png',dpi=300)
plt.close('all')

zz = np.arange(1,5,0.01)
ff = np.log10( 1 + zz ) + 3

## Plot 2
plt.figure()
plt.errorbar(z[ind_r],tau_o_r[ind_r],yerr=tau_r_err[:,ind_r],fmt='ro')
plt.errorbar(z[ind_g],tau_o_g[ind_g],yerr=tau_g_err[:,ind_g],fmt='go')
plt.errorbar(z[ind_i],tau_o_i[ind_i],yerr=tau_i_err[:,ind_i],fmt='mo')

plt.plot(zz,ff,'b--',lw=4,alpha=0.4)

plt.xlabel(r'Redshift')
plt.ylabel(r'$\tau_{DRW}$')
plt.savefig('zvt.png',dpi=300)
plt.close('all')
"""
"""
Optimize - plan is to replace the asymmetric errors with a skewed gaussian


from scipy.optimize import minimize
import scipy.stats as sp 

skewg = np.zeros((len(DBID),3))
skewr = np.zeros((len(DBID),3))
skewi = np.zeros((len(DBID),3))

def skewfit(x,a16,amed,a84):

	a,loc,scale = x 

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=a,loc=loc,scale=scale)

	vmed = sp.skewnorm.median(a=a,loc=loc,scale=scale)

	return ( (v[0] - a16)**2 + (vmed - amed)**2 + (v[2] - a84)**2 )

## Loop over each quasar

for i in range(len(DBID)): 

	print('Working on quasar ',i)

	# Processing g
	v16 = tau_o_g[i] - tau_o_gl[i]
	vmed = tau_o_g[i] 
	v84 = tau_o_g[i] + tau_o_gu[i]

	print('v: ',v16,vmed,v84)

	x0 = np.array([0.0,vmed,1.0])

	res = minimize(skewfit,x0,method='nelder-mead',
                args=(v16,vmed,v84), options={'xatol': 1e-8, 'disp': True})

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=res.x[0],loc=res.x[1],scale=res.x[2])
	print('a: ',v)

	skewg[i,:] = res.x

	# Processing r
	v16 = tau_o_r[i] - tau_o_rl[i]
	vmed = tau_o_r[i] 
	v84 = tau_o_r[i] + tau_o_ru[i]

	print('v: ',v16,vmed,v84)

	x0 = np.array([0.0,vmed,1.0])

	res = minimize(skewfit,x0,method='nelder-mead',
                args=(v16,vmed,v84), options={'xatol': 1e-8, 'disp': True})

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=res.x[0],loc=res.x[1],scale=res.x[2])
	print('a: ',v)

	skewr[i,:] = res.x

	# Processing i
	v16 = tau_o_i[i] - tau_o_il[i]
	vmed = tau_o_i[i] 
	v84 = tau_o_i[i] + tau_o_iu[i]

	print('v: ',v16,vmed,v84)

	x0 = np.array([0.0,vmed,1.0])

	res = minimize(skewfit,x0,method='nelder-mead',
                args=(v16,vmed,v84), options={'xatol': 1e-8, 'disp': True})

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=res.x[0],loc=res.x[1],scale=res.x[2])
	print('a: ',v)

	skewi[i,:] = res.x

import pickle as pkl

file = open('SkewFit.pkl','wb')
pkl.dump(skewg,file)
pkl.dump(skewr,file)
pkl.dump(skewi,file)
file.close()
"""

import pickle as pkl
import scipy.stats as sp 

file = open('SkewFit.pkl','rb')

skewg = pkl.load(file)
skewr = pkl.load(file)
skewi = pkl.load(file)

file.close()

"""
Diagnostic plots

## g plot

fg1 = plt.figure()
fg2 = plt.figure()
fg3 = plt.figure()

for i in range(len(DBID)):

	v16 = tau_o_g[i] - tau_o_gl[i]
	vmed = tau_o_g[i] 
	v84 = tau_o_g[i] + tau_o_gu[i]

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=skewg[i,0],loc=skewg[i,1],scale=skewg[i,2])

	plt.figure(fg1)
	plt.plot( i , (v16-v[0])*100/v16,'o')

	plt.figure(fg2)
	plt.plot( i , (vmed-v[1])*100/vmed,'o')

	plt.figure(fg3)
	plt.plot( i , (v84-v[2])*100/v84,'o')


plt.figure(fg1)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (lower)')
plt.savefig('g_tau_lower.png',dpi=300)

plt.figure(fg2)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (median)')
plt.savefig('g_tau_median.png',dpi=300)

plt.figure(fg3)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (upper)')
plt.savefig('g_tau_upper.png',dpi=300)

## g plot

fg1 = plt.figure()
fg2 = plt.figure()
fg3 = plt.figure()

for i in range(len(DBID)):

	v16 = tau_o_r[i] - tau_o_rl[i]
	vmed = tau_o_r[i] 
	v84 = tau_o_r[i] + tau_o_ru[i]

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=skewr[i,0],loc=skewr[i,1],scale=skewr[i,2])

	plt.figure(fg1)
	plt.plot( i , (v16-v[0])*100/v16,'o')

	plt.figure(fg2)
	plt.plot( i , (vmed-v[1])*100/vmed,'o')

	plt.figure(fg3)
	plt.plot( i , (v84-v[2])*100/v84,'o')


plt.figure(fg1)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (lower)')
plt.savefig('r_tau_lower.png',dpi=300)

plt.figure(fg2)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (median)')
plt.savefig('r_tau_median.png',dpi=300)

plt.figure(fg3)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (upper)')
plt.savefig('r_tau_upper.png',dpi=300)

## i plot

fg1 = plt.figure()
fg2 = plt.figure()
fg3 = plt.figure()

for i in range(len(DBID)):

	v16 = tau_o_i[i] - tau_o_il[i]
	vmed = tau_o_i[i] 
	v84 = tau_o_i[i] + tau_o_iu[i]

	v = sp.skewnorm.ppf([0.16,0.50,0.84],a=skewi[i,0],loc=skewi[i,1],scale=skewi[i,2])

	plt.figure(fg1)
	plt.plot( i , (v16-v[0])*100/v16,'o')

	plt.figure(fg2)
	plt.plot( i , (vmed-v[1])*100/vmed,'o')

	plt.figure(fg3)
	plt.plot( i , (v84-v[2])*100/v84,'o')


plt.figure(fg1)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (lower)')
plt.savefig('i_tau_lower.png',dpi=300)

plt.figure(fg2)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (median)')
plt.savefig('i_tau_median.png',dpi=300)

plt.figure(fg3)
plt.xlabel('Nspec')
plt.ylabel('Percentage error (upper)')
plt.savefig('i_tau_upper.png',dpi=300)

"""

"""
Emcee 
"""

#import emcee
import scipy.stats as sp 

z_fit1 = np.concatenate((z[ind_g1],z[ind_r1],z[ind_i1]))
t_fit1 = np.concatenate((skewg[ind_g1,:],skewr[ind_r1,:],skewi[ind_i1,:]))

z_fit2 = np.concatenate((z[ind_g2],z[ind_r2],z[ind_i2]))
t_fit2 = np.concatenate((skewg[ind_g2,:],skewr[ind_r2,:],skewi[ind_i2,:]))

z_fit3 = np.concatenate((z[ind_g3],z[ind_r3],z[ind_i3]))
t_fit3 = np.concatenate((skewg[ind_g3,:],skewr[ind_r3,:],skewi[ind_i3,:]))

z_fit4 = np.concatenate((z[ind_g4],z[ind_r4],z[ind_i4]))
t_fit4 = np.concatenate((skewg[ind_g4,:],skewr[ind_r4,:],skewi[ind_i4,:]))

z_fit5 = np.concatenate((z[ind_g5],z[ind_r5],z[ind_i5]))
t_fit5 = np.concatenate((skewg[ind_g5,:],skewr[ind_r5,:],skewi[ind_i5,:]))

z_fit6 = np.concatenate((z[ind_g6],z[ind_r6],z[ind_i6]))
t_fit6 = np.concatenate((skewg[ind_g6,:],skewr[ind_r6,:],skewi[ind_i6,:]))

z_fit7 = np.concatenate((z[ind_g7],z[ind_r7],z[ind_i7]))
t_fit7 = np.concatenate((skewg[ind_g7,:],skewr[ind_r7,:],skewi[ind_i7,:]))

z_fit8 = np.concatenate((z[ind_g8],z[ind_r8],z[ind_i8]))
t_fit8 = np.concatenate((skewg[ind_g8,:],skewr[ind_r8,:],skewi[ind_i8,:]))

z_fit9 = np.concatenate((z[ind_g9],z[ind_r9],z[ind_i9]))
t_fit9 = np.concatenate((skewg[ind_g9,:],skewr[ind_r9,:],skewi[ind_i9,:]))

z_fita = np.concatenate((z[ind_ga],z[ind_ra],z[ind_ia]))
t_fita = np.concatenate((skewg[ind_ga,:],skewr[ind_ra,:],skewi[ind_ia,:]))

z_fitb = np.concatenate((z[ind_gb],z[ind_rb],z[ind_ib]))
t_fitb = np.concatenate((skewg[ind_gb,:],skewr[ind_rb,:],skewi[ind_ib,:]))

z_fitc = np.concatenate((z[ind_gc],z[ind_rc],z[ind_ic]))
t_fitc = np.concatenate((skewg[ind_gc,:],skewr[ind_rc,:],skewi[ind_ic,:]))


def logprob(x):

	n,c1,c2,c3,c4,c5,c6,c7,c8,c9,ca,cb,cc = x 

	model1 = c1 + n*np.log10( 1 + z_fit1 )
	model2 = c2 + n*np.log10( 1 + z_fit2 )
	model3 = c3 + n*np.log10( 1 + z_fit3 )
	model4 = c4 + n*np.log10( 1 + z_fit4 )
	model5 = c5 + n*np.log10( 1 + z_fit5 )
	model6 = c6 + n*np.log10( 1 + z_fit6 )
	model7 = c7 + n*np.log10( 1 + z_fit7 )
	model8 = c8 + n*np.log10( 1 + z_fit8 )
	model9 = c9 + n*np.log10( 1 + z_fit9 )
	modela = ca + n*np.log10( 1 + z_fita )
	modelb = cb + n*np.log10( 1 + z_fitb )
	modelc = cc + n*np.log10( 1 + z_fitc )

	t_pdf1 = np.sum( sp.skewnorm.logpdf(model1,a=t_fit1[:,0],loc=t_fit1[:,1],scale=t_fit1[:,2]) )
	t_pdf2 = np.sum( sp.skewnorm.logpdf(model2,a=t_fit2[:,0],loc=t_fit2[:,1],scale=t_fit2[:,2]) )
	t_pdf3 = np.sum( sp.skewnorm.logpdf(model3,a=t_fit3[:,0],loc=t_fit3[:,1],scale=t_fit3[:,2]) )
	t_pdf4 = np.sum( sp.skewnorm.logpdf(model4,a=t_fit4[:,0],loc=t_fit4[:,1],scale=t_fit4[:,2]) )
	t_pdf5 = np.sum( sp.skewnorm.logpdf(model5,a=t_fit5[:,0],loc=t_fit5[:,1],scale=t_fit5[:,2]) )
	t_pdf6 = np.sum( sp.skewnorm.logpdf(model6,a=t_fit6[:,0],loc=t_fit6[:,1],scale=t_fit6[:,2]) )
	t_pdf7 = np.sum( sp.skewnorm.logpdf(model7,a=t_fit7[:,0],loc=t_fit7[:,1],scale=t_fit7[:,2]) )
	t_pdf8 = np.sum( sp.skewnorm.logpdf(model8,a=t_fit8[:,0],loc=t_fit8[:,1],scale=t_fit8[:,2]) )
	t_pdf9 = np.sum( sp.skewnorm.logpdf(model9,a=t_fit9[:,0],loc=t_fit9[:,1],scale=t_fit9[:,2]) )

	t_pdfa = np.sum( sp.skewnorm.logpdf(modela,a=t_fita[:,0],loc=t_fita[:,1],scale=t_fita[:,2]) )
	t_pdfb = np.sum( sp.skewnorm.logpdf(modelb,a=t_fitb[:,0],loc=t_fitb[:,1],scale=t_fitb[:,2]) )
	t_pdfc = np.sum( sp.skewnorm.logpdf(modelc,a=t_fitc[:,0],loc=t_fitc[:,1],scale=t_fitc[:,2]) )


	return ( t_pdf1+t_pdf2+t_pdf3+t_pdf4+t_pdf5+t_pdf6+t_pdf7+t_pdf8+t_pdf9+t_pdfa+t_pdfb+t_pdfc ) 

"""
ndim, nwalkers = 13 , 100

p0 = np.random.rand(nwalkers,ndim) 

p0 = 5 * p0
"""
"""
sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob )
sampler.run_mcmc(p0, 25000, progress=True)
"""

def prior_transform(u):
    """Transforms the uniform random variable `u ~ Unif[0., 1.)`
    to the parameter of interest `x ~ Unif[-10., 10.)`."""

    x = 1. + 4*u 

    # This is n
    x[0] = 1 + 1e-6*(u[0]-0.5)

    return x

#import dynesty

#ndim = 13

## initialize our nested sampler

#dsampler = dynesty.DynamicNestedSampler(logprob, prior_transform, ndim)
#dsampler.run_nested()
#dres = dsampler.results
import dnest4
import numpy.random as rng

class Model(object):
    """
    DNest4 model.
    """

    def __init__(self):
        self.ndim = 13

    def from_prior(self):
        return rng.rand(self.ndim)

    def perturb(self, coords):
        i = np.random.randint(self.ndim)
        coords[i] += dnest4.randh()
        # Note: use the return value of wrap, unlike in C++
        coords[i] = dnest4.wrap(coords[i], 0.0, 1.0)
        return 0.0

    def log_likelihood(self, coords):
        try:
            logl = logprob(prior_transform(coords))
        except:
            logl = -np.Inf
        return logl

if __name__ == "__main__":

    import dnest4
    model = Model()
    sampler = dnest4.DNest4Sampler(model,
                                   backend=dnest4.backends.CSVBackend(".",
                                                                      sep=" "))
    gen = sampler.sample(100, num_steps=100000, new_level_interval=5000,
                         num_per_step=50, thread_steps=100,
                         num_particles=5, lam=1, beta=100, seed=1234)
    for i, sample in enumerate(gen):
        if (i + 1) % 100 == 0:
            print(f"Sample {i+1}.", flush=True)
            import dnest4.classic as dn4
            dn4.postprocess(plot=False, temperature=1)

            # Convert uniform coordinates to parameter space
            posterior_sample = np.atleast_2d(np.loadtxt("posterior_sample.txt"))
            for i in range(posterior_sample.shape[0]):
                posterior_sample[i, :] = prior_transform(posterior_sample[i, :])
            np.savetxt("posterior_sample.txt", posterior_sample)

