#check te PSF of the combined continuum data

this_script = __file__

import argparse

parser=argparse.ArgumentParser(
    description='''Function to make some diagnostic plots of the input PSF fits file.
    				Plots the PSF image (cropped if specified) and plots slices through the PSF.
    				Required modules: argparse, astropy, matplotlib, numpy, scipy''',
    epilog='''Author: R. C. Levy - Last updated: 2021-02-18''')
parser.add_argument('path_to_file',type=str, help='Path to PSF file')
parser.add_argument('--crop_offset', type=float, help='Crop the output PSF cuts to some offset (arcsec)')

args=parser.parse_args()

import numpy as np
from astropy.io import fits
from astropy.visualization import (ManualInterval, ImageNormalize, LogStretch)
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

# open file
psf_file = args.path_to_file
save_path = psf_file.split('.fits')[0]
psf = fits.open(psf_file)
hdr = psf[0].header
psf = psf[0].data

# get slice along x
xcen = int(hdr['CRPIX1'])
x_off = np.arange(1,psf.shape[0]+1,1)-xcen
x_off_arcsec = x_off*hdr['CDELT2']*3600
psf_x = psf[xcen,:]

# get slice along y
ycen = int(hdr['CRPIX2'])
y_off = np.arange(1,psf.shape[1]+1,1)-ycen
y_off_arcsec = y_off*hdr['CDELT2']*3600
psf_y = psf[:,ycen]

# rotate PSF to get diagonal slices
# This assumes the PSF is properly centered
psf_rot = rotate(psf,45,reshape=False)
psf_xy = psf_rot[:,ycen]
psf_yx = psf_rot[xcen,:]

# get Gaussian beam shape
def gauss(x,beam):
	sig = beam/2.355
	g = np.exp(-x**2/(2*sig**2))
	return g/np.nanmax(g)
beam_pix = hdr['BMAJ']/hdr['CDELT2']

# crop the plotted PSF region in offset from center
if args.crop_offset:
	crop_arcsec = args.crop_offset
	crop_pix = crop_arcsec/3600/hdr['CDELT2']
else:
	crop_arcsec= np.max([x_off_arcsec,y_off_arcsec])
	crop_pix = np.max([x_off,y_off])
crop_suff = '_'+str(np.round(crop_arcsec,1))+'arcsec'


# plot the PSF image
plt.figure(1)
plt.clf()
norm=ImageNormalize(psf,stretch=LogStretch(),interval=ManualInterval(vmin=0,vmax=1))
plt.imshow(psf,origin='lower',norm=norm,cmap='gray_r')
plt.colorbar()
plt.xlim(left=xcen-crop_pix,right=xcen+crop_pix)
plt.ylim(bottom=ycen-crop_pix,top=ycen+crop_pix)
plt.savefig(save_path+crop_suff+'_image.pdf',bbox_inches='tight',metadata={'Creator':this_script})
plt.plot([xcen-crop_pix,xcen+crop_pix],[ycen,ycen],linestyle='-',lw=1.0)
plt.plot([xcen,xcen],[ycen-crop_pix,ycen+crop_pix],linestyle='--',lw=1.0)
plt.plot([xcen-crop_pix,xcen+crop_pix],[ycen-crop_pix,ycen+crop_pix],linestyle='-.',lw=1.0)
plt.plot([xcen-crop_pix,xcen+crop_pix],[ycen+crop_pix,ycen-crop_pix],linestyle=':',lw=1.0)
plt.savefig(save_path+crop_suff+'_image_slices.pdf',bbox_inches='tight',metadata={'Creator':this_script})
plt.close()

#make a second x-axis on the top in pix
def arcsec2pix(x):
	return x/3600/hdr['CDELT2']
def pix2arcsec(x):
	return x*hdr['CDELT2']*3600


# plot the slices
plt.figure(2)
plt.clf()
ax1=plt.gca()
plt.axvline(0,color='gray',lw=0.5)
plt.plot(x_off_arcsec,psf_x,label='PSF along x',lw=1.0)
plt.plot(y_off_arcsec,psf_y,'--',label='PSF along y',lw=1.0)
plt.plot(x_off_arcsec,psf_yx,'-.',label='PSF along y=x',lw=1.0)
plt.plot(y_off_arcsec,psf_xy,':',label='PSF along y=-x',lw=1.0)
plt.plot(x_off_arcsec,gauss(x_off_arcsec,hdr['BMAJ']*3600),'k-',label='Pure Gaussian\nFWHM=BMAJ',lw=1.25)
plt.legend()
plt.xlim(left=-crop_arcsec,right=crop_arcsec)
plt.xlabel('Offset from Center (arcsec)')
plt.minorticks_on()
ax2 = ax1.secondary_xaxis('top',functions=(arcsec2pix,pix2arcsec))
ax2.set_xlabel('Offset from Center (pixels)')
ax2.minorticks_on()
plt.savefig(save_path+crop_suff+'_slice.pdf',bbox_inches='tight',metadata={'Creator':this_script})
ax1.set_yscale('log')
ax1.set_ylim(bottom=1E-2,top=1)
plt.savefig(save_path+crop_suff+'_slice_log.pdf',bbox_inches='tight',metadata={'Creator':this_script})
plt.close()



