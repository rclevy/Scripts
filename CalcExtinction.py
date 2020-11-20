#calculate the extinction between two bands
#from Rieke & Lebofsky 1985 Table 3

import argparse

parser=argparse.ArgumentParser(
	prog = 'CalcExtinction',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description='''Calculate the MIR extinction and flux ratios between two wavelengths or bands based on Rieke & Lebofsky (1985). Values are linearly interpolated between those in Table 3. Inputs will be rounded to 0.1 micron. Band names must match Table 3 of Rieke & Lebofsky (1985).
		Examples:
			>>python CalcExtinction.py J 30 8.0
			>>python CalcExtinction.py J 30 M
			>>python CalcExtinction.py 8.0 2.0 J
			>>python CalcExtinction.py 3.3 2.0 8.0
		Required packages: numpy, scipy''',
	epilog='''Author: R. C. Levy (rlevy.astro@gmail.com) - Last updated: 2020-11-20''')
parser.add_argument('ReferenceWavelength',type=str,help='Wavelength or band to use as the reference with known extinction (microns or band letter)')
parser.add_argument('ReferenceExtinction',type=float,help='Extinction (mag) in the reference band (A_ReferenceBand).')
parser.add_argument('TargetWavelength',type=str,help='Wavelength at which to calculate the extinction and flux ratio (microns or band letter).')
args=parser.parse_args()

import numpy as np
from scipy.interpolate import interp1d

#Rieke & Lebofsky 1985 Table 3
band = ['U','B','V','R','I','J','H','K','L','M','N','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0','12.5','13.0'] #microns or band name
wl = np.array([0.365, 0.445, 0.551, 0.658, 0.806, 1.220, 1.630, 2.190, 3.450, 4.750, 10.50, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13.]) #microns
ElamV_EBV = np.array([1.64, 1.0, 0.0, -0.78, -1.60, -2.22, -2.55, -2.744, -2.91, -3.02, -2.93, -3.03, -2.96, -2.87, -2.83, -2.86, -2.87, -2.91, -2.95, -2.98, -3.00, -3.01])
Alam_Av = np.array([1.531, 1.324, 1.000, 0.748, 0.482, 0.282, 0.175, 0.112, 0.058, 0.023, 0.052, 0.020, 0.043, 0.074, 0.087, 0.083, 0.074, 0.060, 0.047, 0.037, 0.030, 0.027])


#sort the wavelengths and interpolate
idx_sort = np.argsort(wl)
wl_sort = wl[idx_sort]
Alam_Av_sort = Alam_Av[idx_sort]
interp = interp1d(wl_sort,Alam_Av_sort)

#set up grid for interpolated wavelengths
l_step = 0.1
l_start = np.round(wl_sort[0],1)
l_end = np.round(wl_sort[-1],1)
ll = np.round(np.arange(l_start,l_end+l_step,l_step),1) #microns
Alam_Av_interp = interp(ll)


#input reference extinction and band
ref_wl = args.ReferenceWavelength
ref_A = args.ReferenceExtinction #mag
try: #input is a number
	ref_wl = np.round(float(args.ReferenceWavelength),1)
	ref_wl_str = str(ref_wl)
	ref_idx = np.where(ll==ref_wl)[0][0]
	ref_AlamAv = Alam_Av_interp[ref_idx]
except (TypeError,ValueError): #input is a string
	ref_wl_str = ref_wl
	ref_idx = band.index(ref_wl)
	ref_AlamAv = Alam_Av[ref_idx]

Av = ref_A/ref_AlamAv #mag

#band in which to compute the extinction
target_wl = args.TargetWavelength
try: #input is a number
	target_wl = np.round(float(args.TargetWavelength),1)
	target_wl_str = str(target_wl)
	target_idx = np.where(ll==target_wl)[0][0]
	target_AlamAv = Alam_Av_interp[target_idx]
except (TypeError,ValueError): #input is a string
	target_wl_str = target_wl
	target_idx = band.index(target_wl)
	target_AlamAv = Alam_Av[target_idx]

target_A = target_AlamAv*Av #mag


#convert this extinction in magnitudes to difference in flux
flux_ratio = 10**(-target_A/2.5)


#print the results
print('Reference extinction: A_'+ref_wl_str+' = %.1f mag' %(ref_A))
print('Target extinction: A_'+target_wl_str+' = %.1f mag' %(target_A))
print('Flux ratio: F_extincted/F_intrinsic = %.3f at %s um' %(flux_ratio,target_wl_str))


