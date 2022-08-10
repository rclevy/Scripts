import argparse

parser=argparse.ArgumentParser(
	prog = 'calcKroupaIMF',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description='''description:
  Calculate the number of stars in a cluster of a given mass using a Kroupa (2002) IMF and the specified min and max stellar masses. Calculate the number of O stars in the cluster. Calculate the average separation between stars in this cluster given a cluster radius and assuming spherical geometry.
	
examples:
  >>python calcKroupaIMF 100000. 1.
  >>python calcKroupaIMF 100000. 1. --Mmax_Msun 200.
  >>python calcKroupaIMF 100000. 1. --Mmin_Msun 0.02 --Mmax_Msun 200.
  >>python calaKroupaIMF 100000. 1. --Mmin_Ostar_Msun 10.

required packages:
  argparse
  numpy''',
	epilog='''Author: R. C. Levy (rlevy.astro@gmail.com) - Last updated: 2022-08-10''')
parser.add_argument('Mcluster_Msun',type=float,help='stellar mass of the cluster in solar masses')
parser.add_argument('Rcluster_pc',type=float,help='radius of the cluster in parsecs')
parser.add_argument('--Mmin_Msun',type=float,default=0.01,help='minimum stellar mass in solar masses (default: 0.01 Msun)')
parser.add_argument('--Mmax_Msun',type=float,default=100.,help='maximum stellar mass in solar masses (default: 100 Msun)')
parser.add_argument('--Mmin_Ostar_Msun',type=float,default=8.,help='minimum mass of an O star in solar masses (default: 8 Msun)')
args=parser.parse_args()

import numpy as np

def filter_Kroupa_minmaxMass(Mmin,Mmax):
	#fiducial mass breaks for the Kroupa IMF
	K_bounds_min = np.array([0.01,0.08,0.5,1.0])
	K_bounds_max = np.array([0.08,0.5,1.0,100.])

	#check that the bounds are valid
	if Mmin <= 0:
		print('ERROR: Mmin must be > 0')
		raise ValueError

	if Mmax <= 0:
		print('ERROR: Mmax must be > 0')
		raise ValueError

	if Mmin > Mmax:
		print('ERROR: Mmin_Msun > Mmax_Msun')
		raise ValueError

	#set the min and max masses to the min and max bounds
	K_bounds_min[0] = Mmin
	K_bounds_max[-1] = Mmax

	if Mmin > K_bounds_max[0]:
		#don't sum over the lowest mass bin, set min mass of second bin to Mmin
		K_bounds_min[1] = Mmin
		K_bounds_min[0] = np.nan
		K_bounds_max[0] = np.nan
	if Mmin > K_bounds_max[1]:
		#don't sum over the lowest two mass bins, set min mass of third bin to Mmin
		K_bounds_min[2] = Mmin
		K_bounds_min[1] = np.nan
		K_bounds_max[1] = np.nan
	if Mmin > K_bounds_max[2]:
		#don't sum over the lowest three mass bins, set min mass of highest bin to Mmin
		K_bounds_min[3] = Mmin
		K_bounds_min[2] = np.nan
		K_bounds_max[2] = np.nan
	if Mmax < K_bounds_min[3]:
		K_bounds_max[2] = Mmax
		K_bounds_max[3] = np.nan
		K_bounds_min[3] = np.nan
	if Mmax < K_bounds_min[2]:
		K_bounds_max[1] = Mmax
		K_bounds_max[2] = np.nan
		K_bounds_min[2] = np.nan
	if Mmax < K_bounds_min[1]:
		K_bounds_max[0] = Mmax
		K_bounds_max[1] = np.nan
		K_bounds_min[1] = np.nan

	K_bounds = (K_bounds_min,K_bounds_max)
	return K_bounds

def Kroupa_A_from_Mtot(Mtot,K_bounds):
	M_piecewise = np.array([(K_bounds[0][0]**-0.7-K_bounds[1][0]**-0.7)/0.7,
							(K_bounds[0][1]**-0.3-K_bounds[1][1]**-0.3)/0.3,
							(K_bounds[1][2]**0.7-K_bounds[0][2]**0.7)/0.7,
							(K_bounds[1][3]**1.7-K_bounds[0][3]**1.7)/1.7])
	A = Mtot/np.nansum(M_piecewise)
	return A

def Kroupa_Ntot(A,K_bounds):
	N_piecewise = np.array([(K_bounds[0][0]**-1.7-K_bounds[1][0]**-1.7)/1.7,
							(K_bounds[0][1]**-1.3-K_bounds[1][1]**-1.3)/1.3,
							(K_bounds[0][2]**-0.3-K_bounds[1][2]**-0.3)/0.3,
							(K_bounds[1][3]**0.7-K_bounds[0][3]**0.7)/0.7])
	Ntot = A*np.nansum(N_piecewise)
	return Ntot

def ave_stellar_separation(Ntot,r):
	V = 4/3*np.pi*r**3
	n = Ntot/V
	sep_ave_pc = (3/(4*np.pi*n))**(1/3)
	sep_ave_au = sep_ave_pc*206265
	sep_ave_ly = sep_ave_pc*3.26156
	sep_ave_lm = sep_ave_ly*12.175
	return sep_ave_pc,sep_ave_au,sep_ave_ly,sep_ave_lm

#parse input parameters
Mtot = args.Mcluster_Msun
r = args.Rcluster_pc
Mmin = args.Mmin_Msun
Mmax = args.Mmax_Msun
Mmin_O = args.Mmin_Ostar_Msun

print('For a Kroupa IMF, a cluster with input parameters:')
print('\tM_cluster = %.1e Msun' %Mtot)
print('\tr_cluster = %.1f pc' %r)
print('\tM_stars = %g - %g Msun' %(Mmin,Mmax))
print('\tM_min_Ostar = %.1f Msun' %Mmin_O)

#calc total number of stars
K_bounds = filter_Kroupa_minmaxMass(Mmin,Mmax)
A = Kroupa_A_from_Mtot(Mtot,K_bounds)
Ntot = Kroupa_Ntot(A,K_bounds)
sep_ave_pc,sep_ave_au,sep_ave_ly,sep_ave_lm = ave_stellar_separation(Ntot,r)

#calc number of O stars
K_bounds_O = filter_Kroupa_minmaxMass(Mmin_O,Mmax)
N_O = Kroupa_Ntot(A,K_bounds_O)

print('yields:')
print('\tN_stars_tot = %.1e' %Ntot)
print('\tN_Ostars = %.1f' %N_O)
print('\tD_ave_interstellar = %.2e pc = %.1f AU = %.2f ly = %.2f lm' %(sep_ave_pc,sep_ave_au,sep_ave_ly,sep_ave_lm))





