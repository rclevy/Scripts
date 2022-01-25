def mk_great_footprint(center,array_rot=0.,which_config='LFA',color='green'):
	r'''
	Make a DS9 region file of the upGREAT LFA or HFA footprint given some center position.

	Parameters
	----------
	center : SkyCoord
		central RA and Dec passed in as an astropy SkyCoord
	array_rot : float
		rotation angle of the array in degrees, default is 0
	which_config : string
		only the LFA and HFA configurations are supported, default is 'LFA'
	color : string
		color of regions, default is 'green'


	Returns
	-------
	None


	Notes
	-----
	Required modules: astropy, numpy, regions
	Developed using astropy 4.0.2, numpy 1.19.2, regions 0.5 in iPython 7.19.0 and Python 3.8.5
	Author: R. C. Levy (rlevy.astro@gmail.com)
	Last Updated: 2022-01-24
	Change log:
		2022-01-21 : file created
		2022-01-24 : generates the footprints in the correct order as numbered in the handbook, https://www-sofia.atlassian.net/wiki/spaces/OHFC1/pages/1147523/6.+GREAT#Figure6-3


	Examples
	--------
	>> ipython
	>> from astropy.coordinates import SkyCoord
	>> from mkGREATfootprint import mk_great_footprint
	>> center = SkyCoord('13h05m27.2776s','-49d28m5.556s') #NGC4945
	>> mk_great_footprint(center,array_rot=0.,which_config='LFA',color='cyan')
	>> #output region file has the name SOFIA_upGREAT_[which_config]_Cen[center]_Rot[array_rot].reg
	'''

	import numpy as np
	from astropy.coordinates import SkyCoord,Angle,concatenate
	from regions import CircleSkyRegion,Regions


	#get radius of each pixel, from https://www-sofia.atlassian.net/wiki/spaces/OHFC1/pages/1147523/6.+GREAT#Table6-1
	#get linear separation of pixel centers, from #https://www-sofia.atlassian.net/wiki/spaces/OHFC1/pages/1147523/6.+GREAT#Figure6-3
	if which_config == 'LFA':
		radius = Angle(14.1/2,'arcsec')
		separation = Angle(31.7,'arcsec') 
	elif which_config == 'HFA':
		radius = Angle(6.3/2,'arcsec')
		separation = Angle(13.8,'arcsec')
	else:
		print("Error: Only LFA and HFA configurations are supported.\nPlease try again with either mk_great_footprint(center,which_config='LFA') or mk_great_footprint(center,which_config='HFA').")
		return

	#get coordinates of the surrounding pixels
	#account for any rotation of the array too
	angles = np.radians(np.arange(0,360,60)-150+array_rot)
	other_pixels = center.directional_offset_by(angles,separation)
	pixel_coords = concatenate((center,other_pixels))

	regions = ['']*len(pixel_coords)
	for i in range(len(pixel_coords)):
		regions[i] = CircleSkyRegion(pixel_coords[i],radius)
		regions[i].visual['color']=color
	reg = Regions(regions)
	#write regions to a file
	fname = 'SOFIA_upGREAT_'+which_config+'_regions_Cen'+center.to_string('hmsdms').replace(' ','')+'_Rot'+str(array_rot)+'.reg'
	reg.write(fname,format='ds9',overwrite=True)

	return