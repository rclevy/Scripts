def mk_great_footprint(center,array_rot=0.,which_config='LFA',color='green'):
	#make a DS9 region file of the SOFIA upGREAT footprint
	#given some center (passed in as a SkyCoord) and array rotation angle (degrees)
	#specify which configuration (LFA or HFA)
	#can optically specify region color (default green)
	#
	#requires numpy, astropy, and regions
	#
	#written by Rebecca Levy (rlevy.astro@gmail.com)
	#last updated 2022-01-21


	from astropy.coordinates import SkyCoord,Angle,concatenate,SkyOffsetFrame
	from regions import CircleSkyRegion,Regions
	import numpy as np
	import astropy.units as u

	def rotate(origin, point, angle):
		"""
		Rotate a point counterclockwise by a given angle around a given origin.

		The angle should be given in degrees.

		Based on https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
		"""
		ox = origin.ra
		oy = origin.dec
		px = point.ra
		py = point.dec

		qx = ox+np.cos(np.radians(angle))*(px-ox)-np.sin(np.radians(angle))*(py-oy)
		qy = oy+np.sin(np.radians(angle))*(px-ox)+np.cos(np.radians(angle))*(py-oy)
		
		point_rot = SkyCoord(qx,qy)
		return point_rot

	#get radius of each pixel, from https://www-sofia.atlassian.net/wiki/spaces/OHFC1/pages/1147523/6.+GREAT#Table6-1
	#get linear separation of pixel centers, from #https://www-sofia.atlassian.net/wiki/spaces/OHFC1/pages/1147523/6.+GREAT#Figure6-3
	if which_config == 'LFA':
		#radius = Angle(14.1/2,'arcsec') 
		radius = Angle(15.1/2,'arcsec')
		separation = Angle(31.7,'arcsec') 
	elif which_config == 'HFA':
		radius = Angle(6.3/2,'arcsec')
		separation = Angle(13.8,'arcsec')
	else:
		print('Error: Only LFA and HFA configurations are supported.\nPlease try again with either mk_great_footprint(center,array_rot,"LFA") or mk_great_footprint(center,array_rot,"HFA").')
		return

	#get coordinates of the surrounding pixels
	#account for any rotation of the array too
	angles = np.radians(np.arange(0,360,60)+30+array_rot)
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