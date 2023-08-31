import numpy as np
import pandas as pd
from astropy.io import fits, ascii
import os


def read_in_img_info(imgnum):
	'''
	This function gets the x and y dimensions of the input image
	and returns the smaller value of the two for later use.

	-----

	Parameters:

	imgnum (string): the name of the image used, without the 
	file extension.

	-----

	Outputs:

	imgsize (integer): the smaller value of the x and y 
	dimensions of the image.
	'''

	image = fits.open(imgnum+'.fits')
	img_x_size = image[0].data.shape[1] #getting the x and y dimensions of the image
	img_y_size = image[0].data.shape[0]
	image.close()
	imgsize = min(img_x_size, img_y_size) #choose the smaller value

	return imgsize


def read_in_sex_output(imgnum):
	'''
	This function reads in the necessary information from the
	Source Extractor output catalog.

	-----

	Parameters:

	imgnum (string): the name of the image used, without the 
	file extension.

	-----

	Outputs:

	sexout (dictionary): the dictionary containing all the 
	important Source Extractor output.
	'''

	#This is the section that takes the catalog and extracts the needed values
	#data = np.loadtxt("output_cold_image" + imgnum + "_" + convfil + ".cat") #uncomment if convolution filter names become an issue
	data = ascii.read("output_cold_" + imgnum + ".cat") #get the source extractor output catalog
	real_x = np.array(data['X_IMAGE']) #the detected x coordinates of the sources
	real_y = np.array(data['Y_IMAGE']) #y coordinates
	fwhm = data['FWHM_IMAGE'] #full width half max-es
	hlr_list = data['FLUX_RADIUS'] #half light radii
	mag_auto = data['MAG_AUTO'] #mag auto-s
	scobj_num = data['NUMBER'] #the id numbers assigned to the objects by source extractor
	
	sexout = {'scobj_num': scobj_num, "real_x": real_x, 'real_y': real_y, 'fwhm': fwhm, 
	   'hlr_list': hlr_list, 'mag_auto': mag_auto}

	return sexout


def read_in_gal_info(imgnum):
	'''
	This function reads in the necessary information about the inserted 
	ultra diffuse galaxies from the .csv files generated when they were
	inserted. 

	-----

	Parameters:

	imgnum (string): the name of the image used, without the 
	file extension.
	
	-----

	Outputs:

	values_add (dictionary): the dictionary containing all the 
	important information about the UDGs.
	'''

	#reads in the neessary csv file containing the actual info about the added galaxies
	obj_d = pd.read_csv(imgnum + ".csv")

	#xadd = where we're going to put the actual x coordinates of the galaxies
	#yadd = y coordinates
	#fwhmadd = full width half max-es
	#sizeadd = sizes
	values_add = {'xadd': [], 'yadd': [], 'fwhmadd': [], 'sizeadd': []}
	keys = list(values_add.keys()) 
	dfkeys = ['x value', 'y value', 'FWHM (px)', 'size (px)'] #the keys for the 

	for i in range(4):
		for a in obj_d[dfkeys[i]]: #iterating through every important column in the obj_d dataframe
			values_add[keys[i]].append(a) #adding the values in it to the correct list in values_add

	return values_add



def find_matches(imgdir, imgnum):

	'''
	This function compares the information generated about the added
	ultra diffuse galaxies in an image at the insertion step to the 
	Source Extractor detections in that same image and identifies 
	the most likely matches.

	-----

	Parameters:

	imgdir (string): the path to the directory where all of the necessary
	files about the image in question are.

	imgnum (string): the filename of the image.
	'''

	imgnum = imgnum.strip('.fits') #gets rid of the .fits extension, just in case
	os.chdir(imgdir) #goes to wherever your target files are

	#reading in the necessary info from:
	imgsize = read_in_img_info(imgnum) #the image
	sexout = read_in_sex_output(imgnum) #Source Extractor
	values_add = read_in_gal_info(imgnum) #the .csv files generated at the galaxy insertion step

	#iterating through all the known added galaxies
	for i in range(len(values_add['xadd'])): #doesn't have to be this list, just needed one with the right length

		match_info = {'SC#':[], 'X':[], 'Y':[], 'Distance':[], 'FWHM':[], 'MAG_AUTO':[], 'HLR':[], 
		'FHWM out of range?':[], 'Matches':[]} #setting up the output for the found matches
		n = imgsize/4 + np.round(values_add['sizeadd'][i]) #threshold for how far the extracted position can be from the actual position
		ul = imgsize*.25 #setting the upper limit 
		ll = 0.5 #and lower limit percentage for a fwhm match
		matches = 0 #keeps track of matches
		
		dist = np.sqrt((sexout['real_x']-values_add['xadd'][i])**2 
		 + (sexout['real_y']-values_add['yadd'][i])**2) #array of the distance of each detection from the actual source
		locs = np.where(dist<n)[0] #gets the indices of the ones within the set threshold

		if locs.shape[0] > 0: #if there are detections close enough, finds the one that's the best match

			lrgstfwhm = 0
			matchloc = 0

			for loc in locs: #iterates through the detections

				#finds the one with the largest fhwm and sets that as the best match so far
				currentfwhm = sexout['fwhm'][loc]
				if currentfwhm > lrgstfwhm:
					lrgstfwhm = currentfwhm
					matchloc = loc

			#recording all the important information about the best match
			match_info['SC#'].append(sexout['scobj_num'][matchloc])
			match_info['X'].append(sexout['real_x'][matchloc])
			match_info['Y'].append(sexout['real_y'][matchloc])
			match_info['Distance'].append(dist[matchloc])
			match_info['FWHM'].append(sexout['fwhm'][matchloc])
			match_info['MAG_AUTO'].append(sexout['mag_auto'][matchloc])
			match_info['HLR'].append(sexout['hlr_list'][matchloc])

			#testing if the fwhm is within the specified range; if so, counts it as a 'real' match
			if sexout['fwhm'][matchloc] > (ll*values_add['fwhmadd'][i]) and sexout['fwhm'][matchloc] <= ul:
				
				matches += 1
				match_info['FHWM out of range?'].append('No')
			
			else:

				match_info['FHWM out of range?'].append('Yes')

		match_info['Matches'].append(matches) #records the amount of matches found

		#this is just to make sure all the lists in the matches dictionary are the same length
		#otherwise pandas can't turn it into a dataframe
		longest = 0
		#identifying the longest list in the dictionary
		for key in list(match_info.keys()):
			if len(match_info[key]) > longest:
				longest = len(match_info[key])
		#filling out the other lists with a bunch of NaN values to make them the same length
		for key in list(match_info.keys()):
			n = longest - len(match_info[key])
			match_info[key].extend([np.nan] * n)
		
		#saving the match info as a .csv file
		match_df = pd.DataFrame(data=match_info)
		match_df.to_csv('galaxy'+str(i+1)+'_matches.csv')
		