import os
import gal_inserter as gi 
import finder
import detection_rates as dr

def extract_params(paramfile):
    '''
    This function extracts the parameters the user has set for
    the program from a text file.

    ------

    Parameters:

    paramfile (string): The path to the text file containing the 
    user's preferred parameters.

    ------

    Outputs:

    dir (string): the path to the directory all of the program's 
    output will be sent to.

    image (string): the path to the image to be used as a 
    base for inserting fake galaxies.

    startsize (integer): the smallest size of galaxy to be
    generated, in arcseconds.

    endsize (integer): the largest size of galaxy to be
    generated, in arcseconds.

    sizestep (integer): the step size for the program to
    use as it increments through galaxy sizes, in 
    arcseconds.

    startmag (integer): the lowest magnitude of galaxy to
    be generated.

    endmag (integer): the highest magnitude of galaxy to be
    generated.

    magstep (integer): the step size for the program to
    use as it increments through galaxy magnitudes.

    sexparams (string): the path to the Source Extractor
    parameter file to be used.
    '''

    params = [0, 0, 0, 0, 0, 0, 0, 0, 0] #setting up the list of extracted parameters

    #reading the text file and extracting the parameters
    txt = open(paramfile, 'r')
    lines = txt.readlines()
    txt.close()
    for i in range(2, len(lines)): #we skip the first two lines because they're just a header
        line = lines[i].strip('\n').split(' ')
        params[i-2] = line[-1] 
        if i > 3 and i < 10: #converting the necessary parameters to integers
            params[i-2] = int(params[i-2])

    #assigning the parameters to the corresponding variables
    dir = params[0]
    image = params[1]
    startsize = params[2]
    endsize = params[3]
    sizestep = params[4]
    startmag = params[5]
    endmag = params[6]
    magstep = params[7]
    sexparams = params[8]
    
    return dir, image, startsize, endsize, sizestep, startmag, endmag, magstep, sexparams


def run_sextr(dir, sexparams):

    '''
    This function runs the images with inserted galaxies through Source
    Extractor.

    -----

    Parameters:

    dir (string): the path to the directory that image to be used is
    located in.

    -----

    Outputs:

    None
    '''

    for img in os.listdir(dir): #runs through every file in the directory
        if img.endswith('.fits'): #and checks if it's a .fits file
            print(img) #tells the user which images is being run 
            imgname = os.path.splitext(os.path.basename(img))[0] #removes the .fits from the filename
            catname = dir + '/output_cold_' + imgname + '.cat' #the name for the SEx catalog output file
            checkname = dir + '/check_' + imgname + '.fits' #the name for the SEx check image

            #Assembling the command to be passed to Source Extractor
            basecmd = 'sex ' + dir+'/'+img + ' -c '+sexparams
            catcmd = '-catalog_name ' + catname 
            checkcmd = '-checkimage_name ' + checkname
            sexcmd = basecmd + ' ' + catcmd + ' ' + checkcmd
            
            print(sexcmd)
            os.system(sexcmd) #passing the command to the operating system 


def run_pipeline(paramfile):
    '''
    This function runs the whole pipeline.

    -----

    Parameters:

    paramfile (string): The path to the text file containing the 
    user's preferred parameters.

    -----

    Outputs:

    None
    '''

    #extracting the parameters
    dir, image, startsize, endsize, sizestep, startmag, endmag, magstep, sexparams = extract_params(paramfile)
    print('Parameters read')

    if ~os.path.exists(dir): #makes the directory if it doesn't exist
        os.mkdir(dir)
    os.chdir(dir) #switches to it

    #inserting galaxies
    print('Starting galaxy insertion')
    gi.run_range(image, startsize, endsize, sizestep, startmag, endmag, magstep)

    #iterating through the subdirectories to find images 
    for container in os.listdir(dir):
        for subcont in os.listdir(dir + '/' + container):
            print("Running SEx for:") #running Source Extractor on the image
            run_sextr(dir + "/" + container + "/" + subcont, sexparams)
            for img in os.listdir(dir + '/' + container + '/' + subcont):
                if img.endswith('.fits') and not img.startswith('check'):
                    print('Finding matches for', img)
                    finder.find_matches(dir + '/' + container + '/' + subcont, img) #find matches for the image
    
    #get and record detection rates for the images
    print("Getting detection rates")
    rates_list, sizemag_pairs = dr.get_rates(dir)
    dr.record_rates(rates_list, sizemag_pairs, dir, startsize, endsize, sizestep, startmag, endmag, magstep)
    print('Finished')