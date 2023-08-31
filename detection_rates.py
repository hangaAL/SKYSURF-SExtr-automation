import pandas as pd
import numpy as np
import os

def get_total_rate(container):
    '''
    This function gets the the total detection rate for a
    specific size/magnitude pair.

    -----

    Parameters:

    container (string): the path to the directory where all
    of the folders containing the images with inserted galaxies
    and the corresponding output files are stored.

    -----

    Outputs:

    detect_rate (integer): the percentage of the inserted 
    galaxies that were detected for that size/magnitude pair.
    '''

    gals_detected = [] #stores the amounts of matches per galaxy
    for imgdir in os.listdir(container): #iterating through every folder
        for file in os.listdir(container + '/' + imgdir): #and file
            if file.endswith('_matches.csv'): #finds the matches file
                #gets the amount of matches
                data = pd.read_csv(file) 
                detects = data['Matches'][0]
                gals_detected.append(detects)
    
    detect_rate = (np.sum(gals_detected) / 15)*100 #calculates the rate in percentage form
    return detect_rate


def get_rates(dir):
    '''
    This function gets every size-magnitude pair and their
    corresponding detection.

    -----

    Parameters:

    dir (string): the path to the directory where all of 
    size-magnitude pair subfolders are located.

    -----

    Outputs:

    rates_list (list): the list of all the detection rates
    for every size-magnitude pair.

    sizemag_pairs (list): the list of every size-magnitude
    pair.
    '''

    rates_list = []
    sizemag_pairs = []

    for container in os.listdir(dir): #goes through every size-magnitude pair subfolder
        rates_list.append(get_total_rate(dir + '/' + container)) #gets the detection rates
        sizemag_pairs.append(container[-5:]) #gets the size-magnitude pair from the folder name

    return rates_list, sizemag_pairs

def record_rates(rates_list, sizemag_pairs, dir, startsize, endsize, sizestep, startmag, endmag, magstep):
    '''
    This function makes a table of of all the detection rates for
    each size-magnitude pair and outputs it as a .csv file.

    -----

    Parameters: 

    rates_list (list): the list of all the detection rates
    for every size-magnitude pair.

    sizemag_pairs (list): the list of every size-magnitude
    pair.

    dir (string): the path to the directory where all of 
    size-magnitude pair subfolders are located.

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

    -----

    Outputs:

    None
    '''

    sizes = [size for size in range(startsize, endsize+1, sizestep)] #sets up a list of all the sizes to have in a separate column
    dic = {'Sizes;Magnitudes':sizes} #setting up the dictionary to use for making the table
    for mag in range(startmag, endmag+1, magstep): #setting up a dictionary key for each magnitude
        dic[str(mag)] = [] 
    
    #matching the appropriate detection rates to each size-magnitude pair and putting them in the appropriate key
    for size in sizes:
        for i in range(len(rates_list)):
            if int(sizemag_pairs[i][:2]) == size:
                mag = sizemag_pairs[i][-2:]
                dic[mag].append(int(rates_list[i]))

    #converting the dictionary to a .csv          
    df = pd.DataFrame(data=dic)
    df.to_csv(dir+'/detection_rates.csv')