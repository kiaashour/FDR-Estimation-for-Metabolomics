#Methods for generating decoys from a target database. 
#The first is a naive approach and the second is a spectrum-based approach, 
#where #added fragment ions are based on fragment ions previously added.



from matchms import Spikes, Spectrum
import numpy as np

#Required function for decoy generation

#Check to see if two mz values are within a p.p.m. tolerance of each other

def within_ppm(mz_1, mz_2, ppm_tolerance):
    
    tolerance = (mz_1*ppm_tolerance)/(10**6)
    within = 1 if (abs(mz_1-mz_2) <= tolerance) else 0
    return within

#Check to see if any mzs in mz_list_2 are within a p.p.m. tolerance of any mzs in mz_list_1

def within_ppm_list(mz_list_1, mz_list_2, ppm_tolerance):
    
    for mz_1 in mz_list_1:
        
        for mz_2 in mz_list_2:
        
            if within_ppm(mz_1, mz_2, ppm_tolerance):
                return True 
    return False

#Checks to see if fragment candidate is NOT within p.p.m. tolerance of list of mzs 
#and does not have a greater mz than the precursor_mz of the spectrum being decoyed.

def mz_check(candidate_mz, mzs, spectrum, ppm_tolerance):
    
    precursor_mz = spectrum.get("precursor_mz")
    if precursor_mz < candidate_mz or within_ppm_list(mzs, [candidate_mz], ppm_tolerance):
        return False
    return True

#Searches for spectra with containing a certian precursor peak within a p.p.m. tolerance.
#If any such spectra are found, returns the first. If none are found, returns None.

#precursor mzs: Sorted precursor mzs of target spectra
#precursor_fragments: Sorted array of precursor fragments (refer to library matching notebook for form of fragments)

def precursor_search(precursor_mz, precursor_mzs, precursor_fragments, ppm_tolerance):
    
    tolerance = (precursor_mz*ppm_tolerance)/(10**6)
    lower_mz = precursor_mz - tolerance
    upper_mz = precursor_mz + tolerance
    precursors_in_tolerance = np.where((precursor_mzs >= lower_mz) & (precursor_mzs <= upper_mz))[0]
    precursor_indices = precursor_fragments[:,2]
    precursors_in_tolerance_indices = np.intersect1d(precursor_indices, precursors_in_tolerance, assume_unique = True, return_indices = True)[1]
    if len(precursors_in_tolerance_indices) > 0:
        
        index = precursors_in_tolerance_indices[0]
        precursor_fragment = precursor_fragments[index,:].tolist()
        return precursor_fragment
    return None
    
#Function to get the initial (precursor) fragment to be added to the decoy spectrum.
#Add precursor fragment, if present. 
#If not present, use precursor mz to find appropriate precursor fragment in target spectra .
#In the case where a precursor fragment isn't found, a random fragment from the spectrum is added as the initial fragment.
#Returns initial fragment as numpy array.

#precursor mzs: Sorted precursor mzs of target spectra
#precursor_fragments: Sorted array of precursor fragments (refer to library matching notebook for form of fragments)


def get_initial_fragment(spectrum, precursor_mzs, precursor_fragments, ppm_tolerance):
    
    intensities = spectrum.peaks.intensities
    mzs = spectrum.peaks.mz
    no_peaks = len(intensities)
    precursor_mz = spectrum.get("precursor_mz")                   #Precursor mz is always non-trivial in dataset
    tolerance = (precursor_mz*ppm_tolerance)/(10**6)
    lower_mz = precursor_mz - tolerance
    upper_mz = precursor_mz + tolerance
    mzs_in_tolerance = np.where((mzs >= lower_mz) & (mzs <= upper_mz))[0]
    if len(mzs_in_tolerance) > 0:
    
        intensities_for_mzs = intensities[mzs_in_tolerance]
        max_intensity_index = np.argmax(intensities_for_mzs)
        fragment_index = mzs_in_tolerance[max_intensity_index]
        initial_fragment = [mzs[fragment_index], intensities[fragment_index], np.nan, 1]
    else:

        new_precursor_fragment = precursor_search(precursor_mz, precursor_mzs, precursor_fragments, ppm_tolerance)
        if new_precursor_fragment:

            initial_fragment = new_precursor_fragment
        else:

            random_draw = np.random.randint(0, no_peaks)                
            initial_fragment = [mzs[random_draw], intensities[random_draw], np.nan, 0]   
    return np.array(initial_fragment) 



#Decoy generation functions

#Function that generates decoy spectra from target spectra based on the naive method
#Returns list of decoy spectrums.

#initial_fragments: Array of the initial fragments to be added to the decoys.
#fragments: Array of all the fragments from target spectra (not sorted)
#no_peaks_for_spectra: List contianing no. of peaks of target spectra.
#Some arguments may be rmeoved in the future

def naive_decoys(spectra, initial_fragments, fragments, no_peaks_for_spectra):

    decoy_spectra = []
    no_fragments = len(fragments)
    for i in range(len(spectra)):             
        
        spectrum = spectra[i]
        decoy_intensities = []
        decoy_mzs = []

#Adding initial (precursor) peak.

        initial_fragment = initial_fragments[i,:]
        decoy_intensities.append(initial_fragment[1])
        decoy_mzs.append(initial_fragment[0])   

#Add remaining fragment ions randomly

        desired_no_peaks = no_peaks_for_spectra[i]
        draws = np.random.randint(no_fragments, size = desired_no_peaks - 1)
        for draw in draws:
            
            fragment = fragments[draw,:]
            mz = fragment[0]
            if mz_check(mz, decoy_mzs, spectrum, 5):             

                decoy_intensities.append(fragment[1])
                decoy_mzs.append(mz)
            
        order = np.argsort(decoy_mzs)
        decoy_spectrum = Spectrum(intensities = (np.array(decoy_intensities))[order], mz = (np.array((decoy_mzs)))[order])   
        parentmass = spectrum.get("parent_mass") if spectrum.get("parent_mass") else get_parent_mass(spectrum)          
        decoy_spectrum.set("parent_mass", parentmass)
        decoy_spectra.append(decoy_spectrum)
        
    return decoy_spectra

#Function that generates decoy spectra from target spectra based on the spectrum-based method
#Return list of decoy spectra

##initial_fragments: Array of the initial fragments to be added to the decoys.
#fragments: Array of all the fragments from target spectra (not sorted)
#sorted_fragments: fragments sorted by mz
#first_fragments: List containing the indices of the first fragment of spectra in fragments
#no_peaks_for_spectra: List contianing no. of peaks of target spectra.
#Some arguments may be rmeoved in the future

def spectrum_based_decoys(spectra, initial_fragments, fragments, sorted_fragments, first_fragments, no_peaks_for_spectra):
    
    decoy_spectra = []
    fragment_mzs = sorted_fragments[:,0]
    for i in range(len(spectra)):
        
        spectrum = spectra[i]
        spectrum_peaks = np.arange(first_fragments[i], first_fragments[i+1])
        desired_no_peaks = no_peaks_for_spectra[i]
        decoy_intensities = []
        decoy_mzs = []

#Adding initial (precursor) peak.

        initial_fragment = initial_fragments[i,:]
        decoy_intensities.append(initial_fragment[1])
        decoy_mzs.append(initial_fragment[0])   

#Adding peaks based on peaks already added to the decoy spectrum. 

        fragment_candidates = []
        j = 1
        while j != desired_no_peaks:                                        

#Find fragments within mz range
        
            mz = decoy_mzs[j-1]
            tolerance = (mz*5)/(10**6)
            lower_mz = mz - tolerance
            upper_mz = mz + tolerance
            lower_index = np.searchsorted(fragment_mzs, lower_mz, side="left")
            upper_index = np.searchsorted(fragment_mzs, upper_mz, side="right")
            indices = np.arange(lower_index, upper_index)
        
#Adding 5 or less fragments ions to fragment candidates
#Fragments are stored by their indices in fragments

#First get spectra which have previously added fragment within our p.p.m. tolerance (5)
#Then randomly add 5 fragments from these spectra
#When there are less than 5 suitable spectra, randomly add fragments from suitable spectra

        
            suitable_spectra = np.unique(sorted_fragments[indices, 2])           #Get unique suitable spectra indices   
            if len(suitable_spectra) >= 5:
            
                spectra_draws = np.random.choice(suitable_spectra, 5, replace = False)
                for index in spectra_draws:
                
                    index = int(index)
                    peak_draw = np.random.choice(no_peaks_for_spectra[index])
                    fragment_draw = first_fragments[index] + peak_draw
                    fragment_candidates.append(fragment_draw)
            else:

                spectra_draws = np.random.choice(suitable_spectra, 5, replace = True)   #Note that each spectra must have >=5 peaks    
                for index in spectra_draws:

                    index = int(index)
                    peak_draw = np.random.choice(no_peaks_for_spectra[index])
                    fragment_draw = first_fragments[index] + peak_draw
                    fragment_candidates.append(fragment_draw) 
            
#Add fragment from fragment candidates for decoy spectrum.
#Only add fragment if it is not within 5 p.p.m. of previously added fragments
#and not greater than the precursor mz of the spectrum being decoyed.

#If chosen fragment is not suitable, remove fragment from fragment_candidates_clone 
#and search again.
    
            search_attempts = 0
            fragment_found = False
            no_candidates = len(fragment_candidates)
            fragment_candidates_clone = fragment_candidates
            while search_attempts != no_candidates and not fragment_found:

                choice = np.random.choice(fragment_candidates_clone)
                fragment_candidate = fragments[choice,:]
                candidate_mz = fragment_candidate[0]
                if mz_check(candidate_mz, decoy_mzs, spectrum, 5):    

                    decoy_intensities.append(fragment_candidate[1])
                    decoy_mzs.append(candidate_mz)
                    fragment_found = True
                    j += 1
                else:
                
                    fragment_candidates_clone.remove(choice)

                search_attempts += 1

#If fragment search is over and a suitable fragment candidate isn't found, then 
#add random fragment from spectrum being decoyed
        
            if not fragment_found:

                choice = np.random.choice(spectrum_peaks)
                fragment_candidate = fragments[choice,:]
                decoy_intensities.append(fragment_candidate[1])
                decoy_mzs.append(fragment_candidate[0])
                j += 1
                print("Suitable fragment candidate not found. Random fragment from real spectrum added.")

        
        order = np.argsort(decoy_mzs)
        decoy_spectrum = Spectrum(intensities = (np.array(decoy_intensities))[order], mz = (np.array((decoy_mzs)))[order]) 
        parentmass = spectrum.get("parent_mass") if spectrum.get("parent_mass") else get_parent_mass(spectrum)          
        decoy_spectrum.set("parent_mass", parentmass)
        decoy_spectra.append(decoy_spectrum)

    return decoy_spectra


#Note that these custom functions will be changed in future work for simpler usage.