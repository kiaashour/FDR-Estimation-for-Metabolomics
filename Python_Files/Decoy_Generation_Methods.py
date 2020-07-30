#!/usr/bin/env python
# coding: utf-8

# # Decoy Generation Methods 

# Methods for generating decoys from a target database. The first is a naive approach and the second is a spectrum-based approach, where added fragment ions are based on fragment ions previously added.
# 
# Both methods can be read about [here](file:///C:/Users/kia/Downloads/Metabolomics%20Internship/Papers/Converting%20from%20proteomics%20to%20metabolomics.pdf).

# ## Setup

# In[4]:


from matchms import Spikes, Spectrum
import numpy as np


# ## Preliminary Functions

# In[22]:


#Preliminary functions for the decoy generation methods.



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

#Checks to see if fragment candidate is NOT within p.p.m. tolerance of added fragments and does not have a greater mz than 
#the precursor_mz of the spectrum being decoyed.

def mz_check(fragment_candidate, added_fragments, spectrum, ppm_tolerance):
    
    candidate_mz = fragment_candidate[1]
    precursor_mz = spectrum.get("precursor_mz") if spectrum.get("precursor_mz") else get_precursor_mz(spectrum)
    mzs = []
    for fragment in added_fragments:
        mzs += [fragment[1]]
    if precursor_mz:
        
        if precursor_mz < candidate_mz or within_ppm_list(mzs, [candidate_mz], ppm_tolerance):
            return False
    
    else:
        
        if within_ppm_list(mzs, [candidate_mz], ppm_tolerance):
            return False
            
    return True

#Searches for spectra with a certain precursor mz and its associated precursor peak.
#If any such spectra are found, return the first. If none are found, return None.

def precursor_search(precursor_mz, target_database):
    
    for spectrum in target_database:
        
        mzs = spectrum.peaks.mz  
        new_precursor_mz = spectrum.get("precursor_mz") if spectrum.get("precursor_mz") else get_precursor_mz(spectrum)
        if (precursor_mz in mzs) and (new_precursor_mz == precursor_mz):
            
            return spectrum
    return None
    
#Gets missing precursor_mz via method used in Florian's spec2vec matching experiments.
#Set 0.0 if can't find

def get_precursor_mz(spectrum):
    
    fix_mass = 0.0
    try:
        for history in spectrum.metadata['annotation_history']:
            fix_mass_test = float(history['Precursor_MZ'])
            fix_mass = max(fix_mass, fix_mass_test)
    except KeyError:
        return fix_mass
    
    return fix_mass

#Gets missing parent_mass via method used in Florian's spec2vec matching experiments.
#Set 0.0 if can't find

def get_parent_mass(spectrum):

    fix_mass = 0.0
    parent_mass = 0.0
    try:
        for history in spec.metadata['annotation_history']:
            fix_mass_test = float(history['Precursor_MZ'])
            fix_mass = max(fix_mass, fix_mass_test)
        charge = spec.get("charge")
        protons_mass = 1.00727645199076 * charge
        precursor_mass = fix_mass * abs(charge)
        parent_mass = precursor_mass - protons_mass
    except KeyError:
        return parent_mass
    
    return parent_mass

#Function to get the initial (precursor) fragment to be added to the decoy spectrum.
#Add precursor fragment, if present. 
#If not present, use selected ion mass to find appropriate precursor fragment in (target) spectra .
#In the case where a precursor fragment isn't found, a random fragment from the spectrum is added as the initial fragment.
#Returns initial fragment to be added to the decoy spectrum.

def get_initial_fragment(spectrum, spectra):
    
    intensities = spectrum.peaks.intensities
    mzs = spectrum.peaks.mz
    no_peaks = len(intensities)
    precursor_mz = spectrum.get("precursor_mz")
    if precursor_mz:
        
        if precursor_mz in mzs:
        
            index = np.where(mzs == precursor_mz)[0][0]
            initial_fragment = [intensities[index], precursor_mz]
        else:
            
            new_spectrum = precursor_search(precursor_mz, spectra)
            if new_spectrum:
                
                new_mzs = new_spectrum.peaks.mz
                new_intensities = new_spectrum.peaks.intensities
                index = np.where(new_mzs == precursor_mz)[0][0]
                initial_fragment = [new_intensities[index], precursor_mz]
            else:
                
                random_draw = np.random.randint(0, no_peaks)                
                initial_fragment = [intensities[random_draw], mzs[random_draw]]   

#The case where new_spectrum is None is accounted for by adding a random initial fragment from the target spectrum.
#The case where new_precursor_mz is 0 is accounted for as can't have mz being 0 for a fragment. 

    else:
        
        new_precursor_mz = get_precursor_mz(spectrum)
        if new_precursor_mz in mzs:
        
            index = np.where(mzs == new_precursor_mz)[0][0]
            initial_fragment = [intensities[index], new_precursor_mz]
        else:
            
            new_spectrum = precursor_search(new_precursor_mz, spectra)
            if new_spectrum:
                
                new_mzs = new_spectrum.peaks.mz
                new_intensities = new_spectrum.peaks.intensities
                index = np.where(new_mzs == new_precursor_mz)[0][0]
                initial_fragment = [new_intensities[index], new_precursor_mz]
            else:
                
                random_draw = np.random.randint(0, no_peaks)                
                initial_fragment = [intensities[random_draw], mzs[random_draw]]
                
    return initial_fragment
    
#Add fragments (from a spectrum) that are not already in a list of fragments to the list.
#Returns the original list of fragment  with the new fragments added

def add_unique_fragments(fragment_set, spectrum):
    
    intensities = spectrum.peaks.intensities
    mzs = spectrum.peaks.mz
    for i in range(0, len(intensities)):
        
        fragment = [intensities[i], mzs[i]]
        if fragment not in fragment_set:
            
            fragment_set += [fragment]
    return fragment_set
    


# ## Naive Method of Generating Decoys

# In[10]:


#Generating decoy spectra from target spectra by the naive method



def naive_decoys(spectra):

#Getting all unique fragment ions from target database

    fragment_ions = []
    for spectrum in spectra:                              

        fragment_ions = add_unique_fragments(fragment_ions, spectrum)

#Creating decoy database by adding random fragment ions to empty spectra until desired no. of fragments is reached
#We allow for fragment ions with the same mz value to be added here (but not identical fragment ions)

    decoy_spectra = []
    no_fragment_ions = len(fragment_ions)
    np.random.seed(0)                          #to make it reproducible
    for spectrum in spectra:             
    
        intensities = spectrum.peaks.intensities            
        mzs = spectrum.peaks.mz  
        decoy_intensities = []
        decoy_mzs = []
        added_fragments = []
        no_peaks_added = 0 

#Adding initial (precursor) peak.

        initial_fragment = get_initial_fragment(spectrum, spectra)
        decoy_intensities += [initial_fragment[0]]
        decoy_mzs += [initial_fragment[1]]
        added_fragments += [initial_fragment]
        no_peaks_added += 1
   

#Add remaining fragment ions randomly
    
        desired_no_peaks = len(intensities) 
        precursor_mz = spectrum.get("precursor_mz") if spectrum.get("precursor_mz") else get_precursor_mz(spectrum)          
        while no_peaks_added != desired_no_peaks:

#Only add fragment ions if mz less than precursor mz. 
#Accounts for case when precursor_mz is 0.

            random_draw = np.random.randint(0, no_fragment_ions)
            random_fragment = fragment_ions[random_draw]
            if mz_check(random_fragment, added_fragments, spectrum, 5):             
                
                decoy_intensities += [random_fragment[0]]
                decoy_mzs += [random_fragment[1]]
                added_fragments += [random_fragment]
                no_peaks_added += 1

        order = np.argsort(decoy_mzs)
        decoy_spectrum = Spectrum(intensities = (np.array(decoy_intensities))[order], mz = (np.array((decoy_mzs)))[order])   
        parentmass = spectrum.get("parent_mass") if spectrum.get("parent_mass") else get_parent_mass(spectrum)          
        decoy_spectrum.set("parent_mass", parentmass)
        decoy_spectra += [decoy_spectrum]
        
    return decoy_spectra


# ## Spectrum-Based Method of Generating Decoys

# In[23]:


# Generating decoy spectra by the spectrum-based method from target spectra



def spectrum_based_decoys(spectra):
    
    decoy_spectra = []
    np.random.seed(0)                           #to make it reproducible
    for spectrum in spectra:

        desired_no_peaks = len(spectrum.peaks)
        decoy_intensities = []
        decoy_mzs = []
        added_fragments = []

#Adding initial (precursor) peak.

        initial_fragment = get_initial_fragment(spectrum, spectra)
        decoy_intensities += [initial_fragment[0]]
        decoy_mzs += [initial_fragment[1]]
        added_fragments += [initial_fragment]

#Adding peaks based on peaks already added to the decoy spectrum. 

        fragment_candidates = []
        suitable_fragments_steps = []                     #suitable fragments for fragment candidates for each peak
        i = 1
        while i != desired_no_peaks:                                        

            mz = decoy_mzs[i-1]
            suitable_fragments = []
            for spectrum_2 in spectra:                                      

#Get suitable fragments for fragment candidates list

                mzs_2 = spectrum_2.peaks.mz                      
                if within_ppm_list([mz], mzs_2, 5):  

                    suitable_fragments = add_unique_fragments(suitable_fragments, spectrum_2)

#Adding 5 or less fragment ions (uniformly) to fragment candidates.
#When there are no suitable fragments, look at previously added fragment candidates.  

            suitable_fragments_steps += [suitable_fragments]                           
            search_attempts = 0
            fragment_found = False
            while search_attempts != i and not fragment_found:

#Note that suitable_fragments is non-empty as "spectrum_2" can be spectrum

                suitable_fragments = suitable_fragments_steps[i-search_attempts-1]          
                no_suitable = len(suitable_fragments)                                  
                if no_suitable >= 5:

                    random_draws = np.random.choice(len(suitable_fragments), 5, replace=False)
                    draws = [suitable_fragments[random_draw] for random_draw in random_draws]
                    for fragment_candidate in draws:                                          

                        if fragment_candidate not in fragment_candidates:             

                            fragment_candidates += [fragment_candidate]
                elif no_suitable > 0:

                    for k in range(0, no_suitable):

                        fragment_candidate = suitable_fragments[k]
                        if fragment_candidate not in fragment_candidates:

                            fragment_candidates += [fragment_candidate]


                random_choice = np.random.choice(len(fragment_candidates))
                fragment_candidate = fragment_candidates[random_choice]                
                if mz_check(fragment_candidate, added_fragments, spectrum, 5):    

                    added_fragments += [fragment_candidate]
                    decoy_intensities += [fragment_candidate[0]]
                    decoy_mzs += [fragment_candidate[1]]
                    fragment_found = True
                    i += 1

                search_attempts += 1

#If fragment search is over, a suitable fragment candidate isn't found and decoy_mzs is not empty, then 
#exit out of the searching for more decoy peaks.

            if not fragment_found:

                print("Suitable fragment candidate not found. Non-empty decoy spectrum is returned without further search.")
                break

#Note that decoy spectrum will always be non-empty because of the initial added fragment.

        order = np.argsort(decoy_mzs)
        decoy_spectrum = Spectrum(intensities = (np.array(decoy_intensities))[order], mz = (np.array((decoy_mzs)))[order]) 
        parentmass = spectrum.get("parent_mass") if spectrum.get("parent_mass") else get_parent_mass(spectrum)          
        decoy_spectrum.set("parent_mass", parentmass)
        decoy_spectra += [decoy_spectrum]
        
    return decoy_spectra

