#!/usr/bin/env python
# coding: utf-8

# # Additional Post-Processing

# Additional post-processing functions for the GNPS dataset used in the decoy generation experiment.

# ## Setup

# In[1]:


import numpy as np
from matchms import Spikes
from matchms.filtering import select_by_relative_intensity




# ## Post-Processing Functions

# In[6]:


#Post-processing functions




from Decoy_Generation_Methods import get_precursor_mz

#Remove peaks based on peaks within a mz_window Da range of the top K peaks
#Keep peak only if it is in the top K peaks within a mz_window Da range, mz_window,, K =  50 Da, K 10 default

def window_process(spectrum, k = 6, mz_window = 50):
    
    mzs = spectrum.peaks.mz
    intensities = spectrum.peaks.intensities
    top_k = np.argsort(intensities)[::-1][0:k]
    K_ordered_mzs = mzs[top_k]
    indices = [i for i in range(len(mzs)) if i not in top_k]
    new_mzs = mzs
    new_intensities = intensities
    for i in indices:
        
        compare = abs(mzs[i]-K_ordered_mzs) <= mz_window
        if True not in compare:
            new_mzs[i] = np.nan
            new_intensities[i] = np.nan
    
    nans = np.isnan(new_mzs)
    new_mzs = new_mzs[~nans]
    new_intensities = new_intensities[~nans]
    spectrum.peaks = Spikes(mz=new_mzs, intensities=new_intensities)
    
    return spectrum

#Keep only spectra with precursor_mz below max_mz, default 1000 Da.

def remove_spectra_below(spectrum, max_mz = 1000):
    
    precursor_mz = spectrum.get("precursor_mz") if spectrum.get("precursor_mz") else get_precursor_mz(spectrum)
    if precursor_mz >= max_mz:
        return None
    
    return spectrum
   
#Remove everything apart from the precursor that is within mz_tolerance (Da) of the precursor, 17 Da default

def remove_spectra_within_tolerance(spectrum, mz_tolerance = 17):
    
    precursor_mz = spectrum.get("precursor_mz") if spectrum.get("precursor_mz") else get_precursor_mz(spectrum)
    mzs = spectrum.peaks.mz
    intensities = spectrum.peaks.intensities
    new_mzs = mzs
    new_intensities = intensities
    if precursor_mz:
        for i in range(len(mzs)):
            
            if abs(precursor_mz-mzs[i]) <= mz_tolerance and mzs[i] != precursor_mz:
                new_mzs[i] = np.nan
                new_intensities[i] = np.nan
                
        nans = np.isnan(new_mzs)
        new_mzs = new_mzs[~nans]
        new_intensities = new_intensities[~nans]
        spectrum.peaks = Spikes(mz=new_mzs, intensities=new_intensities)

    return spectrum

#Discard spectra with <no_peaks peaks with relative intensity above intensity_percent%, no_peaks, intensity_percent = 5, 2 default

def remove_spectra_above_intensity(spectrum, no_peaks = 5, intensity_percent = 2):

    intensities = spectrum.peaks.intensities
    intensities_above_p = select_by_relative_intensity(spectrum, intensity_from=intensity_percent/100, intensity_to=1.0)
    if len(intensities_above_p.peaks) < no_peaks:
        return None
    
    return spectrum

