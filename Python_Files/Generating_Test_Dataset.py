#!/usr/bin/env python
# coding: utf-8

# # Generating Test Dataset

# Function for generating random spectra.

# ## Setup

# In[ ]:


import numpy as np
from matchms import Spectrum


# In[ ]:


#Testing



#Function to generate random spectra, with a certain mz and intensity range, and a no. of peaks range.


def generate_random_spectra(n, id_key, type_spectra, intensity_range = (0.1, 0.7), mz_range = (30, 300), peaks_range = (3, 11)):
    
    spectra = []
    for i in range(n):
        
        no_peaks = np.random.randint(peaks_range[0], peaks_range[1])
        intensities = np.random.uniform(intensity_range[0], intensity_range[1], no_peaks)
        mzs = np.random.uniform(mz_range[0], mz_range[1], no_peaks)
        mzs = np.sort(mzs)
        inchikey = np.random.randint(0, 11)
        spectrum = Spectrum(intensities=intensities, mz=mzs, metadata={'id': f'spectrum {i}{id_key}'})
        if type_spectra:
            
            inchikey = np.random.randint(0, 16)
            spectrum.set("inchikey", f'{inchikey}')
        spectra += [spectrum]
    return spectra

