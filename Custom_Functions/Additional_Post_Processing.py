#Additional post-processing functions for the GNPS dataset used in the decoy generation experiment.

from matchms.filtering import select_by_relative_intensity



#Keep only spectra with precursor_mz below max_mz, default 1000 Da.
    
def require_precursor_below_mz(spectrum_in, max_mz = 1000):

    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    precursor_mz = spectrum.get("precursor_mz", None)
    assert precursor_mz is not None, "Precursor mz absent."
    assert isinstance(precursor_mz, (float, int)), ("Expected 'precursor_mz' to be a scalar number.",
                                                    "Consider applying 'add_precursor_mz' filter first.")
    assert max_mz >= 0, "max_mz must be a positive scalar."
    if precursor_mz >= max_mz:
        return None

    return spectrum
   


#Discard spectra with < no_peaks peaks with relative intensity above intensity_percent
#no_peaks, intensity_percent = 5, 2 default

def require_minimum_of_high_peaks(spectrum_in, no_peaks = 5,
                                  intensity_percent = 2.0):


    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    assert no_peaks >= 1, "no_peaks must be a positive nonzero integer."
    assert 0 <= intensity_percent <= 100, "intensity_percent must be a scalar between 0-100."
    intensities_above_p = select_by_relative_intensity(spectrum, intensity_from=intensity_percent/100, intensity_to=1.0)
    if len(intensities_above_p.peaks) < no_peaks:
        return None

    return spectrum