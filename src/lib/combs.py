def high_frequency_spectrum(f, a, f0, fs, fS, laser_wl=1.645560e-6, c=3.0e8):
    '''
    Performs an approximate transformation of the low frequency comb into the high frequency comb,
    which contains the frequencies that actually interact with and are absorbed by the gas. This is 
    an approximation because it assumes a specific wavelength for the laser, which in real life can
    vary slightly accross measurements due to temperature/current fluctuations and other 
    perturbations.

    Args:
        f (array-like): frequency array of the low-frequency comb
        a (array-like): amplitude array of the low-frequency comb
        f0 (float): center frequency of the low-frequency comb
        fs (float): frequency spacing of the low-frequency comb 
        fS (float): frequency spacing of the high-frequency comb
        laser_wl (float): wavelength of the laser, should be as close as possible to the actual 
            laser wavelength used in the experiment
        c (float): speed of light
    '''
    laser_f = c / laser_wl
    return laser_f + (f - f0)*fS/fs, a


def to_frequency(wl, amplitude, c=3.0e8):
    '''
    Converts a wavelength spectrum to a frequency spectrum. The input wavelength array is assumed to
    be in nm.

    Args:
        wl (array-like): wavelength array in nm
        amplitude (array-like): amplitude array
        c (float): speed of light
    '''
    import numpy as np
    return c * 1e9 / np.array(wl[::-1]), amplitude[::-1]


def to_wavelength(f, amplitude):
    '''
    Converts a frequency spectrum to a wavelength spectrum. The output wavelength array is in nm.

    Args:
        f (array-like): frequency array
        amplitude (array-like): amplitude array
    '''
    return to_frequency(f, amplitude)
