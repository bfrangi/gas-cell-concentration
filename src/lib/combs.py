########################################## Comb Functions ##########################################


def comb_frequencies(f, teeth_separation, N):
    '''
    Returns an array of frequencies corresponding to the frequencies of the teeth in the comb.

    Args:
        f (float): frequency of first tooth
        teeth_separation (float): separation between teeth
        N (int): number of teeth
    '''
    import numpy as np
    return np.array([f + i*teeth_separation for i in range(N)])


def comb_sines(f, teeth_separation, N, t):
    import numpy as np
    '''
    Returns an array of sines with frequencies corresponding to the frequencies of the teeth in the
    comb.

    Args:
        f (float): frequency of first tooth
        teeth_separation (float): separation between teeth
        N (int): number of teeth
        t (array-like): time array
    '''
    frequencies = comb_frequencies(f, teeth_separation, N)
    return np.array([np.sin(2*np.pi*freq*t) for freq in frequencies])


def comb(f, teeth_separation, N, t=5, fs=240):
    import numpy as np
    '''
    Returns an array of sines with frequencies corresponding to the frequencies of the teeth in the
    comb and a time array for those sines.

    Args:
        f (float): frequency of first tooth
        teeth_separation (float): separation between teeth
        N (int): number of teeth
        t (float): time in seconds to create the comb for
        fs (float): sampling frequency (number of samples per second)
    '''
    time_array = np.linspace(0, t, int(t*fs), endpoint=False)
    return comb_sines(f, teeth_separation, N, time_array), time_array


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


def normalise_transmission(f, a, replace_outliers=False):
    '''
    Normalises the transmission spectrum by dividing it by its maximum value, after having removed
    any outliers that could distort the data.

    Args:
        f (array-like): frequency array
        a (array-like): amplitude array
    '''
    from .statistics import zscore
    z, avg, std, m = zscore(a, window=3, thresh=1.3, return_all=True)
    max_avg = max(avg)
    a /= max_avg

    if replace_outliers:
        a[~m] = avg[~m] / max_avg

    return f, a
