from .filters import limit_filter, threshold_filter


######################## Simple Absorption and Transmission Curve Analysis #########################


def absorption_curve(x, y_sample, y_reference, threshold=1.5e-5, f_threshold=20e3, f_limit=60e3):
    """
    Returns the absorption curve of a sample and reference spectrum by directly comparing their 
    amplitudes. The absorption curve is calculated as the percentage difference between the sample 
    and reference amplitudes.

    Args:
        x (array-like): The frequency array.
        y_sample (array-like): The amplitude of the sample spectrum.
        y_reference (array-like): The amplitude of the reference spectrum.
        threshold (float): The threshold for the sample amplitude. All amplitudes below this value 
            are removed from the result.
        f_threshold (float): The threshold for the frequency array. All frequencies below this value
            are removed from the result.
        f_limit (float): The limit for the frequency array. All frequencies above this value are
            removed from the result.
    """
    x, y_sample, y_reference = limit_filter(x, f_limit, y_sample, y_reference)
    x, y_sample, y_reference = threshold_filter(
        x, f_threshold, y_sample, y_reference)
    y_sample, x, y_reference = threshold_filter(
        y_sample, threshold, x, y_reference)

    diff = -(y_sample - y_reference) / y_reference * 100
    diff, x = threshold_filter(diff, 0, x)

    return x, diff


def transmission_curve(x, y_sample, y_reference, threshold=1.5e-5, f_threshold=20e3, f_limit=60e3):
    """
    Returns the transmission curve of a sample and reference spectrum by directly comparing their
    amplitudes. The transmission curve is calculated as the percentage difference between the
    reference and sample amplitudes.

    Args:
        x (array-like): The frequency array.
        y_sample (array-like): The amplitude of the sample spectrum.
        y_reference (array-like): The amplitude of the reference spectrum.
        threshold (float): The threshold for the sample amplitude. All amplitudes below this value
            are removed from the result.
        f_threshold (float): The threshold for the frequency array. All frequencies below this value
            are removed from the result.
        f_limit (float): The limit for the frequency array. All frequencies above this value are
            removed from the result.
    """
    x, diff = absorption_curve(
        x, y_sample, y_reference, threshold, f_threshold, f_limit)
    diff = 100 - diff
    return x, diff

###################################### Comb Spectrum Analysis ######################################


class CombAnalyser:
    """Generates the absorption and transmission curves of a sample and reference spectrum.
    Frequency arrays must be the same for both sample and reference data.
    
    Args:
        f_sample (array-like): The frequency array of the sample spectrum.
        a_sample (array-like): The amplitude of the sample spectrum.
        f_reference (array-like): The frequency array of the reference spectrum.
        a_reference (array-like): The amplitude of the reference spectrum."""

    def __init__(self, f_sample, a_sample, f_reference, a_reference):
        if f_sample.any() != f_reference.any():
            raise ValueError(
                'The frequencies of the sample and reference data must be the same.')
        self.f_sample = f_sample
        self.a_sample = a_sample
        self.f_reference = f_reference
        self.a_reference = a_reference

        self.absorption_amp = None
        self.absorption_freq = None

        self.scaling_factor = 1

    def get_absorption_curve(self):
        self.absorption_amp = (self.a_reference - self.a_sample) / \
            self.a_reference * self.scaling_factor
        self.absorption_freq = self.f_sample
        return self.absorption_freq, self.absorption_amp

    def get_transmission_curve(self):
        if self.absorption_amp is None or self.absorption_freq is None:
            self.get_absorption_curve()
        return self.absorption_freq, self.scaling_factor - self.absorption_amp

    def generate_absorption_plot(self, interp=False):
        if self.absorption_freq is None or self.absorption_amp is None:
            self.get_absorption_curve()

        x = self.absorption_freq
        y = self.absorption_amp

        from matplotlib import pyplot as plt
        plt.scatter(x, y)
        if interp:
            from lib.fitting import loose_interpolation
            x, y = loose_interpolation(x, y)
            plt.plot(x, y)
        plt.title('Absorption Spectrum')
        return plt

    def show_absorption_plot(self, interp=False):
        return self.generate_absorption_plot(interp=interp).show()

    def generate_transmission_plot(self, interp=False):
        if self.absorption_freq is None or self.absorption_amp is None:
            self.get_absorption_curve()

        x = self.absorption_freq
        y = self.absorption_amp

        from matplotlib import pyplot as plt
        plt.scatter(x, self.scaling_factor - y)
        if interp:
            from lib.fitting import loose_interpolation
            x, y = loose_interpolation(x, y)
            plt.plot(x, self.scaling_factor - y)
        plt.title('Transmission Spectrum')
        return plt

    def show_transmission_plot(self, interp=False):
        return self.generate_transmission_plot(interp=interp).show()

    def show_reference_spectrum(self):
        from matplotlib import pyplot as plt
        plt.stem(self.f_reference, self.a_reference)
        plt.title('Reference Spectrum')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()

    def show_sample_spectrum(self):
        from matplotlib import pyplot as plt
        plt.stem(self.f_sample, self.a_sample)
        plt.title('Sample Spectrum')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.show()


class CombSpectrumAnalyser:
    """Obtains the comb spectrum of a given frequency and amplitude data, selecting the specific
    frequencies that correspond to the comb teeth."""
    def __init__(self, freq, amp, center_freq=None, freq_spacing=None, number_of_teeth=None):
        self.freq = freq
        self.amp = amp

        self.center_freq = 40000.0 if center_freq is None else center_freq
        self.freq_spacing = 200.0 if freq_spacing is None else freq_spacing
        self.number_of_teeth = 38 if number_of_teeth is None else number_of_teeth

        self.extracted_freq = None
        self.extracted_amp = None

    def get_comb_frequencies(self):
        center_freq = self.center_freq
        rep_rate = self.freq_spacing
        n_teeth = min(int(center_freq // rep_rate * 2), self.number_of_teeth)
        return [float(center_freq + (i - n_teeth//2) * rep_rate) for i in range(n_teeth)]

    def closest_frequency(self, f):
        from numpy import argmin
        return self.freq[argmin(abs(self.freq - f))]

    def extract_teeth(self):
        from numpy import nonzero, array

        comb_frequencies = array(self.get_comb_frequencies())
        approx_comb_frequencies = array(
            [self.closest_frequency(f) for f in comb_frequencies])
        indices = nonzero(approx_comb_frequencies[:, None] == self.freq)[1]

        y = self.amp[indices]

        self.extracted_freq = approx_comb_frequencies
        self.extracted_amp = y
        return self.extracted_freq, self.extracted_amp

    # Plots

    def generate_raw_spectrum_plot(self):
        from matplotlib import pyplot as plt
        plt.stem(self.freq, self.amp)
        return plt

    def show_raw_spectrum_plot(self):
        self.generate_raw_spectrum_plot().show()

    def generate_spectrum_plot(self):
        from matplotlib import pyplot as plt

        rep = self.freq_spacing
        f0 = self.center_freq

        plt.stem(self.extracted_freq, self.extracted_amp)
        plt.title(
            f'Comb Spectrum showing extracted teeth and their amplitudes.\nCentral frequency ${f0}$ Hz and repetition rate ${rep}$ Hz.')
        return plt

    def show_spectrum_plot(self):
        self.generate_spectrum_plot().show()

class AbsorptionAnalyser:
    def __init__(self, filename_sample, filename_reference, center_freq=40000, freq_spacing=200, 
                 number_of_teeth=38, high_freq_modulation=1e9, laser_wavelength=1645.56e-9):

        from .combs import high_frequency_spectrum, normalise_transmission
        from .fft import fft_plot_data
        from .labview import Reader
        import sys

        reader_sample = Reader(filename_sample)
        sample_data = reader_sample.extract_data()

        reader_reference = Reader(filename_reference)
        reference_data = reader_reference.extract_data()

        t_sample = sample_data['time']
        amplitude_sample = sample_data['amplitude']

        t_reference = reference_data['time']
        amplitude_reference = reference_data['amplitude']


        # Check that time arrays are the same in both time series
        try:
            if not (t_sample == t_reference).all():
                raise ValueError
            t = t_sample
        except ValueError:
            print('Time arrays are not the same in both time series')
            sys.exit(1)

        x_sample, y_sample = fft_plot_data(t, amplitude_sample)
        x_reference, y_reference = fft_plot_data(t, amplitude_reference)

        # Analyse the comb spectrum

        sample_analyzer = CombSpectrumAnalyser(x_sample, y_sample, center_freq, freq_spacing, 
                                               number_of_teeth)
        sample_freq, sample_amp = sample_analyzer.extract_teeth()

        reference_analyzer = CombSpectrumAnalyser(x_reference, y_reference, center_freq, 
                                                  freq_spacing, number_of_teeth)
        ref_freq, ref_amp = reference_analyzer.extract_teeth()

        self.comb_analyser = CombAnalyser(sample_freq, sample_amp, ref_freq, ref_amp)
        self.f_radio_frequency, self.a_radio_frequency = self.comb_analyser.get_transmission_curve()

        # Transform the low frequency comb spectrum to the high frequency comb spectrum
        f, a = high_frequency_spectrum(self.f_radio_frequency, self.a_radio_frequency, 
                                                     f0=center_freq, fs=freq_spacing, 
                                                     fS=high_freq_modulation, laser_wl=laser_wavelength)

        # Normalise the transmission curve
        self.f, self.a = normalise_transmission(f, a, replace_outliers=False)
        
    def generate_absorption_plot(self, interp=False):
        x, y = self.f, self.a

        from matplotlib import pyplot as plt
        plt.scatter(x, y)
        if interp:
            from lib.fitting import loose_interpolation
            x, y = loose_interpolation(x, y)
            plt.plot(x, y)
        plt.title('Absorption Spectrum')
        return plt

    def show_absorption_plot(self, interp=False):
        return self.generate_absorption_plot(interp=interp).show()

    def generate_transmission_plot(self, interp=False):
        x, y = self.f, self.a

        from matplotlib import pyplot as plt
        plt.scatter(x, 1 - y)
        if interp:
            from lib.fitting import loose_interpolation
            x, y = loose_interpolation(x, y)
            plt.plot(x, 1 - y)
        plt.title('Transmission Spectrum')
        return plt

    def show_transmission_plot(self, interp=False):
        return self.generate_transmission_plot(interp=interp).show()


class ConcentrationFitter:
    def __init__(self, filename_sample, filename_reference, temperature=296, length=0.1,
                 center_freq=40000, freq_spacing=200, number_of_teeth=38, high_freq_modulation=1e9,
                 laser_wavelength=1645.56e-9, correction=None):
        self.filename_sample = filename_sample
        self.filename_reference = filename_reference
        self.temperature = temperature
        self.length = length
        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth
        self.high_freq_modulation = high_freq_modulation
        self.laser_wavelength = laser_wavelength

        self.a_transmission_curve_measurement = None
        self.f_transmission_curve_measurement = None
        self.a_transmission_curve_spectralcalc = None
        self.f_transmission_curve_spectralcalc = None

        self.absorption_analyser = None

        self.concentration = None

        self.correction = correction

    def get_concentration(self):
        self.adjust_measurement_transmission_curve()
        return self.concentration

    def get_measurement_transmission_curve(self):
        aa = AbsorptionAnalyser(self.filename_sample, self.filename_reference,
                                center_freq=self.center_freq, freq_spacing=self.freq_spacing,
                                number_of_teeth=self.number_of_teeth, 
                                high_freq_modulation=self.high_freq_modulation, 
                                laser_wavelength=self.laser_wavelength)
        self.absorption_analyser = aa
        if self.correction:
            import numpy as np
            self.a_pre_correction = np.array([i for i in aa.a])
            aa.a /= self.correction(aa.f)

        self.f_transmission_curve_measurement, self.a_transmission_curve_measurement = aa.f, aa.a
        return self.f_transmission_curve_measurement, self.a_transmission_curve_measurement

    def get_spectralcalc_transmissison_curve(self):
        from lib.fitting import intersect_onto_common_grid
        from lib.combs import to_frequency
        from lib.spectralcalc import interpolated_transmission_curve
        import numpy as np
        from scipy.optimize import minimize

        if any(x is None for x in [self.f_transmission_curve_measurement, self.a_transmission_curve_measurement]):
            self.get_measurement_transmission_curve()

        # Fit concentration
        def f(conc, f_sample, a_sample):
            # Get the spectralcalc curve
            wl_ref, a_ref = interpolated_transmission_curve(
                conc, self.temperature, self.length)

            # Transform the spectralcalc data to frequency
            f_ref, a_ref = to_frequency(wl_ref, a_ref)

            # Shift the sample spectrum to overlap with the spectralcalc data as much as possible
            # f_sample = overlap_data(f_sample, a_sample, f_ref, a_ref)

            # Interpolate the sample data onto the common grid
            f_sample_com, a_sample_com, f_ref_com, a_ref_com = intersect_onto_common_grid(
                f_sample, a_sample, f_ref, a_ref)
            return np.sum((a_sample_com - a_ref_com)**2)

        result = minimize(f, 0.001, args=(self.f_transmission_curve_measurement,
                          self.a_transmission_curve_measurement), bounds=[(0, 1)])
        concentration = result.x[0]
        self.concentration = concentration

        # Get the spectralcalc curve
        wl_spectralcalc, amplitude_spectralcalc = interpolated_transmission_curve(
            concentration, self.temperature, self.length)

        # Transform the spectralcalc data to frequency
        self.f_transmission_curve_spectralcalc, self.a_transmission_curve_spectralcalc = to_frequency(
            wl_spectralcalc, amplitude_spectralcalc)

        return self.f_transmission_curve_spectralcalc, self.a_transmission_curve_spectralcalc

    def adjust_measurement_transmission_curve(self):
        from lib.fitting import overlap_data

        if any(x is None for x in [self.f_transmission_curve_measurement, self.a_transmission_curve_measurement, self.f_transmission_curve_spectralcalc, self.a_transmission_curve_spectralcalc]):
            self.get_spectralcalc_transmissison_curve()

        # # Shift the sample spectrum to overlap with the spectralcalc data as much as possible
        # f_measurement = overlap_data(self.f_transmission_curve_measurement, self.a_transmission_curve_measurement,
        #                              self.f_transmission_curve_spectralcalc, self.a_transmission_curve_spectralcalc)
        # self.f_transmission_curve_measurement = f_measurement
        return self.f_transmission_curve_measurement

    def get_transmission_curves_plot(self):
        from lib.combs import to_wavelength
        import matplotlib.pyplot as plt

        if self.concentration is None:
            self.get_concentration()

        # Transform everything back to wavelength
        wl_sample, a_sample = to_wavelength(
            self.f_transmission_curve_measurement, self.a_transmission_curve_measurement)
        wl_spectralcalc, a_spectralcalc = to_wavelength(
            self.f_transmission_curve_spectralcalc, self.a_transmission_curve_spectralcalc)
        
        if self.correction:
            wl, a = to_wavelength(self.absorption_analyser.f, self.correction(self.absorption_analyser.f))
            wl_pre, a_pre = to_wavelength(self.absorption_analyser.f, self.a_pre_correction)

        # Plot the transmission curve
        plt.plot(wl_spectralcalc, a_spectralcalc, 'b-', label=f'SpectralCalc for c = {self.concentration:.4f} VMR')
        plt.plot(wl_sample, a_sample, 'ro-', label='Corrected Experiment' if self.correction else 'Experiment')
        if self.correction:
            plt.plot(wl_pre, a_pre, 'o-', label='Non Corrected Experiment')
            plt.plot(wl, a, 'g-', label='Correction')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Transmission')
        plt.title(f'Plot of the transmission as a function of wavelength ({self.temperature} K, {self.concentration:.4f} VMR)')
        plt.legend()
        return plt

    def show_transmission_curves_plot(self):
        self.get_transmission_curves_plot().show()
