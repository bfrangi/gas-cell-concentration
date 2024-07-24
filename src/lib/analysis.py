class CombAnalyser:
    """Generates the absorption and transmission curves of a sample and reference spectrum.
    Frequency arrays must be the same for both sample and reference data.

    Args:
        f_sample (array-like): The frequency array of the sample spectrum.
        a_sample (array-like): The amplitude of the sample spectrum.
        f_reference (array-like): The frequency array of the reference spectrum.
        a_reference (array-like): The amplitude of the reference spectrum."""

    def __init__(self, f_sample, a_sample, f_reference, a_reference):
        if not (f_sample == f_reference).all():
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

    # Plots

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

    def __init__(self, freq, amp, center_freq=40000.0, freq_spacing=200.0, number_of_teeth=38):
        self.freq = freq
        self.amp = amp

        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth

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

class SpectralCalcFitter:
    """Fits the experimental data to a given spectralcalc data, removing the etalon effect if 
    requested."""

    def __init__(self, filename_sample, filename_reference, filename_spectralcalc,
                 scaling_factor=1.0, center_freq=40000, freq_spacing=200, number_of_teeth=38,
                 high_freq_modulation=1e9, etalon_removal=()):
        self.filename_sample = filename_sample
        self.filename_reference = filename_reference
        self.filename_spectralcalc = filename_spectralcalc
        self.scaling_factor = scaling_factor
        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth
        self.high_freq_modulation = high_freq_modulation
        self.etalon_removal = etalon_removal

        self.t = None
        self.amplitude_sample = None
        self.amplitude_reference = None

        self.wl_spectralcalc = None
        self.power_spectralcalc = None
        self.pressure = None
        self.temperature = None
        self.concentration = None
        self.length = None

        self.fft_sample_x = None
        self.fft_sample_y = None
        self.fft_reference_x = None
        self.fft_reference_y = None

        self.comb_analyser = None

        self.wl_sample = None
        self.power_sample = None
        self.wl_reference = None
        self.power_reference = None

        self.wl_sample_corrected = None
        self.power_sample_corrected = None
        self.x_sine = None
        self.y_sine = None

    def read_data(self):
        from lib.readers import LVMReader, SpectralcalcReader

        reader_sample = LVMReader(self.filename_sample)
        sample_data = reader_sample.extract_data()

        reader_reference = LVMReader(self.filename_reference)
        reference_data = reader_reference.extract_data()

        reader_spectralcalc = SpectralcalcReader(self.filename_spectralcalc)
        spectralcalc_data = reader_spectralcalc.extract_data()

        t_sample = sample_data['time']
        self.amplitude_sample = sample_data['amplitude']

        t_reference = reference_data['time']
        self.amplitude_reference = reference_data['amplitude']

        # Check that time arrays are the same in both time series
        if not (t_sample == t_reference).all():
            raise ValueError(
                'Time arrays are not the same in both time series')
        self.t = t_sample

        self.wl_spectralcalc = spectralcalc_data['wavelength']
        self.power_spectralcalc = spectralcalc_data['amplitude']
        self.pressure = spectralcalc_data['pressure']
        self.temperature = spectralcalc_data['temperature']
        self.concentration = spectralcalc_data['concentration']
        self.length = spectralcalc_data['length']

    def compute_fft(self):
        if any(x is None for x in (self.t, self.amplitude_sample, self.amplitude_reference)):
            self.read_data()

        from lib.fft import fft_plot_data

        self.fft_sample_x, self.fft_sample_y = fft_plot_data(
            self.t, self.amplitude_sample)
        self.fft_reference_x, self.fft_reference_y = fft_plot_data(
            self.t, self.amplitude_reference)

        return self.fft_sample_x, self.fft_sample_y, self.fft_reference_x, self.fft_reference_y

    def create_comb_analyser(self):
        if any(x is None for x in (self.fft_sample_x, self.fft_sample_y, self.fft_reference_x, self.fft_reference_y)):
            self.compute_fft()

        sample_analyzer = CombSpectrumAnalyser(
            self.fft_sample_x, self.fft_sample_y, self.center_freq, self.freq_spacing,
            self.number_of_teeth)
        sample_freq, sample_amp = sample_analyzer.extract_teeth()

        reference_analyzer = CombSpectrumAnalyser(
            self.fft_reference_x, self.fft_reference_y, self.center_freq, self.freq_spacing,
            self.number_of_teeth)
        ref_freq, ref_amp = reference_analyzer.extract_teeth()

        comb_analyser = CombAnalyser(
            sample_freq, sample_amp, ref_freq, ref_amp)
        comb_analyser.scaling_factor = self.scaling_factor

        self.comb_analyser = comb_analyser
        return self.comb_analyser

    def get_transmission_curves(self):
        if self.comb_analyser is None:
            self.create_comb_analyser()

        from lib.combs import (high_frequency_spectrum, to_frequency,
                               to_wavelength)
        from lib.fitting import overlap_data

        # Get the transmission curve
        f_sample, a_sample = self.comb_analyser.get_transmission_curve()

        # Transform the spectralcalc data to frequency
        f_spectralcalc, a_spectralcalc = to_frequency(
            self.wl_spectralcalc, self.power_spectralcalc)

        # Transform the low frequency comb spectrum to the high frequency comb spectrum
        f_sample, a_sample = high_frequency_spectrum(f_sample, a_sample, f0=self.center_freq,
                                                     fs=self.freq_spacing, fS=self.high_freq_modulation)

        # Shift the sample spectrum to overlap with the spectralcalc data as much as possible
        f_sample = overlap_data(
            f_sample, a_sample, f_spectralcalc, a_spectralcalc)

        # Transform everything back to wavelength
        self.wl_sample, self.power_sample = to_wavelength(f_sample, a_sample)
        self.wl_spectralcalc, self.power_spectralcalc = to_wavelength(
            f_spectralcalc, a_spectralcalc)

        if self.etalon_removal:
            from lib.fitting import try_remove_etalon
            wl_sample_corrected, self.power_sample_corrected, self.x_sine, self.y_sine = try_remove_etalon(
                self.wl_sample, self.power_sample, ignore_regions=[self.etalon_removal])
            self.wl_sample_corrected = overlap_data(wl_sample_corrected, self.power_sample_corrected,
                                                    self.wl_spectralcalc, self.power_spectralcalc)

        return self.wl_sample, self.power_sample, self.wl_spectralcalc, self.power_spectralcalc

    # Plots

    def time_series_plot(self):
        if any(x is None for x in (self.t, self.amplitude_sample, self.amplitude_reference)):
            self.read_data()

        import matplotlib.pyplot as plt

        # Subplot the time series of the signal that has travelled through the sample
        plt.subplot(2, 1, 1)
        plt.plot(self.t, self.amplitude_sample, 'g-', linewidth=2)
        plt.xlabel('Time (t)')
        plt.ylabel('Amplitude')
        plt.title('Time series of the sample')

        # Subplot the time series of the reference signal
        plt.subplot(2, 1, 2)
        plt.plot(self.t, self.amplitude_reference, 'b-', linewidth=2)
        plt.xlabel('Time (t)')
        plt.ylabel('Amplitude')
        plt.title('Time series of the reference')

        return plt

    def show_time_series_plot(self):
        self.time_series_plot().show()

    def fft_plot(self):
        if any(x is None for x in (self.fft_sample_x, self.fft_sample_y, self.fft_reference_x, self.fft_reference_y)):
            self.compute_fft()

        import matplotlib.pyplot as plt

        # Subplot the FFT of the signal that has travelled through the sample
        plt.subplot(2, 1, 1)
        plt.plot(self.fft_sample_x[1::],
                 self.fft_sample_y[1::], 'g-', linewidth=2)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.title('FFT of the sample')

        # Subplot the FFT of the reference signal
        plt.subplot(2, 1, 2)
        plt.plot(self.fft_reference_x[1::],
                 self.fft_reference_y[1::], 'b-', linewidth=2)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        plt.title('FFT of the reference')

        return plt

    def show_fft_plot(self):
        self.fft_plot().show()

    def transmission_plot(self, save_figure=False):
        if any(x is None for x in (self.wl_sample, self.power_sample, self.wl_spectralcalc, self.power_spectralcalc)):
            self.get_transmission_curves()

        import matplotlib.pyplot as plt

        # Plot the transmission curve
        plt.plot(self.wl_sample, self.power_sample, 'ro-', label='Experiment')

        if self.etalon_removal:
            plt.plot(self.x_sine, self.y_sine, 'yo-', label='Sine fit')
            plt.plot(self.wl_sample_corrected, self.power_sample_corrected,
                     'go-', label='Corrected Experiment')

        plt.plot(self.wl_spectralcalc, self.power_spectralcalc,
                 'b-', label='Spectralcalc Simulation')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Transmission (%)')
        title = f'Plot of the percentage transmission ({self.temperature}, {self.concentration})'
        plt.title(title)
        plt.legend()
        plt.xlim(self.wl_sample[0] - 0.05, self.wl_sample[-1] + 0.05)
        plt.gcf().set_size_inches(10, 6)

        if save_figure:
            import os

            sample = os.path.splitext(self.filename_sample)[0]
            reference = os.path.splitext(self.filename_reference)[0]

            filename = f"{sample} and {reference}"
            if self.etalon_removal:
                filename += " - ETALON Removed"
            figure_filename = f'../figures/{filename}.pdf'

            plt.savefig(figure_filename)
            f = os.path.dirname(os.path.realpath(__file__)) + figure_filename
            print(f"Saved figure to `{f}`")

        return plt

    def show_transmission_plot(self, save_figure=False):
        self.transmission_plot(save_figure=save_figure).show()