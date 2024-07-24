import os
import numpy as np
from .fitting import interpolate_onto_common_grid


class SpectralcalcReader:
    def __init__(self, filename):
        self.filename = filename
        script_dir = os.path.dirname(__file__)  # Script directory
        self.full_path = os.path.join(
            script_dir, f'../../spectralcalc-simulations/{self.filename}')
        self.raw_data = None
        self.data = None
        self.pressure = None
        self.temperature = None
        self.concentration = None
        self.length = None

    def extract_metadata(self, metadata):
        import re
        for m in re.finditer(r'Wavelength \(((?:(?!microns).)+)\)', metadata):
            raise ValueError(
                f'Wavelength units must be microns: `{m.group(1)}` found.')
        length_match = re.search(r'Length \(([^\)]*)\) = ([^\n]*)', metadata)
        self.length = length_match.group(2) + ' ' + length_match.group(1)
        pressure_match = re.search(
            r'Pressure \(([^\)]*)\) = ([^\n]*)', metadata)
        self.pressure = pressure_match.group(2) + ' ' + pressure_match.group(1)
        temperature_match = re.search(
            r'Temperature \(([^\)]*)\) = ([^\n]*)', metadata)
        self.temperature = temperature_match.group(
            2) + ' ' + temperature_match.group(1)
        concentration_match = re.search(r'\s([\.0-9]+)\s+\n', metadata)
        self.concentration = concentration_match.group(1) + ' VMR'

    def read(self):
        with open(self.full_path, 'r') as f:
            metadata = ''
            line = "#"
            while line.startswith("#"):
                line = f.readline()
                if line.startswith("#"):
                    metadata += line
            self.extract_metadata(metadata)
            lines = [line] + f.readlines()
            self.raw_data = lines[1::]
        return self.raw_data

    def clean_data(self):
        data = self.raw_data if self.raw_data is not None else self.read()
        wavelength_data = []
        amplitude_data = []
        for line in data:
            wl, amp = line.strip().split()
            wavelength_data.append(float(wl)*1e3)
            amplitude_data.append(float(amp))

        self.raw_data = {
            'wavelength': wavelength_data,
            'amplitude': amplitude_data,
            'pressure': self.pressure,
            'temperature': self.temperature,
            'concentration': self.concentration,
            'length': self.length
        }
        return self.raw_data

    def extract_data(self):
        self.read()
        return self.clean_data()


class SpectralCalcFitter:

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
        from lib.labview import Reader

        reader_sample = Reader(self.filename_sample)
        sample_data = reader_sample.extract_data()

        reader_reference = Reader(self.filename_reference)
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

    def compute_fft(self):
        if any(x is None for x in (self.t, self.amplitude_sample, self.amplitude_reference)):
            self.read_data()

        from lib.fft import fft_plot_data

        self.fft_sample_x, self.fft_sample_y = fft_plot_data(
            self.t, self.amplitude_sample)
        self.fft_reference_x, self.fft_reference_y = fft_plot_data(
            self.t, self.amplitude_reference)

        return self.fft_sample_x, self.fft_sample_y, self.fft_reference_x, self.fft_reference_y

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

    def create_comb_analyser(self):
        if any(x is None for x in (self.fft_sample_x, self.fft_sample_y, self.fft_reference_x, self.fft_reference_y)):
            self.compute_fft()

        from lib.analysis import CombAnalyser, CombSpectrumAnalyser

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

        from lib.combs import high_frequency_spectrum, to_frequency, to_wavelength
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


def get_length(length_str):
    if "cm" in length_str:
        return float(length_str.replace("cm", "").strip()) * 1e-2
    raise ValueError(f"Length unit not recognized: {length_str}")


def get_concentration(concentration_str):
    if "VMR" in concentration_str:
        return float(concentration_str.replace("VMR", "").strip())
    raise ValueError(f"Concentration unit not recognized: {concentration_str}")


def get_temperature(temperature_str):
    if "K" in temperature_str:
        temp = float(temperature_str.replace("K", "").strip())
        if temp % 1 == 0:
            return int(temp)
        return temp
    raise ValueError(f"Temperature unit not recognized: {temperature_str}")


def real_molar_attenuation_coefficient(filename):
    reader = SpectralcalcReader(filename)
    data = reader.extract_data()
    wl = data['wavelength']
    amplitude = np.array(data['amplitude'])
    length = get_length(data['length'])
    concentration = get_concentration(data['concentration'])
    return wl, -np.log10(amplitude) / (length*concentration)


def interpolated_molar_attenuation_coefficient(concentration, temperature=296):
    lower_reference_filename = f"CH4_{temperature}K_0.022VMR.txt"
    higher_reference_filename = f"CH4_{temperature}K_1VMR.txt"

    wl_lower, mac_lower = real_molar_attenuation_coefficient(
        lower_reference_filename)
    wl_higher, mac_higher = real_molar_attenuation_coefficient(
        higher_reference_filename)

    wl_lower, mac_lower, wl_higher, mac_higher = interpolate_onto_common_grid(
        wl_lower, mac_lower, wl_higher, mac_higher, remove_nan=False)

    mac = mac_lower + (mac_higher - mac_lower) * \
        (concentration - 0.022) / (1 - 0.022)

    mask = np.isnan(mac)
    return wl_lower[~mask], mac[~mask]


def interpolated_transmission_curve(concentration, temperature=296, length=0.1):
    wl, mac = interpolated_molar_attenuation_coefficient(
        concentration, temperature)
    return wl, 10**(-mac*length*concentration)
