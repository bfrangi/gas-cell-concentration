import sys

import matplotlib.pyplot as plt

from lib.analysis import CombAnalyser, CombSpectrumAnalyser
from lib.combs import high_frequency_spectrum, to_frequency, to_wavelength
from lib.fft import fft_plot_data
from lib.fitting import overlap_data
from lib.labview import Reader
from lib.spectralcalc import SpectralcalcReader


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filename_sample', nargs='?', default=None)
    parser.add_argument('filename_reference', nargs='?', default=None)
    parser.add_argument('filename_spectralcalc', nargs='?', default=None)
    parser.add_argument('--center-freq', '-cf', 
                        help='Center frequency in the low-frequency domain comb.', type=float, default=40000)
    parser.add_argument('--freq-spacing', '-fs', 
                        help='Frequency spacing in the low-frequency domain comb.', type=float, default=200)
    parser.add_argument('--number-of-teeth', '-not', 
                        help="Number of teeth in the comb.", type=int, default=38)
    parser.add_argument('--high-freq-modulation', '-hfm', 
                        help='Frequency spacing in the high-frequency domain comb.', type=float, default=1e9)
    parser.add_argument('--scaling-factor', '-sf', 
                        help='Amplitude scale factor applied to the measured spectrum.', type=float, default=1.0)
    parser.add_argument('--remove-etalon', '-re', 
                        help='Perform a sine fitting to remove ETALON effects.', nargs=2, type=float, default=())
    parser.add_argument('--save-figure', '-s', 
                        help='Save the figure to PDF.', action=argparse.BooleanOptionalAction)

    # Parse args
    args = sys.argv[1::]
    args = parser.parse_args(args)
    filename_sample = args.filename_sample
    filename_reference = args.filename_reference
    filename_spectralcalc = args.filename_spectralcalc
    center_freq = args.center_freq
    freq_spacing = args.freq_spacing
    number_of_teeth = args.number_of_teeth
    high_freq_modulation = args.high_freq_modulation
    etalon_removal = tuple(args.remove_etalon)
    scaling_factor = args.scaling_factor
    save_figure = args.save_figure

    if not filename_sample or not filename_reference or not filename_spectralcalc:
        parser.print_help()
        exit()
    return filename_reference, filename_sample, filename_spectralcalc, center_freq, \
        freq_spacing, number_of_teeth, high_freq_modulation, etalon_removal, scaling_factor, save_figure


filename_reference, filename_sample, filename_spectralcalc, center_freq, freq_spacing, \
    number_of_teeth, high_freq_modulation, etalon_removal, scaling_factor, save_figure = parse_arguments()


########################################### Reading Data ###########################################

# Read data from files
reader_sample = Reader(filename_sample)
sample_data = reader_sample.extract_data()

reader_reference = Reader(filename_reference)
reference_data = reader_reference.extract_data()

reader_spectralcalc = SpectralcalcReader(filename_spectralcalc)
spectralcalc_data = reader_spectralcalc.extract_data()


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

wl_spectralcalc = spectralcalc_data['wavelength']
amplitude_spectralcalc = spectralcalc_data['amplitude']
pressure = spectralcalc_data['pressure']
temperature = spectralcalc_data['temperature']
concentration = spectralcalc_data['concentration']
length = spectralcalc_data['length']

############################ Computing FFT of Sample and Reference Data ############################

x_sample, y_sample = fft_plot_data(t, amplitude_sample)
x_reference, y_reference = fft_plot_data(t, amplitude_reference)

if False:
    # Subplot the FFT of the signal that has travelled through the sample
    plt.subplot(2, 1, 1)
    plt.plot(x_sample[1::], y_sample[1::], 'g-', linewidth=2)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Plot of the FFT of the signal that has travelled through the sample')

    # Subplot the FFT of the reference signal
    plt.subplot(2, 1, 2)
    plt.plot(x_reference[1::], y_reference[1::], 'b-', linewidth=2)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Plot of the FFT of the reference signal')

    # Show the plots
    plt.show()


####################################### Get absorption data ########################################

# Analyse the comb spectrum

sample_analyzer = CombSpectrumAnalyser(x_sample, y_sample, center_freq, freq_spacing, number_of_teeth)
sample_freq, sample_amp = sample_analyzer.extract_teeth()

reference_analyzer = CombSpectrumAnalyser(x_reference, y_reference, center_freq, freq_spacing, number_of_teeth)
ref_freq, ref_amp = reference_analyzer.extract_teeth()

caa = CombAnalyser(sample_freq, sample_amp, ref_freq, ref_amp)
caa.scaling_factor = scaling_factor

if False:
    caa.show_reference_spectrum()
    caa.show_sample_spectrum()
    caa.show_transmission_plot(interp=True)

################################# Plot alongside Spectralcalc data #################################

# Get the transmission curve and interpolate it to get a smooth curve
f_sample, a_sample = caa.get_transmission_curve()
# f_sample, a_sample = loose_interpolation(f_sample, a_sample)

# Transform the spectralcalc data to frequency
f_spectralcalc, a_spectralcalc = to_frequency(wl_spectralcalc, amplitude_spectralcalc)

# Transform the low frequency comb spectrum to the high frequency comb spectrum
f_sample, a_sample = high_frequency_spectrum(f_sample, a_sample, f0=center_freq, fs=freq_spacing, 
                                       fS=high_freq_modulation)

# Shift the sample spectrum to overlap with the spectralcalc data as much as possible
f_sample = overlap_data(f_sample, a_sample, f_spectralcalc, a_spectralcalc)

# Transform everything back to wavelength
wl_sample, a_sample = to_wavelength(f_sample, a_sample)
wl_spectralcalc, a_spectralcalc = to_wavelength(f_spectralcalc, a_spectralcalc)

# Plot the transmission curve
plt.plot(wl_sample, a_sample, 'ro-', label='Experiment')

if etalon_removal:
    from lib.fitting import remove_etalon

    ignore_regions = [etalon_removal]
    initial_guess = [1.5, 5, 3.14, 37]

    error = True
    while error:
        error = False
        try:
            x_sine, y_sine, wl_sample_corrected, a_sample_corrected = remove_etalon(
                wl_sample, a_sample, ignore_regions=ignore_regions, 
                initial_guess=initial_guess, verbose=True)
        except Exception as e:
            if initial_guess[2] < 9.42:
                error = True
            else:
                raise e
        finally:
            if error:
                print("Error in fitting, trying again with different initial guess")
                initial_guess[2] += 0.05
    
    wl_sample_corrected = overlap_data(wl_sample_corrected, a_sample_corrected, wl_spectralcalc, a_spectralcalc)

    plt.plot(x_sine, y_sine, 'yo-', label='Sine fit')
    plt.plot(wl_sample_corrected, a_sample_corrected, 'go-', label='Corrected Experiment')

plt.plot(wl_spectralcalc, a_spectralcalc, 'b-', label='Spectralcalc Simulation')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (%)')
title = f'Plot of the percentage transmission ({temperature}, {concentration})'
plt.title(title)
plt.legend()
plt.xlim(wl_sample[0] - 0.05, wl_sample[-1] + 0.05)
plt.gcf().set_size_inches(10, 6)

if save_figure:
    import os

    sample = os.path.splitext(filename_sample)[0]
    reference = os.path.splitext(filename_reference)[0]

    filename = f"{sample} and {reference}"
    if etalon_removal:
        filename += " - ETALON Removed"
    figure_filename = f'../figures/{filename}.pdf'

    plt.savefig(figure_filename)
    f = os.path.dirname(os.path.realpath(__file__)) + figure_filename
    print(f"Saved figure to `{f}`")

plt.show()
