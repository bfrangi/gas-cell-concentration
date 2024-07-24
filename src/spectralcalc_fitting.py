from lib.parser import spectralcalc_fitting_parse_arguments
from lib.spectralcalc import SpectralCalcFitter

filename_reference, filename_sample, filename_spectralcalc, center_freq, freq_spacing, \
    number_of_teeth, high_freq_modulation, etalon_removal, scaling_factor, save_figure \
    = spectralcalc_fitting_parse_arguments()


fitter = SpectralCalcFitter(filename_sample, filename_reference, filename_spectralcalc,
                            scaling_factor=scaling_factor, center_freq=center_freq,
                            freq_spacing=freq_spacing, number_of_teeth=number_of_teeth,
                            high_freq_modulation=high_freq_modulation,
                            etalon_removal=etalon_removal)
if False:
    fitter.show_time_series_plot()

if False:
    fitter.show_fft_plot()

if False:
    caa = fitter.create_comb_analyser()
    caa.show_reference_spectrum()
    caa.show_sample_spectrum()
    caa.show_transmission_plot(interp=True)

fitter.show_transmission_plot(save_figure=save_figure)
