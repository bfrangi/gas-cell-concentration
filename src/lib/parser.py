def spectralcalc_fitting_parse_arguments():
    import argparse
    import sys

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
