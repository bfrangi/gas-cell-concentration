# Dual-Comb Spectroscopy Measurements of CH4 Concentration in a Glass Cell

This repository contains the code to process the data obtained from a dual-comb spectroscopy
measurement of CH4 concentration in a glass cell. 

## Setup

To install all requirements, make sure you have Python 3.6 or later installed and create a virtual
environment to work in (this last step is optional).

## Usage

Run the script with:
```bash
pip install -r requirements.txt
python3 spectralcalc_fitting.py filename_sample filename_reference filename_spectralcalc [-h] [--center-freq CENTER_FREQ] [--freq-spacing FREQ_SPACING] [--number-of-teeth NUMBER_OF_TEETH] [--high-freq-modulation HIGH_FREQ_MODULATION] [--scaling-factor SCALING_FACTOR] [--remove-etalon REMOVE_ETALON REMOVE_ETALON] [--save-figure | --no-save-figure | -s] 
```

The three positional arguments are the filenames of the sample (`filename_sample`), reference 
(`filename_reference`), and spectralcalc files (`filename_spectralcalc`), respectively. The sample
and reference files are the measured spectra, while the spectralcalc file is the calculated 
spectrum, obtained from the Hitran database.

The optional arguments are:
* `--center-freq` (`-cf`): Center frequency in the low-frequency domain comb. Default is `40000`.
* `--freq-spacing` (`-fs`): Frequency spacing in the low-frequency domain comb. Default is `200`.
* `--number-of-teeth` (`-not`): Number of teeth in the comb. Default is `38`.
* `--high-freq-modulation` (`-hfm`): Frequency spacing in the high-frequency domain comb. Default is
  `1e9`.
* `--scaling-factor` (`-sf`): Amplitude scale factor applied to the measured spectrum. Default is 
  `1.0`.
* `--remove-etalon` (`-re`): Tuple with the start and end wavelengths in nm that make up the range
  of the spectrum where the absorption lines are present, so that they can be ignored when 
  performing the sine fitting to remove ETALON effects. Default is `()` (no ETALON removal).
* `--save-figure` (`-s`): Save the figure to PDF when the option is given. Default is `False`.

## Data format

The measurement data files should be in the following format, with a `.lvm` file extension:

```
f=400000
	0,039023
	0,036448
	0,036448
	0,037735
```

The first line is the sampling frequency, and the following lines are the measured intensities.

The spectralcalc file should be in the default `.txt` format provided by the Spectral Calculator 
from G & A Technical Software, Inc. The following is an example:

```
# G & A Technical Software, Inc.  (GATS) Spectral Calculator
# http://www.spectralcalc.com
# Date: Fri Feb  2 11:31:09 2024
# Line list: hitran2020
#
# Cell Number: 1
# Pressure (mb) = 1013.25
# Temperature (K) = 288
# Length (cm) = 10
# Gas               Isotope Id            Volume Mixing Ratio   
#  CH4              0  All Isotopes       0.23                 
#
# Bandpass Starting Wavelength (microns) = 1.64500101564336
# Bandpass Ending Wavelength (microns) = 1.64599898312268
# Number of Spectral Points = 492
# Mean Transmittance = 0.936559267036084
#
#   Wavelength (microns)   Transmittance
1.645001015643356e+00 9.973234640772504e-01
1.645003046933873e+00 9.973718428892576e-01
1.645005078229408e+00 9.974192038622194e-01
1.645007109529959e+00 9.974627353508370e-01
1.645009140835526e+00 9.975007223918826e-01
1.645011172146111e+00 9.975322882814269e-01
```


## Examples

There are several measurements provided in this repository. The correct configuration of the script 
for each of these measurements is:

1. `sample_1.lvm`, `reference_1.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_1.lvm reference_1.lvm CH4_296K_0.15VMR.txt -cf 40000 -fs 200 -not 35 -hfm 1000000000.0 -sf 1.04
```

2. `sample_1.lvm`, `reference_1.lvm` with ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_1.lvm reference_1.lvm CH4_296K_0.08VMR.txt -cf 40000 -fs 200 -not 35 -hfm 1000000000.0 -sf 1.03 -re 1645.44 1645.65
```

3. `sample_2.lvm`, `reference_2.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_2.lvm reference_2.lvm CH4_296K_0.31VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 1.11
```

4. `sample_2.lvm`, `reference_2.lvm` with ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_2.lvm reference_2.lvm CH4_296K_0.28VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 1.075 -re 1645.39 1645.7
```

5. `sample_3.lvm`, `reference_3.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_3.lvm reference_3.lvm CH4_296K_0.29VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 1.1
```

6. `sample_3.lvm`, `reference_3.lvm` with ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_3.lvm reference_3.lvm CH4_296K_0.32VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 0.95 -re 1645.45 1645.75
```

7. `sample_4.lvm`, `reference_4.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_4.lvm reference_4.lvm CH4_296K_0.29VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 1.1
```

8. `sample_4.lvm`, `reference_4.lvm` with ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_4.lvm reference_4.lvm CH4_296K_0.3VMR.txt -cf 40000 -fs 100 -not 34 -hfm 1500000000.0 -sf 0.95 -re 1645.48 1645.75
```

9. `sample_5.lvm`, `reference_5.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_5.lvm reference_5.lvm CH4_296K_0.04VMR.txt -cf 40000 -fs 200 -not 38 -hfm 1500000000.0 -sf 0.96
```

10. `sample_6.lvm`, `reference_6.lvm` without ETALON removal:
```bash
python3 spectralcalc_fitting.py sample_6.lvm reference_6.lvm CH4_296K_0.043VMR.txt -cf 40000 -fs 200 -not 36 -hfm 1500000000.0 -sf 0.98
```