import os
import numpy as np
from .fitting import interpolate_onto_common_grid


class SpectralcalcReader:
    def __init__(self, filename):
        self.filename = filename
        script_dir = os.path.dirname(__file__)  # Script directory
        self.full_path = os.path.join(script_dir, f'../../spectralcalc-simulations/{self.filename}')
        self.raw_data = None
        self.data = None
        self.pressure = None
        self.temperature = None
        self.concentration = None
        self.length = None

    def extract_metadata(self, metadata):
        import re
        for m in re.finditer(r'Wavelength \(((?:(?!microns).)+)\)', metadata):
            raise ValueError(f'Wavelength units must be microns: `{m.group(1)}` found.')
        length_match = re.search(r'Length \(([^\)]*)\) = ([^\n]*)', metadata)
        self.length = length_match.group(2) + ' ' + length_match.group(1)
        pressure_match = re.search(r'Pressure \(([^\)]*)\) = ([^\n]*)', metadata)
        self.pressure = pressure_match.group(2) + ' ' + pressure_match.group(1)
        temperature_match = re.search(r'Temperature \(([^\)]*)\) = ([^\n]*)', metadata)
        self.temperature = temperature_match.group(2) + ' ' + temperature_match.group(1)
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
        temp =  float(temperature_str.replace("K", "").strip())
        if temp%1 == 0:
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

    wl_lower, mac_lower = real_molar_attenuation_coefficient(lower_reference_filename)
    wl_higher, mac_higher = real_molar_attenuation_coefficient(higher_reference_filename)

    wl_lower, mac_lower, wl_higher, mac_higher = interpolate_onto_common_grid(
        wl_lower, mac_lower, wl_higher, mac_higher, remove_nan=False)

    mac = mac_lower + (mac_higher - mac_lower) * (concentration - 0.022) / (1 - 0.022)

    mask = np.isnan(mac)
    return wl_lower[~mask], mac[~mask]


def interpolated_transmission_curve(concentration, temperature=296, length=0.1):
    wl, mac = interpolated_molar_attenuation_coefficient(concentration, temperature)
    return wl, 10**(-mac*length*concentration)
