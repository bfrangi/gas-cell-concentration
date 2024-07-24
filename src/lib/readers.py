class LVMReader:
    """Reads a LabView Measurement file (.lvm) and extracts the data."""
    def __init__(self, filename):
        import os

        self.filename = filename
        script_dir = os.path.dirname(__file__)  # Script directory
        self.full_path = os.path.join(
            script_dir, f'../../measurements/{self.filename}')
        self.acq_freq = None
        self.raw_data = None
        self.data = None

    def read(self):
        with open(self.full_path, 'r') as f:
            lines = f.readlines()
            self.acq_freq = float(lines[0].split('=')[1].strip())
            if not self.acq_freq:
                raise ValueError(
                    f'Acquisition frequency not found in file {self.filename}. Must be specified in the first line as `f=<value in Hz>`.')
            self.raw_data = lines[1::]
        return self.raw_data

    def clean_data(self):
        data = self.raw_data if self.raw_data is not None else self.read()
        data = [d.replace(',', '.') for d in data]
        data = [float(d.strip()) for d in data]
        self.raw_data = data
        return self.raw_data

    def extract_data(self):
        import numpy as np

        self.read()
        self.clean_data()

        time = np.linspace(0, len(self.raw_data) /
                           self.acq_freq, len(self.raw_data))
        amplitude = self.raw_data

        time_delta = time[1] - time[0]
        time = time_delta * np.arange(len(time))

        data = {'time': time, 'amplitude': amplitude}
        self.data = data
        return self.data


class SpectralcalcReader:
    """Reads a Spectralcalc Simulation file and extracts the data."""
    def __init__(self, filename):
        import os

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
