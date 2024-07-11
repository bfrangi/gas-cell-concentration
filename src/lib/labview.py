import pandas as pd
import os 
import numpy as np

class CSVReader:
    def __init__(self, filename):
        self.filename = filename
        script_dir = os.path.dirname(__file__)  # Script directory
        self.full_path = os.path.join(script_dir, f'../../measurements/{self.filename}')
        self.raw_data = None
        self.data = None

    def read(self):
        self.raw_data = pd.read_csv(self.full_path)
        return self.raw_data
    
    def clean_data(self):
        data = self.raw_data if self.raw_data is not None else self.read()
        data.replace(r'u', r'e-6', regex=True, inplace=True)
        data.replace(r'm', r'e-3', regex=True, inplace=True)
        data.replace(r',', r'.', regex=True, inplace=True)
        data = data.astype(float)
        self.raw_data = data
        return self.raw_data
    
    def extract_data(self):
        self.read()
        data = self.clean_data()

        time = data['Time - Voltage'].to_numpy()
        amplitude = data['Amplitude - Voltage'].to_numpy()

        time_delta = time[1] - time[0]
        time = time_delta * np.arange(len(time))

        data = {'time': time, 'amplitude': amplitude}
        self.data = data
        return self.data

class LVMReader:
    def __init__(self, filename):
        self.filename = filename
        script_dir = os.path.dirname(__file__)  # Script directory
        self.full_path = os.path.join(script_dir, f'../../measurements/{self.filename}')
        self.acq_freq = None
        self.raw_data = None
        self.data = None

    def read(self):
        with open(self.full_path, 'r') as f:
            lines = f.readlines()
            self.acq_freq = float(lines[0].split('=')[1].strip())
            if not self.acq_freq:
                raise ValueError(f'Acquisition frequency not found in file {self.filename}. Must be specified in the first line as `f=<value in Hz>`.')
            self.raw_data = lines[1::]
        return self.raw_data
    
    def clean_data(self):
        data = self.raw_data if self.raw_data is not None else self.read()
        data = [d.replace(',', '.') for d in data]
        data = [float(d.strip()) for d in data]
        self.raw_data = data
        return self.raw_data
    
    def extract_data(self):
        self.read()
        self.clean_data()

        time = np.linspace(0, len(self.raw_data) / self.acq_freq, len(self.raw_data))
        amplitude = self.raw_data

        time_delta = time[1] - time[0]
        time = time_delta * np.arange(len(time))

        data = {'time': time, 'amplitude': amplitude}
        self.data = data
        return self.data

class Reader:
    def __init__(self, filename):
        import sys, os
        name, extension = os.path.splitext(filename)
        self.reader = None
        if extension == '.csv':
            self.reader = CSVReader(filename=filename)
        elif extension == '.lvm':
            self.reader = LVMReader(filename=filename)
        else:
            print(f'File extension `{extension}` not supported')
            sys.exit(1)
    
    def extract_data(self):
        return self.reader.extract_data()
