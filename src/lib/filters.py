########################## Filter Functions ##########################

# Bandpass Filter


def bandpass_params(cutoff, fs, order=5):
    from scipy.signal import butter
    return butter(order, cutoff, fs=fs, btype='bandpass', analog=False)


def bandpass_filter(data, cutoff, fs, order=5):
    from scipy.signal import lfilter
    b, a = bandpass_params(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


# Lowpass Filter


def lowpass_params(cutoff, fs, order=5):
    from scipy.signal import butter
    nyq = 0.5 * fs  # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    return butter(order, normal_cutoff, btype='low', analog=False)


def lowpass_filter(data, cutoff, fs, order=5):
    from scipy.signal import filtfilt
    b, a = lowpass_params(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

####################### Filter Response Plots ########################

# Lowpass Filter Parameters


def lowpass_response_plot(cutoff, order=5, fs=240.0):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.signal import freqz
    b, a = lowpass_params(cutoff, fs, order)   # Filter coefficients
    w, h = freqz(b, a, fs=fs, worN=8000)
    plt.plot(w, np.abs(h), 'b')
    plt.plot([cutoff], [0.5*np.sqrt(2)], 'ko')
    plt.axvline(cutoff, color='k')
    plt.xlim(0, 0.5*fs)
    plt.title("Lowpass Filter Frequency Response")
    plt.xlabel('Frequency [Hz]')
    plt.grid()
    return plt

######################## Threshold and Limit #########################

# Threshold Filter

def threshold_filter(arr, threshold, *other_arrs):
    import numpy as np
    arr = np.array(arr)
    indices_to_delete = np.where(arr < threshold)
    arr =  np.delete(arr, indices_to_delete)
    other_arrs = [np.delete(other_arr, indices_to_delete) for other_arr in other_arrs]
    return arr, *other_arrs

# Limit Filter

def limit_filter(arr, limit, *other_arrs):
    import numpy as np
    arr = np.array(arr)
    indices_to_delete = np.where(arr > limit)
    arr =  np.delete(arr, indices_to_delete)
    other_arrs = [np.delete(other_arr, indices_to_delete) for other_arr in other_arrs]
    return arr, *other_arrs