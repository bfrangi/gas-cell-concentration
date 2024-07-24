def fft(x, y):
    from scipy.fft import fft as scipy_fft
    from scipy.fft import fftfreq
    N = len(x)  # Number of sample points
    T = x[1] - x[0]  # Sample spacing
    yf = scipy_fft(y)
    xf = fftfreq(N, T)[:N//2]
    return xf, yf

def fft_plot_data(x, y):
    import numpy as np
    xf, yf = fft(x, y)
    N = len(x)  # Number of sample points
    x_data = xf
    y_data = 2.0/N * np.abs(yf[0:N//2])
    return x_data, y_data
