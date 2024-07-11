####################################### Linear interpolation #######################################


def loose_interpolation(x, y, smoothness=100):
    """
    Use this function to obtain a loose interpolation of a set of data points.
    The interpolation will not necessarily go through all the data points, but
    will try to go as close to all as possible, while maintainig smoothness.

    Args:
        x (array-like): The x-values of the data points.
        y (array-like): The y-values of the data points.
        smoothness (int, optional): The smoothness of the interpolation. Defaults to 100.

    Example:
        ```
        import matplotlib.pyplot as plt

        x = [1, 2, 3, 4, 5, 6]
        y = [0.13, 0.27, 2.10, 0.52, 0.23, 0.13]

        plt.plot(x, y, 'o', label='Data points')
        x_plt, y_plt = loose_interpolation(x, y)
        plt.plot(x_plt, y_plt, '-', label='Loose interpolation')

        plt.show()
        ```
    """
    from numpy import linspace
    from csaps import csaps
    xs = linspace(x[0], x[-1], len(x)*smoothness, endpoint=True)
    return xs, csaps(x, y, xs, smooth=0.85)


######################################## Sinusoidal fitting ########################################


def sine_function(x, amplitude, frequency, phase, offset):
    '''
    Returns a sine function of the form:
    y = amplitude * sin(2 * pi * frequency * x + phase) + offset

    Args:
        x (array-like): The x-values.
        amplitude (float): The amplitude of the sine wave.
        frequency (float): The frequency of the sine wave.
        phase (float): The phase of the sine wave.
        offset (float): The offset of the sine wave.
    '''
    import numpy as np
    return amplitude * np.sin(2 * np.pi * frequency * x + phase) + offset


def fit_sine_to_data(x, y, initial_guess, ignore_regions=[], replace_with="NaN", nan_policy='omit'):
    '''
    Fits a sine function to the input data and returns the fitting parameters.

    Args:
        x (array-like): The x-values of the data.
        y (array-like): The y-values of the data.
        initial_guess (tuple): The initial guess for the fitting parameters. Should be of the form
            (amplitude, frequency, phase, offset).
        ignore_regions (list, optional): A list of tuples, each containing two values (start and 
            end) that define a region to ignore during fitting. Defaults to [].
        nan_policy (str, optional): How to handle NaN values. Can be 'omit', 'raise' or None.
            Defaults to 'omit'.
    '''
    from scipy.optimize import curve_fit
    import numpy as np
    y = y.copy()

    # Set unwanted x regions to NaN
    ignore_indices = []
    for region in ignore_regions:
        indices = np.where(np.logical_and(x > region[0], x < region[1]))
        ignore_indices.extend(indices[0])

    ignore_indices = np.array(ignore_indices)

    if replace_with == "NaN":
        y[ignore_indices] = np.nan
    elif replace_with == "mean":
        y[ignore_indices] = np.mean(y[~ignore_indices])
    else:
        y[ignore_indices] = replace_with

    # Define the model
    def model(x, amplitude, frequency, phase, offset):
        return sine_function(x, amplitude, frequency, phase, offset)

    # Perform curve fitting
    popt, pcov = curve_fit(model, x, y, p0=initial_guess,
                           nan_policy=nan_policy, maxfev=10000)

    return popt, pcov

########################################## Linear fitting ##########################################

def fit_line_to_data(x, y):
    """Find the slope and intercept of a line that fits the data.
    
    Args:
        x (array-like): The x-values of the data.
        y (array-like): The y-values of the data.
        
    Returns:
        tuple: A tuple containing the slope and the intercept of the line."""
    from scipy.optimize import curve_fit

    def f(x, A, B):
        return A*x + B

    popt, pcov = curve_fit(f, x, y)

    return popt[0], popt[1]


######################################## Offset calculation ########################################


def interpolate_onto_common_grid(x1, y1, x2, y2, step=1000j, remove_nan=True):
    """
    Interpolates the data sets y1 and y2 onto a common grid. This is useful when the data sets have
    different x-values, and we want to compare them.

    Args:
        x1 (array-like): The x-values of the first data set.
        y1 (array-like): The y-values of the first data set.
        x2 (array-like): The x-values of the second data set.
        y2 (array-like): The y-values of the second data set.
        step (complex, optional): The step size of the grid. Defaults to 1000j.
        remove_nan (bool, optional): Whether to remove NaN values from the interpolated data. 
            Defaults to True.

    Returns:
        tuple: A tuple containing the x-values and y-values of the interpolated data sets. Tuple
            is of the form (new_x1, new_y1, new_x2, new_y2).
    """
    import numpy as np
    from scipy.interpolate import griddata
    x1, y1, x2, y2 = np.array(x1), np.array(y1), np.array(x2), np.array(y2)

    uneven_grid_x = np.unique(np.concatenate((x1, x2)))
    grid_x = np.mgrid[uneven_grid_x[0]:uneven_grid_x[-1]:step]

    new_y1 = griddata(x1, y1, grid_x, method='cubic')
    new_y2 = griddata(x2, y2, grid_x, method='cubic')

    if not remove_nan:
        return grid_x, new_y1, grid_x, new_y2

    is_nan = np.isnan(new_y1)
    new_y1 = new_y1[~is_nan]
    new_x1 = grid_x[~is_nan]

    is_nan = np.isnan(new_y2)
    new_y2 = new_y2[~is_nan]
    new_x2 = grid_x[~is_nan]

    return new_x1, new_y1, new_x2, new_y2


def intersect_onto_common_grid(x1, y1, x2, y2, step=1000j):
    """
    Interpolates the data sets y1 and y2 onto a common grid only keeping the values where both data
    sets are not NaN.

    Args:
        x1 (array-like): The x-values of the first data set.
        y1 (array-like): The y-values of the first data set.
        x2 (array-like): The x-values of the second data set.
        y2 (array-like): The y-values of the second data set.
        step (complex, optional): The step size of the grid. Defaults to 1000j.

    Returns:
        tuple: A tuple containing the x-values and y-values of the interpolated data sets. Tuple
            is of the form (new_x1, new_y1, new_x2, new_y2).
    """
    import numpy as np
    new_x1, new_y1, new_x2, new_y2 = interpolate_onto_common_grid(x1, y1, x2, y2, step=step, remove_nan=False)

    is_nan_1 = np.isnan(new_y1)
    is_nan_2 = np.isnan(new_y2)

    is_nan = np.logical_or(is_nan_1, is_nan_2)

    new_y1 = new_y1[~is_nan]
    new_x1 = new_x1[~is_nan]
    new_y2 = new_y2[~is_nan]
    new_x2 = new_x2[~is_nan]

    return new_x1, new_y1, new_x2, new_y2


def mean_square_difference(y1, y2):
    """
    Computes the mean square difference between y1 and y2.

    Args:
        y1 (array): The first array.
        y2 (array): The second array.
    """
    if len(y1) != len(y2):
        raise ValueError('The two arrays must have the same length')
    import numpy as np
    return np.sum((y1 - y2)**2) / len(y1)


def get_offset(x1, y1, x2, y2, step=1000j):
    """
    Returns the horzontal shift between the data y1 and y2. This shift is taken as the distance that
    y1 should be shifted to the right in order for the two data sets to overlap as much as possible.
    Note that for this function to work, one of the data sets must be sufficiently shorter than the
    other, so that it can be compared to as many overlapping sections of the longer data set as
    possible.

    Args:
        x1 (array-like): The x-values of the first data set.
        y1 (array-like): The y-values of the first data set.
        x2 (array-like): The x-values of the second data set.
        y2 (array-like): The y-values of the second data set.
        step (complex, optional): The step size of the common grid built for the two data sets. 
            Defaults to 1000j.
    """
    x1, y1, x2, y2 = interpolate_onto_common_grid(x1, y1, x2, y2, step=step)
    step = x1[1] - x1[0]

    longest = y1 if len(y1) >= len(y2) else y2
    shortest = y2 if len(y1) >= len(y2) else y1

    def target(i):
        return mean_square_difference(longest[i:i+len(shortest)], shortest)

    min_distance = None
    min_i = 0
    for i in range(len(longest) - len(shortest)):
        res = target(i)
        if min_distance is None or res < min_distance:
            min_distance = res
            min_i = i

    result = step * min_i

    if len(y1) >= len(y2):
        return -x1[0] + x2[0] - result
    return -x1[0] + x2[0] + result


def overlap_data(x_to_fit, y_to_fit, x_reference, y_reference, step=1000j):
    """
    Overlaps the data y_to_fit with the reference data y_reference. This is done by shifting 
    y_to_fit horizontally so that it overlaps with y_reference as much as possible. The function
    returns the x-values of the shifted data y_to_fit (the y-values remain the same).

    Args:
        x_to_fit (array-like): The x-values of the data to be shifted.
        y_to_fit (array-like): The y-values of the data to be shifted.
        x_reference (array-like): The x-values of the reference data.
        y_reference (array-like): The y-values of the reference data.
        step (complex, optional): The step size of the common grid built for the two data sets. 
            Defaults to 1000j.

    Example:

        ```
        import matplotlib.pyplot as plt

        x1 = [1, 1.5, 2, 3, 4, 4.5, 5, 6, 7, 8, 9]
        y1 = [-10, -4, 1, 2, 3, 3.8, 4, 5, 4, 3, 2]
        x2 = [7, 8.5, 9, 10]
        y2 = [1, 2, 3, 4 ]

        x1 = overlap_data(x1, y1, x2, y2)

        plt.plot(x2, y2, label=f'shifted y2')
        plt.plot(x1, y1, label='y1')
        plt.legend()
        plt.show()
        ```
    """
    offset = get_offset(x_to_fit, y_to_fit, x_reference, y_reference, step=step)
    return x_to_fit + offset


########################################## Remove Etalon ###########################################

def remove_etalon(x, y, ignore_regions=[], initial_guess=[5, 5, 3.14, 37], verbose=False):
    """
    Removes the etalon signal from the input data by fitting a sine function to the data and
    subtracting it from the original data.
    
    Args:
        x (array-like): The x-values of the data.
        y (array-like): The y-values of the data.
        ignore_regions (list, optional): A list of tuples, each containing two values (start and 
            end) that define a region to ignore during fitting. Defaults to [].
        initial_guess (list, optional): The initial guess for the fitting parameters. Should be of 
            the form (amplitude, frequency, phase, offset). Defaults to [5, 5, 3.14, 37].
        verbose (bool, optional): Whether to print the fitting parameters. Defaults to False.
    """
    import numpy as np
    sine_parameters, covariance = fit_sine_to_data(
        x, y, initial_guess, replace_with="NaN", ignore_regions=ignore_regions)

    x_sine = x
    y_sine = sine_function(x_sine, *sine_parameters)

    if verbose:
        print("Sine fitting parameters:")
        print(f"Amplitude: {sine_parameters[0]}")
        print(f"Frequency: {sine_parameters[1]}")
        print(f"Phase: {sine_parameters[2]}")
        print(f"Offset: {sine_parameters[3]}")

    return x_sine, y_sine, x, y - (y_sine - np.max(y_sine))
