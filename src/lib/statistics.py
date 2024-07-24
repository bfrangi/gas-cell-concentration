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