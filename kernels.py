import numpy as np

def boxcar(x, h):
    """
    Univriate boxcar kernel smoothing function.

    Returns 1 if x is in (-h, h), zero otherwise.

    Parameters:
        - x (array-like): A value or array of values to apply the smoothing 
        kernel to. 
        - h (float): Kernel bandwidth. 

    Returns:
        (array-like): A value or array matching the size of x containing the 
        kernel applied to x. 
    """
    x = np.divide(x, h)
    return np.where(abs(x) < 1, 1, 0)


def epanechnikov(x, h):
    """
    Univriate epanechnikov kernel smoothing function.

    Returns (3/4)*(1 - x^2) if x is in (-h, h), zero otherwise. 

    Parameters:
        - x (array-like): A value or array of values to apply the smoothing 
        kernel to. 
        - h (float): Kernel bandwidth. 

    Returns:
        (array-like): A value or array matching the size of x containing the 
        kernel applied to x.
    """
    x = np.divide(x, h)
    return np.where(abs(x) < 1, 0.75*(1 - np.power(x, 2)), 0)