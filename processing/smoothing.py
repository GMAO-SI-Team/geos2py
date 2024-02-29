import numpy as np
import scipy.signal
from scipy.fft import fftshift, ifftshift, fftn, ifftn

def savitzky_golay2d(z: np.ndarray, window_size: int, order: int, derivative=None) -> np.ndarray:
    """
    A low pass filter for smoothing data implementing the Savitzky Golay algorithm
    """
    # number of terms in the polynomial expression
    n_terms = (order + 1) * (order + 2) / 2.0

    if window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size ** 2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [(k - n, n) for k in range(order + 1) for n in range(k + 1)]

    # coordinates of points
    ind = np.arange(-half_size, half_size + 1, dtype=np.float64)
    dx = np.repeat(ind, window_size)
    dy = np.tile(ind, [window_size, 1]).reshape(window_size ** 2, )

    # Build matrix of the system of equations
    A = np.empty((window_size ** 2, len(exps)))
    for i, exp in enumerate(exps):
        A[:, i] = (dx ** exp[0]) * (dy ** exp[1])

    # Pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2 * half_size, z.shape[1] + 2 * half_size
    Z = np.zeros((new_shape))
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs(np.flipud(z[1:half_size + 1, :]) - band)
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs(np.flipud(z[-half_size - 1:-1, :]) - band)
    # left band
    band = np.tile(z[:, 0].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs(np.fliplr(z[:, 1:half_size + 1]) - band)
    # right band
    band = np.tile(z[:, -1].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, -half_size:] = band + np.abs(np.fliplr(z[:, -half_size - 1:-1]) - band)
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs(np.flipud(np.fliplr(z[1:half_size + 1, 1:half_size + 1])) - band)
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs(np.flipud(np.fliplr(z[-half_size - 1:-1, -half_size - 1:-1])) - band)

    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs(np.flipud(Z[half_size + 1:2 * half_size + 1, -half_size:]) - band)
    # bottom left corner
    band = Z[-half_size:, half_size].reshape(-1, 1)
    Z[-half_size:, :half_size] = band - np.abs(np.fliplr(Z[-half_size:, half_size + 1:2 * half_size + 1]) - band)

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')

def dist(n: list, m=None) -> np.ndarray:
    """
    Return a rectangular array in which each pixel represents the Euclidean distance from the origin.

    Parameters:
    - n (list or tuple): A list or tuple containing a single integer representing the size of the array.
    - m (list or tuple, optional): A list or tuple containing a single integer representing the size of the array 
                                   along the second dimension. If not provided, defaults to the same value as n.

    Returns:
    - numpy.ndarray: A 2D array representing the Euclidean distances.

    Example:
    >>> n_val = 10
    >>> m_val = 6
    >>> result = dist([n_val], [m_val])
    >>> print(result)
    [[0.         1.          2.          3.          4.          5.          4.          3.          2.          1.        ]
     [1.         1.41421356  2.23606798  3.16227766  4.12310563  5.09901951  4.12310563  3.16227766  2.23606798  1.41421356]
     [2.         2.23606798  2.82842712  3.60555128  4.47213595  5.38516481  4.47213595  3.60555128  2.82842712  2.23606798]
     [3.         3.16227766  3.60555128  4.24264069  5.          5.83095189  5.          4.24264069  3.60555128  3.16227766]
     [2.         2.23606798  2.82842712  3.60555128  4.47213595  5.38516481  4.47213595  3.60555128  2.82842712  2.23606798]
     [1.         1.41421356  2.23606798  3.16227766  4.12310563  5.09901951  4.12310563  3.16227766  2.23606798  1.41421356]]
    """

    n1 = n[0]
    m1 = n1 if m is None or len(m) <= 0 else m[0]
    x = np.arange(n1)  # Make a row
    x = np.concatenate([x[x < (10 - x)], (10 - x)[(10 - x) <= x]]) ** 2 # Squared rows

    a = np.zeros((m1, n1))  # Make array

    for i in range(m1 // 2 + 1):  # Row loop
        y = np.sqrt(x + i**2)  # Euclidean distance
        a[i, :] = y  # Insert the row
        if i != 0:
            a[m1 - i, :] = y  # Symmetrical
    return a

def bandpass_filter(input_array: np.ndarray, low_freq: float, high_freq: float, ideal=False, butterworth=None, gaussian=False) -> np.ndarray:
    """
    This implementation is derived from bandpass_filter.pro in IDL lib subdirectory.

    Apply a bandpass filter to the input array in the frequency domain.

    Parameters:
    - input_array (ndarray): Input array to be filtered.
    - low_freq (float): Lower cutoff frequency of the bandpass filter (normalized frequency, [0, 1]).
    - high_freq (float): Upper cutoff frequency of the bandpass filter (normalized frequency, [0, 1]).
    - ideal (bool, optional): If True, use an ideal bandpass filter. Default is False.
    - butterworth (int, optional): If provided, use a Butterworth filter of the given order. Default is None.
    - gaussian (bool, optional): If True, use a Gaussian bandpass filter. Default is False.

    Returns:
    - ndarray: Filtered array.
    """

    # Handle all the flags
    if sum([ideal, butterworth is not None, gaussian]) != 1:
        raise ValueError('Only set one of IDEAL, BUTTERWORTH, or GAUSSIAN.')
    
    if len(input_array.shape) != 2:
        raise ValueError('Input must be a two-dimensional array.')
    
    if butterworth is None:
        butterworth = 1
    
    if low_freq is None or high_freq is None:
        raise ValueError('Flow and Fhigh must be supplied.')
    
    if low_freq > 1 or low_freq < 0 or high_freq > 1 or high_freq < 0:
        raise ValueError('Frequency is out of range ([0,1]).')
    
    if high_freq < low_freq:
        raise ValueError('Fhigh must be greater than Flow.')
    
    if butterworth <= 0:
        raise ValueError('Butterworth dimension must be a positive value.')
    
    # Perform Fourier Transform
    fourier_transform = fftn(input_array, axes=(0, 1))
    fourier_transform = fftshift(fourier_transform)
    
    # Compute distance matrix
    dimensions = np.array(input_array.shape)
    dist_matrix = dist([dimensions[1]], [dimensions[0]])
    dist_matrix = np.roll(dist_matrix, dimensions[0] // 2 + 1, axis=1)
    dist_matrix = np.roll(dist_matrix, dimensions[1] // 2 + 1, axis=0)
    dist_matrix /= np.max(dist_matrix)
    
    # Compute filter function based on flags
    D = dist_matrix
    W = high_freq - low_freq
    D0 = (high_freq + low_freq) / 2.0
    
    if ideal:
        H = np.where((low_freq > D) | (high_freq < D), 0, 1)
    elif butterworth is not None:
        if low_freq != 0 and high_freq != 1:
            H = 1.0 - 1.0 / (1 + ((D * W) / (D ** 2 - D0 ** 2)) ** (2 * butterworth))
        elif low_freq == 0:
            H = 1.0 / (1 + (D / high_freq) ** (2 * butterworth))
        else:
            H = 1.0 / (1 + (low_freq / D) ** (2 * butterworth))
    else:  # gaussian
        if low_freq != 0 and high_freq != 1:
            H = np.exp(-(((D ** 2 - D0 ** 2) / (D * W)) ** 2))
        elif low_freq == 0:
            H = np.exp(-(D ** 2 / (2 * high_freq ** 2)))
        else:
            H = 1.0 - np.exp(-(D ** 2 / (2 * low_freq ** 2)))
    
    result_fourier = H * fourier_transform
    
    # Inverse Fourier Transform
    result = ifftn(ifftshift(result_fourier), axes=(0, 1)).real
    return result