import numpy as np
import scipy.signal
import scipy.ndimage
import scipy.interpolate
from scipy.fft import fft2, ifft2, fftshift, ifftshift

def savitzky_golay2d(z, window_size, order, derivative=None):
    """
    A low pass filter for smoothing data
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
    
# def ideal_bandpass_filter(signal, lowcut, highcut):
#     fft_signal = np.fft.fft(signal)
#     frequencies = np.fft.fftfreq(len(signal))

#     # Create an ideal bandpass filter
#     ideal_filter = np.zeros_like(frequencies)
#     ideal_filter[(frequencies >= lowcut) & (frequencies <= highcut)] = 1

#     # Reshape filter to column vector for broadcasting
#     ideal_filter = ideal_filter.reshape(-1, 1)

#     # Apply the filter to the frequency domain representation of the signal
#     filtered_signal = fft_signal * ideal_filter

#     # Inverse Fourier transform to get back to the time domain
#     filtered_signal = np.fft.ifft(filtered_signal)

#     return np.real(filtered_signal)

def ideal_bandpass_filter(img, low_cutoff, high_cutoff):
    rows, cols = img.shape
    crow, ccol = rows // 2, cols // 2

    # Ensure integer values for slicing
    crow, ccol = int(crow), int(ccol)

    # Fourier transform
    f_transform = fft2(img)
    f_transform_shifted = fftshift(f_transform)

    # Ideal bandpass filter
    mask = np.zeros_like(img)
    
    if isinstance(high_cutoff, int):
        high_cutoff = int(high_cutoff)
    elif isinstance(high_cutoff, float):
        high_cutoff = int(high_cutoff * min(crow, ccol))  # Scale based on image dimensions

    mask[crow - high_cutoff:crow + high_cutoff, ccol - high_cutoff:ccol + high_cutoff] = 1
    mask[crow - int(low_cutoff):crow + int(low_cutoff), ccol - int(low_cutoff):ccol + int(low_cutoff)] = 0

    # Apply the filter
    f_transform_shifted_filtered = f_transform_shifted * mask

    # Inverse Fourier transform
    img_filtered = ifft2(ifftshift(f_transform_shifted_filtered)).real

    return img_filtered