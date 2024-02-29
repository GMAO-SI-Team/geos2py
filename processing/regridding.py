import subprocess
import numpy as np

from scipy.io import FortranFile
from sunpy.image import resample
from numba import jit, float32, float64, int64, int32

def congrid(data: np.ndarray, shape: tuple, center=True, method='linear') -> np.ndarray:
    """
    Resamples data to a specified shape.

    Parameters:
    - data (np.ndarray): Input data to be resampled.
    - shape (tuple): Desired shape of the output data.
    - center (bool, optional): If True, the coordinates are corrected to represent the center of the bins. Defaults to True.
    - method (str, optional): Method of interpolation. Defaults to 'linear'. Supports
    """
    return resample.resample(data, shape, center=center, method=method)

def regrid(data: np.ndarray, method='conservative', gridspec=None, undef=1e15) -> np.ndarray:
    """
    Parameters:
    - method (str): Method of regridding. Supported values are 'conservative' and 'bilinear'. Defaults to 'conservative'.
    """
    if method == 'conservative':
        ny, nx = data.shape[1:] if len(data.shape) == 3 else data.shape
        n_lons = nx * 4
        n_lats = nx * 2 + 1
        nt, i_out, j_out, w_out, i_in, j_in, w_in = gridspec
        data_map = np.zeros((n_lats, n_lons))
        ff = np.zeros((n_lats, n_lons))
        shell = np.zeros((2, n_lats, n_lons), dtype=np.float64)
        regridded_data, ff = conservative_regrid(shell, data, data_map, nt[0], undef, i_out, i_in, j_out, j_in, w_out, w_in, ff)
        regridded_data[np.where(ff != 0.0)] = regridded_data[np.where(ff != 0.0)] / ff[np.where(ff != 0.0)]
    return regridded_data

@jit(float64[:, :, :](float64[:, :, :], float64[:, :], float64[:, :], int64, float64, float32[:], float32[:], float32[:], float32[:], float32[:], float32[:], float64[:, :]), nopython=True)
def conservative_regrid(shell, cube_data, cube_data_0, nt, undef, i_out, i_in, j_out, j_in, w_out, w_in, ff):
    """
    Conservative regrid algorithm to interpolate data onto another grid.

    Parameters:
    - shell (np.ndarray): Shell to store intermediate and final results.
    - cube_data (np.ndarray): Input data to be interpolated.
    - cube_data_0 (np.ndarray): Interpolated data.
    - nt (int): Total number of iterations.
    - undef (float): Value representing undefined or missing data.
    - i_out (np.ndarray): Indices for output grid in x-direction.
    - i_in (np.ndarray): Indices for input grid in x-direction.
    - j_out (np.ndarray): Indices for output grid in y-direction.
    - j_in (np.ndarray): Indices for input grid in y-direction.
    - w_out (np.ndarray): Weighting factors for output grid.
    - w_in (np.ndarray): Weighting factors for input grid.
    - ff (np.ndarray): Fraction of total accumulation.

    Returns:
    - np.ndarray: Shell containing the final interpolated data and accumulation fractions.
    """
    for n in range(nt):
        i_out_sample = np.int32(i_out[n] - 1)
        j_out_sample = np.int32(j_out[n] - 1)
        i_in_sample = np.int32(i_in[n] - 1)
        j_in_sample = np.int32(j_in[n] - 1)
        validator = cube_data[j_in_sample, i_in_sample]

        if validator != undef:
            w_out_sample = w_out[n]
            cube_data_0[j_out_sample, i_out_sample] = cube_data_0[j_out_sample, i_out_sample] + w_out_sample * validator
            ff[j_out_sample, i_out_sample] = ff[j_out_sample, i_out_sample] + w_out_sample
            
        if n % np.int32((nt - 1) * 5e-2) == 0:
            print(str(np.int32(100 * n / (nt - 1))).strip() + '% Complete')

    print('100% Complete\n')
    shell[0, :, :] = cube_data_0
    shell[1, :, :] = ff
    return shell

def read_nt(tile_file: str, include_n_grids=False) -> int:
    """
    Reads the number of timesteps from a tile file.

    Parameters:
    - tile_file (str): Path to the tile file.
    - include_n_grids (bool, optional): If True, also includes the number of grids. Defaults to False.

    Returns:
    - int or tuple: Number of timesteps or a tuple containing the number of timesteps and the number of grids, if include_n_grids is True.
    """
    print(tile_file)
    subprocess.run(['ls', '-l', tile_file])
    with FortranFile(tile_file, 'r') as f:
        nt = f.read_ints('i4')
        if include_n_grids:
            n_grids = f.read_ints('i4')
            return nt, n_grids
    return nt

def read_tile_file(tile_file: str) -> float:
    """
    Reads data from a tile file.

    Parameters:
    - tile_file (str): Path to the tile file.

    Returns:
    - tuple: Tuple containing data read from the tile file.
    """
    with FortranFile(tile_file, 'r') as f:
        nt = f.read_ints('i4')
        n_grids = f.read_ints('i4')
        for grid in range(n_grids[0]):
            grid_name = f.read_ints('i4').tobytes()
            i_m = f.read_ints('i4')
            j_m = f.read_ints('i4')
            print(grid_name.decode().strip(), i_m[0], j_m[0])
            
        mask = f.read_reals('f4')
        x = f.read_reals('f4')
        y = f.read_reals('f4')
        
        i_out = f.read_reals('f4')
        j_out = f.read_reals('f4')
        w_out = f.read_reals('f4')
        i_in = f.read_reals('f4')
        j_in = f.read_reals('f4')
        w_in = f.read_reals('f4')

    return nt, i_out, j_out, w_out, i_in, j_in, w_in