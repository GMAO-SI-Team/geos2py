import subprocess
import numpy as np

from scipy.io import FortranFile
from sunpy.image import resample
from numba import jit, float32, float64, int64, int32

def congrid(data, shape, center=True):
    return resample.resample(data, shape, center=center)

def regrid(data, method='conservative', gridspec=None, undef=1e15):
    """
    method: str (default: 'conservative')
        Supported values are 'conservative' and 'bilinear'.
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

def read_nt(tile_file: str, include_n_grids=False):
    print(tile_file)
    subprocess.run(['ls', '-l', tile_file])
    with FortranFile(tile_file, 'r') as f:
        nt = f.read_ints('i4')
        if include_n_grids:
            n_grids = f.read_ints('i4')
            return nt, n_grids
    return nt

def read_tile_file(tile_file: str):
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