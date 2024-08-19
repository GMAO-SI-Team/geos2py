import os
import sys
import netCDF4
import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date
from datetime import datetime

class Cube(object):
    """
    Represents a cubed sphere dataset.

    Attributes:
    - data (np.ndarray): The cube data.
    - is1d (bool): True if the cube is 1-dimensional, False otherwise.
    - is2d (bool): True if the cube is 2-dimensional, False otherwise.
    - is3d (bool): True if the cube is 3-dimensional, False otherwise.
    - settings (dict): Settings for loading the cube.
    """
    def __init__(self, data: np.ndarray, is1d: int, is2d: int, is3d: int, settings: dict):
        self.data = data
        self.is1d = is1d
        self.is2d = is2d
        self.is3d = is3d
        self.settings = settings

def load_cube(infile: str, search_variable: str, undef: float, **kwargs) -> Cube:
    """
    Loads cube data from a NetCDF file.

    Parameters:
    - infile (str): Path to the NetCDF file.
    - search_variable (str): Variable to search in the NetCDF file.
    - undef (float): Undefined value for cleaning cube data.
    - **kwargs: Additional keyword arguments.

    Returns:
    - class Cube: An instance of the Cube class containing loaded data.
    """
    print(infile)
    settings = extract_load_settings(kwargs)

    level = settings['level']
    time = 0 if not settings['time'] else settings['time']
    do_multiply, do_subtract = [0, 0]
    search_variables = np.array(search_variable.split('*'))
    shape = search_variables.shape
    if shape[0] == 1:
        search_variables = np.array(search_variable.split('-'))
        do_subtract = 1 if search_variables.shape[0] == 2 else do_subtract
    else:
        do_multiply = 1

    if os.path.exists(infile):
        dataset = Dataset(infile)
        nc_variables = dataset.variables
        data = get_nc_variable_data(dataset, search_variable)

        nx, ny, nf, nt, nz = [0] * 5

        dims = dataset.dimensions
        keys = list(dims.keys())

        nt, ny, nx, nz = [dims[keys[i]].size for i in range(4)] if 'ecmwf.inst3_3d_wxm_Np' in infile else [0] * 4

        nx, ny, nt = [dims[keys[i]].size for i in list(range(2, 4)) + [5]] if 'NOinst_30mn_met_c0720sfc' in infile else [0] * 3
        nz = 6 if 'NOinst_30mn_met_c0720sfc' in infile else nz

        if 'NOinst_30mn_met_c0720sfc' not in infile and 'ecmwf.inst3_3d_wxm_Np' not in infile:
            nx, ny, nf, nt = get_xy_dimension_sizes(dims, nx, ny, nf, nt)

            if nx == 0 and ny == 0:
                if len(dims) >= 2:
                    for i in range(2):
                        dim_size = dims[keys[i]].size
                        nx = dim_size if keys[i].strip() == 'lon' else nx
                        ny = dim_size if keys[i].strip() == 'lat' else ny
                nz = 0
                nt = 1
                
                if len(dims) >= 3:
                    dim_size = dims[keys[2]].size
                    nt = dim_size if keys[2] == 'time' else nt
                    nz = dim_size if keys[2] != 'time' else nz

                nt = dims[keys[3]].size if len(dims) >= 4 else nt

                if len(dims) >= 8:
                    nx = dims[keys[1]].size
                    ny = dims[keys[3]].size
                    nz = dims[keys[4]].size
                    nf = dims[keys[6]].size

                nt = dims[keys[8]].size if len(dims) >= 9 else nt

        is_cube_old_ensemble_dim = 1 if nx == 6 else 0
        is_cube_ensemble_dim = 1 if nz == 6 else 0

        if not undef:
            print('KeyError: Missing keyword \'undef\'. Cannot clean cube data.')
            sys.exit()

        if nf == 6 or is_cube_ensemble_dim or is_cube_old_ensemble_dim:
            cube_data_c = get_nc_variable_data(dataset, search_variables[0])[:, :, :]
            cube_data_c = clean_data(cube_data_c, undef=undef)
            cube_data_m = np.zeros((ny * 6, ny))
            for c in range(6):
                j_0 = ny * c
                j_1 = ny * c + ny
                cube_data_m[j_0:j_1, :] = cube_data_c[c, :, :]
            nx = ny
            ny = ny * 6
            nz = 1
        else:
            if 'ecmwf.inst3_3d_wxm_Np' in infile:
                levels = nc_variables['isobaric']
                settings['flip'] = True
            else:
                levels = nc_variables['lev'] if level and nz >=1 else None

                if search_variable == 'WSPD50M':
                    u = nc_variables['U50M']
                    v = nc_variables['V50M']
                    cube_data_m = np.sqrt(u**2 + v**2)

                else:
                    if shape[0] == 1:
                        if level:
                            cube_data_m = curate_levels(nc_variables, search_variable, levels, level, nx, ny, time=time)
                        else:
                            cube_data_m = get_nc_variable_data(dataset, search_variable)
                    else:
                        if level:
                            cube_data_a = curate_levels(nc_variables, search_variables[0], levels, level)
                            cube_data_a = clean_data(cube_data_a, undef=undef, fill='zero')
                            cube_data_b = curate_levels(nc_variables, search_variables[1], levels, level)
                            cube_data_b = clean_data(cube_data_b, undef=undef, fill='zero')
                            cube_data_m = cube_data_a * cube_data_b if do_multiply else cube_data_m
                            cube_data_m = cube_data_a - cube_data_b if do_subtract else cube_data_m
                            cube_data_a, cube_data_b = (0, 0)

        cube_data_m = clean_data(cube_data_m, undef=undef, fill='zero')

        if settings['scale_factor']:
            scale_factor = settings['scale_factor']
            print(f'SCALING by: {scale_factor}')
            cube_data_m = cube_data_m * scale_factor

        if dataset.isopen():
            dataset.close()

        data = cube_data_m

    else:
        print(search_variable)
        print(infile)
        sys.exit()

    is1d, is2d, is3d = [len(data.shape) == dim for dim in range(1, 4)]
    return Cube(data, is1d, is2d, is3d, settings)

def load_lcc(infile, search_variable):
    dataset = Dataset(infile)
    dims = dataset.dimensions
    nx, ny, nf, nt, nz = [0] * 5
    nx, ny, nf, nt = get_xy_dimension_sizes(dims, nx, ny, nf, nt)
    data = get_nc_variable_data(dataset, search_variable)
    if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
        grid = Dataset('/discover/nobackup/projects/gmao/osse2//stage/BCS_FILES/lambert_grid.nc4')
    else:
        grid = Dataset('data/lambert_grid.nc4')
    lons = get_nc_variable_data(grid, 'lons')
    lats = get_nc_variable_data(grid, 'lats')
    return data, lons, lats

def extract_load_settings(params: dict, param_type='cube') -> dict:
    """
    Extracts and defines default optional parameters for loading cube or GOES data.

    Parameters:
    - params (dict): Optional parameters passed to load_cube or load_goes.
    - param_type (str, optional): Type of parameters ('cube' or 'goes'). Defaults to 'cube'.

    Returns:
    - dict: A dictionary with all optional parameters defined.
    """
    if param_type not in ['cube', 'goes']:
        raise Exception('ValueError: Invalid param_type. Only \'cube\' and \'goes\' params are supported.')
    
    if param_type == 'cube':
        optional_params = {
            'triangles': False,
            'scale_factor': None,
            'smooth': False,
            'pixel': False,
            'flip': False,
            'regrid_cube_to_latlon': False,
            'time': None,
            'level': False,
            'data_levels': None,
            'range': False,
            'regrid_resolution': False,
            'log': False,
            'n_colors': False,
            'project': False,
            'detail': False, 
            'stretch': False,
            'bilinear': False,
            'conservative': False,
            'no_project': False,
            'no_map_set': False
        }
        
    if param_type == 'goes':
        optional_params = {
            'olr': False, 
            'flip': False, 
            'pixel': False, 
            'scale_factor': None, 
            'no_project': False
        }
        
    for optional_param in optional_params:
        params[optional_param] = params[optional_param] if optional_param in params else optional_params[optional_param]

    return params

def clean_data(data: np.ndarray, undef=None, fill='nan', data_type='cube') -> np.ndarray:
    """
    Cleans cube or GOES data by replacing undefined values with np.nan or zero.

    Parameters:
    - data (np.ndarray): Data to be cleaned.
    - undef (float): Undefined value for cleaning cube or GOES data.
    - fill (str): Fill value ('nan' or 'zero').
    - data_type (str): Type of data ('cube' or 'goes').

    Returns:
    - np.ndarray: Cleaned data.
    """
    if data_type not in ['cube', 'goes']:
        raise Exception('TypeError: Invalid data_type parameter. Only \'cube\' and \'goes\' are supported.')
    if undef:
        fill = {'nan': np.nan, 'zero': 0.0}[fill]
        is_undef = np.where(data == undef) if data_type == 'cube' else np.where(data >= undef)
        if len(is_undef) == len(data):
            sys.exit()
        if len(is_undef[0]) == 0:                 
            data[is_undef] = fill
            return data
        else:
            mask = data == undef
            return np.ma.masked_where(mask, data)
        
def curate_levels(nc_variables, search_variable, levels, level: int, nx: int, ny: int, time=None) -> np.ndarray:
    """
    Curates levels for cube data based on the specified level.

    Parameters:
    - nc_variables: NetCDF variables.
    - search_variable (str): Variable to search in the NetCDF file.
    - levels: Levels in the NetCDF file.
    - level (int): Specified level.
    - nx (int): Size of x-dimension.
    - ny (int): Size of y-dimension.
    - time: Specified time.

    Returns:
    - np.ndarray: Curated cube data.
    """
    if level:
        print(f'Reading LEVEL: {level}')
        if time:
            print(f'Reading TIME: {time}')
        same_levels = np.where(levels == level)
        if len(same_levels[0]) == 0:
            cube_data_m = nc_variables[search_variable][0].data # figure out count and offset
        else:
            print(f'ERROR: bad level: {level}')
            print(levels)
            sys.exit()

        cube_data = np.zeros((ny, nx))
        cube_data[:, :] = cube_data_m[same_levels[0], :, :]
    else:
        cube_data = nc_variables[search_variable][0].data
    return cube_data

def get_xy_dimension_sizes(dims: list, nx: int, ny: int, nf: int, nt:int) -> tuple:
    """
    Retrieves x and y dimension sizes from NetCDF dimensions.

    Parameters:
    - dims (list): NetCDF dimensions.
    - nx (int): Size of x-dimension.
    - ny (int): Size of y-dimension.
    - nf (int): Size of nf-dimension.
    - nt (int): Size of time-dimension.

    Returns:
    - tuple: Sizes of x and y dimensions.
    """
    print(len(dims))
    print('----------')
    for dim in dims:
        dim_size = dims[dim].size
        print(f'{dim} {dim_size}')
        nx = dim_size if dim.strip() in ['lon', 'Xdim'] else nx
        ny = dim_size if dim.strip() in ['lat', 'Ydim'] else ny
        nf = dim_size if dim.strip() == 'nf' else nf
        nt = dim_size if dim.strip() == 'time' else nt
    print('----------')
    print(len(dims))
    return nx, ny, nf, nt

def get_nc_variable_data(dataset: netCDF4.Dataset, search_variable: str) -> np.ndarray:
    """
    Retrieves data for a specific variable from a NetCDF dataset.

    Parameters:
    - dataset (netCDF4.Dataset): NetCDF dataset.
    - search_variable (str): Variable to search in the NetCDF file.

    Returns:
    - np.ndarray: Data array for the specified variable.
    """

    return dataset.variables[search_variable][0].data


def get_unlimited_dimensions(dataset) -> list:
    """
    Retrieves unlimited dimensions from a NetCDF dataset.

    Parameters:
    - dataset (netCDF4.Dataset): NetCDF dataset.

    Returns:
    - list: List of unlimited dimensions.
    """
    return [dataset.dimensions[dim] for dim in dataset.dimensions if dataset.dimensions[dim].isunlimited()]


def get_nc_timestamp(datadir: str) -> datetime: 
    """
    Retrieves timestamp from a NetCDF dataset.

    Parameters:
    - string (str): NetCDF filename.

    Returns:
    - datetime (datetime.datetime): datetime retrieved from NetCDF
    """


    dataset = Dataset(datadir)
    time_var = dataset.variables['time']
    return num2date(time_var[:],time_var.units)[0]
