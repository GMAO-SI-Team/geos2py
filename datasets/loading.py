import os
import sys
import numpy as np

from netCDF4 import Dataset

class Cube(object):
    def __init__(self, data, is1d, is2d, is3d, settings):
        self.data = data
        self.is1d = is1d
        self.is2d = is2d
        self.is3d = is3d
        self.settings = settings

def load_cube(infile: str, search_variable: str, undef, **kwargs):
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

        unlimited_dims = get_unlimited_dimensions(dataset)

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
                cube_data = curate_levels(nc_variables, search_variable, levels, level, nx, ny, time=time)
                settings['flip'] = True
                # cube_data_m = images.flip_hemispheres(cube_data[:, ::-1, :])
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

        # if not settings['no_project']:
        #     payload = arrays.ready_cube_interpolation(
        #         cube_data_m, g_lons, g_lats, lon_beg, lon_end, lat_beg, lat_end, undef, file_tag, settings
        #     )
        #     data = arrays.interpolate_cube(payload, settings)
        # else:
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
    Extracts defined optional parameters and defines default values for undefined optional parameters.
    Supports load_cube and load_goes.

    param: params - the optional parameters passed to load_cube or load_goes (arbitrary types and length)
    param: param_type - 'cube' (default) or 'goes'

    Returns a dictonary with all optional parameters defined.
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

def clean_data(data: np.ndarray, undef=None, fill='nan', data_type='cube'):
    """
    param: fill - np.nan (default) or 'zero' (0.0)
    param: data_type - 'cube' (default) or 'goes'
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
        
def curate_levels(nc_variables, search_variable, levels, level, nx, ny, time=None):
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

def get_xy_dimension_sizes(dims, nx, ny, nf, nt) -> tuple:
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

def get_nc_variable_data(dataset, search_variable: str) -> np.ndarray:
    return dataset.variables[search_variable][0].data

def get_unlimited_dimensions(dataset) -> list:
    return [dataset.dimensions[dim] for dim in dataset.dimensions if dataset.dimensions[dim].isunlimited()]