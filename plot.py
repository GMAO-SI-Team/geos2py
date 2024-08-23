import os
import csv
import json
import time
import pickle
import calendar
import argparse
import numpy as np
import cartopy.crs as ccrs
import matplotlib.cm as cm

from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import congrid, regrid, read_tile_file
from processing.smoothing import savitzky_golay2d, bandpass_filter
from processing.scaling import bytscl
from datasets import loading
from epilogue.annotation import annotate
from scipy.interpolate import RegularGridInterpolator

# Start timer
t0 = time.time()

# Parse command line arguments to specify data source and plot type

# Argument structure:
# python plot.py file_location --region=region_code plot_type

try:
    parser = argparse.ArgumentParser(description='Generate plots for given NetCDF file')

    parser.add_argument('data_dir', help="File location of NetCDF to be plotted")

    parser.add_argument('plot_type', choices=['aerosols', 'cape', 'ir8', 'precrain', 'precsnow', 't2m', 'tpw', 'vort500mb', 'winds10m', 'radar', 'wxtype', 'slp'], 
                        help="Variable to be plotted. A full list of names and their corresponding variables can be found in the ReadMe.")
    
    parser.add_argument('--region', help="Region code. A full list of region codes can be found in the ReadMe (default: all regions)")

    parser.add_argument('--cache_dir', help="Location of cache directory if not located in folder (default: /cache)")

    parser.add_argument('--results_dir', help="Location of directory in which to save images (default: /results)")

    args = parser.parse_args()
except argparse.ArgumentError as e:
    print(f"Error parsing arguments: {e}")
    exit(1)

try:
    data_dir = args.data_dir
    plot_type = args.plot_type
except AttributeError as e:
    print(f"Error accessing arguments: {e}")
    exit(1)

try:
    timestamp = loading.get_nc_timestamp(data_dir)
except Exception as e:
    print(f"Error loading timestamp from NetCDF file: {e}")
    exit(1)

try:
    year = timestamp.year
    month = timestamp.month
    day = timestamp.day
    hour = timestamp.hour
    minute = timestamp.minute
except AttributeError as e:
    print(f"Error accessing timestamp attributes: {e}")
    exit(1)

region = args.region if args.region else '-1'
cache_dir = args.cache_dir if args.cache_dir else 'cache'
results_dir = args.results_dir if args.results_dir else 'results'
f_date = f'{year}{month}{day}_{hour}z'
s_tag =  f'f5295_fp-{f_date}'

try:
    # Load regions JSON file
    with open('regions.json', 'r') as infile:
        region_info = json.load(infile)
except FileNotFoundError as e:
    print(f"Error: regions.json file not found: {e}")
    exit(1)
except json.JSONDecodeError as e:
    print(f"Error decoding regions.json file: {e}")
    exit(1)
except Exception as e:
    print(f"Error loading regions.json file: {e}")
    exit(1)

try:
    regions = region_info[region] if region in ('0', '-1') else [region]
except KeyError as e:
    print(f"Error: region '{region}' not found in regions.json: {e}")
    exit(1)

times = {}

def plot_cape():
    """
    Plotting routine for CAPE data
    """
    try:
        # Define the CAPE colormap and normalization
        cape_cmap = Colormap('plotall_cape', 'CAPE')
        cmap = cape_cmap.cmap
        norm = cape_cmap.norm
    except AttributeError as e:
        print(f"Error initializing CAPE colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Load and extract the CAPE data
        cube = loading.load_cube(data_dir, 'CAPE', 1e15, pixel=True, no_project=True)
        data = loading.clean_data(cube.data, undef=1e15)
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing CAPE data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or cleaning CAPE data: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        # Mask data values below 100 (the lower bound)
        mask = data < 100
        data = np.ma.masked_where(mask, data)
        data = np.ma.masked_where(data > 10000, data)
    except Exception as e:
        print(f"Error during data masking: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_cape')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_ir8():
    """
    Plotting routine for IR8 data
    """
    try:
        # Define the infrared colormap and normalization
        ir8_cmap = Colormap('plotall_ir8', 'TBISCCP')
        cmap = ir8_cmap.cmap
        norm = ir8_cmap.norm
    except AttributeError as e:
        print(f"Error initializing infrared colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read and extract the data
        cube = loading.load_cube(data_dir, 'TBISCCP', 1e15, no_map_set=True)
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing infrared data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading infrared data: {e}")
        return

    try:
        if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
            tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
        else:
            tile_file = f'data/DC2880xPC1441_CF0720x6C.bin'
        
        gridspec = read_tile_file(tile_file)
    except FileNotFoundError as e:
        print(f"Tile file not found: {e}")
        return
    except Exception as e:
        print(f"Unexpected error reading tile file: {e}")
        return

    try:
        # Convert Kelvin to Celsius
        data = cube.data - 273.15
    except AttributeError as e:
        print(f"Error accessing data for conversion: {e}")
        return
    except Exception as e:
        print(f"Unexpected error converting temperature data: {e}")
        return

    # First order conservative regridding using 4320x720 to 2880x1441 gridspec
    if cube.is2d:
        data = congrid(data, (2760, 1441), center=True)
    else:
        data = regrid(data, method='conservative', gridspec=gridspec, undef=1e15)

    
       

    try:
        # Resample data to 5760x2760 shape
        data = loading.clean_data(data, undef=1e15)
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_ir8')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_aerosols():
    """
    Plotting routine for Aerosol data
    """
    data_dir = args.data_dir

    try:
        # Load aerosol data
        ss = loading.load_cube(data_dir, 'SSEXTTAU', 1e15, no_map_set=True).data
        du = loading.load_cube(data_dir, 'DUEXTTAU', 1e15, no_map_set=True).data
        oc = loading.load_cube(data_dir, 'OCEXTTAU', 1e15, no_map_set=True).data
        bc = loading.load_cube(data_dir, 'BCEXTTAU', 1e15, no_map_set=True).data
        su = loading.load_cube(data_dir, 'SUEXTTAU', 1e15, no_map_set=True).data
        ni = loading.load_cube(data_dir, 'NIEXTTAU', 1e15, no_map_set=True).data
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing aerosol data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading aerosol data: {e}")
        return

    try:
        # Resample the aerosols data to 5760x2760 shape
        aerosols = [congrid(aerosol, (2760, 5760), method='nearest', center=True) for aerosol in [ss, du, oc, bc, su, ni]]
        print(f'Finished regridding: {time.time() - t0} seconds')
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    # Define colormap and normalization for each aerosol
    aerosol_types = ['SSEXTTAU', 'DUEXTTAU', 'OCEXTTAU', 'BCEXTTAU', 'SUEXTTAU', 'NIEXTTAU']
    cmaps = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).cmap for aerosol in aerosol_types]
    norms = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).norm for aerosol in aerosol_types]

    try:
        # Convert scalar mappable objects to RGBA arrays then scale alpha channels for blending
        data = []
        for aerosol, cmap, norm in zip(aerosols, cmaps, norms):
            sm = cm.ScalarMappable(norm=norm, cmap=cmap)
            blend = sm.to_rgba(aerosol)
            blend[:, :, 3] = bytscl(aerosol, low=0, high=0.25)
            data.append(blend)
    except Exception as e:
        print(f"Error during RGBA conversion or alpha channel scaling: {e}")
        return

    plot_data(data, cmaps, norms, 'plotall_aerosols')
    


def plot_precrain():
    """
    Plotting routine for Total Precipitable Rain data
    """
    try:
        # Define zero array for summing accumulated data
        saved_data = f'tmp/raindata-{f_date}.pkl'
        if os.path.exists(saved_data):
            with open(saved_data, 'rb') as pkl:
                acc_data = pickle.load(pkl)
        else:
            acc_data = np.zeros((2760, 5760))
    except FileNotFoundError as e:
        print(f"Error: Saved data file not found: {e}")
        return
    except pickle.UnpicklingError as e:
        print(f"Error loading pickled data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data initialization: {e}")
        return

    try:
        # Define the rain colormap and normalization
        rain_cmap = Colormap('plotall_precrain', 'PRECTOT')
        cmap = rain_cmap.cmap
        norm = rain_cmap.norm
    except AttributeError as e:
        print(f"Error initializing rain colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Load and extract the rain data
        cube = loading.load_cube(data_dir, 'PRECTOT', 1e15, no_map_set=True).data
        acc_data += congrid(cube, (2760, 5760), center=True)        
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing rain data: {e}")
        return
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or processing rain data: {e}")
        return
    
    acc_data *= 1440

    try:
        plot_data(acc_data, cmap, norm, 'plotall_precrain')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_precsnow():
    """
    Plotting routine for Total Precipitable Snow data
    """
    try:
        # Define zero array for summing accumulated data
        saved_data = f'tmp/snowdata-{f_date}.pkl'
        if os.path.exists(saved_data):
            with open(saved_data, 'rb') as pkl:
                acc_data = pickle.load(pkl)
        else:
            acc_data = np.zeros((2760, 5760))
    except FileNotFoundError as e:
        print(f"Error: Saved data file not found: {e}")
        return
    except pickle.UnpicklingError as e:
        print(f"Error loading pickled data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data initialization: {e}")
        return

    try:
        # Define the snow colormap and normalization
        snow_cmap = Colormap('plotall_precsnow', 'PRECSNO')
        cmap = snow_cmap.cmap
        norm = snow_cmap.norm
    except AttributeError as e:
        print(f"Error initializing snow colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Load and extract the data
        cube = loading.load_cube(data_dir, 'PRECSNO', 1e15, no_map_set=True, pixel=True, scale_factor=3600/25.4)
        data = loading.clean_data(cube.data, undef=1e15)
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing snow data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or cleaning snow data: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        # Mask data values below 0.1 (the lower bound)
        data = np.ma.masked_where(data < 0.1, data)
    except Exception as e:
        print(f"Error during data masking: {e}")
        return

    try:
        # Sum historic data and current data
        acc_data = acc_data + data
        acc_data *= 1440

    except Exception as e:
        print(f"Error summing data: {e}")
        return

    try:
        plot_data(acc_data, cmap, norm, 'plotall_precsnow')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_t2m():
    """
    Plotting routine for T2M data
    """
  

    try:
        # Define the temperature colormaps and normalization
        temp_cmap = Colormap('plotall_t2m', 'T2M')
        cmap = temp_cmap.cmap
        norm = temp_cmap.norm
    except AttributeError as e:
        print(f"Error initializing temperature colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read and extract the data
        data = loading.load_cube(data_dir, 'T2M', 1e15, no_map_set=True).data
        data = (data - 273.15) * 1.8 + 32  # convert Kelvin to Fahrenheit
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing temperature data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or processing temperature data: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_t2m')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_tpw():
    """
    Plotting routine for TPW data
    """
    try:
        # Define the water colormap and normalization
        tpw_cmap = Colormap('plotall_tpw', 'TQV')
        cmap = tpw_cmap.cmap
        norm = tpw_cmap.norm
    except AttributeError as e:
        print(f"Error initializing TPW colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read and extract the data
        cube = loading.load_cube(data_dir, 'TQV', 1e15, no_map_set=True, pixel=True)
        data = loading.clean_data(cube.data, undef=1e15)
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing TPW data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or processing TPW data: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_tpw')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_vort500mb():
    """
    Plotting routine for VORT500MB data
    """
    try:
        data_dir_vort = data_dir
    except Exception as e:
        print(f"Error in setting data directory: {e}")
        return

    try:
        # Define the vorticity colormap and normalization
        vorticity_cmap = Colormap('plotall_vort500mb', 'vort')
        cmap = vorticity_cmap.cmap
        norm = vorticity_cmap.norm
    except AttributeError as e:
        print(f"Error initializing vorticity colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read and extract the data
        heights = loading.load_cube(data_dir_vort, 'H500', 1e15, no_map_set=True).data
        u = loading.load_cube(data_dir_vort, 'U500', 1e15, no_map_set=True).data
        v = loading.load_cube(data_dir_vort, 'V500', 1e15, no_map_set=True).data
    except FileNotFoundError as e:
        print(f"NetCDF file not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing vorticity data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading or processing vorticity data: {e}")
        return

    try:
        # Compute forces and derivatives
        rad = np.pi / 180.0
        re = 6371220.0
        lats = 2.0 * (np.arange(720) / (720 - 1) - 0.5) * 90 * rad
        deg_m = 2 * re * np.pi / 720
        dlon_m = deg_m * np.cos(lats * np.pi / 180)
        dlat_m = deg_m
        dvdx = np.zeros((720, 720))
        dudy = np.zeros((720, 720))
        data = np.zeros((720, 720))

        for j in range(720):
            for i in range(720):
                ip1 = i + 1
                if i == 719:
                    ip1 = 0
                dvdx[j, i] = (v[j, ip1] - v[j, i]) / dlon_m[j]
                jp1 = j + 1
                if j == 719:
                    jp1 = j
                dudy[j, i] = (u[jp1, i] - u[j, i]) / dlat_m
                data[j, i] = (dvdx[j, i] - dudy[j, i])

        data = data * 1.e5

        for j in range(720):
            if lats[j] <= 0:
                data[j, :] = -1 * data[j, :]
    except Exception as e:
        print(f"Error in computing forces and derivatives: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        heights = congrid(heights, (2760, 5760), center=True)
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        # Smooth heights
        heights = savitzky_golay2d(heights, int(5760 * 0.025) + 1, 2)
    except Exception as e:
        print(f"Error during height smoothing: {e}")
        return

    try:
        print('data', data.min(), data.max())
        print('heights', heights.min(), heights.max())
    except Exception as e:
        print(f"Error printing data/height statistics: {e}")
        return

    try:
        # Mask data values < 2.5 (lower bound) and > 60 (upper bound)
        mask = np.logical_or(data < 2.5, data > 60)
        data = np.ma.masked_array(data, mask)
        print('data', data.min(), data.max())

        data = [data, heights]
    except Exception as e:
        print(f"Error in masking data: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_vort500mb')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_winds10m():
    """
    Plotting routine for WINDS10M data
    """
    try:
        # Define the wind colormap and normalization
        wind_cmap = Colormap('plotall_winds10m', 'spd')
        cmap = wind_cmap.cmap
        norm = wind_cmap.norm
    except AttributeError as e:
        print(f"Error initializing wind colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read, extract, and scale the data
        vectors = {
            'ULML': {
                'data_file': data_dir,
                'scale_factor': 2.23694
            },
            'VLML': {
                'data_file': data_dir,
                'scale_factor': 2.23694
            }
        }
        wind_data = []
        for vector in vectors:
            try:
                print(f'Reading {vector.lower()}')
                cube = loading.load_cube(
                    vectors[vector]['data_file'], vector, 1e15, pixel=True, scale_factor=vectors[vector]['scale_factor']
                )
                data = loading.clean_data(cube.data, undef=1e15)
                wind_data.append(data)
                print()
            except FileNotFoundError as e:
                print(f"NetCDF file for {vector} not found: {e}")
                return
            except AttributeError as e:
                print(f"Error accessing {vector} data: {e}")
                return
            except Exception as e:
                print(f"Unexpected error loading or processing {vector} data: {e}")
                return
    except Exception as e:
        print(f"Error in reading or extracting wind data: {e}")
        return

    try:
        # Compute vector magnitude
        ulml, vlml = wind_data
        data = np.sqrt(ulml ** 2 + vlml ** 2)
    except Exception as e:
        print(f"Error in computing vector magnitude: {e}")
        return

    try:
        # Resample data to 5760x2760 shape
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_winds10m')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return


def plot_radar():
    """
    Plotting routine for RADAR data
    """
    try:
        data_dir_radar = data_dir
    except NameError as e:
        print(f"Error: data_dir is not defined: {e}")
        return

    try:
        # Define the radar colormap and normalization
        radar_cmap = Colormap('plotall_radar', 'DBZ_MAX')
        cmap = radar_cmap.cmap
        norm = radar_cmap.norm
    except AttributeError as e:
        print(f"Error initializing radar colormap: {e}")
        return
    except Exception as e:
        print(f"Unexpected error in colormap setup: {e}")
        return

    try:
        # Read and extract the data
        cube = loading.load_cube(data_dir_radar, 'DBZ_MAX', 1e15, no_map_set=True)
    except FileNotFoundError as e:
        print(f"NetCDF file for radar data not found: {e}")
        return
    except AttributeError as e:
        print(f"Error accessing radar data: {e}")
        return
    except Exception as e:
        print(f"Unexpected error loading radar data: {e}")
        return

    try:
        # Determine the appropriate tile file path
        if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
            tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
        else:
            tile_file = f'data/DC2880xPC1441_CF0720x6C.bin'
    except Exception as e:
        print(f"Error determining tile file path: {e}")
        return

    try:
        gridspec = read_tile_file(tile_file)
    except FileNotFoundError as e:
        print(f"Tile file not found: {e}")
        return
    except Exception as e:
        print(f"Unexpected error reading tile file: {e}")
        return

    try:
         # First order conservative regridding using 4320x720 to 2880x1441 gridspec
        if cube.is2d:
            data = congrid(data, (2760, 1441), center=True)
        else:
            data = regrid(data, method='conservative', gridspec=gridspec, undef=1e15)
    except Exception as e:
        print(f"Error during regridding: {e}")
        return

    try:
        data = loading.clean_data(data, undef=1e15)
    except Exception as e:
        print(f"Error during data cleaning: {e}")
        return

    try:
        data = congrid(data, (2760, 5760), center=True)
    except ValueError as e:
        print(f"Error in data resampling: {e}")
        return
    except Exception as e:
        print(f"Unexpected error during data resampling: {e}")
        return

    try:
        # Mask data values below 5 (the lower bound)
        mask = data < 5
        data = np.ma.masked_where(mask, data)
    except Exception as e:
        print(f"Error masking data: {e}")
        return

    try:
        plot_data(data, cmap, norm, 'plotall_radar')
    except Exception as e:
        print(f"Error during plotting: {e}")
        return



# These functions are attempting to locate sea level pressure minima
# Not yet complete
# Does not identify the appropriate values
def find_slp_mins(slp, dlons, dlats):
    slp_threshold = 1004.0  # Mb
    slp_smoothed = bandpass_filter(slp, 0, 0.8, ideal=True)

    nlats, nlons = slp.shape

    storm = np.zeros((nlats, nlons))

    slp_minima = find_storm_minima_mask(slp_smoothed, 1e-8)

    ifind = np.where(slp_smoothed > slp_threshold)
    if len(ifind[0]):
        slp_minima[ifind] = 0
    ifind = np.where(slp_minima == 1)
    print(np.min(slp_smoothed[ifind]))
    iict1 = np.size(ifind[0])
    storm[ifind] = 1
    print("There are", iict1, "points with slp_minima=1")

    ifind = np.where(storm == 1)
    sorting = np.argsort(slp[ifind])
    ifind = ifind[0][sorting]
    iict4 = 0

    dx = 10 / dlons
    dy = 10 / dlats
    nstorms = np.size(ifind)
    for n in range(nstorms):
        i = ifind[n] % nlons
        j = ifind[n] // nlons
        i_start = int(max(i - dx, 0))
        i_end = int(min(i + dx, nlons - 1))
        j_start = int(max(j - dy, 0))
        j_end = int(min(j + dy, nlats - 1))
        storm_local = storm[j_start:j_end+1, i_start:i_end+1]
        slp_local = slp[j_start:j_end+1, i_start:i_end+1]
        icount = np.where(storm_local == 4)
        if np.size(icount) == 0:
            imark = np.where(storm_local == 1)
            if np.size(imark) != 0:
                sorting = np.argsort(slp_local[imark])
                imark = imark[0][sorting]
                storm_local[imark] = 0
                storm_local[imark[0]] = 4
                iict4 += 1
            else:
                storm_local = 0
            storm[j_start:j_end+1, i_start:i_end+1] = storm_local
        
    count4 = np.count_nonzero(storm == 4)

    print("There are", count4, "points with storm=4")

    ilocs = np.zeros((count4, 2), dtype=int)
    ij = 0
    for j in range(nlats):
        for i in range(nlons):
            if storm[j, i] == 4:
                ilocs[ij, 0] = i
                ilocs[ij, 1] = j
                ij += 1

    return ilocs

def find_storm_minima_mask(slp_smoothed, fuzz_threshold):
    mask = np.logical_and.reduce((
        slp_smoothed < np.roll(slp_smoothed, 1, axis=1),
        slp_smoothed < np.roll(slp_smoothed, -1, axis=1),
        slp_smoothed < np.roll(slp_smoothed, 1, axis=0),
        slp_smoothed < np.roll(slp_smoothed, -1, axis=0)
    ))

    fuzz = slp_smoothed[mask]
    fuzz_find = np.where(fuzz < fuzz_threshold)

    if len(fuzz_find[0]):
        mask[np.where(mask)[fuzz_find]] = False
    
    return mask
 
def plot_slp():
    # Define the sea level pressure colormap and normalization
    pressure_cmap = Colormap('plotall_slp', 'SLP')
    cmap = pressure_cmap.cmap
    norm = pressure_cmap.norm

    # Read and extract the data
    data = loading.load_cube(data_dir, 'SLP', 1e15, no_map_set=True, scale_factor=0.01).data

    # Resample to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)

    # Open regions JSON
    with open('regions.json', 'r') as infile:
        region_info = json.load(infile)

    # Set region(s)
    regions = region_info[region] if region in ('0', '-1') else [region]

    times = {}

    for region in regions:
        # Define plot parameters for specified region
        region = str(region)
        file_tag = region_info[region]['file_tag']
        lon_cen, lat_cen = region_info[region]['center']
        lon_beg, lon_end, lat_beg, lat_end = region_info[region]['extent'] if 'extent' in region_info[region] else (-180, 180, -90, 90)
        proj_name = region_info[region]['proj'] if file_tag not in ['australia_mapset', 'southamerica_mapset'] else 'sub'
        projs = {
            'sub': ccrs.PlateCarree(),
            'ortho': ccrs.Orthographic(lon_cen, lat_cen),
            'laea': ccrs.LambertAzimuthalEqualArea(lon_cen, lat_cen),
            'geos': ccrs.Geostationary(lon_cen),
            'nsper': ccrs.NearsidePerspective(lon_cen, lat_cen)
        }
        target_proj = projs[proj_name]
        for proj in projs:
            times[proj] = [] if proj not in times else times[proj]

        # Locate sea level pressure minima
        x_window = (lon_end - lon_beg) / 5760
        y_window = (lat_end - lat_beg) / (2760 - 1)
        # slp_mins = find_storm_minima(data, x_window, y_window)
        slp_mins = find_slp_mins(data, x_window, y_window)
        print(f'SLP minima array shape: {slp_mins.shape}')
        
        n_mins = slp_mins.shape[0]
        slp_min_locs = np.zeros((n_mins, 3))
        for s in range(n_mins):
            slp_min_locs[s, 0] = slp_mins[s, 0] * x_window + 0.5 * x_window - 180
            slp_min_locs[s, 1] = slp_mins[s, 1] * y_window - 90
            slp_min_locs[s, 2] = data[slp_mins[s, 1], slp_mins[s, 1]]

        # Start regional timer
        dt0 = time.time()

        # Initialize sea level pressure plotter and plot data
        pressure_plotter = Plotter('plotall_slp', region, file_tag, target_proj, proj_name, cache_dir, label_coords=(slp_min_locs, x_window))
        pressure_plotter.render(data, cmap, norm)
        
        # Set image annotation text
        forecast_hours = f'000 Forecast Hours'  
        forecast_hours = f'{forecast_hours}\n INIT: {f_date}'
        forecast_p_tag = f'GEOS {s_tag.split("-")[0]}'
        forecast_str = f'{forecast_hours}\n {forecast_p_tag}'
        date_index = f'{year}-{month}-{day} {hour}:{minute}Z'
        date_str = f'{year} {calendar.month_name[int(month)]} {day}'
        time_index = f'{hour}:{minute}am EDT {calendar.day_name[calendar.weekday(int(year), int(month), int(day))]}'
        date_index = f'{date_index}\n{date_str}\n{time_index}'

        satellite = ['geos', 'nsper']
        mode = 'dark' if proj_name in satellite else 'light'

        # Annotate final image
        annotate(f'tmp/{proj_name}-slp-{file_tag}.png', 'plotall_slp', results_dir, mode=mode, forecast=forecast_str, date=date_index)
        print(f'{file_tag} saved successfully')

        # Record region time
        times[proj_name].append(time.time() - dt0)

    t = time.time() - t0



def plot_data(data, cmap, norm, plot_tag):
    try:
        for region in regions:
            try:
                # Define plot parameters for specified region
                region = str(region)
                file_tag = region_info[region]['file_tag']
                lon_cen, lat_cen = region_info[region]['center']
                proj_name = region_info[region]['proj'] if file_tag not in ['australia_mapset', 'southamerica_mapset'] else 'sub'
                projs = {
                    'sub': ccrs.PlateCarree(),
                    'ortho': ccrs.Orthographic(lon_cen, lat_cen),
                    'laea': ccrs.LambertAzimuthalEqualArea(lon_cen, lat_cen),
                    'geos': ccrs.Geostationary(lon_cen),
                    'nsper': ccrs.NearsidePerspective(lon_cen, lat_cen)
                }
                target_proj = projs[proj_name]
                for proj in projs:
                    times[proj] = [] if proj not in times else times[proj]
            except KeyError as e:
                print(f"Error: Region information missing or incorrect: {e}")
                continue
            except Exception as e:
                print(f"Unexpected error in setting up region plot parameters: {e}")
                continue

            try:
                # Start regional timer
                dt0 = time.time()

                    # Read city coordinates
                cities = 'cities/all_cities_md.txt' if file_tag == 'maryland_mapset' else 'cities/all_cities.txt'
                city_coords = []
                with open(cities, 'r') as csvfile:
                    reader = csv.reader(csvfile)
                    for row in reader:
                        row = [item for item in row if item]
                        if len(row) == 4 or 'world' in cities:
                            city_lat, city_lon = [float(coord.strip().replace('+', '')) for coord in row[2:4]]
                            city_coords.append((city_lat, city_lon))

                # Define linear interpolator to locate cities by pixel
                lons = np.linspace(-180, 180, 5760)
                lats = np.linspace(-90, 90, 2760)
                city_interpolator = RegularGridInterpolator((lats, lons), data, method='linear')


                # Initialize plotter with tag and plot data
                plotter = Plotter(plot_tag, region, file_tag, target_proj, proj_name, cache_dir, label_coords=city_coords, interpolator=city_interpolator)
                plotter.render(data, cmap, norm)
            except Exception as e:
                print(f"Error during plotting for region {region}: {e}")
                continue

            try:
                # Set image annotation text
                forecast_hours = f'000 Forecast Hours'
                forecast_hours = f'{forecast_hours}\n INIT: {f_date}'
                forecast_p_tag = f'GEOS {s_tag.split("-")[0]}'
                forecast_str = f'{forecast_hours}\n {forecast_p_tag}'
                date_index = f'{year}-{month}-{day} {hour}:{minute}Z'
                date_str = f'{year} {calendar.month_name[int(month)]} {day}'
                time_index = f'{hour}:{minute}am EDT {calendar.day_name[calendar.weekday(int(year), int(month), int(day))]}'
                date_index = f'{date_index}\n{date_str}\n{time_index}'

                satellite = ['geos', 'nsper']
                mode = 'dark' if proj_name in satellite else 'light'

                # Annotate final image
                annotate(f'tmp/{proj_name}-{plot_type}-{file_tag}.png', plot_tag, results_dir, mode=mode, forecast=forecast_str, date=date_index)
                print(f'{file_tag} saved successfully')
            except Exception as e:
                print(f"Error during image annotation for region {region}: {e}")
                continue

            try:
                # Record region time
                times[proj_name].append(time.time() - dt0)
            except Exception as e:
                print(f"Error recording time for region {region}: {e}")
                continue

        try:
            t = time.time() - t0

            # Evaluate total compute time and average plotting time by projection type
            print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
            for proj in times:
                if times[proj]:
                    t_avg = sum(times[proj]) / len(times[proj])
                    print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')
        except Exception as e:
            print(f"Error in computing total time and averages: {e}")

    except Exception as e:
        print(f"Unexpected error in plot_data function: {e}")



# Main routine to call appropriate plotting function based on specified plot type
case = {
    'cape': plot_cape,
    'ir8': plot_ir8,
    'aerosols': plot_aerosols,
    'precrain': plot_precrain,
    'precsnow': plot_precsnow,
    't2m': plot_t2m,
    'tpw': plot_tpw,
    'vort500mb': plot_vort500mb,
    'winds10m': plot_winds10m,
    'radar': plot_radar,
    'slp': plot_slp
}

# Call the function using the key 'plot_type'
routine = case.get(plot_type)
if routine:
    routine()
else:
    print(f"No plotting routine found for: {plot_type}")


# End timer
t1 = time.time()
total = t1 - t0
print(f'Total runtime: {total} seconds')
