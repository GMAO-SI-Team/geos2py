import os
import json
import time
import calendar
import argparse
import numpy as np
import cartopy.crs as ccrs
import matplotlib.cm as cm

from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import congrid, regrid, read_tile_file
from processing.scaling import bytscl
from datasets import loading
from epilogue.annotation import annotate

# Start timer
t0 = time.time()

# Parse command line arguments to specify data source and plot type
parser = argparse.ArgumentParser(description='Process date specified data and plot type')
parser.add_argument('year')
parser.add_argument('month')
parser.add_argument('day')
parser.add_argument('hour')
parser.add_argument('minute')
parser.add_argument('tag')
parser.add_argument('s_tag')
parser.add_argument('data_dir')
parser.add_argument('stream')
parser.add_argument('plot_type', choices=['aerosols', 'cape', 'ir8', 'precrain'])
parser.add_argument('--region')
parser.add_argument('--f_date')
args = parser.parse_args()

year = args.year
month = args.month
day = args.day
hour = args.hour
minute = args.minute

region = args.region if args.region else '-1'
tag = args.tag
f_date = args.f_date if args.f_date else f'{year}{month}{day}_{hour}z'
s_tag = args.s_tag if args.s_tag else f'f5295_fp-{f_date}'
stream = args.stream
data_dir = args.data_dir
plot_type = args.plot_type

# Load regions JSON file
with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

regions = region_info[region] if region in ('0', '-1') else [region]

times = {}

def plot_aerosols():
    # File path for aerosols data
    data_dir_aerosols = f'{data_dir}/GEOS.fp.fcst.inst1_2d_hwl_Nx.{f_date[:-1]}+20231013_0000.V01.nc4'
    
    # Load aerosol data
    ss = loading.load_cube(data_dir_aerosols, 'SSEXTTAU', 1e15, no_map_set=True).data
    du = loading.load_cube(data_dir_aerosols, 'DUEXTTAU', 1e15, no_map_set=True).data
    oc = loading.load_cube(data_dir_aerosols, 'OCEXTTAU', 1e15, no_map_set=True).data
    bc = loading.load_cube(data_dir_aerosols, 'BCEXTTAU', 1e15, no_map_set=True).data
    su = loading.load_cube(data_dir_aerosols, 'SUEXTTAU', 1e15, no_map_set=True).data
    ni = loading.load_cube(data_dir_aerosols, 'NIEXTTAU', 1e15, no_map_set=True).data
    
    # Resample the aerosols data to 5760x2760 shape
    aerosols = [congrid(aerosol, (2760, 5760), method='nearest', center=True) for aerosol in [ss, du, oc, bc, su, ni]]
    print(f'finished regridding: {time.time() - t0} seconds')

    # Define colormap and normalization for each aerosol
    aerosol_types = ['SSEXTTAU', 'DUEXTTAU', 'OCEXTTAU', 'BCEXTTAU', 'SUEXTTAU', 'NIEXTTAU']
    cmaps = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).cmap for aerosol in aerosol_types]
    norms = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).norm for aerosol in aerosol_types]
    
    # Convert scalar mappable objects to rgba arrays then scale alpha channels for blending
    data = []
    for aerosol, cmap, norm in zip(aerosols, cmaps, norms):
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        blend = sm.to_rgba(aerosol)
        blend[:, :, 3] = bytscl(aerosol, low=0, high=0.25)
        data.append(blend)
    
    plot_data(data, cmaps, norms, 'plotall_aerosols')

def plot_cape():
    # File path for CAPE data
    data_dir_cape = f'{data_dir}/GEOS.fp.fcst.inst3_2d_met_Nx.{f_date[:-1]}+20231224_0900.V01.nc4'
    
    # Define the CAPE colormap and normalization
    cape_cmap = Colormap('plotall_cape', 'CAPE')
    cmap = cape_cmap.cmap
    norm = cape_cmap.norm

    # Load and extract the CAPE data
    cube = loading.load_cube(data_dir_cape, 'CAPE', 1e15, pixel=True, no_project=True)
    data = loading.clean_data(cube.data, undef=1e15)

    # Resample data to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)

    # Mask data values below 100 (the lower bound)
    mask = data < 100
    data = np.ma.masked_where(mask, data)
    data = np.ma.masked_where(data > 10000, data)
    
    plot_data([data], [cmap], [norm], 'plotall_cape')

def plot_ir8():
    
    #
    data_dir = f'C:/Users/emonc/Documents/internship2024/stockv12-2024Jul04-1day-c24.geosgcm_surf.20000415_1930z.nc4'

    # Define the infrared colormap and normalization
    ir8_cmap = Colormap('plotall_ir8', 'TBISCCP')
    cmap = ir8_cmap.cmap
    norm = ir8_cmap.norm

    # Read and extract the data
    cube = loading.load_cube(data_dir, 'TBISCCP', 1e15, no_map_set=True)


    #
    '''
    if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
        tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
    else:
        tile_file = 'tiles/DC2880xPC1441_CF0720x6C.bin'
    gridspec = read_tile_file(tile_file)
    '''

    # Convert Kelvin to Celsius
    data = cube.data - 273.15

    # First order conservative regridding using 4320x720 to 2880x1441 gridspec
#    data = regrid(data, method='conservative', gridspec=gridspec, undef=1e15)

    # Resample data to 5760x2760 shape
    data = loading.clean_data(data, undef=1e15)
    data = congrid(data, (2760, 5760), center=True)

    # Open regions JSON and set region(s)
    with open('regions.json', 'r') as infile:
        region_info = json.load(infile)

    # Set region(s)
    regions = region_info['-1']

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

        # Start regional timer
        dt0 = time.time()

        # Initialize infrared plotter and plot data
        ir8_plotter = Plotter('plotall_ir8', region, file_tag, target_proj, proj_name)
        ir8_plotter.render(data, cmap, norm)
        
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
        annotate(f'tmp/{proj_name}-ir8-{file_tag}.png', 'plotall_ir8', mode=mode, forecast=forecast_str, date=date_index)
        print(f'{file_tag} saved successfully')

        # Record region time
        times[proj_name].append(time.time() - dt0)

    t = time.time() - t0

    # Evalute total compute time and average plotting time by projection type
    print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
    for proj in times:
        if times[proj]:
            t_avg = sum(times[proj]) / len(times[proj])
            print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')



def plot_precrain():
    # File path for rain data
    data_dir_precrain = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'
    
    # Define zero array for summing accumulated data
    saved_data = f'tmp/raindata-{f_date}.pkl'
    if os.path.exists(saved_data):
        with open(saved_data, 'rb') as pkl:
            acc_data = pickle.load(pkl)
    else:
        acc_data = np.zeros((2760, 5760))

    # Define the rain colormap and normalization
    rain_cmap = Colormap('plotall_precrain', 'PRECTOT')
    cmap = rain_cmap.cmap
    norm = rain_cmap.norm

    # Load and extract the rain data
    cube = loading.load_cube(data_dir_precrain, 'PRECTOT', 1e15, no_map_set=True).data
    acc_data += congrid(cube, (2760, 5760), center=True)
    
    plot_data([acc_data], [cmap], [norm], 'plotall_precrain')

def plot_data(data, cmap, norm, plot_tag):
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
            'rotpole': ccrs.RotatedPole(lon_cen, lat_cen)
        }
        target_proj = projs[proj_name]

        proj_extent = [lon_beg, lon_end, lat_beg, lat_end]

         # Initialize infrared plotter and plot data
        ir8_plotter = Plotter(plot_tag, region, file_tag, target_proj, proj_name)
        ir8_plotter.render(data, cmap, norm)
    
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
        annotate(f'tmp/{proj_name}-ir8-{file_tag}.png', 'plotall_ir8', mode=mode, forecast=forecast_str, date=date_index)
        print(f'{file_tag} saved successfully')

        # Record region time
#        times[proj_name].append(time.time() - dt0)

# Main routine to call appropriate plotting function based on specified plot type
if plot_type == 'aerosols':
    plot_aerosols()
elif plot_type == 'cape':
    plot_cape()
elif plot_type == 'ir8':
    plot_ir8()
elif plot_type == 'precrain':
    plot_precrain()

# End timer
t1 = time.time()
total = t1 - t0
print(f'Total runtime: {total} seconds')