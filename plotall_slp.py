import os
import csv
import json
import time
import calendar
import argparse
import numpy as np
import PIL.Image as image
import cartopy.crs as ccrs

from scipy.interpolate import RegularGridInterpolator
from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import regrid, congrid, read_tile_file
from processing.smoothing import bandpass_filter
from datasets import loading
from epilogue.annotation import annotate

# Allows for large images
image.MAX_IMAGE_PIXELS = None

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

# Start timer
t0 = time.time()

# Parse command line arguments to specify data source
parser = argparse.ArgumentParser(description='Process date specified data')
parser.add_argument('year')
parser.add_argument('month')
parser.add_argument('day')
parser.add_argument('hour')
parser.add_argument('minute')
parser.add_argument('tag')
parser.add_argument('s_tag')
parser.add_argument('data_dir')
parser.add_argument('stream')
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

data_dir = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

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
    pressure_plotter = Plotter('plotall_slp', region, file_tag, target_proj, proj_name, label_coords=(slp_min_locs, x_window))
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
    annotate(f'tmp/{proj_name}-slp-{file_tag}.png', 'plotall_slp', mode=mode, forecast=forecast_str, date=date_index)
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