import os
import json
import time
import argparse
import numpy as np
import cartopy.crs as ccrs

from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import congrid
from processing.smoothing import savitzky_golay2d
from datasets import loading
from epilogue.annotation import annotate

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

# Define the vorticity colormap and normalization
vorticity_cmap = Colormap('plotall_vort500mb', 'vort')
cmap = vorticity_cmap.cmap
norm = vorticity_cmap.norm

# Read and extract the data
heights = loading.load_cube(data_dir, 'H500', 1e15, no_map_set=True).data
u = loading.load_cube(data_dir, 'U500', 1e15, no_map_set=True).data
v = loading.load_cube(data_dir, 'V500', 1e15, no_map_set=True).data

# Compute forces and derivatives
rad = np.pi / 180.0
re = 6371220.0
rr = re * rad
lats = 2.0 * (np.arange(720) / (720 - 1) - 0.5) * 90 * rad
omeg = 7.0721e-5
coriolis = 2 * omeg * np.sin(lats * np.pi / 180)
deg_m = 2 * re * np.pi / 720
dlon_m = deg_m * np.cos(lats * np.pi / 180)
dlat_m = deg_m
dvdx = np.zeros((720, 720))
dudy = np.zeros((720, 720))
data = np.zeros((720, 720))

for j in range(720):
    for i in range(720):
        ip1 = i + 1
        if (i == 719):
            ip1 = 0
        dvdx[j, i] = (v[j, ip1] - v[j, i]) / dlon_m[j]
        jp1 = j + 1
        if (j == 719):
            jp1 = j
        dudy[j, i] = (u[jp1, i] - u[j, i]) / dlat_m
        data[j, i] = (dvdx[j, i] - dudy[j, i])

data = data * 1.e5

for j in range(720):
    if lats[j] <= 0:
        data[j, :] = -1 * data[j, :]

# Resample data to 5760x2760 shape
heights = congrid(heights, (2760, 5760), center=True)
data = congrid(data, (2760, 5760), center=True)

# Smooth heights
heights = savitzky_golay2d(heights, int(5760 * 0.025) + 1, 2)

print('data', data.min(), data.max())
print('heights', heights.min(), heights.max())

# Mask data values < 2.5 (lower bound) and > 60 (upper bound)
mask = np.logical_or(data < 2.5, data > 60)
data = np.ma.masked_array(data, mask)
print('data', data.min(), data.max())

data = [data, heights]

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
    proj_name = region_info[region]['proj'] if not file_tag in ['australia_mapset', 'southamerica_mapset'] else 'sub'
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

    # Initialize vorticity plotter and plot data
    vorticity_plotter = Plotter('plotall_vort500mb', region, file_tag, target_proj, proj_name)
    vorticity_plotter.render(data, cmap, norm)
    
    # Set image annotation text
    forecast_hours = f'000 Forecast Hours'
    forecast_hours = f'{forecast_hours}\n INIT: 20230529_12z'
    forecast_p_tag = 'GEOS ' + 'f5295_fp'
    forecast_str = f'{forecast_hours}\n {forecast_p_tag}'
    date_index = f'{2023}-{"05"}-{29} {12}:{"00"}Z'
    date_str = '2023 May 29'
    time_index = '08:00am EDT Monday'
    date_index = f'{date_index}\n{date_str}\n{time_index}'

    satellite = ['geos', 'nsper']
    mode = 'dark' if proj_name in satellite else 'light'

    # Annotate final image
    annotate(f'tmp/{proj_name}-vort500mb-{file_tag}.png', 'plotall_vort500mb', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')

    # Record region time
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

# Evalute total compute time and average plotting time by projection type
print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    if times[proj]:
        t_avg = sum(times[proj]) / len(times[proj]) if len(times[proj]) else 0
        print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')