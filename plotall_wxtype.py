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

# Define the weather colormaps and normalization
weather = ['snow', 'ice', 'frzr', 'rain']
weather_cmaps = [Colormap('plotall_wxtype', item) for item in weather]
cmap = [weather_cmap.cmap for weather_cmap in weather_cmaps]
norm = [weather_cmap.norm for weather_cmap in weather_cmaps]

asmfile = f'{data_dir}/GEOS.fp.asm.inst3_3d_asm_Np.{f_date[:-1]}00.V01.nc4'
flxfile = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'
slvfile = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

# Read, extract, and scale the data
wxtypes = {
    'PHIS': {
        'data_file': asmfile,
        'scale_factor': 1 / 9.81
    },
    'PRECSNO': {
        'data_file': flxfile,
        'scale_factor': 3600 / 25.4
    },
    'PRECTOT': {
        'data_file': flxfile,
        'scale_factor': 3600 / 25.4
    },
    'T2M': {
        'data_file': slvfile,
        'scale_factor': 1
    },
    'H1000': {
        'data_file': slvfile,
        'scale_factor': 1
    },
    'H500': {
        'data_file': slvfile,
        'scale_factor': 1
    },
    'SLP': {
        'data_file': slvfile,
        'scale_factor': 1 / 100
    }
}
data_wxtypes = []
for wxtype in wxtypes:
    print(f'Reading {wxtype}')
    cube = loading.load_cube(
        wxtypes[wxtype]['data_file'], wxtype, 1e15, no_map_set=True, scale_factor=wxtypes[wxtype]['scale_factor']
    )
    data = loading.clean_data(cube.data, undef=1e15)
    data_wxtypes.append(data)
    print()

# Segment weather types
phis, snow, rain, t2m, h1000, h500, slp = data_wxtypes
ice = snow * 0
frzr = snow * 0
thick = h500 - h1000
elev_factor = (phis - 305) / 915
elev_factor[np.where(elev_factor < 0)] = 0.0
elev_factor[np.where(elev_factor > 1)] = 1.0
elev_factor = elev_factor * 100
low = 5400 + elev_factor
ice[np.where((thick.any() >= low.any()) and (snow > 0))] = snow[np.where((thick.any() >= low.any()) and (snow > 0))]
snow[np.where((thick.any() >= low.any()) and (snow > 0))] = 0.0
rain[np.where(snow >= 0.1)] = 0.0
rain[np.where(ice >= 0.1)] = 0.0

# Resample data to 5760x2760 shape
snow = congrid(snow, (2760, 5760), center=True)
ice = congrid(ice, (2760, 5760), center=True)
frzr = congrid(frzr, (2760, 5760), center=True)
rain = congrid(rain, (2760, 5760), center=True)
slp = congrid(slp, (2760, 5760), center=True)
t2m = congrid(t2m, (2760, 5760), center=True)

# Smooth sea level pressure and 2-meter temperature
t2m = savitzky_golay2d(t2m, int(5760 * 0.025) + 1, 2)
slp = savitzky_golay2d(slp, int(5760 * 0.005) + 1, 2)

# Mask data values below 0.01 (the lower bound)
data = [np.ma.masked_where(wxtype < 0.01, wxtype) for wxtype in (snow, ice, frzr, rain, t2m, slp)]

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

    # Initialize weather plotter and plot data
    wxtype_plotter = Plotter('plotall_wxtype', region, file_tag, target_proj, proj_name)
    wxtype_plotter.render(data, cmap, norm)
    
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
    annotate(f'tmp/{proj_name}-wxtype-{file_tag}.png', 'plotall_wxtype', mode=mode, forecast=forecast_str, date=date_index)
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