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
from datasets import loading
from epilogue.annotation import annotate

# Allows for large images
image.MAX_IMAGE_PIXELS = None

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

data_dir = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

# Define the wind colormap and normalization
wind_cmap = Colormap('plotall_winds10m', 'spd')
cmap = wind_cmap.cmap
norm = wind_cmap.norm

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
    # pixel = wxtype in ['PHIS', 'PRECSNO', 'PRECTOT']
    print(f'Reading {vector.lower()}')
    cube = loading.load_cube(
        vectors[vector]['data_file'], vector, 1e15, pixel=True, scale_factor=vectors[vector]['scale_factor']
    )
    data = loading.clean_data(cube.data, undef=1e15)
    wind_data.append(data)
    print()

# Compute vector magnitude
ulml, vlml = wind_data
data = np.sqrt(ulml ** 2 + vlml ** 2)

# Resample data to 5760x2760 shape
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

    # cities = 'cities/all_cities.txt' if file_tag in ('midatlantic_mapset', 'maryland_mapset') else 'cities/world_cities.csv'
    
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

    # Start regional timer
    dt0 = time.time()

    # Initialize wind plotter and plot data
    wind_plotter = Plotter('plotall_winds10m', region, file_tag, target_proj, proj_name, label_coords=city_coords, interpolator=city_interpolator)
    wind_plotter.render(data, cmap, norm)
    
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
    annotate(f'tmp/{proj_name}-winds10m-{file_tag}.png', 'plotall_winds10m', mode=mode, forecast=forecast_str, date=date_index)
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