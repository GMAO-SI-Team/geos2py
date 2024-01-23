import os
import json
import time
import calendar
import argparse
import numpy as np
import cartopy.crs as ccrs

from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import regrid, congrid, read_tile_file
from datasets import loading
from epilogue.annotation import annotate
from utils.directories.tiles import find_tile_file

n_dates = 1

t0 = time.time()
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

# if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
#     data_dir = f'{data_dir}/Y{year}/M{month}/D{day}/H{hour}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'
# else:
#     data_dir = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

data_dir = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

acc_data = np.zeros((2760, 5760))

snow_cmap = Colormap('plotall_precsnow', 'PRECSNO')
cmap = snow_cmap.cmap
norm = snow_cmap.norm

with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

for _ in range(n_dates):
    cube = loading.load_cube(data_dir, 'PRECSNO', 1e15, no_map_set=True, pixel=True, scale_factor=3600/25.4)
    data = loading.clean_data(cube.data, undef=1e15)
    data = congrid(data, (2760, 5760), center=True)
    data = np.ma.masked_where(data < 0.1, data)
    acc_data = acc_data + data

    regions = region_info[region] if region in ('0', '-1') else [region]

    times = {}

    for region in regions:
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

        dt0 = time.time()
        snow_plotter = Plotter('plotall_precsnow', region, file_tag, target_proj, proj_name)
        snow_plotter.render(acc_data, cmap, norm)
        # plot(region, file_tag, target_proj, proj_name, data, cmap, norm)
        # test
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

        annotate(f'tmp/{proj_name}-precsnow-{file_tag}.png', 'plotall_precsnow', mode=mode, forecast=forecast_str, date=date_index)
        print(f'{file_tag} saved successfully')
        times[proj_name].append(time.time() - dt0)

    t = time.time() - t0

    print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
    for proj in times:
        if times[proj]:
            t_avg = sum(times[proj]) / len(times[proj])
            print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')