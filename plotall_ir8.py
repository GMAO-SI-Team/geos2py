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
# parser.add_argument('--use_test_dirs')
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
# override = True if args.use_test_dirs else False
# data_dir = 'data/GEOS.fp.fcst.inst_30mn_met_c0720sfc.20230529_12+20230529_1200.V01.nc4' #if override else args.data_dir
data_dir = f'{data_dir}/GEOS.fp.fcst.inst_30mn_met_c0720sfc.{f_date[:-1]}+{f_date[:-1]}{minute}.V01.nc4'

ir8_cmap = Colormap('plotall_ir8', 'TBISCCP')
cmap = ir8_cmap.cmap
norm = ir8_cmap.norm

cube = loading.load_cube(data_dir, 'TBISCCP', 1e15, no_map_set=True)
# tile_file = find_tile_file(720, 720) if not override else 'data/DC2880xPC1441_CF0720x6C.bin'
if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
    tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
else:
    tile_file = 'tiles/DC2880xPC1441_CF0720x6C.bin'
gridspec = read_tile_file(tile_file)
data = regrid(cube.data - 273.15, method='conservative', gridspec=gridspec, undef=1e15)
data = loading.clean_data(data, undef=1e15)
data = congrid(data, (2760, 5760), center=True)

with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

regions = region_info[region] if region in ('0', '-1') else [region]
# regions = [50, 89, 53, 51, 49, 65, 66, 67, 68, 69, 70, 73, 74, 75, 34, 35, 42, 43, 44, 46, 47, 71, 72]
# regions = [34, 35, 42, 43, 44, 46, 47, 71, 72]

# import matplotlib.pyplot as plt

# fig, axes = plt.subplots(nrows=2, subplot_kw={'projection': ccrs.PlateCarree()})
# lons = np.linspace(-180, 180, 5760)
# lats = np.linspace(-90, 90, 2760)
# lons, lats = np.meshgrid(lons, lats)
# axes[0].imshow(data, cmap=cmap, norm=norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
# axes[1].pcolormesh(lons, lats, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
# plt.show()

# import sys

# sys.exit()

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
    ir8_plotter = Plotter('plotall_ir8', region, file_tag, target_proj, proj_name)
    ir8_plotter.render(data, cmap, norm)
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

    annotate(f'tmp/{proj_name}-ir8-{file_tag}.png', 'plotall_ir8', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    if times[proj]:
        t_avg = sum(times[proj]) / len(times[proj])
        print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')