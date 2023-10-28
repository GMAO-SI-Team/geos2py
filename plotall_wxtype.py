import os
import json
import time
import argparse
import numpy as np
import cartopy.crs as ccrs
import matplotlib.colors as colors

from plotting.weather import get_contour_colors, get_contour_levels, plot, annotate
from processing.regridding import congrid
from processing.smoothing import savitzky_golay2d
from datasets import loading
from utils.directories.tiles import find_tile_file

t0 = time.time()

# parser = argparse.ArgumentParser(description='Process date specified data')
# parser.add_argument('year')
# parser.add_argument('month')
# parser.add_argument('day')
# parser.add_argument('hour')
# parser.add_argument('minute')
# parser.add_argument('tag')
# parser.add_argument('s_tag')
# parser.add_argument('data_path')
# parser.add_argument('stream')
# parser.add_argument('--region')
# parser.add_argument('--f_date')
# parser.add_argument('--use_test_dirs')
# args = parser.parse_args()

# year = args.year
# month = args.month
# day = args.day
# hour = args.hour
# minute = args.minute

# region = int(args.region) if args.region else -1
# tag = args.tag if args.tag else os.getenv('MAIN_TAG')
# f_date = args.f_date if args.f_date else f'{year}{month}{day}_{hour}z'
# s_tag = args.s_tag if args.s_tag else f'f5295_fp-{f_date}'
# stream = args.stream if args.stream else os.getenv('STREAM')
# override = True if args.use_test_dirs else False

levels, slp_levels, t2m_levels = get_contour_levels()
rgb_snow, rgb_ice, rgb_frzr, rgb_rain = get_contour_colors('plotall_wxtype')

asmfile = 'data/GEOS.fp.asm.inst3_3d_asm_Np.20230529_0000.V01.nc4'
flxfile = 'data/GEOS.fp.asm.tavg1_2d_flx_Nx.20230529_0030.V01.nc4'
slvfile = 'data/GEOS.fp.asm.tavg1_2d_slv_Nx.20230529_0030.V01.nc4'
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
    # pixel = wxtype in ['PHIS', 'PRECSNO', 'PRECTOT']
    print(f'Reading {wxtype}')
    cube = loading.load_cube(
        wxtypes[wxtype]['data_file'], wxtype, 1e15, no_map_set=True, scale_factor=wxtypes[wxtype]['scale_factor']
    )
    data = loading.clean_data(cube.data, undef=1e15)
    data_wxtypes.append(data)
    print()

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

snow_cmap = colors.ListedColormap(rgb_snow)
snow_norm = colors.BoundaryNorm(levels, snow_cmap.N)
ice_cmap = colors.ListedColormap(rgb_ice)
ice_norm = colors.BoundaryNorm(levels, ice_cmap.N)
frzr_cmap = colors.ListedColormap(rgb_frzr)
frzr_norm = colors.BoundaryNorm(levels, frzr_cmap.N)
rain_cmap = colors.ListedColormap(rgb_rain)
rain_norm = colors.BoundaryNorm(levels, rain_cmap.N)
cmaps = [snow_cmap, ice_cmap, frzr_cmap, rain_cmap]
norms = [snow_norm, ice_norm, frzr_norm, rain_norm]

snow = congrid(snow, (2760, 5760), center=True)
ice = congrid(ice, (2760, 5760), center=True)
frzr = congrid(frzr, (2760, 5760), center=True)
rain = congrid(rain, (2760, 5760), center=True)
slp = congrid(slp, (2760, 5760), center=True)
t2m = congrid(t2m, (2760, 5760), center=True)

t2m = savitzky_golay2d(t2m, int(5760 * 0.025) + 1, 2)
slp = savitzky_golay2d(slp, int(5760 * 0.005) + 1, 2)

data = [snow, ice, frzr, rain, t2m, slp]
levels = [levels, slp_levels, t2m_levels]

with open('regions.json', 'r') as infile:
    regions_dict = json.load(infile)

# regions = regions_dict[region] if region in (0, -1) else [region]
regions = [50, 89, 53, 51, 49, 65, 66, 67, 68, 69, 70, 73, 74, 75,34, 35, 42, 43, 44, 46, 47, 71, 72]

times = {}

for region in regions:
    region = str(region)
    file_tag = regions_dict[region]['file_tag']
    lon_cen, lat_cen = regions_dict[region]['center']
    lon_beg, lon_end, lat_beg, lat_end = regions_dict[region]['extent'] if 'extent' in regions_dict[region] else (-180, 180, -90, 90)
    proj_name = regions_dict[region]['proj'] if not file_tag in ['australia_mapset', 'southamerica_mapset'] else 'sub'
    projs = {
        'subp': ccrs.PlateCarree(),
        'ortho': ccrs.Orthographic(lon_cen, lat_cen),
        'laea': ccrs.LambertAzimuthalEqualArea(lon_cen, lat_cen),
        'geos': ccrs.Geostationary(lon_cen),
        'nsper': ccrs.NearsidePerspective(lon_cen, lat_cen)
    }
    target_proj = projs[proj_name]
    for proj in projs:
        times[proj] = [] if proj not in times else times[proj]

    dt0 = time.time()
    plot(region, file_tag, target_proj, proj_name, data, cmaps, norms, levels)
    # test
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

    annotate(f'tmp/{proj_name}-wxtype-{file_tag}.png', 'wxtype', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    t_avg = sum(times[proj]) / len(times[proj])
    print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')