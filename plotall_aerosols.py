import time
import json
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import matplotlib.colors as colors

from datasets import loading
from processing.regridding import congrid
from plotting.aerosols import get_contour_levels, get_colormap, plot
from epilogue.annotation import annotate

t0 = time.time()

aerfile = 'data/GEOS.fp.fcst.inst1_2d_hwl_Nx.20231011_12+20231013_0000.V01.nc4'

aerosol_types = ['SSEXTTAU', 'DUEXTTAU', 'OCEXTTAU', 'BCEXTTAU', 'SUEXTTAU', 'NIEXTTAU']
aerosols_data = []
for aerosol in aerosol_types:
    print(f'Reading {aerosol}')
    cube = loading.load_cube(aerfile, aerosol, 1e15, no_map_set=True)
    data = loading.clean_data(cube.data, undef=1e15)
    print(data.shape)
    data = congrid(data, (2880, 1441), center=True)
    data = congrid(data, (2760, 5760), center=True)
    data = np.ma.masked_where(data < 0.05, data)
    aerosols_data.append(data)
    print()

clevs, alevs = get_contour_levels()
cmaps = []
norms = []
for aerosol in aerosol_types:
    color = get_colormap(aerosol)
    cmap = mpl.colormaps[color] if type(color) == str else colors.ListedColormap(color)
    norm = colors.BoundaryNorm(clevs, cmap.N)
    cmaps.append(cmap)
    norms.append(norm)

with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

regions = [50, 89, 53, 51, 49, 65, 66, 67, 68, 69, 70, 73, 74, 75, 34, 35, 42, 43, 44, 46, 47, 71, 72]

times = {}

for region in regions:
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

    dt0 = time.time()
    plot(region, file_tag, target_proj, proj_name, aerosols_data, cmaps, norms)
    # test
    forecast_hours = f'000 Forecast Hours'
    forecast_hours = f'{forecast_hours}\n INIT: 20231013_12z'
    forecast_p_tag = 'GEOS ' + 'f5295_fp'
    forecast_str = f'{forecast_hours}\n {forecast_p_tag}'
    date_index = f'{2023}-{10}-{13} {12}:{"00"}Z'
    date_str = '2023 Oct 11'
    time_index = '12:00am EDT Monday'
    date_index = f'{date_index}\n{date_str}\n{time_index}'

    satellite = ['geos', 'nsper']
    mode = 'dark' if proj_name in satellite else 'light'

    annotate(f'tmp/{proj_name}-aerosols-{file_tag}.png', 'aerosols', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    t_avg = sum(times[proj]) / len(times[proj]) if len(times[proj]) else 0
    print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')