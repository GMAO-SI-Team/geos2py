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
from processing.smoothing import ideal_bandpass_filter
from datasets import loading
from epilogue.annotation import annotate

image.MAX_IMAGE_PIXELS = None

def find_storm_minima(slp, x_window, y_window, threshold=1004, smoothing_band=0.8, fuzz_threshold=1e-8):
    print(f'SLP (max, min): {slp.max()} {slp.min()}')

    slp_smoothed = ideal_bandpass_filter(slp, 0, smoothing_band)
    print(f'SLP smoothed (max, min): {slp_smoothed.max()} {slp_smoothed.min()}')

    n_lats, n_lons = slp.shape
    storm = np.zeros_like(slp)
    storm_minima_mask = find_storm_minima_mask(slp_smoothed, threshold, fuzz_threshold)

    count = len(np.where(storm_minima_mask)[0])
    storm[storm_minima_mask] = 1
    print(f'Storm max: {storm.max()}')
    print(f'There are {count} points with storm minima (slpSMin)=1')

    count_storm4 = label_storm_minima(storm, slp, slp_smoothed, x_window, y_window)

    ilocs = get_storm_minima_locations(storm, n_lons, n_lats, count_storm4)
    return ilocs

def find_storm_minima_mask(slp_smoothed, threshold, fuzz_threshold):
    mask = np.logical_and.reduce((
        slp_smoothed < np.roll(slp_smoothed, 1, axis=1),
        slp_smoothed < np.roll(slp_smoothed, -1, axis=1),
        slp_smoothed < np.roll(slp_smoothed, 1, axis=0),
        slp_smoothed < np.roll(slp_smoothed, -1, axis=0)
    ))

    fzz = slp_smoothed[mask]
    wfzz = np.where(fzz < fuzz_threshold)

    if len(wfzz[0]):
        mask[np.where(mask)[wfzz]] = False
    
    return mask

def label_storm_minima(storm, slp, slp_smoothed, x_window, y_window):
    count_storm4 = 0
    storm_mask = np.where(storm == 1)
    slp_idx = np.argsort(slp[storm_mask])
    storm_mask = storm_mask[0][slp_idx]
    n_lats, n_lons = slp_smoothed.shape

    print(f'storm mask (max, min): {storm_mask.max()} {storm_mask.min()}')

    if len(storm_mask):
        x_dist = 10 / x_window
        y_dist = 10 / y_window
        for n in range(len(storm_mask)):
            i = int(storm_mask[n] % n_lons)
            j = int(storm_mask[n] / n_lons)
            istart = int(max(i - x_dist, 0))
            iend = int(min(i + x_dist, n_lons - 1))
            jstart = int(max(j - y_dist, 0))
            jend = int(min(j + y_dist, n_lats - 1))

            storm_local = storm[jstart:jend, istart:iend]
            slp_local = slp[jstart:jend, istart:iend]
            icount = np.where(storm_local == 4)

            if n < 15:
                print(f'storm mask location (n={n}): {storm_mask[n]}')

            if len(icount[0]) == 0:
                imark = np.where(storm_local == 1)
                if len(imark[0]):
                    local_idx = slp_local[imark].argsort()
                    imark = imark[0][local_idx]
                    storm_local[imark] = 4
                    count_storm4 += 1
                else:
                    storm_local = 0
                storm[jstart:jend, istart:iend] = storm_local
    
    print(f'There are {count_storm4} points with storm=4')

    return count_storm4

def get_storm_minima_locations(storm, n_lons, n_lats, count_storm4):
    ilocs = np.zeros((count_storm4, 2))
    ij = 0

    for j in range(n_lats):
        for i in range(n_lons):
            if storm[j, i] == 4:
                ilocs[ij, 0] = i
                ilocs[ij, 1] = j
                ij += 1

    return ilocs

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

data_dir = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

pressure_cmap = Colormap('plotall_slp', 'SLP')
cmap = pressure_cmap.cmap
norm = pressure_cmap.norm

data = loading.load_cube(data_dir, 'SLP', 1e15, no_map_set=True, scale_factor=0.01).data
# data = (data - 273.15) * 1.8 + 32 
data = congrid(data, (2760, 5760), center=True)

# mask = data < 10
# data = np.ma.masked_where(mask, data)

with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

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

    x_window = (lon_end - lon_beg) / 5760
    y_window = (lat_end - lat_beg) / (2760 - 1)
    slp_mins = find_storm_minima(data, x_window, y_window)
    print(f'SLP minima array shape: {slp_mins.shape}')
    # slp_maxs = find_slp_maxs(data, data_lons, data_lats)
    slp_min_locs = slp_max_locs = np.array((2760, 3))
    for s in range(2760):
        slp_min_locs[s, 0] = slp_mins[s, 0] * x_window + 0.5 * x_window - 180
        slp_min_locs[s, 1] = slp_mins[s, 1] * y_window - 90
        slp_min_locs[s, 2] = data[slp_mins[s, 1], slp_mins[s, 1]]
        print(slp_min_locs[s, :])
        # slp_max_locs[s, 0] = slp_maxs[s, 0] * data_lons + 0.5 * data_lons - 180
        # slp_max_locs[s, 1] = slp_maxs[s, 1] * data_lats - 90
        # slp_max_locs[s, 2] = data[slp_maxs[s, 1], slp_maxs[s, 1]]
        # print(slp_max_locs[s, :])

    dt0 = time.time()
    pressure_plotter = Plotter('plotall_slp', region, file_tag, target_proj, proj_name, label_coords=(slp_min_locs, data_lons))
    pressure_plotter.render(data, cmap, norm)
    
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

    annotate(f'tmp/{proj_name}-slp-{file_tag}.png', 'plotall_slp', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    if times[proj]:
        t_avg = sum(times[proj]) / len(times[proj])
        print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')