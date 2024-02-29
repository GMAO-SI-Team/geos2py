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
from processing.regridding import congrid
from processing.scaling import bytscl
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

data_dir = f'{data_dir}/GEOS.fp.fcst.inst1_2d_hwl_Nx.{f_date[:-1]}+20231013_0000.V01.nc4' # 20231011_12

# Read and extract the data
ss = loading.load_cube(data_dir, 'SSEXTTAU', 1e15, no_map_set=True).data
du = loading.load_cube(data_dir, 'DUEXTTAU', 1e15, no_map_set=True).data
oc = loading.load_cube(data_dir, 'OCEXTTAU', 1e15, no_map_set=True).data
bc = loading.load_cube(data_dir, 'BCEXTTAU', 1e15, no_map_set=True).data
su = loading.load_cube(data_dir, 'SUEXTTAU', 1e15, no_map_set=True).data
ni = loading.load_cube(data_dir, 'NIEXTTAU', 1e15, no_map_set=True).data

# Resample the aerosols data to 5760x2760 shape
aerosols = [congrid(aerosol, (2760, 5760), method='nearest', center=True) for aerosol in [ss, du, oc, bc, su, ni]]
print(f'finished regridding: {time.time() - t0} seconds')

# Define the colormap and normalization for each aerosol
ss_cmap = Colormap('plotall_aerosols', 'SSEXTTAU', data_min=0, data_max=0.5).cmap
ss_norm = Colormap('plotall_aerosols', 'SSEXTTAU', data_min=0, data_max=0.5).norm
du_cmap = Colormap('plotall_aerosols', 'DUEXTTAU', data_min=0, data_max=0.5).cmap
du_norm = Colormap('plotall_aerosols', 'DUEXTTAU', data_min=0, data_max=0.5).norm
oc_cmap = Colormap('plotall_aerosols', 'OCEXTTAU', data_min=0, data_max=0.5).cmap
oc_norm = Colormap('plotall_aerosols', 'OCEXTTAU', data_min=0, data_max=0.5).norm
bc_cmap = Colormap('plotall_aerosols', 'BCEXTTAU', data_min=0, data_max=0.5).cmap
bc_norm = Colormap('plotall_aerosols', 'BCEXTTAU', data_min=0, data_max=0.5).norm
su_cmap = Colormap('plotall_aerosols', 'SUEXTTAU', data_min=0, data_max=0.5).cmap
su_norm = Colormap('plotall_aerosols', 'SUEXTTAU', data_min=0, data_max=0.5).norm
ni_cmap = Colormap('plotall_aerosols', 'NIEXTTAU', data_min=0, data_max=0.5).cmap
ni_norm = Colormap('plotall_aerosols', 'NIEXTTAU', data_min=0, data_max=0.5).norm

ss, du, oc, bc, su, ni = aerosols

# Define scalar mappable objects
ss_sm = cm.ScalarMappable(norm=ss_norm, cmap=ss_cmap)
du_sm = cm.ScalarMappable(norm=du_norm, cmap=du_cmap)
oc_sm = cm.ScalarMappable(norm=oc_norm, cmap=oc_cmap)
bc_sm = cm.ScalarMappable(norm=bc_norm, cmap=bc_cmap)
su_sm = cm.ScalarMappable(norm=su_norm, cmap=su_cmap)
ni_sm = cm.ScalarMappable(norm=ni_norm, cmap=ni_cmap)

# Convert scalar mappable objects to rgba arrays then scale alpha channels for blending
ss_blend = ss_sm.to_rgba(ss)
ss_blend[:, :, 3] = bytscl(ss, low=0, high=0.25)
du_blend = du_sm.to_rgba(du)
du_blend[:, :, 3] = bytscl(du, low=0, high=0.25)
oc_blend = oc_sm.to_rgba(oc)
oc_blend[:, :, 3] = bytscl(oc, low=0, high=0.25)
bc_blend = bc_sm.to_rgba(bc)
bc_blend[:, :, 3] = bytscl(bc, low=0, high=0.25)
su_blend = su_sm.to_rgba(su)
su_blend[:, :, 3] = bytscl(su, low=0, high=0.25)
ni_blend = ni_sm.to_rgba(ni)
ni_blend[:, :, 3] = bytscl(ni, low=0, high=0.25)

data = [ss_blend, du_blend, oc_blend, bc_blend, su_blend, ni_blend]
cmaps = [ss_cmap, du_cmap, oc_cmap, bc_cmap, su_cmap, ni_cmap]
norms = [ss_norm, du_norm, oc_norm, bc_norm, su_norm, ni_norm]

# Open regions JSON and set region(s)
with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

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

    # Set regional timer
    dt0 = time.time()

    # Initialize aerosols plotter and plot data
    aerosol_plotter = Plotter('plotall_aerosols', region, file_tag, target_proj, proj_name)
    aerosol_plotter.render(data, cmaps, norms)

    # Set image annotation text
    forecast_hours = f'000 Forecast Hours'
    forecast_hours = f'{forecast_hours}\n INIT: 20231013_12z'
    forecast_p_tag = 'GEOS ' + 'f5295_fp'
    forecast_str = f'{forecast_hours}\n {forecast_p_tag}'
    date_index = f'{2023}-{10}-{13} {12}:{"00"}Z'
    date_str = '2023 Oct 11'
    time_index = '12:00pm EDT Monday'
    date_index = f'{date_index}\n{date_str}\n{time_index}'

    satellite = ['geos', 'nsper']
    mode = 'dark' if proj_name in satellite else 'light'

    # Annotate final image
    annotate(f'tmp/{proj_name}-aerosols-{file_tag}.png', 'plotall_aerosols', mode=mode, forecast=forecast_str, date=date_index)
    print(f'{file_tag} saved successfully')

    # Record region time
    times[proj_name].append(time.time() - dt0)

t = time.time() - t0

# Evalute total compute time and average plotting time by projection type
print(f'total run time: {t // 3600} hours {t % 3600 // 60} minutes {t % 3600 % 60} seconds')
for proj in times:
    t_avg = sum(times[proj]) / len(times[proj]) if len(times[proj]) else 0
    print(f'average time \'{proj}\': {t_avg // 3600} hours {t_avg % 3600 // 60} minutes {t_avg % 3600 % 60} seconds')