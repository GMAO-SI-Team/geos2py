import os
import json
import time
import pickle
import calendar
import argparse
import numpy as np
import cartopy.crs as ccrs
import matplotlib.cm as cm

from plotting.colormaps import Colormap
from plotting.plots import Plotter
from processing.regridding import congrid, regrid, read_tile_file
from processing.smoothing import savitzky_golay2d
from processing.scaling import bytscl
from datasets import loading
from epilogue.annotation import annotate

# Start timer
t0 = time.time()

# Parse command line arguments to specify data source and plot type
parser = argparse.ArgumentParser(description='Process date specified data and plot type')
parser.add_argument('year')
parser.add_argument('month')
parser.add_argument('day')
parser.add_argument('hour')
parser.add_argument('minute')
parser.add_argument('tag')
parser.add_argument('s_tag')
parser.add_argument('data_dir')
parser.add_argument('stream')
parser.add_argument('plot_type', choices=['aerosols', 'cape', 'ir8', 'precrain'])
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
plot_type = args.plot_type

# Load regions JSON file
with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

regions = region_info[region] if region in ('0', '-1') else [region]

times = {}



def plot_aerosols():
    # File path for aerosols data
    data_dir_aerosols = f'{data_dir}/GEOS.fp.fcst.inst1_2d_hwl_Nx.{f_date[:-1]}+20231013_0000.V01.nc4'
    
    # Load aerosol data
    ss = loading.load_cube(data_dir_aerosols, 'SSEXTTAU', 1e15, no_map_set=True).data
    du = loading.load_cube(data_dir_aerosols, 'DUEXTTAU', 1e15, no_map_set=True).data
    oc = loading.load_cube(data_dir_aerosols, 'OCEXTTAU', 1e15, no_map_set=True).data
    bc = loading.load_cube(data_dir_aerosols, 'BCEXTTAU', 1e15, no_map_set=True).data
    su = loading.load_cube(data_dir_aerosols, 'SUEXTTAU', 1e15, no_map_set=True).data
    ni = loading.load_cube(data_dir_aerosols, 'NIEXTTAU', 1e15, no_map_set=True).data
    
    # Resample the aerosols data to 5760x2760 shape
    aerosols = [congrid(aerosol, (2760, 5760), method='nearest', center=True) for aerosol in [ss, du, oc, bc, su, ni]]
    print(f'finished regridding: {time.time() - t0} seconds')

    # Define colormap and normalization for each aerosol
    aerosol_types = ['SSEXTTAU', 'DUEXTTAU', 'OCEXTTAU', 'BCEXTTAU', 'SUEXTTAU', 'NIEXTTAU']
    cmaps = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).cmap for aerosol in aerosol_types]
    norms = [Colormap('plotall_aerosols', aerosol, data_min=0, data_max=0.5).norm for aerosol in aerosol_types]
    
    # Convert scalar mappable objects to rgba arrays then scale alpha channels for blending
    data = []
    for aerosol, cmap, norm in zip(aerosols, cmaps, norms):
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        blend = sm.to_rgba(aerosol)
        blend[:, :, 3] = bytscl(aerosol, low=0, high=0.25)
        data.append(blend)
    
    plot_data(data, cmaps, norms, 'plotall_aerosols')

    

def plot_cape():
    # File path for CAPE data
    data_dir_cape = f'{data_dir}/GEOS.fp.fcst.inst3_2d_met_Nx.{f_date[:-1]}+20231224_0900.V01.nc4'
    
    # Define the CAPE colormap and normalization
    cape_cmap = Colormap('plotall_cape', 'CAPE')
    cmap = cape_cmap.cmap
    norm = cape_cmap.norm

    # Load and extract the CAPE data
    cube = loading.load_cube(data_dir_cape, 'CAPE', 1e15, pixel=True, no_project=True)
    data = loading.clean_data(cube.data, undef=1e15)

    # Resample data to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)

    # Mask data values below 100 (the lower bound)
    mask = data < 100
    data = np.ma.masked_where(mask, data)
    data = np.ma.masked_where(data > 10000, data)
    
    plot_data([data], [cmap], [norm], 'plotall_cape')


def plot_ir8():
    
    #
    data_dir_ir8 = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

    # Define the infrared colormap and normalization
    ir8_cmap = Colormap('plotall_ir8', 'TBISCCP')
    cmap = ir8_cmap.cmap
    norm = ir8_cmap.norm

    # Read and extract the data
    cube = loading.load_cube(data_dir_ir8, 'TBISCCP', 1e15, no_map_set=True)


    
    if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
        tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
    else:
        tile_file = 'tiles/DC2880xPC1441_CF0720x6C.bin'
    gridspec = read_tile_file(tile_file)

    # Convert Kelvin to Celsius
    data = cube.data - 273.15

    # First order conservative regridding using 4320x720 to 2880x1441 gridspec
    data = regrid(data, method='conservative', gridspec=gridspec, undef=1e15)

    # Resample data to 5760x2760 shape
    data = loading.clean_data(data, undef=1e15)
    data = congrid(data, (2760, 5760), center=True) 

    plot_data(data, cmap, norm, 'plotall_ir8')


def plot_precrain():
    # File path for rain data
    data_dir_precrain = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'
    
    # Define zero array for summing accumulated data
    saved_data = f'tmp/raindata-{f_date}.pkl'
    if os.path.exists(saved_data):
        with open(saved_data, 'rb') as pkl:
            acc_data = pickle.load(pkl)
    else:
        acc_data = np.zeros((2760, 5760))

    # Define the rain colormap and normalization
    rain_cmap = Colormap('plotall_precrain', 'PRECTOT')
    cmap = rain_cmap.cmap
    norm = rain_cmap.norm

    # Load and extract the rain data
    cube = loading.load_cube(data_dir_precrain, 'PRECTOT', 1e15, no_map_set=True).data
    acc_data += congrid(cube, (2760, 5760), center=True)
    
    plot_data(acc_data, cmap, norm, 'plotall_precrain')


def plot_precsnow():
    # File path for snow data
    data_dir_precsnow = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

    # Define zero array for summing accumulated data
    saved_data = f'tmp/snowdata-{f_date}.pkl'
    if os.path.exists(saved_data):
        with open(saved_data, 'rb') as pkl:
            acc_data = pickle.load(pkl)
    else:
        acc_data = np.zeros((2760, 5760))

    # Define the snow colormap and normalization
    snow_cmap = Colormap('plotall_precsnow', 'PRECSNO')
    cmap = snow_cmap.cmap
    norm = snow_cmap.norm

    # Load and extract the data
    cube = loading.load_cube(data_dir_precsnow, 'PRECSNO', 1e15, no_map_set=True, pixel=True, scale_factor=3600/25.4)
    data = loading.clean_data(cube.data, undef=1e15)

    # Resample data to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)

    # Mask data values below 0.1 (the lower bound)
    data = np.ma.masked_where(data < 0.1, data)

    # Sum historic data and current data
    acc_data = acc_data + data
    plot_data(acc_data, cmap, norm, 'plotall_precsnow')


def plot_t2m():
    data_dir_t2m = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

    # Define the temperature colormaps and normalization
    temp_cmap = Colormap('plotall_t2m', 'T2M')
    cmap = temp_cmap.cmap
    norm = temp_cmap.norm

    # Read and extract the data
    data = loading.load_cube(data_dir_t2m, 'T2M', 1e15, no_map_set=True).data
    data = (data - 273.15) * 1.8 + 32 # convert Kelvin to Farenheit

    # Resample data to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)
    plot_data(data, cmap, norm, 'plotall_t2m')


def plot_tpw():
    data_dir_tpw = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

    # Define the water colormap and normalization
    tpw_cmap = Colormap('plotall_tpw', 'TQV')
    cmap = tpw_cmap.cmap
    norm = tpw_cmap.norm

    # Read and extract the data
    cube = loading.load_cube(data_dir_tpw, 'TQV', 1e15, no_map_set=True, pixel=True)
    data = loading.clean_data(cube.data, undef=1e15)

    # Resample data to 5760x2760 shape
    data = congrid(data, (2760, 5760), center=True)
    plot_data(data, cmap, norm, 'plotall_tpw')


def plot_vort500mb():
    data_dir_vort = f'{data_dir}/GEOS.fp.asm.tavg1_2d_slv_Nx.{f_date[:-1]}30.V01.nc4'

    # Define the vorticity colormap and normalization
    vorticity_cmap = Colormap('plotall_vort500mb', 'vort')
    cmap = vorticity_cmap.cmap
    norm = vorticity_cmap.norm

    # Read and extract the data
    heights = loading.load_cube(data_dir_vort, 'H500', 1e15, no_map_set=True).data
    u = loading.load_cube(data_dir_vort, 'U500', 1e15, no_map_set=True).data
    v = loading.load_cube(data_dir_vort, 'V500', 1e15, no_map_set=True).data

    # Compute forces and derivatives
    rad = np.pi / 180.0
    re = 6371220.0
    lats = 2.0 * (np.arange(720) / (720 - 1) - 0.5) * 90 * rad
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

    plot_data(data, cmap, norm, 'plotall_vort500mb')


def plot_winds10m():
    data_dir_winds = f'{data_dir}/GEOS.fp.asm.tavg1_2d_flx_Nx.{f_date[:-1]}30.V01.nc4'

    # Define the wind colormap and normalization
    wind_cmap = Colormap('plotall_winds10m', 'spd')
    cmap = wind_cmap.cmap
    norm = wind_cmap.norm

    # Read, extract, and scale the data
    vectors = {
        'ULML': {
            'data_file': data_dir_winds,
            'scale_factor': 2.23694
        },
        'VLML': {
            'data_file': data_dir_winds,
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
    plot_data(data, cmap, norm, 'plotall_winds10m')


def plot_radar():
    data_dir_radar = f'{data_dir}/GEOS.fp.fcst.inst_30mn_met_c0720sfc.{f_date[:-1]}+{f_date[:-1]}{minute}.V01.nc4'

    # Define the radar colormap and normalization
    radar_cmap = Colormap('plotall_radar', 'DBZ_MAX')
    cmap = radar_cmap.cmap
    norm = radar_cmap.norm

    # Read and extract the data
    cube = loading.load_cube(data_dir_radar, 'DBZ_MAX', 1e15, no_map_set=True)
    if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
        tile_file = '/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC2880xPC1441_CF0720x6C.bin'
    else:
        tile_file = 'tiles/DC2880xPC1441_CF0720x6C.bin'
    gridspec = read_tile_file(tile_file)

    # First order conservative regridding using 4320x720 to 2880x1441 gridspec
    data = regrid(cube.data, method='conservative', gridspec=gridspec, undef=1e15)
    data = loading.clean_data(data, undef=1e15)
    data = congrid(data, (2760, 5760), center=True)

    # Mask data values below 5 (the lower bound)
    mask = data < 5
    data = np.ma.masked_where(mask, data)

    plot_data(data, cmap, norm, 'plotall_radar')


def plot_wxtype():
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

    plot_data(data, cmap, norm, 'plotall_wxtype')


def plot_slp():
    print("Not yet implemented")



def plot_data(data, cmap, norm, plot_tag):

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

        # Start regional timer
        dt0 = time.time()

        # Initialize plotter with tag and plot data
        plotter = Plotter(plot_tag, region, file_tag, target_proj, proj_name)
        plotter.render(data, cmap, norm)
        
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
        annotate(f'tmp/{proj_name}-ir8-{file_tag}.png', 'plotall_ir8', mode=mode, forecast=forecast_str, date=date_index)
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






# Main routine to call appropriate plotting function based on specified plot type
match plot_type:
    case 'aerosols':
        plot_aerosols()
    case 'cape':
        plot_cape()
    case 'ir8':
        plot_ir8()
    case 'precrain':
        plot_precrain()
    case 'precsnow':
        plot_precsnow()
    case 't2m':
        plot_t2m()
    case 'tpw':
        plot_tpw()
    case 'vort500mb':
        plot_vort500mb()
    case 'winds10m':
        plot_winds10m()
    case 'radar':
        plot_radar()
    case 'wxtype':
        plot_wxtype()
    case 'slp':
        plot_slp()
    case _:
        print("Invalid plot type")





# End timer
t1 = time.time()
total = t1 - t0
print(f'Total runtime: {total} seconds')
