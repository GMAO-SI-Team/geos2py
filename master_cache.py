import os
import time
import json
import pickle
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import PIL.Image as image
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from pyogrio import read_dataframe

# Parse command line arguments
# plot_type: specifies plots to cache (ex: plotall_ir8)
# region: specifies region to cache for specified plot type (ex: usa_mapset)
# element (optional): accepts base_image (natural Earth backgrounds) and features (coastlines, borders, etc.) 
parser = argparse.ArgumentParser(description='Get base image information')
parser.add_argument('plot_type')
parser.add_argument('region')
parser.add_argument('--element')
args = parser.parse_args()

plot_type = args.plot_type
region = args.region
element = args.element if args.element and args.element in ('features', 'base_image') else None

regions = [str(x) for x in ['globe', 50, 89, 53, 51, 49, 65, 66, 67, 68, 69, 70, 73, 74, 75, 34, 35, 42, 43, 44, 46, 47, 71, 72]] if region == 'all' else [str(region)]
plot_types = ['plotall_ir8', 'plotall_wxtype', 'plotall_radar', 'plotall_aerosols', 'plotall_precrain', 'plotall_precsnow', 'plotall_slp', 'plotall_t2m', 'plotall_tpw', 'plotall_cape', 'plotall_vort500mb', 'plotall_winds10m'] if plot_type == 'all' else [plot_type]

# Allows for large images
image.MAX_IMAGE_PIXELS = None

with open('regions.json', 'r') as infile:
    region_info = json.load(infile)

# Designates paths for the base images and shapefiles
shp_path = '/discover/nobackup/qcambrel/gmaopy/SHAPE_FILES'
base_images = '/discover/nobackup/qcambrel/gmaopy/base_images'

def get_coastlines() -> gpd.GeoDataFrame:
    """
    Get coastlines geometries.

    Returns:
    - geopandas.GeoDataFrame: Coastlines geometries.
    """
    if os.path.exists('cache/coastlines'):
            with open('cache/coastlines', 'rb') as picklefile:
                coastlines = pickle.load(picklefile)
    else:
        gshhs_paths = [f'{shp_path}/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L{i}.shp' for i in range(1, 7)]
        shapefile_l1 = read_dataframe(gshhs_paths[0])
        shapefile_l2 = read_dataframe(gshhs_paths[1])
        shapefile_l3 = read_dataframe(gshhs_paths[2])
        shapefile_l4 = read_dataframe(gshhs_paths[3])
        shapefile_l5 = read_dataframe(gshhs_paths[4])
        shapefile_l6 = read_dataframe(gshhs_paths[5])
        shapefiles = [shapefile_l1, shapefile_l2, shapefile_l3, shapefile_l4, shapefile_l5, shapefile_l6]
        coastlines = gpd.GeoDataFrame(pd.concat(shapefiles, ignore_index=True))
        with open('cache/coastlines', 'wb') as picklefile:
            print('serializing coastlines')
            pickle.dump(coastlines, picklefile)
    return coastlines

def get_states() -> gpd.GeoDataFrame:
    """
    Get states geometries.

    Returns:
    - geopandas.GeoDataFrame: States geometries.
    """
    if os.path.exists('cache/states'):
        with open('cache/states', 'rb') as picklefile:
            states = pickle.load(picklefile)
    else:
        states = read_dataframe(f'{shp_path}/state_bounds.shp')
        with open('cache/states', 'wb') as picklefile:
            print('serializing states')
            pickle.dump(states, picklefile)
    return states

def get_countries() -> gpd.GeoDataFrame:
    """
    Get countries geometries.

    Returns:
    - geopandas.GeoDataFrame: Countries geometries.
    """
    if os.path.exists('cache/countries'):
        with open('cache/countries', 'rb') as picklefile:
            countries = pickle.load(picklefile)
    else:
        countries = read_dataframe(f'{shp_path}//WB_Adm0_boundary_lines_10m/WB_Adm0_boundary_lines_10m.shp')
        with open('cache/countries', 'wb') as picklefile:
            print('serializing countries')
            pickle.dump(countries, picklefile)
    return countries

def get_roads() -> gpd.GeoDataFrame:
    """
    Get roads geometries.

    Returns:
    - geopandas.GeoDataFrame: Roads geometries.
    """
    if os.path.exists('cache/roads'):
        with open('cache/roads', 'rb') as picklefile:
            roads = pickle.load(picklefile)
    else:
        roads = read_dataframe(f'{shp_path}/intrstat.shp')
        with open('cache/roads', 'wb') as picklefile:
            print('serializing roads')
            pickle.dump(roads, picklefile)
    return roads

def get_counties() -> gpd.GeoDataFrame:
    """
    Get counties geometries.

    Returns:
    - geopandas.GeoDataFrame: Counties geometries.
    """
    if os.path.exists('cache/counties'):
        with open('cache/counties', 'rb') as picklefile:
            counties = pickle.load(picklefile)
    else:
        counties = read_dataframe(f'{shp_path}/UScounties.shp')
        with open('cache/counties', 'wb') as picklefile:
            print('serializing counties')
            pickle.dump(counties, picklefile)
    return counties

def generate_white_blue_bg():
    """
    Generate white-blue background image for plotall_precrain and plotall_precsnow.
    """
    fig = plt.figure(dpi=1500)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ocean_color = mpl.colors.rgb2hex(np.array([190,232,255]) / 255)
    coastlines = get_coastlines()
    shore = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='white', edgecolor='black', linewidth=0)
    ax.set_facecolor(ocean_color)
    ax.add_feature(shore)
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.savefig(f'{base_images}/geodetic_white_blue_bg.png', bbox_inches='tight', pad_inches=0)
    plt.close()
    print('white blue background saved successfully')

def cache_plotall_ir8(region: str, file_tag: str, element=None):
    """
    Cache plotall_ir8 image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        plt.axis('off')
        plt.savefig(f'cache/cb_ir_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    
def cache_plotall_wxtype(region, file_tag, element=None):
    """
    Cache plotall_wxtype image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        plt.axis('off')
        plt.savefig(f'cache/cb_weather_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/eo_base_2020_clean_geo.3x21600x10800.jpg'
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_auto_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')
    
def cache_plotall_radar(region, file_tag, element=None):
    """
    Cache plotall_radar image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('usa_mapset', 'midatlantic_mapset'):
            ax.add_feature(roads)
        plt.axis('off')
        plt.savefig(f'cache/cb_radar_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/eo_base_2020_clean_geo.3x21600x10800.jpg'
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_auto_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')

def cache_plotall_aerosols(region, file_tag, element=None):
    """
    Cache plotall_aerosols image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        map_color = np.array([50, 50, 50]) / 255
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor=map_color, linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=map_color, linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=map_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        plt.axis('off')
        plt.savefig(f'cache/cb_aero_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/natural_earth_grey_16200x8100.jpeg'
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        img = plt.imread(natural_earth)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), cmap='gray', vmin=0, vmax=255)
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760), cmap='gray', vmin=0, vmax=255)
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_grey_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')

def cache_plotall_precrain(region, file_tag, element=None):
    """
    Cache plotall_precrain image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('usa_mapset', 'midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(roads)
        plt.axis('off')
        plt.savefig(f'cache/cb_rain_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/geodetic_white_blue_bg.png'
        if not os.path.exists(natural_earth):
            generate_white_blue_bg()
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_white_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')
    
def cache_plotall_precsnow(region, file_tag, element=None):
    """
    Cache plotall_precsnow image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('usa_mapset', 'midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(roads)
        plt.axis('off')
        plt.savefig(f'cache/cb_snow_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/geodetic_white_blue_bg.png'
        if not os.path.exists(natural_earth):
            generate_white_blue_bg()
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_white_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')

def cache_plotall_slp(region, file_tag, element=None):
    """
    Cache plotall_slp image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='white', linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        plt.axis('off')
        plt.savefig(f'cache/cb_slp_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')

def cache_plotall_t2m(region, file_tag, element=None):
    """
    Cache plotall_t2m image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        print('getting counties')
        counties = get_counties()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        counties = cfeature.ShapelyFeature(counties.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255 if file_tag == 'maryland_mapset' else 'black' 
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(counties)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset', 'usa_mapset'):
            ax.add_feature(roads)
        plt.axis('off')
        plt.savefig(f'cache/cb_t2m_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')

def cache_plotall_tpw(region, file_tag, element=None):
    """
    Cache plotall_tpw image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        print('getting counties')
        counties = get_counties()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        counties = cfeature.ShapelyFeature(counties.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255 # needs to be changed
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset', 'usa_mapset'):
            ax.add_feature(roads)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(counties)
        plt.axis('off')
        plt.savefig(f'cache/cb_tpw_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')

def cache_plotall_cape(region, file_tag, element=None):
    """
    Cache plotall_cape image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255 # need to change
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('usa_mapset', 'midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(roads)
        plt.axis('off')
        plt.savefig(f'cache/cb_cape_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/geodetic_white_blue_bg.png'
        if not os.path.exists(natural_earth):
            generate_white_blue_bg()
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        plt.savefig(f'cache/natural_earth_white_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')

def cache_plotall_vort500mb(region, file_tag, element=None):
    """
    Cache plotall_vort500mb image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        plt.axis('off')
        plt.savefig(f'cache/cb_vort_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')
    if element in ('base_image', None):
        print(f'starting {file_tag} natural earth base image')
        natural_earth = f'{base_images}/eo_base_2020_clean_geo.3x21600x10800.jpg'
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        img = plt.imread(natural_earth)
        if file_tag in ('australia_mapset', 'southamerica_mapset', 'globe'):
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        else:
            ax.imshow(img, extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree(), regrid_shape=(2760, 5760))
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        plt.axis('off')
        plt.savefig(f'cache/natural_earth_auto_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} base image saved successfully')

def cache_plotall_winds10m(region, file_tag, element=None):
    """
    Cache plotall_winds10m image features for a specified region.

    Parameters:
    - region (str): Region identifier.
    - file_tag (str): File tag.
    - element (str, optional): Element to cache. Defaults to None.
    """
    proj = region_info[region]['proj']
    center = region_info[region]['center']
    if proj == 'sub':
        target_proj = ccrs.PlateCarree()
    if proj == 'ortho':
        target_proj = ccrs.Orthographic(center[0], center[1])
    if proj == 'laea':
        file_tag = region_info[region]['file_tag']
        target_proj = ccrs.PlateCarree() if file_tag in ('australia_mapset', 'southamerica_mapset') else ccrs.LambertAzimuthalEqualArea(center[0], center[1])
    if proj == 'geos':
        target_proj = ccrs.Geostationary(center[0])
    if proj == 'nsper':
        target_proj = ccrs.NearsidePerspective(center[0], center[1])
    if element in ('features', None):
        print('getting coastlines')
        coastlines = get_coastlines()
        print('getting states')
        states = get_states()
        print('getting countries')
        countries = get_countries()
        print('getting roads')
        roads = get_roads()
        print('getting counties')
        counties = get_counties()
        gshhs = cfeature.ShapelyFeature(coastlines.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        states = cfeature.ShapelyFeature(states['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        countries = cfeature.ShapelyFeature(countries['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        counties = cfeature.ShapelyFeature(counties.boundary, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.125)
        roads_color = np.array([125, 15, 15]) / 255 # needs to be changed
        roads = cfeature.ShapelyFeature(roads['geometry'], ccrs.PlateCarree(), facecolor='none', edgecolor=roads_color, linewidth=0.125)
        print('plotting')
        fig = plt.figure(dpi=1500)
        ax = plt.axes(projection=target_proj)
        if proj in ('ortho', 'laea'):
            ax.set_extent(region_info[region]['extent'], ccrs.PlateCarree())
        ax.add_feature(gshhs)
        ax.add_feature(states)
        ax.add_feature(countries)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset', 'usa_mapset'):
            ax.add_feature(roads)
        if file_tag in ('midatlantic_mapset', 'maryland_mapset'):
            ax.add_feature(counties)
        plt.axis('off')
        plt.savefig(f'cache/cb_wind_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        print(f'{file_tag} features saved successfully')

# Dictionary to map caching functions to plot type
cache_dict = {
    'plotall_ir8': cache_plotall_ir8,
    'plotall_wxtype': cache_plotall_wxtype,
    'plotall_radar': cache_plotall_radar,
    'plotall_aerosols': cache_plotall_aerosols,
    'plotall_precrain': cache_plotall_precrain,
    'plotall_precsnow': cache_plotall_precsnow,
    'plotall_slp': cache_plotall_slp,
    'plotall_t2m': cache_plotall_t2m,
    'plotall_tpw': cache_plotall_tpw,
    'plotall_cape': cache_plotall_cape,
    'plotall_vort500mb': cache_plotall_vort500mb,
    'plotall_winds10m': cache_plotall_winds10m
}

# Start timer
t0 = time.time()

for plot_type in plot_types:
    for region in regions:
        file_tag = region_info[region]['file_tag']
        print(f'caching {plot_type} assets for {file_tag}')
        cache_dict[plot_type](region, file_tag, element)

print('caching complete')

# Evalute total run time
runtime = time.time() - t0
print(f'took {runtime // 3600} hours {runtime % 3600 // 60} minutes {runtime % 3600 % 60} seconds')