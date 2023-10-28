import os
import csv
import json
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.interpolate as interp

from PIL import Image, ImageDraw, ImageFont

def get_contour_levels():
    levels = [0.01,0.02,0.03,0.05,0.10,0.25,0.5,1.0]
    levs1 = 884 + np.arange(12) * 8
    levs2 = 980 + np.arange(5) * 4
    levs3 = 1000 + np.arange(31) * 2
    slp_levels = np.append(np.append(levs1, levs2), levs3)
    t2m_levels = [273.15]
    return levels, slp_levels, t2m_levels

def get_contour_colors(plot_type):
    rgb_snow = [
        [172,196,225],
        [118,172,204],
        [ 67,139,186],
        [ 41,101,162],
        [ 24, 68,145],
        [ 12, 44,113],
        [ 12, 27, 64],
        [ 12, 27, 64]
    ]
    rgb_ice = [
        [242,179,251],
        [228,144,249],
        [214,107,243],
        [211, 68,248],
        [205, 40,247],
        [161, 12,209],
        [122,  7,156],
        [122,  7,156]
    ]
    rgb_frzr= [
        [253,169,195],
        [252,119,154],
        [250, 59, 95],
        [200, 16, 57],
        [161,  8, 30],
        [105,  3, 23],
        [ 58,  1, 16],
        [ 58,  1, 16]
    ]
    rgb_rain = [
        [114,198,114],
        [ 38,136, 67],
        [ 26, 81, 35],
        [ 12, 32, 18],
        [236,246, 80],
        [253,179, 57],
        [252, 98, 33],
        [252, 98, 33]
    ]
    rgb_snow = np.array([np.array(x) for x in rgb_snow]) / 255
    rgb_ice = np.array([np.array(x) for x in rgb_ice]) / 255
    rgb_frzr = np.array([np.array(x) for x in rgb_frzr]) / 255
    rgb_rain = np.array([np.array(x) for x in rgb_rain]) / 255
    return rgb_snow, rgb_ice, rgb_frzr, rgb_rain

def plot(region, file_tag, target_proj, proj_name, data, cmaps, norms, levels):
    snow, ice, frzr, rain, t2m, slp = data
    snow_cmap, ice_cmap, frzr_cmap, rain_cmap = cmaps
    snow_norm, ice_norm, frzr_norm, rain_norm = norms
    prec_levels, slp_levels, t2m_levels = levels
    natural_earth = f'cache/natural_earth_{file_tag}.png'
    fig = plt.figure(dpi=1500)
    ax = plt.axes(projection=target_proj)
    lons = np.linspace(-180, 180, 5760)
    lats = np.linspace(-90, 90, 2760)
    lons, lats = np.meshgrid(lons, lats)
    img = plt.imread(natural_earth)
    with open('extents.json', 'r') as infile:
        extents = json.load(infile)
    ax.imshow(img, extent=extents[file_tag], transform=target_proj)
    if proj_name in ['sub', 'ortho', 'laea']:
        with open('regions.json', 'r') as infile:
            regions_dict = json.load(infile)
        ax.set_extent(regions_dict[region]['extent'], ccrs.PlateCarree())
    if proj_name == 'laea':
        transform = target_proj.transform_points(ccrs.PlateCarree(), lons, lats, rain)
        transform[np.where(np.isnan(transform))] = 0
        lons = transform[:, :, 0]
        lats = transform[:, :, 1]
        ax.contourf(lons, lats, snow, transform=target_proj, levels=prec_levels, cmap=snow_cmap, norm=snow_norm)
        ax.contourf(lons, lats, rain, transform=target_proj, levels=prec_levels, cmap=rain_cmap, norm=rain_norm)
        ax.contourf(lons, lats, frzr, transform=target_proj, levels=prec_levels, cmap=frzr_cmap, norm=frzr_norm)
        ax.contourf(lons, lats, ice, transform=target_proj, levels=prec_levels, cmap=ice_cmap, norm=ice_norm)
        slpc = ax.contour(lons, lats, slp, transform=target_proj, colors='black', alpha=0.7, levels=slp_levels, linewidths=0.125)
        t2mc = ax.contour(lons, lats, t2m, transform=target_proj, colors='black', alpha=0.7, levels=t2m_levels, linewidths=0.125)
        ax.clabel(slpc, inline=True, fontsize=4)
        ax.clabel(t2mc, inline=True, fontsize=4)
    else:
        ax.contourf(lons, lats, snow, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=snow_cmap, norm=snow_norm)
        ax.contourf(lons, lats, rain, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=rain_cmap, norm=rain_norm)
        ax.contourf(lons, lats, frzr, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=frzr_cmap, norm=frzr_norm)
        ax.contourf(lons, lats, ice, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=ice_cmap, norm=ice_norm)
        slpc = ax.contour(lons, lats, slp, transform=ccrs.PlateCarree(), colors='black', alpha=0.7, levels=slp_levels, linewidths=0.125)
        t2mc = ax.contour(lons, lats, t2m, transform=ccrs.PlateCarree(), colors='black', alpha=0.7, levels=t2m_levels, linewidths=0.125)
        ax.clabel(slpc, inline=True, fontsize=4)
        ax.clabel(t2mc, inline=True, fontsize=4)
    plt.axis('off')
    if not os.path.exists('tmp/'):
        print('creating temp folder')
        os.mkdir('tmp')
    plt.savefig(f'tmp/tmp_weather_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()
    img = Image.open(f'tmp/tmp_weather_{file_tag}.png')
    img2 = Image.open(f'cache/cb_weather_{file_tag}.png')
    img.paste(img2, mask=img2)
    img.save(f'tmp/{proj_name}-wxtype-{file_tag}.png')

def annotate(path: str, plot_type: str, mode='light', **kwargs):
    """
    mode (optional, str) - 'light' or 'dark' 
    """
    themes = {
        'light': 'black',
        'dark': 'white'
    }
    mapset = Image.open(path)
    if plot_type == 'ir8':
        colorbar = Image.open(f'assets/plot_ir_nesdis-cbar-{themes[mode]}.png')
    if plot_type == 'wxtype':
        colorbar = Image.open(f'assets/plot_wxtype4all-cbar-{themes[mode]}.png')
    gmao_logo = Image.open(f'assets/GMAO_logo_{themes[mode]}_628x180.png')
    nasa_logo = Image.open('assets/nasa-logo_transparent-210x180.png')

    bgcolor = (0, 0, 0) if mode == 'dark' else (255, 255, 255)
    dst = Image.new('RGB', (5760, 3240), bgcolor)
    map_corner = 80
    if 'geos' in path or 'nsper' in path:
        mapset = mapset.resize((2760, 2760))
        dst.paste(mapset, (int(dst.width / 2) - int(mapset.width / 2), 180))
    else:
        mapset = mapset.resize((5600, 2760))
        dst.paste(mapset, (map_corner, 180))

    dst.paste(colorbar, (1320, 2940), mask=colorbar)
    dst.paste(nasa_logo, (0, 0), mask=nasa_logo)
    dst.paste(gmao_logo, (dst.width - gmao_logo.width, 0))

    canvas = ImageDraw.Draw(dst)
    font = ImageFont.truetype('fonts/Open-Sans-Bold/OpenSans-Bold.ttf', 72)
    fill = (0, 0, 0) if mode == 'light' else (255, 255, 255)
    canvas.text((215, 75), 'Global Modeling and Assimilation Office', spacing=0, font=font, fill=fill)

    if 'forecast' in kwargs and 'date' in kwargs:
        canvas.text((0, 2940), kwargs['date'], spacing=0, font=font, fill=fill)
        canvas.text((dst.width - 715, 2940), kwargs['forecast'], align='right', spacing=0, font=font, fill=fill)

    base_image = path.split('-')[-1]
    if 'discover' in os.getcwd() or 'gpfsm' in os.getcwd():
        dst.save(f'/discover/nobackup/projects/gmao/g6dev/pub/qcambrel/plotall_{plot_type}/plotall_{plot_type}-{base_image}')
    else:
        new_path = f'results/plotall_{plot_type}/plotall_{plot_type}-{base_image}'
        dst.save(new_path)