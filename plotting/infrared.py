import os
import csv
import json
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.interpolate as interp

from PIL import Image, ImageDraw, ImageFont

def get_contour_levels(ticks: list):
    ticks = np.array(ticks) # default: [-110, -59, -20, 6, 31, 57]
    levels = np.array(ticks, dtype=np.float64)
    function = interp.interp1d(np.arange(len(levels)), levels)
    target = 5 * np.arange(256) / 255.0
    return function(np.linspace(0.0, len(levels) - 1, len(target)))

def get_contour_colors(plot_type):
    color_table = 'NESDIS_IR_10p3micron.txt' if plot_type == 'plotall_ir8' else None # infrared color table
    rgb_values = []
    with open(f'cmaps/{color_table}', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            rgb_value = [item.strip() for item in row[0].split('     ')[1:]]
            rgb_values.append(np.array(rgb_value))
    rgb_values = np.array(rgb_values).astype(int)
    return rgb_values / 255

def plot(region, file_tag, target_proj, proj_name, data, cmap, norm):
    fig = plt.figure(dpi=1500)
    ax = plt.axes(projection=target_proj)
    lons = np.linspace(-180, 180, 5760)
    lats = np.linspace(-90, 90, 2760)
    lons, lats = np.meshgrid(lons, lats)
    if proj_name in ['sub', 'ortho', 'laea']:
        with open('regions.json', 'r') as infile:
            regions_dict = json.load(infile)
        ax.set_extent(regions_dict[region]['extent'], ccrs.PlateCarree())
    ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    plt.axis('off')
    if not os.path.exists('tmp/'):
        print('creating temp folder')
        os.mkdir('tmp')
    plt.savefig(f'tmp/tmp_ir_{file_tag}.png', bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close()
    img = Image.open(f'tmp/tmp_ir_{file_tag}.png')
    img2 = Image.open(f'cache/cb_ir_{file_tag}.png')
    img.paste(img2, mask=img2)
    img.save(f'tmp/{proj_name}-ir8-{file_tag}.png')

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
