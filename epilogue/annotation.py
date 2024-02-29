import os
from PIL import Image, ImageDraw, ImageFont

def annotate(path: str, plot_type: str, mode='light', **kwargs):
    """
    Annotates map images with colorbars, logos, and additional information.

    Parameters:
    - path (str): Path to the base map image.
    - plot_type (str): Type of plot for which annotation is being done.
    - mode (str, optional): Color mode of the annotation (either 'light' or 'dark'). Defaults to 'light'.
    - **kwargs: Additional keyword arguments including forecast information and date.

    Returns:
    - None
    """
    themes = {
        'light': 'black',
        'dark': 'white'
    }
    colorbars = {
        'plotall_ir8': f'assets/plot_ir_nesdis-cbar-{themes[mode]}.png',
        'plotall_wxtype': f'assets/plot_wxtype4all-cbar-{themes[mode]}.png',
        'plotall_aerosols': f'assets/aerosols3_cbar-{themes[mode]}.png',
        'plotall_radar': f'assets/plot_radar-cbar-{themes[mode]}.png',
        'plotall_precrain': f'assets/precrain-cbar-{themes[mode]}.png',
        'plotall_precsnow': f'assets/precsnow10t1pivotal-cbar-{themes[mode]}.png',
        'plotall_slp': f'assets/plot_slp-cbar-{themes[mode]}.png',
        'plotall_t2m': f'assets/plot_t2m-cbar-{themes[mode]}.png',
        'plotall_tpw': f'assets/plot_tpw2-cbar_{themes[mode]}.png',
        'plotall_cape': f'assets/plot_cape-cbar-{themes[mode]}.png',
        'plotall_vort500mb': f'assets/plot_vort500mb-cbar-{themes[mode]}.png',
        'plotall_winds10m': f'assets/plot_winds10m-cbar-{themes[mode]}.png'
    }
    mapset = Image.open(path)
    colorbar = Image.open(colorbars[plot_type])
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

    if plot_type == 'plotall_aerosols':
        dst.paste(colorbar, (820, 2960), mask=colorbar)
    else:
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
        dst.save(f'/discover/nobackup/projects/gmao/g6dev/pub/qcambrel/{plot_type}/{plot_type}-{base_image}')
    else:
        new_path = f'results/{plot_type}/{plot_type}-{base_image}'
        dst.save(new_path)