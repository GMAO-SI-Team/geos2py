import os
import json
import numpy as np
import PIL.Image as image
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

class Plotter(object):
    """
    Plotter instance to render GEOS forecasts
    """
    def __init__(self, plot_type: str, region: str, file_tag: str, target_proj: ccrs.Projection, proj_name: str, label_coords=None, interpolator=None):
        """
        Initialize Plotter object.

        Parameters:
        - plot_type (str): Type of plot.
        - region (str): Region for the plot.
        - file_tag (str): Tag for the file.
        - target_proj (cartopy.crs.Projection): Target projection.
        - proj_name (str): Name of the projection.
        - label_coords (list of tuples, optional): Coordinates for labels. Defaults to None.
        - interpolator (callable, optional): Interpolator function. Defaults to None.
        """
        self.plot_type = plot_type
        self.region = region
        self.file_tag = file_tag
        lons = np.linspace(-180, 180, 5760)
        lats = np.linspace(-90, 90, 2760)
        self.lons, self.lats = np.meshgrid(lons, lats)
        self.target_proj = target_proj
        self.proj_name = proj_name
        self.label_coords = label_coords
        self.interpolator = interpolator
        self.limit_extent = proj_name in ('sub', 'ortho', 'laea')
        if self.limit_extent:
            with open('regions.json', 'r') as infile:
                region_info = json.load(infile)
            self.extent = region_info[self.region]['extent']
        self._set_paths()
        self._set_plotters()

    def _set_paths(self):
        """
        Set paths for different elements of the plot.
        """
        plot_abbr = {
            'plotall_ir8': 'ir',
            'plotall_wxtype': 'weather',
            'plotall_radar': 'radar',
            'plotall_aerosols': 'aero',
            'plotall_precrain': 'rain',
            'plotall_precsnow': 'rain',
            'plotall_slp': 'slp',
            'plotall_t2m': 't2m',
            'plotall_tpw': 'tpw',
            'plotall_cape': 'cape',
            'plotall_vort500mb': 'vort',
            'plotall_winds10m': 'wind'
        }
        self.contour_path = f'tmp/tmp_{plot_abbr[self.plot_type]}_{self.file_tag}.png'
        self.feature_path = f'cache/cb_{plot_abbr[self.plot_type]}_{self.file_tag}.png'
        self.completed_path = f'tmp/{self.proj_name}-{self.plot_type.split("_")[1]}-{self.file_tag}.png'
        
        grayscale = {
            'plotall_wxtype': 'auto',
            'plotall_radar': 'auto',
            'plotall_aerosols': 'grey',
            'plotall_precrain': 'white',
            'plotall_precsnow': 'white',
            'plotall_cape': 'white',
            'plotall_vort500mb': 'auto'
        }
        self.natural_earth_path = f'cache/natural_earth_{grayscale[self.plot_type]}_{self.file_tag}.png' if self.plot_type in grayscale else None

    def _set_plotters(self):
        """
        Set plotter functions based on plot type.
        """
        self.plotters = {
            'plotall_ir8': self._plot_infrared,
            'plotall_wxtype': self._plot_weather,
            'plotall_radar': self._plot_radar,
            'plotall_aerosols': self._plot_aerosols,
            'plotall_precrain': self._plot_rain,
            'plotall_precsnow': self._plot_snow,
            'plotall_slp': self._plot_pressure,
            'plotall_t2m': self._plot_temp,
            'plotall_tpw': self._plot_water,
            'plotall_cape': self._plot_energy,
            'plotall_vort500mb': self._plot_vorticity,
            'plotall_winds10m': self._plot_wind
        }

    def render(self, data, cmap, norm, save=True):
        """
        Render the plot.

        Parameters:
        - data (np.ndarray or list of nd.arrays): Data for the plot.
        - cmap (mpl.colors.ListedColormap): Colormap for the plot.
        - norm (mpl.colors.Normalize): Normalization for the plot.
        - save (bool, optional): Whether to save the plot. Defaults to True.
        """
        self.data = data
        self.cmap = cmap
        self.norm = norm
        self.fig = plt.figure(dpi=1500)
        self.ax = plt.axes(projection=self.target_proj)
        if self.limit_extent:
            self.ax.set_extent(self.extent, ccrs.PlateCarree())
        self.plotters[self.plot_type]()
        plt.axis('off')
        if not os.path.exists('tmp/'):
            print('creating temp folder')
            os.mkdir('tmp')
        plt.savefig(self.contour_path, bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()
        if self.natural_earth_path:
            img = image.open(self.natural_earth_path)
            img2 = image.open(self.contour_path)
            img3 = image.open(self.feature_path)
            img.paste(img2, mask=img2)
            img.paste(img3, mask=img3)
            img.save(self.completed_path)
        else:
            img = image.open(self.contour_path)
            img2 = image.open(self.feature_path)
            img.paste(img2, mask=img2)
            img.save(self.completed_path)

    def _plot_infrared(self):
        """
        Plot clean longwave infrared data.
        """
        # self.ax.pcolormesh(self.lons, self.lats, self.data, transform=ccrs.PlateCarree(), cmap=self.cmap, norm=self.norm)
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())

    def _plot_weather(self):
        """
        Plot weather data. Includes snow, ice, freezing rain, rain, 2-meter temperature, and sea level pressure.
        """
        snow, ice, frzr, rain, t2m, slp = self.data
        snow_cmap, ice_cmap, frzr_cmap, rain_cmap = self.cmap
        snow_norm, ice_norm, frzr_norm, rain_norm = self.norm
        prec_levels = [0.01,0.02,0.03,0.05,0.10,0.25,0.5,1.0]
        levs1 = 884 + np.arange(12) * 8
        levs2 = 980 + np.arange(5) * 4
        levs3 = 1000 + np.arange(31) * 2
        slp_levels = np.append(np.append(levs1, levs2), levs3)
        t2m_levels = [273.15]

        snow_mask = snow < 0.01
        ice_mask = ice < 0.01
        frzr_mask = frzr < 0.01
        rain_mask = rain < 0.01
        snow_lons = np.ma.masked_where(snow_mask, self.lons)
        snow_lats = np.ma.masked_where(snow_lats, self.lats)
        ice_lons = np.ma.masked_where(ice_mask, self.lons)
        ice_lats = np.ma.masked_where(ice_lats, self.lats)
        frzr_lons = np.ma.masked_where(frzr_mask, self.lons)
        frzr_lats = np.ma.masked_where(frzr_lats, self.lats)
        rain_lons = np.ma.masked_where(rain_mask, self.lons)
        rain_lats = np.ma.masked_where(rain_lats, self.lats)
        self.ax.contourf(snow_lons, snow_lats, snow, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=snow_cmap, norm=snow_norm)
        self.ax.contourf(rain_lons, rain_lats, rain, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=rain_cmap, norm=rain_norm)
        self.ax.contourf(frzr_lons, frzr_lats, frzr, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=frzr_cmap, norm=frzr_norm)
        self.ax.contourf(ice_lons, ice_lats, ice, transform=ccrs.PlateCarree(), levels=prec_levels, cmap=ice_cmap, norm=ice_norm)
        slpc = self.ax.contour(self.lons, self.lats, slp, transform=ccrs.PlateCarree(), colors='black', alpha=0.7, levels=slp_levels, linewidths=0.125)
        t2mc = self.ax.contour(self.lons, self.lats, t2m, transform=ccrs.PlateCarree(), colors='black', alpha=0.7, levels=t2m_levels, linewidths=0.125)
        self.ax.clabel(slpc, inline=True, fontsize=4)
        self.ax.clabel(t2mc, inline=True, fontsize=4)

    def _plot_aerosols(self):
        """
        Plot aerosols data. Includes seasalt, dust, organic carbon, black carbon, sulfate and nitrate.
        """
        ss, du, oc, bc, su, ni = self.data
        self.ax.imshow(ss, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')
        self.ax.imshow(du, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')
        self.ax.imshow(ni, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')
        self.ax.imshow(su, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')
        self.ax.imshow(oc, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')
        self.ax.imshow(bc, interpolation='nearest', transform=ccrs.PlateCarree(), extent=(-180, 180, -90, 90), origin='lower')

    def _plot_radar(self):
        """
        Plot simulated radar sensitivity data.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())

    def _plot_rain(self):
        """
        Plot accumulated rain data.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        for label_coord in self.label_coords:
            y_coord, x_coord = label_coord
            interpolated_value = self.interpolator((y_coord, x_coord))
            if self.limit_extent:
                x_0, x_1, y_0, y_1 = self.extent
                if x_0 <= x_coord <= x_1 and y_0 <= y_coord <= y_1:
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='white', clip_on=True, fontsize=1, transform=ccrs.PlateCarree())
            else:
                point = self.target_proj.transform_point(x_coord, y_coord, ccrs.PlateCarree())
                if True not in np.isnan(point):
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='white', fontsize=1, transform=ccrs.PlateCarree())

    def _plot_snow(self):
        """
        Plot accumulated snow data.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        for label_coord in self.label_coords:
            y_coord, x_coord = label_coord
            interpolated_value = self.interpolator((y_coord, x_coord))
            if self.limit_extent:
                x_0, x_1, y_0, y_1 = self.extent
                if x_0 <= x_coord <= x_1 and y_0 <= y_coord <= y_1:
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='white', clip_on=True, fontsize=1, transform=ccrs.PlateCarree())
            else:
                point = self.target_proj.transform_point(x_coord, y_coord, ccrs.PlateCarree())
                if True not in np.isnan(point):
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='white', fontsize=1, transform=ccrs.PlateCarree())

    def _plot_water(self):
        """
        Plot total precipitable water data.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())

    def _plot_pressure(self):
        """
        Plot sea level pressure data. Labels pressure minima.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        slp_min_locs, x_window = self.label_coords
        for slp_min_loc in slp_min_locs:
            x_coord, y_coord, pressure = slp_min_loc
            label_size = max(3, 2.5 / (360 * 111.11111 / (360 / x_window)))
            if self.limit_extent:
                x_0, x_1, y_0, y_1 = self.extent
                if x_0 <= x_coord <= x_1 and y_0 <= y_coord <= y_1:
                    self.ax.text(x_coord, y_coord, f'{np.round(pressure).astype(int)}', color='black', clip_on=True, fontsize=label_size, transform=ccrs.PlateCarree())
            else:
                point = self.target_proj.transform_point(x_coord, y_coord, ccrs.PlateCarree())
                if True not in np.isnan(point):
                    label_size *= 0.5
                    self.ax.text(x_coord, y_coord, f'{np.round(pressure).astype(int)}', color='black', clip_on=True, fontsize=label_size, transform=ccrs.PlateCarree())

    def _plot_temp(self):
        """
        Plot 2-meter temperature data. Labels temp by city.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        for label_coord in self.label_coords:
            y_coord, x_coord = label_coord
            interpolated_value = self.interpolator((y_coord, x_coord))
            label_by_temp = {
                40 <= interpolated_value < 70: ('black', 0.9),
                -20 <= interpolated_value < 40: ('white', 0.95),
                -20 > interpolated_value: ('magenta', 1),
                70 <= interpolated_value < 80: (np.array([255, 255, 0]) / 255, 0.9),
                80 <= interpolated_value < 90: (np.array([255, 255, 175]) / 255, 0.95),
                90 <= interpolated_value < 100: (np.array([255, 255, 255]) / 255, 1),
                100 <= interpolated_value < 105: (np.array([255, 75, 210]) / 255, 1.05),
                105 <= interpolated_value < 110: (np.array([255, 0, 255]) / 255, 1.1),
                110 <= interpolated_value: (np.array([225, 0, 225]) / 255, 1.15)
            }
            label_color = label_by_temp[True][0]
            label_size = 3 * label_by_temp[True][1]
            if self.limit_extent:
                x_0, x_1, y_0, y_1 = self.extent
                if x_0 <= x_coord <= x_1 and y_0 <= y_coord <= y_1:
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color=label_color, clip_on=True, fontsize=label_size, transform=ccrs.PlateCarree())
            else:
                point = self.target_proj.transform_point(x_coord, y_coord, ccrs.PlateCarree())
                if True not in np.isnan(point):
                    self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color=label_color, fontsize=label_size, transform=ccrs.PlateCarree())

    def _plot_energy(self):
        """
        Plot surface convective available potential energy (CAPE) data.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())                
    
    def _plot_wind(self):
        """
        Plot 10-meter winds data. Labels wind vectors by city.
        """
        self.ax.imshow(self.data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        for label_coord in self.label_coords:
            y_coord, x_coord = label_coord
            interpolated_value = self.interpolator((y_coord, x_coord))
            if self.limit_extent:
                x_0, x_1, y_0, y_1 = self.extent
                if x_0 <= x_coord <= x_1 and y_0 <= y_coord <= y_1:
                    if interpolated_value >= 10:
                        self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='black', clip_on=True, fontsize=3, transform=ccrs.PlateCarree())
            else:
                point = self.target_proj.transform_point(x_coord, y_coord, ccrs.PlateCarree())
                if True not in np.isnan(point):
                    if interpolated_value >= 10:
                        self.ax.text(x_coord, y_coord, f'{np.round(interpolated_value).astype(int)}', color='black', fontsize=3, transform=ccrs.PlateCarree())

    def _plot_vorticity(self):
        """
        Plot 500 meter vorticity data. Includes vorticity and heights.
        """
        data, heights = self.data
        levels = np.arange(60) * 30 + 4500
        self.ax.imshow(data, cmap=self.cmap, norm=self.norm, origin='lower', interpolation='nearest', extent=(-180, 180, -90, 90), transform=ccrs.PlateCarree())
        heightsc = self.ax.contour(self.lons, self.lats, heights, transform=ccrs.PlateCarree(), levels=levels, colors='black', linewidths=0.125)
        self.ax.clabel(heightsc, inline=True, fontsize=4)