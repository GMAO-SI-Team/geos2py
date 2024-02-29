import csv
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interp

class Colormap(object):
    """
    Colormap instance based on raw NetCDF variable or processed variable
    """
    def __init__(self, plot_type: str, variable: str, data_min=None, data_max=None):
        """
        Initializes Colormap object with plot type, variable, and optional data range.
        
        Parameters:
        - plot_type (str): Type of plot.
        - variable (str): Type of variable.
        - data_min (float, optional): Minimum data value. Defaults to None.
        - data_max (float, optional): Maximum data value. Defaults to None.
        """
        self.plot_type = plot_type
        self.variable = variable
        self.data_min = data_min
        self.data_max = data_max
        self.set_cmap()
    
    def __repr__(self):
        """
        Returns a string representation of the Colormap object.
        """
        state = self.variable in ('snow', 'rain', 'frzr', 'ice', 'vort', 'spd')
        state_dict = {
            True: 'processed',
            False: 'raw'
        }
        return f'Colormap [{self.plot_type}] -> {state_dict[state]} variable: {self.variable}'

    # cmap_dict: double nested dictionary
    # outer keys -> plot types
    # inner keys -> raw or processed variable
    # inner values -> associated cmap definition function
    # function returns dict(rgb/cmap, levels)
    def set_cmap(self):
        """
        Sets the colormap based on plot type and variable type.
        """
        self.cmap_dict = {
            # Dictionary mapping plot types to variable types and associated functions
            'plotall_ir8': {
                'TBISCCP': self._set_TBISCCP
            },
            'plotall_wxtype': {
                'snow': self._set_snow,
                'rain': self._set_rain,
                'frzr': self._set_frzr,
                'ice': self._set_ice
            },
            'plotall_radar': {
                'DBZ_MAX': self._set_DBZ_MAX
            },
            'plotall_aerosols': {
                'SSEXTTAU': self._set_SSEXTTAU,
                'DUEXTTAU': self._set_DUEXTTAU,
                'OCEXTTAU': self._set_OCEXTTAU,
                'BCEXTTAU': self._set_BCEXTTAU,
                'SUEXTTAU': self._set_SUEXTTAU,
                'NIEXTTAU': self._set_NIEXTTAU
            },
            'plotall_precrain': {
                'PRECTOT': self._set_PRECTOT
            },
            'plotall_precsnow': {
                'PRECSNO': self._set_PRECSNO
            },
            'plotall_slp': {
                'SLP': self._set_SLP
            },
            'plotall_t2m': {
                'T2M': self._set_T2M
            },
            'plotall_tpw': {
                'TQV': self._set_TQV
            },
            'plotall_cape': {
                'CAPE': self._set_CAPE
            },
            'plotall_vort500mb': {
                'vort': self._set_vort
            },
            'plotall_winds10m': {
                'spd': self._set_spd
            }
        }
        try:
            cmap_data = self.cmap_dict[self.plot_type][self.variable]()
        except KeyError:
            self._handle_error()
        rgb = cmap_data['rgb']
        levels = cmap_data['levels']
        cmap = plt.colormaps[rgb] if type(rgb) == str else mpl.colors.ListedColormap(rgb)
        if self.plot_type == 'plotall_aerosols':
            norm = mpl.colors.Normalize(vmin=self.data_min, vmax=self.data_max)
        else:
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
        self.cmap = cmap
        self.norm = norm
        self.levels = levels

    def _set_TBISCCP(self):
        """
        Sets colormap for TBISCCP plot type.
        """
        color_table = 'NESDIS_IR_10p3micron.txt'
        rgb = []
        with open(f'colortables/{color_table}', 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                rgb_row = [item.strip() for item in row[0].split('     ')[1:]]
                rgb.append(rgb_row)
        rgb = np.array(rgb).astype(int) / 255
        ticks = np.array([-110, -59, -20, 6, 31, 57])
        levels = np.array(ticks, dtype=np.float64)
        function = interp.interp1d(np.arange(len(levels)), levels)
        target = 5 * np.arange(256) / 255.0
        levels = function(np.linspace(0.0, len(levels) - 1, len(target)))
        return {'rgb': rgb, 'levels': levels}
    
    @staticmethod
    def _set_levels_wxtype(func):
        """
        Decorator to set levels for weather type colormaps.
        """
        def setter(*args):
            prec_dict = func(*args)
            levels = np.array([0.01, 0.02, 0.03, 0.05, 0.10, 0.25, 0.5, 1.0])
            prec_dict['levels'] = levels
            return prec_dict
        return setter

    @_set_levels_wxtype
    def _set_snow(self):
        """
        Sets colormap for snow data.
        """
        rgb_snow = np.array([
            [172,196,225],
            [118,172,204],
            [ 67,139,186],
            [ 41,101,162],
            [ 24, 68,145],
            [ 12, 44,113],
            [ 12, 27, 64],
            [ 12, 27, 64]
        ]) / 255
        return {'rgb': rgb_snow}

    @_set_levels_wxtype
    def _set_rain(self):
        """
        Sets colormap for rain data.
        """
        rgb_rain = np.array([
            [114,198,114],
            [ 38,136, 67],
            [ 26, 81, 35],
            [ 12, 32, 18],
            [236,246, 80],
            [253,179, 57],
            [252, 98, 33],
            [252, 98, 33]
        ]) / 255
        return {'rgb': rgb_rain}

    @_set_levels_wxtype
    def _set_frzr(self):
        """
        Sets colormap for freezing rain data.
        """
        rgb_frzr= np.array([
            [253,169,195],
            [252,119,154],
            [250, 59, 95],
            [200, 16, 57],
            [161,  8, 30],
            [105,  3, 23],
            [ 58,  1, 16],
            [ 58,  1, 16]
        ]) / 255
        return {'rgb': rgb_frzr}

    @_set_levels_wxtype
    def _set_ice(self):
        """
        Sets colormap for ice data.
        """
        rgb_ice = np.array([
            [242,179,251],
            [228,144,249],
            [214,107,243],
            [211, 68,248],
            [205, 40,247],
            [161, 12,209],
            [122,  7,156],
            [122,  7,156]
        ]) / 255
        return {'rgb': rgb_ice}

    def _set_DBZ_MAX(self):
        """
        Sets colormap for radar data.
        """
        rgb = np.array([
            [255, 255, 255],
            [0, 224, 227],
            [0, 141, 243],
            [0, 12, 243],
            [0, 239, 8],
            [0, 183, 0],
            [0, 123, 0],
            [255, 246, 0],
            [228, 173, 0],
            [255, 129, 0],
            [255, 0, 0],
            [209, 0, 0],
            [180, 0, 0],
            [249, 7, 253],
            [133, 67, 186],
            [245, 245, 245]
        ]) / 255
        return {'rgb': rgb, 'levels': np.arange(16) * 5.0}

    @staticmethod
    def _set_levels_aerosols(func):
        """
        Decorator to set levels for aerosol colormaps.
        """
        def setter(*args):
            prec_dict = func(*args)
            levels = np.array([0, 0.5])
            prec_dict['levels'] = levels
            return prec_dict
        return setter

    @_set_levels_aerosols
    def _set_SSEXTTAU(self):
        """
        Sets colormap for seasalt data.
        """
        return {'rgb': 'Blues'}
    
    @_set_levels_aerosols
    def _set_DUEXTTAU(self):
        """
        Sets colormap for dust data.
        """
        return {'rgb': 'YlOrBr'}

    @_set_levels_aerosols
    def _set_OCEXTTAU(self):
        """
        Sets colormap for organic carbon data.
        """
        return {'rgb': 'Greens'}
    
    @_set_levels_aerosols
    def _set_BCEXTTAU(self):
        """
        Sets colormap for black carbon data.
        """
        color = np.array([75, 5, 200]) / 255
        r, g, b = color
        rgb = np.zeros((256, 3))
        rgb[:,0] = ((1 - r) * (np.arange(256) / 255.0) + r)[::-1]
        rgb[:,1] = ((1 - g) * (np.arange(256) / 255.0) + g)[::-1]
        rgb[:,2] = ((1 - b) * (np.arange(256) / 255.0) + b)[::-1]
        return {'rgb': rgb}
    
    @_set_levels_aerosols
    def _set_SUEXTTAU(self):
        """
        Sets colormap for sulfate data.
        """
        return {'rgb': 'RdPu'}
    
    @_set_levels_aerosols
    def _set_NIEXTTAU(self):
        """
        Sets colormap for nitrate data.
        """
        color = np.array([255, 255, 50]) / 255
        r, g, b = color
        rgb = np.zeros((256, 3))
        rgb[:,0] = ((1 - r) * (np.arange(256) / 255.0) + r)[::-1]
        rgb[:,1] = ((1 - g) * (np.arange(256) / 255.0) + g)[::-1]
        rgb[:,2] = ((1 - b) * (np.arange(256) / 255.0) + b)[::-1]
        return {'rgb': rgb}
    
    def _set_PRECTOT(self):
        """
        Sets colormap for accumulated rain data.
        """
        rgb = np.array([
            [255, 255, 255],
            [175, 210, 235],
            [130, 160, 230],
            [90, 105, 220],
            [65, 175, 45],
            [95, 215, 75],
            [140, 240, 130],
            [175, 255, 165],
            [250, 250, 35],
            [250, 215, 30],
            [250, 190, 30],
            [245, 150, 30],
            [240, 110, 30],
            [230, 30, 30],
            [220, 0, 25],
            [200, 0, 25],
            [180, 0, 15],
            [140, 0, 10],
            [170, 50, 60],
            [205, 130, 130],
            [230, 185, 185],
            [245, 225, 225],
            [215, 200, 230],
            [185, 165, 210],
            [155, 125, 185],
            [140, 100, 170],
            [120, 70, 160]
        ]) / 255
        levels = np.array([
            0.01, 0.1, 0.25, 0.50, 0.75, 
            1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
            4.0, 4.5, 5.0, 6.0, 7.0, 8.0,
            9.0, 10, 12, 14, 16, 18, 20, 
            22, 24
        ])
        return {'rgb': rgb, 'levels': levels}
    
    def _set_PRECSNO(self):
        """
        Sets colormap for accumulated snow data.
        """
        rgb = np.array([
            [221,221,221], 
            [193,193,193], 
            [165,165,165], 
            [139,139,139], 
            [160,227,239], 
            [104,183,207], 
            [ 53,141,176], 
            [ 11,101,147], 
            [  0, 61,162], 
            [ 35, 89,175], 
            [ 68,116,189], 
            [105,149,202], 
            [144,182,216], 
            [185,215,231], 
            [193,164,212], 
            [184,138,201], 
            [176,115,191], 
            [166, 92,181], 
            [157, 70,171], 
            [149, 49,162], 
            [123, 16, 62], 
            [139, 28, 79], 
            [156, 41, 97], 
            [172, 56,116], 
            [190, 72,136], 
            [207, 88,156], 
            [231,157,174], 
            [227,145,151], 
            [225,132,128], 
            [222,119,107], 
            [220,107, 86], 
            [216, 96, 67], 
            [214,118, 76], 
            [219,137, 96], 
            [226,159,119], 
            [231,182,145], 
            [238,203,168], 
            [243,226,194], 
            [250,248,219]
        ]) / 255
        levels = np.array([
            0.1, 0.25, 0.5, 0.75,
            1, 1.5, 2, 2.5, 3,
            3.5, 4, 4.5, 5, 5.5,
            6, 7, 8, 9, 10, 11,
            12, 14, 16, 18, 20,
            22, 24, 26, 28, 30,
            32, 34, 36, 40, 44,
            48, 52, 56, 60
        ])
        return {'rgb': rgb, 'levels': levels}
    
    def _set_SLP(self):
        """
        Sets colormap for sea level pressure data.
        """
        rgb = np.array([
            [250,250,250], 
            [225,225,225], 
            [200,200,200], 
            [175,175,175], 
            [150,150,150], 
            [125,125,125], 
            [110,110,110], 
            [ 90, 90, 90], 
            [ 70, 70, 70], 
            [ 50, 50, 50], 
            [ 30, 30, 30], 
            [ 10,  0, 10], 
            [ 75,  0,  0], 
            [148,  0,  0], 
            [178,  0,  0], 
            [217,  0,  0], 
            [255, 20,  0], 
            [255, 72,  0], 
            [255,143,  0], 
            [255,181, 31], 
            [255,231, 93], 
            [255,252,148], 
            [218,255,255], 
            [131,198,251], 
            [ 62,143,245], 
            [ 18, 66,240], 
            [ 27, 10,204], 
            [ 30, 11, 60], 
            [ 58, 30,102], 
            [142, 47,150], 
            [199, 95,201], 
            [224,152,238], 
            [218,205,218], 
            [175,154,175], 
            [225,200,200]
        ]) / 255
        levels = np.array([
            880, 885, 890, 895, 900, 905, 910, 915,
            920, 925, 930, 935, 940, 945, 950, 955,
            960, 965, 970, 975, 980, 985, 990, 995,
            1000, 1004, 1006, 1008, 1010, 1012, 1014,
            1016, 1020, 1028
        ])
        return {'rgb': rgb, 'levels': levels}
    
    def _set_T2M(self):
        """
        Sets colormap for 2-meter temperature data.
        """
        rgb = np.array([
            [ 30, 30, 75], 
            [ 60, 60,125], 
            [105,120,190], 
            [125,150,206], 
            [153,190,230], 
            [200,225,240], 
            [235,226,235], 
            [203,139,190], 
            [180, 56,152], 
            [137,  0,123], 
            [ 90,  0,129], 
            [ 60,  0,131], 
            [ 21,  0,135], 
            [ 14,  0,216], 
            [ 26, 51,241], 
            [ 43,129,241], 
            [ 68,213,240], 
            [ 90,240,200], 
            [127,239,136], 
            [186,239, 54], 
            [235,235,  0], 
            [233,175,  0], 
            [232,123,  0], 
            [231, 46,  0], 
            [218,  0,  0], 
            [160,  0,  0], 
            [110,  0,  0], 
            [ 80,  0,  0], 
            [ 65,  0,  0], 
            [ 50,  0,  0], 
            [ 75, 65, 65], 
            [110,110,110], 
            [180,180,180]
        ]) / 255
        levels = np.array([
            -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 
            10, 15, 20 ,25, 30, 35, 40, 45, 50, 55, 60, 65, 
            70, 75, 80, 85, 90, 95, 100, 105, 110, 115
        ])
        return {'rgb': rgb, 'levels': levels}

    def _set_TQV(self):
        """
        Sets colormap for total precipitable water data.
        """
        with open('colortables/MPL_terrain.rgb', 'r') as infile:
            rgb = infile.readlines()[2:]
        rgb = [line.strip() for line in rgb]
        rgb = [line.split() for line in rgb]
        rgb = rgb[::-1]
        rgb = np.array(rgb).astype(float)
        levels = np.arange(128) / 2 + 5
        return {'rgb': rgb, 'levels': levels}

    def _set_CAPE(self):
        """
        Sets colormap for surface CAPE data.
        """
        rgb = np.array([
            [200,200,200], 
            [160,160,160], 
            [125,125,125], 
            [ 91, 91, 91], 
            [116,170,255], 
            [ 83,125,226], 
            [ 52, 84,197], 
            [ 23, 43,158], 
            [117,255,117], 
            [ 79,197, 79], 
            [ 45,142, 45], 
            [ 16, 91, 16], 
            [255,255, 91], 
            [221,170, 60], 
            [188, 91, 31], 
            [156, 24,  0], 
            [255,142,255], 
            [212,106,212], 
            [170, 70,170], 
            [129, 39,129], 
            [ 90, 10, 90]
        ]) / 255
        levels = np.array([
            100, 325, 550, 775,
            1000, 1250, 1500, 1750,
            2000, 2250, 2500, 2750,
            3000, 3500, 4000, 4500,
            5000, 6000, 7000, 8000,
            9000, 10000
        ])
        return {'rgb': rgb, 'levels': levels}

    def _set_vort(self):
        """
        Sets colormap for 500mb vorticity data.
        """
        with open('colortables/005_STD_GAMMA-II.dat', 'r') as infile:
            rgb = infile.readlines()
        rgb = [line.strip() for line in rgb]
        rgb = [line.split(',') for line in rgb]
        rgb = np.array(rgb).astype(float)[::-1]
        return {'rgb': rgb, 'levels': np.arange(61)}
    
    def _set_spd(self):
        """
        Sets colormap for 10-meter winds data.
        """
        rgb = np.array([
            [240,240,240], 
            [210,210,210], 
            [187,187,187], 
            [150,150,150], 
            [ 20,100,210], 
            [ 40,130,240], 
            [ 82,165,243], 
            [ 53,213, 51], 
            [ 80,242, 79], 
            [147,246,137], 
            [203,253,194], 
            [255,246,169], 
            [255,233,124], 
            [253,193, 63], 
            [255,161,  0], 
            [255, 95,  4], 
            [255, 50,  0], 
            [225, 20,  0], 
            [197,  0,  0], 
            [163,  0,  0], 
            [118, 82, 70], 
            [138,102, 90], 
            [177,143,133], 
            [225,191,182], 
            [238,220,210], 
            [255,200,200], 
            [245,160,160], 
            [225,136,130], 
            [232,108,100], 
            [229, 96, 87]
        ]) / 255
        levels = np.array([
            0, 3, 6, 9, 12, 15, 18, 20, 22, 24, 27,
            30, 35, 40, 45, 50, 55, 60, 65, 70, 75,
            80, 85, 90, 100, 105, 110, 120, 135
        ])
        return {'rgb': rgb, 'levels': levels}
    
    def _handle_error(self):
        error_type = 0
        if self.plot_type not in self.cmap_dict:
            error_type = 1
        if True not in [self.variable in val for val in self.cmap_dict.values()]:
            error_type = 2
        if not error_type and self.variable not in self.cmap_dict[self.plot_type]:
            error_type = 3
        message = {
            1: 'Invalid or unsupported plot type provided',
            2: 'Invalid or unsupported variable provided',
            3: 'Incorrect variable for provided plot type'
        }
        print(message[error_type])
        sys.exit()
    

if __name__ == '__main__':
    import os
    import io
    import unittest
    from contextlib import redirect_stdout

    def test_tbisccp_ctable_exists():
        pass

    def test_vort_ctable_exists():
        pass

    def test_tqv_ctable_exists():
        pass

    def test_wxtypes_levels():
        pass

    def test_aerosols_levels():
        pass

    def test_method_privacy():
        pass

    def test_bad_plot_type():
        f = io.StringIO()
        with redirect_stdout(f):
            test_cmap = Colormap('plotall_bad', 'TBISCCP')
        out = f.getvalue()
        return out

    def test_bad_variable():
        f = io.StringIO()
        with redirect_stdout(f):
            test_cmap = Colormap('plotall_ir8', 'TBISCUP')
        out = f.getvalue()
        return out

    def test_plot_type_variable_mismatch():
        f = io.StringIO()
        with redirect_stdout(f):
            test_cmap = Colormap('plotall_ir8', 'DBZ_MAX')
        out = f.getvalue()
        return out
    
    assert test_bad_plot_type() == 'Invalid or unsupported plot type provided'
    assert test_bad_variable() == 'Invalid or unsupported variable provided'
    assert test_plot_type_variable_mismatch() == 'Incorrect variable for provided plot type'