# Experimental GEOS Forecasts with Python
GMAOpy has been deprecated as it has fulfilled its role as a proof of concept. Future development efforts are now concentrated on GEOSpy. The entire codebase has been consolidated into a single script, plot.py, which now handles all data processing, visualization, and utility functions. This unification simplifies workflow execution and enhances performance while maintaining the modular architecture that allows for flexibility and efficiency.

## Table of Contents
1.	[Environment Setup](#environment-setup)
2.	[Performance](#performance)
3.	[Master Cache Script](#master-cache-script)
4.	[Base Images and Shapefiles](#base-images-and-shapefiles)
5.	[Modules Overview](#modules-overview)
6.	[How to Run the Scripts](#how-to-run-the-scripts)
7.	[Plot Types and Region Codes](#plot-types-and-region-codes)

## Environment Setup


Ensure the following dependencies are installed:

*	python >= 3.11.4
*	matplotlib >= 3.9.1
*	cartopy >= 0.23.0
*	pillow >= 10.4.0
*	netcdf4 >= 1.7.1
*	numpy >= 2.0.1
*	scipy >= 1.14.0
*	numba >= 0.60.0
*	pyogrio >= 0.9.0
*	geopandas >= 1.0.1
*	pandas >= 2.2.2
*	contourpy >= 1.2.1
*	sunpy >= 6.0.0

## Performance
The entire project has been refactored into a single script, plot.py, leading to significant performance improvements:

*	Data Processing: Conservative regridding completes in approximately one second or less. Serialization and deserialization processes are highly optimized, while loading base images takes only seconds.
*	Plotting: Some plots now render in as little as 10 seconds, with most graphical bugs resolved. However, there remains a computational issue with SLP minima location in the plotall_slp function.

## Master Cache Script

The master_cache.py script is responsible for managing the caching of
reprojected features and base images. It requires the following command-line
arguments:

*	plot_type: Name of the plot (e.g., plotall_ir8, plotall_radar) or all for all plot types.
*	region: Numeric region code or all for all regions (_mapset and _proj).
Optional argument:
*	--element: Specifies whether to cache "features" (e.g., coastlines, borders) or "base_image" (e.g., background images).

## Base Images and Shapefiles

### Base Images

*	Color: /discover/nobackup/projects/gmao/g6dev/pub/BMNG/New/eo_base_2020_clean_geo.3x21600x10800.jpg
*	Grayscale: /discover/nobackup/projects/gmao/g6dev/pub/BMNG/natural_earth_grey_16200x8100.jpeg

### Shapefiles

*	Used in Code: /discover/nobackup/qcambrel/gmaopy/SHAPE_FILES
*	Original Location: /home/wputman/IDL_ESSENTIALS/SHAPE_FILES

## Modules Overview

All functionalities are managed within the plot.py script, organized into the following modules:

### Plotting

* Plotter Class
  * Utilizes ax.imshow for filled contours, enhancing rendering speed and accuracy.
  * Labels are added using ax.text or ax.clabel with checks for accuracy in limited domains.
* Colormap Class
  * Supports both existing Matplotlib colormaps and custom colormaps from NumPy arrays.

### Processing

* Regridding
  * congrid, regrid, conservative_regrid: Handle various regridding tasks, optimized with Numba.
* Smoothing
  * savitzky_golay2d: Applies Savitzky-Golay smoothing to 2D data.
  * bandpass_filter: Implements various bandpass filters for data denoising.
* Scaling
  * bytscl: Scales images using IDL's bytscl function.

### Epilogue

* Annotate
  *	Finalizes plot images with annotations using Pillow.

### Utilities

* Paths
  * Manages file paths and directories within the project.

## How to Run the Scripts

The main plotting function's arguments follow the following structure:

```
python plot.py file_location plot_type --region_code 'XX' --cache_dir <path> --results_dir <path>
```

## Plot Types and Region Codes

### Plot Types (plot_type)

*	'ir8': Infrared channel 8
*	'radar': Radar
*	'aerosols': Aerosols
*	'precrain': Precipitable rain
*	'precsnow': Precipitable snow
*	'tpw': Total precipitable water
*	't2m': Air temperature 2 meters above sea-level
*	'slp': Surface-level pressure
*	'cape': Convective available potential energy
*	'vort500mb': 500 mb absolute vorticity
*	'winds10m': Surface wind speed

### Region Codes (region_code) 

Default value: -1 (all regions)

These codes are serialized in regions.json.

*	-1: All regions
*	0: All global projections
*	50: USA mapset
*	89: Central USA mapset
*	53: Mid-Atlantic mapset
*	51: Maryland mapset
*	49: North Atlantic mapset
*	65: Europe mapset
*	66: Asia mapset
*	67: Australia mapset
*	68: Africa mapset
*	69: South America mapset
*	70: North America mapset
*	73: East Pacific mapset
*	74: Indian Ocean mapset
*	75: West Atlantic mapset
*	34: North America projection
*	35: West Pacific projection
*	42: GOES East projection
*	43: GOES West projection
*	44: Meteosat-8 projection
*	46: Meteosat-10 projection
*	47: Himawari projection
*	71: Northern Hemisphere projection
*	72: Southern Hemisphere projection

Use region code globe to plot the lat/lon globe projection.

### cache_dir

The default value for this optional argument `data/`. If the cached assets are located elsewhere, use this argument to specify its location.

### results_dir

The default value for this optional arugment is `results/`. To change the location that the plots are saved to, use this argument to specify another folder. 

## Example Plots

PRECRAIN:

![plotall_precrain-usa_mapset](https://github.com/user-attachments/assets/7535f920-719b-4cb4-a006-4958c24327e5)

 
(PRECRAIN Plot for USA)

CAPE:

 ![plotall_cape-northamerica_mapset](https://github.com/user-attachments/assets/21ba3c42-1e67-4337-bf68-6be51b1b56ae)

(CAPE plot for North America)

T2M:
 
 ![plotall_t2m-epacific_mapset](https://github.com/user-attachments/assets/c0edce4b-3eb0-4573-82a3-8e0773d89854)

(T2M plot for Eastern Pacific)

Important Notes

*	A known issue exists with find_slp_mins in plotall_slp, which incorrectly identifies too many minima.
*	Additional code is needed to track forecast hours.
*	Data accumulation for plotall_precsnow and plotall_precrain is managed by pickling the array.
*	The aerosols plotting routine has a known error when generating the plots. 
*	The conservative regridding function has errors when attempting to process 2-dimensional NetCDF files
