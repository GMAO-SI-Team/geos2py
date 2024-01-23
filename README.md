# Experimental GEOS Forecasts with Python

GMAOpy has effectively served its purpose as a proof of concept and is now depreciated. All further development will be focused on GEOSpy. All data processing and visualization functionality has been refactored and adapted into a more modular architecture. Utilities for logging, file management, asset management, and various other supplemental tasks have been provided. Relevant regional and global attributes have been serialized with JSON and YAML for improved readability and workflow. Coastlines, borders, and other shapefile data are now reprojected and cached for faster plotting. Base images have similarly been reprojected and cached for faster plotting. Modules for reproducing these images from scratch have been provided as utilities. Bear in mind that the aforementioned process is time consuming. Data processing procedures are highly optimized. Conservative regridding takes approximately a second or less. Serialization and deserialization is very fast. Loading of base images requires only seconds. Plotting is even further optimized with some average times dipping down to around 10 seconds. Bugs limiting accuracy have also been significantly reduced. The only inaccuracies to resolves are in plotall_precsnow. As for the remaining plots, plotall_slp will be complete as soon as the issue regarding locating sea level pressure minima is resolved and plotall_aerosols will be complete as soon as the blending issue is resolved. Otherwise, the plotting is virtual complete with just some stretch goals remaining. Additionally, the prologue component needs to be updated to begin handling differing dates and directories as well as the logging, which shouldn't be too arduous.  

## Environment Setup
```sh
ml use -a /home/mathomp4/modulefiles-SLES12
ml python/MINIpyD
```

If you are working locally, the dependencies are:
- python==3.11.4
- matplotlib==3.8.0
- cartopy==0.22.0
- pillow==10.1.0
- netcdf4==1.6.5
- numpy==1.24.4
- scipy==1.11.1
- numba==0.57.1
- pyogrio==0.6.0
- geopandas==0.13.2
- pandas==2.0.3
- contourpy=1.1.0

## Performance

#### plotall_ir8
- average time 'sub': 0.0 hours 0.0 minutes 10.232033967971802 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 11.514650583267212 seconds
- average time 'laea': 0.0 hours 0.0 minutes 13.195988631248474 seconds
- average time 'geos': 0.0 hours 0.0 minutes 14.224856567382812 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 13.688894152641296 seconds

#### plotall_wxtype
- average time 'sub': 0.0 hours 0.0 minutes 39.75187683105469 seconds
- average time 'ortho': 0.0 hours 1.0 minutes 47.156792879104614 seconds
- average time 'laea': 0.0 hours 0.0 minutes 48.601584911346436 seconds
- average time 'geos': 0.0 hours 1.0 minutes 46.67049479484558 seconds
- average time 'nsper': 0.0 hours 1.0 minutes 50.50699996948242 seconds

#### plotall_radar
- average time 'sub': 0.0 hours 0.0 minutes 14.13889753818512 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 15.65565550327301 seconds
- average time 'laea': 0.0 hours 0.0 minutes 16.701356053352356 seconds
- average time 'geos': 0.0 hours 0.0 minutes 16.995589017868042 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 17.368183076381683 seconds

#### plotall_precrain
- average time 'sub': 0.0 hours 0.0 minutes 8.57942807674408 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 11.32456088066101 seconds
- average time 'laea': 0.0 hours 0.0 minutes 13.072884774208068 seconds
- average time 'geos': 0.0 hours 0.0 minutes 12.116543436050415 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 12.069870054721832 seconds

#### plotall_precsnow
- average time 'sub': 0.0 hours 0.0 minutes 9.807256579399109 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 10.264758110046387 seconds
- average time 'laea': 0.0 hours 0.0 minutes 11.609730076789855 seconds
- average time 'geos': 0.0 hours 0.0 minutes 11.448821258544921 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 10.888549625873566 seconds

#### plotall_aerosols
**pending**

#### plotall_slp
**pending**

#### plotall_t2m
- average time 'sub': 0.0 hours 0.0 minutes 9.49286139011383 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 12.913065075874329 seconds
- average time 'laea': 0.0 hours 0.0 minutes 13.569564175605773 seconds
- average time 'geos': 0.0 hours 0.0 minutes 13.966146612167359 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 14.035657703876495 seconds

#### plotall_tpw
- average time 'sub': 0.0 hours 0.0 minutes 11.00346064567566 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 12.062952399253845 seconds
- average time 'laea': 0.0 hours 0.0 minutes 14.307450294494629 seconds
- average time 'geos': 0.0 hours 0.0 minutes 13.70802993774414 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 13.780909240245819 seconds

#### plotall_cape
- average time 'sub': 0.0 hours 0.0 minutes 10.370289087295532 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 12.009145975112915 seconds
- average time 'laea': 0.0 hours 0.0 minutes 13.465498185157776 seconds
- average time 'geos': 0.0 hours 0.0 minutes 13.688392782211304 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 12.895038187503815 seconds

#### plotall_vort500mb
- average time 'sub': 0.0 hours 0.0 minutes 24.911166191101074 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 42.29221761226654 seconds
- average time 'laea': 0.0 hours 0.0 minutes 50.836144137382504 seconds
- average time 'geos': 0.0 hours 0.0 minutes 43.15956325531006 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 43.311415791511536 seconds

#### plotall_winds10m
- average time 'sub': 0.0 hours 0.0 minutes 10.466203808784485 seconds
- average time 'ortho': 0.0 hours 0.0 minutes 12.390382051467896 seconds
- average time 'laea': 0.0 hours 0.0 minutes 14.286669969558716 seconds
- average time 'geos': 0.0 hours 0.0 minutes 15.273489904403686 seconds
- average time 'nsper': 0.0 hours 0.0 minutes 14.948720812797546 seconds

## TODO
- Fix storm location bug in plotall_slp
    - Likely due to issues with the ideal bandpass filter implementation
- Fix data inaccuracy bug in plotall_precsnow
    - Likely due to a masking issue
- Get better blending with plotall_aerosols
    - Still not sure how to improve this just yet
- Stay up to date on Cartopy and Matplotlib releases
- Update prologue to handle different dates and directories
- Figure out how the logging should work
- Complete conditional Python or IDL functionality for forecast Shell script
- Continue to add more documentation
- Refactor for better cohesion where worth it

### Caching
A lot of performance improvements are resultant from pre-emptive caching. Chief cached elements include coastlines, borders, natural earth base images, and so on. These assets are pretransformed to suit the corresponding plot and region and are named accordingly. It is strongly encouraged that anyone running these scripts simply copy the cached images from a folder in this repo or on Discover. There is a script responsible for caching the assets, but beware that it takes close to a day to complete caching for all plot types and all regions in one go. The total cache takes up 1.4 GB of disk space.

#### master_cache.py
When running this script, it expects two commmand line arguments: "plot_type" and "region".
- plot_type: name of plot such as plotall_ir8 or plotall_radar or all for all plot_types
- region: numeric region code or all for all regions (_mapset and _proj)

I have added an optional argument: "--element",
- element: expects "features", which just grabs features like coastlines, borders, roads, etc or "base_image", which just grabs the relevant natural Earth background

Generally, I recommended just copying the zip file of all the cached elements.
```sh
cp /discover/nobackup/qcambrel/geospy/cache.zip destination_directory
```

## Modules
- plotting
    - plots.Plotter class
        - Filled contours are depicted with ax.imshow in place of ax.contourf or ax.pcolormesh. This is all around faster and less prone to corrupted or inaccurately reprojected data.
        - Labels are plotted using ax.text based on interpolated coordinates unless the labels belong to contour lines then ax.clabel is used instead. The former is prone to error but using checks for limited domains or null reprojected values serves as an effective defense.   
    - colormaps.Colormap class
        - Colormaps are either a pre-existing Matplotlib colormap or generated from a NumPy array of rgb values using colors.ListedColormap then normalized based on the levels for the respective plot.
- processing
    - regridding
        - congrid: resample data array to specified dimensions
        - regrid: regrid data array based on specified method and gridspec
        - conservative_regrid: implement first order conservative regridding based on gridspec; vectorized by Numba
        - read_nt: returns tile file dimensions
        - read_tile_file: returns tile file to prepare gridspec
    - smoothing
        - savitzky_golay2d: implements the Savitzky Golay 2D algorithm to smooth noisy contours
        - ideal_bandpass_filter: implements an ideal bandpass filter to denoise data for minima search
- epilogue
    - annotate: composes the final plot image using Pillow

## How to run the scripts
```sh
python plot_type.py year month day hour minute tag s_tag data_dir stream --region=region_code --f_date=f_date
```

#### Plot Types
- plotall_ir8
- plotall_radar
- plotall_wxtype
- plotall_aerosols
- plotall_precrain
- plotall_precsnow
- plotall_tpw
- plotall_t2m
- plotall_slp
- plotall_cape
- plotall_vort500mb
- plotall_winds10m

#### Region Codes
These are all serialized in regions.json.
- region code -1: all regions
- region code 0: all global projections
- region code 50: usa_mapset
- region code 89: centralusa_mapset
- region code 53: midatlantic_mapset
- region code 51: maryland_mapset
- region code 49: northatlantic_mapset
- region code 65: europe_mapset
- region code 66: asia_mapset
- region code 67: australia_mapset
- region code 68: africa_mapset
- region code 69: southamerica_mapset
- region code 70: northamerica_mapset
- region code 73: epacific_mapset
- region code 74: indianocean_mapset
- region code 75: westatlantic_mapset
- region code 34: northamerica_proj
- region code 35: westpacific_proj
- region code 42: goeseast_proj
- region code 43: goeswest_proj
- region code 44: meteosat8_proj
- region code 46: meteosat10_proj
- region code 47: himawari_proj
- region code 71: nh_proj
- region code 72: sh_proj

### Using canned values

#### plotall_ir8
test arguments
```sh
python plotall_ir8.py 2023 05 29 12 00 G5GMAO f5295_fp-20230529_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_12z' && python plotall_ir8.py 2023 05 29 12 00 G5GMAO f5295_fp-20230529_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_12z'
```

#### plotall_wxtype
test arguments
```sh
python plotall_wxtype.py 2023 05 29 00 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_wxtype.py 2023 05 29 00 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```
#### plotall_radar
test arguments
```sh
python plotall_radar.py 2023 05 29 12 00 G5GMAO f5295_fp-20230529_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_12z' && python plotall_radar.py 2023 05 29 12 00 G5GMAO f5295_fp-20230529_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_12z'
```

#### plotall_precrain
test arguments
```sh
python plotall_precrain.py 2023 05 29 00 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_precrain.py 2023 05 29 00 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_precsnow
test arguments
```sh
python plotall_precsnow.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_precsnow.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_cape
test arguments
```sh
python plotall_cape.py 2023 12 23 12 00 G5GMAO f5295_fp-20231223_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20231223_12z' && python plotall_tpw.py 2023 12 23 12 00 G5GMAO f5295_fp-20231223_12z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20231223_12z'
```

#### plotall_winds10m
test arguments
```sh
python plotall_winds10m.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_winds10m.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_vort500mb
test arguments
```sh
python plotall_vort500mb.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_vort500mb.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_t2m
test arguments
```sh
python plotall_t2m.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_t2m.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_tpw
test arguments
```sh
python plotall_tpw.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_tpw.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_slp
test arguments
```sh
python plotall_slp.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_slp.py 2023 05 29 00 30 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```

#### plotall_aerosols
test arguments
```sh
python plotall_aerosols.py 2023 10 11 12 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='-1' --f_date='20230529_00z' && python plotall_aerosols.py 2023 10 11 12 00 G5GMAO f5295_fp-20230529_00z /discover/nobackup/qcambrel/geospy/data inst1_2d_asm_Mx --region='0' --f_date='20230529_00z'
```