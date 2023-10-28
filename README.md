# Experimental GEOS Forecasts with Python

GMAOpy has effectively served its purpose as a proof of concept and is now depreciated. All further development will be focused on GEOSpy. All data processing and visualization functionality has been refactored and adapted into a more modular architecture. Utilities for logging, file management, asset management, and various other supplemental tasks have been provided. Relevant regional and global attributes have been serialized with JSON and YAML for improved readability and workflow. Coastlines, borders, and other shapefile data are now reprojected and cached for faster plotting. Base images have similarly been reprojected and cached for faster plotting. Modules for reproducing these images from scratch have been provided as utilities. Bear in mind that the aforementioned process is time consuming. Data processing procedures are highly optimized. Conservative regridding takes approximately a second or less. Serialization and deserialization is very fast. Loading of base images requires only seconds. Less intricate plots have times reduced down to well under a minute. Even less optimized plots take less than two minutes, which is still much improved. On average, speed is equivalent to IDL or better. Accuracy is also equivalent barring one edge case.

#### plotall_ir8
average time 'sub': 0.0 hours 0.0 minutes 37.928170919418335 seconds
average time 'ortho': 0.0 hours 0.0 minutes 27.064877033233643 seconds
average time 'laea': 0.0 hours 0.0 minutes 58.88895893096924 seconds
average time 'geos': 0.0 hours 0.0 minutes 31.202008962631226 seconds
average time 'nsper': 0.0 hours 0.0 minutes 33.32458209991455 seconds

#### plotall_wxtype
average time 'sub': 0.0 hours 0.0 minutes 39.75187683105469 seconds
average time 'ortho': 0.0 hours 1.0 minutes 47.156792879104614 seconds
average time 'laea': 0.0 hours 0.0 minutes 48.601584911346436 seconds
average time 'geos': 0.0 hours 1.0 minutes 46.67049479484558 seconds
average time 'nsper': 0.0 hours 1.0 minutes 50.50699996948242 seconds

## TODO
- Implement other forecasts
- Implement parallelization
- Continue plot optimization
- Stay up to date on Cartopy and Matplotlib releases
- Refactor plot prologue feature
- Implement full directory traversal
- Test forecast Shell script
- Improve code coverage
- Implement command line arguments