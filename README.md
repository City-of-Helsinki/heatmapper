# heatmapper

Depends on a few libraries which need some hoop jumping to get working inside a virtualenv (SciPy is huge and for such a small use, so removing is TODO):
-Installing gdal-python only installs the bindings, but the library is still needed. You need to go  `sudo apt-get install libgdal-dev`.
-SciPy: need to (or at least *I* needed to) `sudo apt-get install libblas-dev liblapack-dev libblas3gf libc6 libgcc1 libgfortran3 liblapack3gf libstdc++6 build-essential gfortran python-all-dev libatlas-base-dev`


Scripts for turning a dataset containing x- and y-coordinates and their respective weights/frequencies into a heatmap and publishing it on a geoserver. Python to handle reading the data & constructing the heatmap GeoTiffs and bash for uploading the file to and communicating with geoserver via it's REST interface.

Use by copying the data-csv into the same directory and then going `export USER=your_geoserver_username; export PASSWORD=your_geoserver_password; sh publisher.sh`.

Original purpose is to visualize population density using a data set where each datum is a building with coordinates and number of residents living in it.

### TODO
- rely on configuration files instead of hardcoding for:
  * names of aggregate fields (and the names of which fields to aggregate)
  * names of the fields to use as x- and y-coordinates
  * epsg code
  * geoserver file structure and authenticaton
  * parameters controlling the heatmapping (radius, binsize, etc.)
- polish colormap (or even make it possible to roll your own?)