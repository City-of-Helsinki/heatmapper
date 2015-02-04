# heatmapper

A set of scripts for turning a dataset containing x- and y-coordinates and their respective weights/frequencies into a heatmap and publishing it on a geoserver. Python to handle reading the data & constructing the heatmap GeoTiffs and bash for uploading the file to communicating with geoserver.

Original purpose is to visualize population density using a data set where each datum is a building with coordinates and number of residents living in it.

### TODO
- rely on configuration files instead of hardcoding for:
  * names of aggregate fields (and the names of which fields to aggregate)
  * names of the fields to use as x- and y-coordinates
  * epsg code
  * geoserver file structure and authenticaton
  * parameters controlling the heatmapping (radius, binsize, etc.)
- polish colormap (or even make it possible to roll your own?)