#!/bin/bash

# # We assume that current folder contains this file, the input .csv-files, layersetter.xml, and styletemplate.xml

EPSG="3879"
BINSIZE="50"
WINDOWSCALE="2"
WORKSPACE="popdensity"

# # Create output folders
mkdir output_tiffs
mkdir figs
mkdir sld

# # Create tiffs with python script from csv-files in current directory
echo "Creating GeoTIFFs from files in current directory"
python heatmapper.py $EPSG $BINSIZE $WINDOWSCALE

cd output_tiffs

# Copy the files to the server, make sure to export variables USER and PASSWORD
echo "Copying them over SSH"
for tiff in *.tif; do
	scp ./$tiff $GEOSERVERUSER@geoserver.hel.fi:/home/$GEOSERVERUSER/geoserver-data/$WORKSPACE/$tiff
done

# Create workspace (not needed probably)
# echo "POST to create workspace $WORKSPACE"
# curl -u $USER:$PASSWORD -v -POST -H "Content-type:text/xml" -d "<workspace><name>$WORKSPACE</name></workspace>" http://geoserver.hel.fi/geoserver/rest/workspaces

# Set up each of the files one at a time as a coverage store on the geoserver
echo "POST and PUT to create new store, specify their data source file and the allowed styles"
for f in *.tif; do
	# new coverage store without any defining features
	curl -u $USER:$PASSWORD -v -POST -H "Content-type:text/xml" -d "<coverageStore><name>$f</name><workspace>$WORKSPACE</workspace><enabled>true</enabled></coverageStore>" http://geoserver.hel.fi/geoserver/rest/workspaces/$WORKSPACE/coveragestores

	# add a file to act as the store for the coverage store we created above, and publish it as a layer as well (??)
	curl -u $USER:$PASSWORD -v -XPUT -H "Content-type:text/plain" -d "file:///srv/geoserver/data/data/$GEOSERVERUSER/$WORKSPACE/$f" http://geoserver.hel.fi/geoserver/rest/workspaces/$WORKSPACE/coveragestores/$f/external.geotiff

	# set the default and allowed styles for this layer by name (all the tiffs have the same default and allowed styles)
	curl -u $USER:$PASSWORD -v -XPUT -H "Content-type:text/xml" -d @../layerstylesetter.xml http://geoserver.hel.fi/geoserver/rest/layers/${f%%.*}
done

# # PUT the modified style files into the server
cd ../sld # we're in output_tiffs folder
echo "PUT to update all the styles"
for f in *.xml; do
	curl -u $USER:$PASSWORD -v -XPUT -H "Content-type:application/vnd.ogc.sld+xml" -d @$f http://geoserver.hel.fi/geoserver/rest/styles/${f%%.*}
done