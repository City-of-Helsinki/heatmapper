import os, sys
import ConfigParser #TODO lol

from osgeo import gdal
from osgeo import osr

import numpy as np
import numpy.lib.recfunctions as rf
from scipy.ndimage.filters import gaussian_filter

def retain_relevant_fields(data):
    #TODO read these from a conf file
    aggregate_fields = {}
    aggregate_fields['aggr0_6'] = data['ika0'] + data['ika1'] + data['ika2'] + data['ika3'] + data['ika4'] + data['ika5'] + data['ika6']
    aggregate_fields['aggr7_12'] = data['ika7'] + data['ika8'] + data['ika9'] + data['ika10'] + data['ika11'] + data['ika12']
    aggregate_fields['aggr13_17'] = data['ika13'] + data['ika14'] + data['ika15'] + data['ika16'] + data['ika17']
    aggregate_fields['aggr18_29'] = data['ika18'] + data['ika19'] + data['ika20'] + data['ika21'] + data['ika22'] + data['ika23'] + data['ika24'] + data['ika25_29']
    aggregate_fields['aggr30_64'] = data['ika30_34'] + data['ika35_39'] + data['ika40_44'] + data['ika45_49'] + data['ika50_54'] + data['ika55_59'] + data['ika60_64']
    aggregate_fields['aggr64_'] = data['ika65_69'] + data['ika70_74'] + data['ika75_79'] + data['ika80_84'] + data['ika85_89'] + data['ika90_94'] + data['ika95_']

    #add the fields to data rec array
    augmented_data = rf.rec_append_fields(data, aggregate_fields.keys(), aggregate_fields.values())
    #...and add these fields later. they are here for their column names, but the line above would cause an exception if we added them before
    aggregate_fields['asyht'] = data['asyht']
    aggregate_fields['ruots'] = data['ruots']
    aggregate_fields['ekoord'] = data['ekoord']
    aggregate_fields['nkoord'] = data['nkoord']

    # drop all fields whose names are not in aggregate_fields
    fields2drop = [d for d in data.dtype.names if d not in aggregate_fields.keys()]
    return rf.rec_drop_fields(augmented_data, fields2drop)

def read_file_prune_fields_clean_values(infile_name, x_name, y_name):
    data = np.recfromcsv(infile_name, delimiter=',')
    data = retain_relevant_fields(data)
    data = data[data[y_name] != -1] #this takes care of the garbage rows
    return data[y_name], data[x_name], rf.rec_drop_fields(data, [y_name, x_name])

def compute_geotransform(x, y, binsize=1):
    #take these into account for lattice transformation
    min_y = np.min(y)
    max_y = np.max(y)
    min_x = np.min(x)
    max_x = np.max(x)
    pextent = max_y - min_y
    iextent = max_x - min_x

    #output raster image dimensions & resolutions
    nrows, ncols = (pextent/binsize, iextent/binsize)
    nres = (max_y-min_y)/float(nrows)
    eres = (max_x-min_x)/float(ncols)
    #...which is basically the geotransform
    geotransform = [min_x, eres, 0, max_y, 0, -nres]

    return ncols, nrows, geotransform

def createRaster(outfilename, nrows, ncols, geotransform, EPSG, n_bands):
    outputraster = gdal.GetDriverByName('GTiff').Create('./output_tiffs/'+outfilename+'.tif', ncols+1, nrows+1, n_bands, gdal.GDT_Float32)
    outputraster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    # ETRS89-GK25FIN maps to 3879
    # before 2012 everything is in KKJ2, maps to 2392
    srs.ImportFromEPSG(int(epsg))
    # source: http://www.maanmittauslaitos.fi/sites/default/files/tiedostolataukset/kartat/koordinaatit/epsg_koodit.pdf
    outputraster.SetProjection( srs.ExportToWkt() )
    return outputraster

def heatmap(x, y, weights, nrows, ncols, cutoff=0, noise=0.0, binsize=1, windowscale=1, windowarea_squareroot=100):
    '''
    Generates a lattice and adds weights[i] to coordinate ((x/binsize)[i], (y/binsize)[i]). Reverts all cell with value < cutoff to 0 and adds noise to the freq at each cell.
    Then generates the heatmap by smoothing with a gaussian kernel which has area windowarea_squareroot(default=100)**2 coordinate units.
    Binsize affects this, as we need the window area to be defined in terms of the original coordinate units.

    E.g. we expect the x and y to correspond to meters and want a kernel with an area of 100m*100m (i.e. 1ha). However, a lattice where each cell is 1m*1m is too fine grained,
    so we can set binsize=10, meaning that each cell represents 10m*10m. This way we will still be able to make sure that the area of the kernel window is windowarea_squareroot**2
    regardless of the "size" of the cells in the lattice where the weights are aggregated.
    '''
    lattice = np.zeros((nrows+1, ncols+1))
    for i,freq in enumerate(weights):
        ycell = (y[i]-np.min(y)) / binsize
        xcell = (x[i]-np.min(x)) / binsize
        # print ycell, xcell
        lattice[ycell][xcell] = freq

    highpassed = np.copy(lattice)
    highpassed[highpassed < cutoff] = 0

    #what a monster. for each element: add or subtract noise that's up to (noise*100)% of the value. 0 adds nothing, 0.1 adds 10%
    noised = np.array([np.floor(datum + np.random.uniform(-(np.ceil(datum*noise)),np.ceil(datum*noise))) for datum in highpassed])

    #FUNCTION GIVES SKEWED RESULTS because of ndi.gaussian_filter's truncation algorithm. it always adds one unit after truncation to make the function nicely trail to 0
    # this skews results quite quickly ESPECIALLY with large bins because the one unit is squared when computing the area
    #...
    # so it would go something like this:
    # (sigma * truncate) + 1 = r_window
    # sigma = (r_window - 1) / truncate

    scale = windowscale*windowarea_squareroot/binsize
    r_window = np.sqrt(scale*scale / np.pi)
    #parameter for controlling when the filter is truncated, and thus effects our window size. default is 4, but just explicating it here
    truncate = 4.0
    # bandwidth_sigma = r_window / truncate
    bandwidth_sigma = (r_window - 1) / truncate

    smoothed = gaussian_filter(noised, bandwidth_sigma, truncate=truncate)

    return smoothed

if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print 'Usage: python heatmapper.py source_srs binsize windowscale'
        sys.exit(0)

    # Check requisite xml existence
    try:
        styletemplate = open('styletemplate.xml').read()
    except e:
        print "need to have the style template named styletemplate.xml in the same folder"
        sys.exit(1)

    epsg = sys.argv[1]
    binsize = int(sys.argv[2])
    windowscale = int(sys.argv[3]) # for visual purposes

    # do the stuff for all csv-files in this folder
    infiles = [f for f in os.listdir('.') if 'csv' in  f.split('.')[-1]]
    for infile_name in infiles[::-1]: #reverse to process latest file last
        print 'processing file '+infile_name
        # discard all those sensitive and useless columns
        yarr, xarr, datarec = read_file_prune_fields_clean_values(infile_name, 'ekoord', 'nkoord')

        # get some simple parameters
        ncols, nrows, geotransform = compute_geotransform(xarr, yarr, binsize)

        outfilename = infile_name.split('.')[0]
        raster = createRaster(outfilename, nrows, ncols, geotransform, epsg, len(datarec.dtype.names))

        # for all columns we have left
        for i,field in enumerate(datarec.dtype.names):
            # compute heatmap and save it
            data = datarec[field]
            heatmap_lattice = heatmap(xarr, yarr, data, nrows, ncols, cutoff=5, noise=0.1, binsize=binsize, windowscale=windowscale, windowarea_squareroot=100) #array like
            raster.GetRasterBand(i+1).WriteArray(np.flipud(heatmap_lattice))

            # save corresponding style file
            with open('./sld/popdensity_%s.xml' % (field), 'w') as f:
                lowlimit1 = np.percentile(heatmap_lattice[heatmap_lattice > 0], 1)
                lowlimit2 = np.percentile(heatmap_lattice[heatmap_lattice > 0], 2)
                lowlimit3 = np.percentile(heatmap_lattice[heatmap_lattice > 0], 3)
                lowlimit4 = np.percentile(heatmap_lattice[heatmap_lattice > 0], 10)
                # lowlimit = np.min(heatmap_lattice[heatmap_lattice > 0])
                highlimit = np.percentile(heatmap_lattice[heatmap_lattice > 0], 95)
                midlimit = lowlimit4 + ((highlimit - lowlimit4) * 0.75) #we want the midlimit to be between min & max.

                print 'column '+field+', low4:', lowlimit4,' mid:', midlimit,' high:', highlimit
                f.write(styletemplate % {'stratum':field, 'band_n':(i+1), 'minlimit':np.min(heatmap_lattice[heatmap_lattice > 0]), 'lowlimit1':lowlimit1, 'lowlimit2':lowlimit2, 'lowlimit3':lowlimit3, 'lowlimit4':lowlimit4, 'midlimit':midlimit, 'highlimit':highlimit})