import os

from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import numpy.lib.recfunctions as rf
from scipy.ndimage.filters import gaussian_filter
import sys

# quite a senseless encapsulation, done originally for parallelization purposes
def data2heatmap2multibandtif(cutoff, infile, epsg='3879', binsize=10, windowscale=1):
    
    try:
        styletemplate = open('styletemplate.xml').read()
    except e:
        print "need to have the style template named styletemplate.xml in the same folder"
        sys.exit(0)

    data = np.recfromcsv(infile, delimiter=',')
    #construct array with only the fields we want
    data = retain_relevant_fields(data)

    #prune invalid values
    data = data[data['nkoord'] != -1]

    #take these into account for lattice transformation
    min_nkoord = np.min(data['nkoord'])
    max_nkoord = np.max(data['nkoord'])
    min_ekoord = np.min(data['ekoord'])
    max_ekoord = np.max(data['ekoord'])
    pextent = np.arange(min_nkoord, max_nkoord)
    iextent = np.arange(min_ekoord, max_ekoord)

    #output image dimensions & resolutions
    nrows, ncols = (pextent.shape[0]/binsize, iextent.shape[0]/binsize)
    nres = (max_nkoord-min_nkoord)/float(nrows)
    eres = (max_ekoord-min_ekoord)/float(ncols)
    #...which is basically the geotransform
    geotransform = [min_ekoord, eres, 0, max_nkoord, 0, -nres]
    print '[min_ekoord, eres, 0, max_nkoord, 0, -nres] -> ', geotransform

    #===============================================
    #GENERATE HEATMAPS AND STORE THEM IN THE OUTFILE
    #===============================================
    # iterate through all strata, compute heatmap and write as a band in the image file (also a pyplot image)

    datakeys = [d for d in data.dtype.names if d not in ['nkoord', 'ekoord']]
    print 'Writing a layer for the following keys: ', datakeys
    
    #create output file and set parameters
    outfile = infile.split('.')[0]

    # outputraster = gdal.GetDriverByName('GTiff').Create('./output_tiffs/cutoff'+str(cutoff)+'_binsize'+str(binsize)+'_windowscale'+str(windowscale)+'_'+outfile+'.tif', ncols+1, nrows+1, len(datakeys), gdal.GDT_Float32)
    outputraster = gdal.GetDriverByName('GTiff').Create('./output_tiffs/'+outfile+'.tif', ncols+1, nrows+1, len(datakeys), gdal.GDT_Float32)
    outputraster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    # ETRS89-GK25FIN maps to 3879
    # before 2012 everything is in KKJ2, maps to 2392
    srs.ImportFromEPSG(int(epsg))
    # source: http://www.maanmittauslaitos.fi/sites/default/files/tiedostolataukset/kartat/koordinaatit/epsg_koodit.pdf
    outputraster.SetProjection( srs.ExportToWkt() )

    # strata_fieldnames = ['asyht']
    for i,field in enumerate(datakeys):

        #build 2d lattice to cumulate the frequencies on (even though every row is unique to that position in the lattice)
        #NB! this is indexed by [0,4xxx] instead of the actual coordinate
        lattice = np.zeros(((pextent.shape[0]/binsize)+1, (iextent.shape[0]/binsize)+1))

        #cumulate
        # print 'cumulating frequencies'
        for building in data:
            lattice[(building['nkoord']-min_nkoord)/binsize][(building['ekoord']-min_ekoord)/binsize] = building[field]
        #now we have a "sparse" (i.e. most values are 0) 2d lattice where each datum (x,y) has value that represents the frequency of people living in that coordinate's building
        
        #cut off elements where population in a house is below threshold
        highpassed = np.copy(lattice)
        highpassed[highpassed < cutoff] = 0

        #what a monster. add up to +-10% of the value (rounded up) as noise to each element, and then round down
        noised = np.array([np.floor(datum + np.random.uniform(-(np.ceil(datum*0.1)),np.ceil(datum*0.1))) for datum in highpassed])

        #FUNCTION GIVES SKEWED RESULTS because of ndi.gaussian_filter's truncation algorithm. it always adds one unit after truncation to make the function nicely trail to 0
        # this skews results quite quickly ESPECIALLY with large bins because the one unit is squared when computing the area
        #...

        # so it would go something like this:
        # (sigma * truncate) + 1 = r_window
        # sigma * truncate = r_window - 1
        # sigma = (r_window - 1) / truncate

        scale = windowscale*100/binsize
        r_window = np.sqrt(scale*scale / np.pi)
        #parameter for controlling when the filter is truncated, and thus effects our window size. default is 4, but just explicating it here
        truncate = 4.0
        # bandwidth_sigma = r_window / truncate
        bandwidth_sigma = (r_window - 1) / truncate
        print 'smoothing: gaussian with sigma ', bandwidth_sigma
        preprocessed = noised #which preprocessing stages were done

        # draw .png pyplot figures for easy previewing of of the geotiff
        # plt.figure()
        smoothed = gaussian_filter(preprocessed, bandwidth_sigma, truncate=truncate)
        # plt.imshow(smoothed, origin='lower')
        # plt.colorbar()
        # fname = '_binsize'+str(binsize)+'_cutoff'+str(cutoff)+'_bw'+str(bandwidth_sigma)+'_windowscale'+str(windowscale)+'_'+field+'.png'
        # plt.savefig('./figs/windowsize_controlled_geotiff'+fname)
        # plt.close('all')

        print "writing " + field + " to band ",i+1
        outputraster.GetRasterBand(i+1).WriteArray(np.flipud(smoothed)) #flip updown

        #write corresponding style file
        with open('./sld/popdensity_%s.xml' % (field), 'w') as f:
            lowlimit = np.percentile(smoothed[smoothed > 0], 2)
            highlimit = np.percentile(smoothed[smoothed > 0], 98)
            midlimit = lowlimit + ((highlimit - lowlimit) * 0.75) #we want the midlimit to be betwee min & max. 
            # midlimit = np.mean(smoothed[smoothed > 0])
            print lowlimit, midlimit, highlimit
            f.write(styletemplate % {'stratum':field, 'band_n':(i+1), 'minlimit':np.min(smoothed[smoothed > 0]), 'lowlimit':lowlimit, 'midlimit':midlimit, 'highlimit':highlimit})


    outputraster = None

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


if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print 'Usage: python heatmapper.py source_srs'
        sys.exit(0)

    # inpath = sys.argv[1]
    # outpath = sys.argv[2]
    epsg = sys.argv[1]
    binsize = int(sys.argv[2])
    windowscale = int(sys.argv[3]) # for visual purposes

    #read data
    infiles = [f for f in os.listdir('.') if 'csv' in  f.split('.')[-1]]

    for infile_name in infiles[::-1]: #reverse to process latest file last
        data2heatmap2multibandtif(5, infile_name, epsg, binsize=binsize, windowscale=windowscale)