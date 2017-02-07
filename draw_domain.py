"""
Creates a plot that shows the domain size of a given SpecifiedDataset
"""

import os
from datetime import datetime as dtime
from datetime import timedelta as tdelta
from ConfigParser import ConfigParser
import logging

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

#from PyNIO import Nio
#import objects as specdata_objects
#import utils as specdata_utils
from objects import SpecifiedForecastDataset
from pycane.postproc.viz.map import bling as map_bling

##
# PARAMETERS
##
# If True, this will plot what the actual grid looks like (i.e. where there
# are no missing values. Otherwise, it will just use the given field's lat
# and lon values, so it will draw a rectangle
plot_actual_grid = True
# Top-level directory containing Input Specification files
INSPEC_TOPDIR = "./conf/inspec"

def _default_log(log2stdout=logging.INFO, log2file=None, name='el_default'):
    global _logger
    if _logger is None:
        _logger = logging.getLogger(name)
        _logger.setLevel(log2stdout)
        msg_str = '%(asctime)s::%(name)s::%(lineno)s::%(levelname)s - %(message)s'
        #msg_str = '%(asctime)s::%(funcName)s::%(filename)s:%(lineno)s::%(levelname)s - %(message)s'
        # TODO : Use this one if running in debug mode
        #msg_str = '%(levelname)s::%(asctime)s::%(funcName)s::%(filename)s:%(lineno)s - %(message)s'
        msg_str = '%(levelname)s :: %(filename)s:%(lineno)s - %(message)s'
        date_format = "%H:%M:%S"
        formatter = logging.Formatter(msg_str, datefmt=date_format)
        if log2file is not None:
            fh = logging.FileHandler('log.txt')
            fh.setLevel(log2file)
            fh.setFormatter(formatter)
            _logger.addHandler(fh)
        if log2stdout is not None:
            ch = logging.StreamHandler()
            ch.setLevel(log2stdout)
            ch.setFormatter(formatter)
            _logger.addHandler(ch)
    return _logger
_logger = None

def create_map(conf, log=None):
    """  
    Create and return a Basemap object set according to options specified
    in ConfigParser object `conf'. Use extents from Nio dataset `dataset' 
    as a backup.
    """
    if log is None: log=_default_log()
    confmap = lambda param: conf.get("map_settings", param)
    if conf.has_option("map_settings", "map_extents"):
        extents = confmap("map_extents")
        extents = extents.replace("'", "").replace('"', '')
        extents_set = False
        try: 
            extents = confmap('map_extents').replace("'","").replace('"','')
            toks = extents.strip().split()
            south_lat = float(toks[0])
            west_lon = float(toks[1])
            north_lat = float(toks[2])
            east_lon = float(toks[3])
            extents_set  = True 
        except:
            log.info("Error parsing 'map_extents' section, will use dataset")
    if not extents_set:
        raise NotImplementedError()
        #(south_lat, west_lon, north_lat, east_lon) = \
        #            _get_extents_from_dataset(dataset, log) 
    padding = [0,0,0,0]
    if conf.has_option("map_settings", "padding"):
        padding = get_list(confmap("padding"), cast=float, minValues=4, maxValues=4)
    projection = confmap("projection")
    resolution = confmap("resolution")
    basemap = Basemap(llcrnrlon=west_lon-padding[1], llcrnrlat=south_lat-padding[2],
                      urcrnrlon=east_lon+padding[3], urcrnrlat=north_lat+padding[0],
                      projection=projection, resolution=resolution)
                  #   area_thresh=area_thresh)
    return basemap

def plot_data(lats, lons, outfile=None):
    conf = ConfigParser()
    conf.readfp(open("main.conf"))
    basemap = create_map(conf)
    map_settings = dict(conf.items("map_settings"))
    (x,y) = basemap(lons, lats)
    #cs = basemap.contourf(x, y, data, np.arange(data.min(),data.max(), 2))
    #cs = basemap.contourf(x, y, data, np.arange(245, 280)) # 500mb Temp
    basemap.plot(x, y, ls='none', marker='.', markersize=1, color='gray', zorder=0)
    map_bling.decorate_map(basemap, **map_settings)
    #cb = basemap.colorbar(cs, "top", pad=padding, ax=cb_ax, size=0.30)
    if outfile:
        plt.savefig(outfile)

if __name__ == "__main__":

    # 3km quasi-basin
    inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_5000x2500/2006090600/postprd"
    # 800x800
    #inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    #topdir =  "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_now/"

    start_date = dtime(year=2006, month=9, day=6, hour=0)
    domain = 1
    fhr = tdelta(0)
    # lat/lon will be obtained for this Field
    #        levType.levVal.levUnits.fieldName.fieldUnits
    field = "isobaric.500.hPa.air_temperature.K"

    src_data = SpecifiedForecastDataset(inspec, topdir, start_date, 
                                       fcst_offset=fhr, domain=domain,
                                       inspecsTopdir=INSPEC_TOPDIR,
                                       #log=log,
                                       #fieldTransforms=
                                       #coordTransforms=
                                       )
    src_data.set_coord_transform("make_2d", True) # for basemap
    src_data.set_coord_transform("pm_relative_longitude", True) # ditto

    
    (lev_type, lev_val, lev_units, field_name, units) = field.split(".")
    lev_val = int(lev_val)

    src_field = src_data.get_field(field_name, units, lev_type)
    src_slice = src_field.get_horz_slice(lev_val, levUnits=lev_units)
    
    src_data = np.ma.masked_equal(src_slice.data[:], src_field.missing_value)
    
    (lats,lons) = src_field.coords.latlons    # can't use this isince they may have been transposed (ie if dim_order was ji)
    if plot_actual_grid:
        lats[src_data.mask] = src_field.missing_value
        lons[src_data.mask] = src_field.missing_value

    plot_data(lats, lons, "foo.png")
    #import pdb ; pdb.set_trace()
    #plot_before_and_after(src_data, latsBefore, lonsBefore, regridded_data, lats, lons, "bar.png")
