"""
NOTES: For saving Regrid (weights) to database, if missingVal are different, we
should re-calculate the weights
 -> would require knowing the missing value points in addition to the coords

TODO
 - database of Regrid weights
 - colormap that's more applicable to diff plots (i.e. center is always 0 and white or light gray)

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
import ESMF
from objects import SpecifiedForecastDataset
from pycane.postproc.viz.map import bling as map_bling


# enable logging
esmpy = ESMF.Manager(debug=True)

plot_it = True

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

def plot_before_and_after(dataBefore, latsBefore, lonsBefore, dataAfter, latsAfter, lonsAfter, outfile):
    plt.subplot(1,2,1)
    plot_data(dataBefore, latsBefore, lonsBefore)
    plt.subplot(1,2,2)
    plot_data(dataAfter, latsAfter, lonsAfter)
    plt.savefig(outfile)

def plot_data(data, lats, lons, outfile=None):
    conf = ConfigParser()
    conf.readfp(open("test.conf"))
    basemap = create_map(conf)
    map_settings = dict(conf.items("map_settings"))
    (x,y) = basemap(lons, lats)
    #cs = basemap.contourf(x, y, data, np.arange(data.min(),data.max(), 2))
    #cs = basemap.contourf(x, y, data, np.arange(245, 280)) # 500mb Temp
    cs = basemap.contourf(x, y, data, np.linspace(245, 280)) # 500mb Temp
    map_bling.decorate_map(basemap, **map_settings)
    #cb = basemap.colorbar(cs, "top", pad=padding, ax=cb_ax, size=0.30)
    cb = basemap.colorbar(cs, "top", pad=0.3, size=0.30)
    if outfile:
        plt.savefig(outfile)

INSPEC_TOPDIR = "./inspec"
if __name__ == "__main__":

    # ** NOTE: _new_beta_pl and gamma have different start dates, so we're comparing
    #         different dates here
    field = "isobaric.500.hPa.air_temperature.K"
    src_inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    src_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_now/"

    #dest_inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    #dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_5000x2500/postprd_now"
    #dest_inspec = ["grib/upp/upp_nmb_grib1.conf"]
    #dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_grib1"
    #dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/_new_beta_pl/800x800/postprd.orig"
    dest_inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    dest_topdir =  "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_now/"

    start_date = dtime(year=2006, month=9, day=6, hour=0)
    domain = 1
    fhr = tdelta(0)

    src_data = SpecifiedForecastDataset(src_inspec, src_topdir, start_date, 
                                       fcst_offset=fhr, domain=domain,
                                       inspecsTopdir=INSPEC_TOPDIR,
                                       #log=log,
                                       #fieldTransforms=
                                       #coordTransforms=
                                       )
    src_data.set_coord_transform("make_2d", True) # ESMF wants 2d, I think

    dest_data = SpecifiedForecastDataset(dest_inspec, dest_topdir, start_date,
                                         inspecsTopdir=INSPEC_TOPDIR)
    dest_data.set_coord_transform("make_2d", True)

    
    (lev_type, lev_val, lev_units, field_name, units) = field.split(".")
    lev_val = int(lev_val)

    src_field = src_data.get_field(field_name, units, lev_type)
    dest_field = dest_data.get_field(field_name, units, lev_type)

    src_slice = src_field.get_horz_slice(lev_val, levUnits=lev_units)
    dest_slice = dest_field.get_horz_slice(lev_val, levUnits=lev_units)
    
    src_esmf_field = src_slice.get_esmf_field()
    dest_esmf_field = dest_slice.get_esmf_field()

    regrid = ESMF.Regrid(src_esmf_field, dest_esmf_field, 
                         src_mask_values=src_field.missing_value, 
                         dst_mask_values=dest_field.missing_value,
                         unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE)
    #jzanote : need the unmapped_action setting or it will fail due to their 
    # being dest. points that do not map to source points
    #import pdb ; pdb.set_trace()
    dest_esmf_field = regrid(src_esmf_field, dest_esmf_field)
    
    #diff = src_esmf_field.data[:] - dest_esmf_field.data[:]
    #print diff.min(), diff.max(), diff.mean()

    src_data = np.ma.masked_equal(src_slice.data[:], src_field.missing_value)
        
    #regridded_data = np.ma.masked_equal(src_esmf_field.data[:], src_field.missing_value)
    regridded_data = np.ma.masked_equal(src_esmf_field.data[:], src_field.missing_value)
    dest_data = np.ma.masked_equal(dest_esmf_field.data[:], dest_field.missing_value)
    diff = dest_data - regridded_data
    print diff.shape, diff.min(), diff.max(), diff.mean()

    #src_data.inputspec
    if plot_it:
        #(lats,lons) = dest_field.coords.latlons    # can't use this isince they may have been transposed (ie if dim_order was ji)
        lats = dest_esmf_field.grid.get_coords(1, dest_slice.staggerloc_2d)
        lons = dest_esmf_field.grid.get_coords(0, dest_slice.staggerloc_2d)
        plot_data(diff, lats, lons, "foo.png")
        #import pdb ; pdb.set_trace()
        (latsBefore, lonsBefore) = src_slice.coords.latlons
        # Transpose src data if necessary
        if src_slice.dim_order == "ji":
            src_data = src_data.T
        if src_slice.coords.lat_dim_order == "ji":
            latsBefore = latsBefore.T
        if src_slice.coords.lon_dim_order == "ji":
            lonsBefore = lonsBefore.T
        # Change to -180->180 if necessary
        if np.max(lonsBefore) > 180.0 or np.min(lonsBefore) < -180.0:
            lonsBefore += 180. 
            lonsBefore %= 360. 
            lonsBefore -= 180. 
            #src_data = np.roll(src_data, 180) 

        plot_before_and_after(src_data, latsBefore, lonsBefore, regridded_data, lats, lons, "bar.png")
