#!/usr/bin/env python
"""
WARNING: This is alpha software. It is provided in the hope that it will be 
useful, but there is no guarantee that it will work, provide correct results, 
(not) destroy your computer, etc. Use at your own risk.

Generate a set of images from one or more input files.
The images to generate are specified in the given <config file(s)>. The input 
file(s) can be specified in the command line or specified via the 
config file. See sample config file for details. It is possible to specify
multiple config files, either by separating them by commas _only_ (no spaces!). e.g. one.conf,two.conf
or as separate arguments with ending in .conf, .cfg, or .config

TODO
 - Read fhr from metadata when possible - add get_fhr() to encapsulate
 - horz interp for vectors ?

"""

import os
import logging 
from ConfigParser import ConfigParser, NoSectionError
from optparse import OptionParser
import sys
from datetime import datetime as dtime
from datetime import timedelta as tdelta
from distutils.util import strtobool
import glob
import importlib 
import copy
import subprocess

import numpy as np
from cfunits import Units # !! Must be loaded before matplotlib or SEGFAULTs !!
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.animation import FFMpegFileWriter,ImageMagickFileWriter
from matplotlib.cm import cmap_d
from PyNIO import Nio
from nwpy.viz.map import bling as map_bling
#from pycane.interpolation.vertical_interpolation.hwrf_vertical_interp import interp_press2press_lin_wrapper
from pycane.postproc.tracker import utils as trkutils
from  nwpy.postproc.io.specdata import objects as specdata_objects
from nwpy.postproc.io.specdata.objects import SpecifiedForecastDataset as SpecData
from nwpy.postproc.io.specdata import utils as specdata_utils

#
# GLOBALS
#
_logger = None
_model2standard_fields = None
# 2-d array containing the pressure levels (2-d) of the input source. 
# the convention is to use 'None' for non-isobaric levels
_model_plevs = None

#
# CONSTANTS
#

##
#	Basic constants
#
# Value used for config settings that should be ignored
IGNORED_OPTION = "-9999"
# how high to make each figure
INCHES_PER_FIGURE = 4.2
# Size of title
FIGURE_TITLE_SIZE = 26 # points
# Width of FIgure
FIGURE_WIDTH = 8.0
# Color for text of all colorbar labels
CBAR_LABEL_COLOR = "#16a085"

# New
INSPEC_PATH = os.path.join(os.getcwd(), "inspec")


#
# FUNCTIONS
#
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

def __parse_args():
    usage="usage: %prog [options] <config file1[,config_file2,...]> "\
          "[config_fileN ... config_fileN+M] [input_file1,...input_fileN]"\
          "".format() + __doc__
    parser = OptionParser(usage=usage)
    #parser.add_option("-f", "--file", dest="input_file", optional=False)
    #parser.add_option("-c", "--config", dest="config_file")
    parser.add_option("-l", "--log-level", dest="log_level", default=logging.INFO,
                      help="(0-100 or constant defined in the logging module")
    (options, args) = parser.parse_args()
    
    input_files = None
    #import pdb ; pdb.set_trace()
    if len(args) < 1:
        print usage
        sys.exit(1)
    
    config_files = args[0].split(",")
    if not os.path.exists(config_files[0]):
        raise Exception("Specified config file '{cf}' does not exist"
                        .format(cf=config_files[0]))

    # Parse the rest of the args. Values ending in .conf, .config, .cfg will
    # be processed as additional config files. All others will be processed
    # as data files to create products from
    if len(args) > 1:
        #input_files = args[1:]
        input_files = []
        other_args = args[1:]
        for arg in other_args:
            ext = os.path.splitext(arg)[1]
            if ext in (".conf", ".config", ".cfg"):
                config_files.append(arg)
            else:
                input_files.append(arg)
    if len(input_files) == 0:
        input_files = None

    # Warn user if there are non-existant config files
    if len(config_files) > 1:
        for cf in config_files[1:]:
            if not os.path.exists(cf):
                print("WARN :: Specified config file '{0}' does not exist".format(cf))
                #log.warn("Specified config file '{0}' does not exist".format(cf))

    try:
        log_level = int(options.log_level)
    except ValueError:
        try:
            log_level = getattr(logging, options.log_level)
        except:
            print 'Unrecognized log level:', options.log_level, '. Not setting.'
    return (input_files, config_files, log_level)

def libstring_to_func(lib_and_func):
    """ 
    Given a string containing the full module path to a function, including
    the function name, return the corresponding function (i.e. callback)
    @param lib_and_func The "fully qualified" function name.
    e.g. foo.bar.baz will return function baz from module foo.bar
    """
    lib_name = lib_and_func[0:lib_and_func.rindex(".")]
    func_name = lib_and_func[lib_and_func.rindex(".")+1:]
    lib = importlib.import_module(lib_name)
    func = getattr(lib, func_name)
    return func

def datestr_to_datetime(startDate):
    """
    Convert a string in `MM-DD-YYYY hh:mm' format to a datetime
    @param startDate The String
    @return the datetime object
    """
    try:
        mdy = startDate.split(" ")[0]
        hm = startDate.split(" ")[1]
        (month, day, year) = [int(x) for x in mdy.split("-")]
        (hour, minute) = [int(x) for x in hm.split(":")]
        #startDate = startDate.replace("'", "").replace('"', '')
        #os.environ['TZ'] = 'UTC'
        #tzset()
        #return time.mktime(time.strptime(startDate, '%m-%d-%Y %H:%M'))
        return dtime(year=year, month=month, day=day,
                                 hour=hour, minute=minute)
    except ValueError:
        print 'Given start date', startDate, 'does not match expected format MM-DD-YYYY hh:mm'
        sys.exit(1)


def _get_files_in_path(inspecs, topdir, startDate, duration, interval, firstFhr=0, domain=1):
    """
    Gets a list of model output files in a given path.
    @param inspecs The input specification files
    @param topdir The path containing the input files from the model/postproc
    @param startDate The first output date to look at. 
    @param duration Duration (as TimeDelta) of files to look at
    @param interval (as TimeDelta) between outputs
    @param firstFhr The first forecast hour to look at
    @param domain The domain number of the file
    """
    # TODO ? This function is not robus since different collections may be
    #        available at different times. It is not really necessary to
    #        get a  list of files since we're looping through dates anyway, but
    #        it does make it more consistent with the case of the user
    #        passing in th from the command line.
    #        The best solution is to guess the dates from the given input
    #        files and eliminate this and loop through dates in main()
    #       as aopposed to idx,infile as is done now
    file_list = []
    conf = ConfigParser()
    conf.optionxform = str
    conf.read(inspecs)
    try:
        fil_pattern = conf.get("BASIC", "file_name")
    except NoSectionError:
        fil_pattern = conf.items("field2filename_mappings")[0][1].format({"adate":startDate, "fhr":fhrInt,"dom":domain, "fdate":currDate})
        log.info("No 'file_name' specified in inputspec::[BASIC]. Will use first entry "
                 "in [field2filename_mappings] ==> {0}."
                 .format(fil_pattern))
    frequency = int(interval.total_seconds())
    duration = int(duration.total_seconds())
    dateRange = range(0, duration+1, frequency)
    all_dates = [startDate + tdelta(seconds=curr) for curr in dateRange]
    fhr = tdelta(hours=firstFhr)
    for currDate in all_dates:
        assert fhr.total_seconds() % 3600 == 0
        fhrInt = int(fhr.total_seconds() / 3600)
        conversion_args = {"fdate":currDate, "adate":startDate, "dom":domain, 
                           "fhr":fhrInt, "fmin":0, "cal_date":currDate}
        #file_path = os.path.join(path, fil_pattern.format(**conversion_args))
        file_path = os.path.join(topdir, fil_pattern)
        file_path = file_path.format(**conversion_args)
        if not os.path.exists(file_path):
            raise Exception("Could not find file '{0}'. Make sure configuration"
                            " settings for 'model', 'path', 'start_date', etc. "
                            " are correct."
                            .format(file_path))
        file_list.append(file_path)
        fhr += interval
    log.debug("Will process the following files: {a}".format(a=file_list))
    return file_list

def rearrange_infiles_by_date(infiles):
    # TODO : Get filepattern from inspec 
    fil_pattern = __inspec.get("BASIC", "file_name")
    #"gfs.t{adate:%H}z.pgrb2f{fhr:02d}"
    #gfs.t12z.pgrb2f00
    lastFhr = -1
    strlen = len(fil_pattern)
    pattern_sIdx = 0
    while pattern_sIdx < strlen and pattern_sIdx != -1:
        pattern_sIdx = fil_pattern.find("{", start=pattern_sIdx)
        pattern_eIdx = fil_pattern.find("}", start=pattern_sIdx)
        curr_attr = fil_pattern[pattern_sIdx+1:pattern_eIdx].split(":")
        if len(curr_attr) == 1:
           pass # TODO 
        elif len(curr_attr) == 2:
           pass # TODO
        else:
            raise Exception("Unexpected attribute in file pattern: {at}"
                            .format(at=curr_attr))
    return infiles, lastFhr

def ___get_extents_from_dataset(dataset, log=None):
    """
    Determine extents of grid from the given Nio file objec
    """
    #TODO : Use some of this for SpecData
    if log is None: log=_default_log()
    if input_file_source is G5NR_NPS_HIRES_NETCDF:
        return (-90., -180., 90., 180.)
    elif input_file_source is NMMB_HIST_NETCDF:
        try:
            attrs = ("CEN_LAT", "CEN_LON", "DX", "DY")
            (cen_lat, cen_lon, dx, dy) = (getattr(dataset, attr)[0] for attr in attrs)
        except:
            log.warn("No CEN_LAT, CEN_LON in dataset. Using (25,-80)")
            cen_lat = 25.
            cen_lon = -80.
            try:
                attrs = ("DX", "DY")
                (dx, dy) = (getattr(dataset, attr)[0] for attr in attrs)
            except:
                log.warn("No 'DX','DY' in Dataset. Using 'DPHD','DLMD'")
                attrs = ("DPHD", "DLMD")
                (dx, dy) = (getattr(dataset, attr)[0] for attr in attrs)

        log.debug("center lat,lon = {},{}".format(cen_lat, cen_lon))
        log.debug("DX,DY = {},{}".format(dx, dy))
        num_lats =  dataset.dimensions["south_north"]
        num_lons =  dataset.dimensions["west_east"]
        log.debug("num_lats,num_lons = {},{}".format(num_lats, num_lons))
        # get map extents
        west_lon = cen_lon - (num_lons*dx/2)
        east_lon = cen_lon + (num_lons*dx/2)
        south_lat = cen_lat - (num_lats*dy/2)
        north_lat = cen_lat + (num_lats*dy/2)
        return (south_lat, west_lon, north_lat, east_lon)
    else:
        raise Exception("Do not know what model was used. Can't continue")

def get_list(list_string, cast=None, minValues=-1, maxValues=-1):
    '''
    Convert space separated list (as string) to a list object, trimming the
    leading or trailing single or double quotes.
    If `cast` is not None: use method named `cast` to convert each element to
     this Type
    If minValues > -1,  raise an exception if we don't have at least minValues values
    If maxValues > -1, raise an exception if we have more than maxValues values
    '''
    list_string = list_string.replace("'", "").replace('"', '')
    l = list_string.strip().split(',')
    if cast is not None:
        l = [ cast(e) for e in l ]
    if minValues > -1:
        if len(l) < minValues: raise Exception('Value given [%s] should be at least %i values long' %(list_string, minValues) )
    if maxValues > -1:
        if len(l) > maxValues: raise Exception('Value given [%s] should be at most %i values long' %(list_string, maxValues) )
    return l

def _interp_filename(filename, pLev='', levIdx='', **conversions):
    """ 
    Interpolate ``filename'' by replacing <plev> with the (passed-in) pressure 
    level and <levIdx> with the passed-in level index
    """
    if pLev is None: pLev = ''
    if levIdx is None: levIdx = ''
    if pLev not in (None,''):
        pLev = str(int(pLev))
    if levIdx not in (None,''):
        levIdx = str(int(levIdx))
    filename = filename.replace("<lev>", pLev).replace("<levIdx>", levIdx)
    for k,v in conversions.iteritems():
        key = "{" + k + "}"
        filename = filename.replace(key, v)
    return filename

def _get_track_fdates_dict(conf, curr_fdate):
    """
    Determine TC vitals path for the storm and get the dictionary mapping
    forecast dates to pycane.postproc.tracker.objects.TrackerEntry objects
    :param conf: ConfigParser object containing the necessary items
    :returns: the dict returned by trkutils.fcst_date_dict()
    """
    if curr_fdate is None:
        raise Exception("curr_fdate must be passed if giving tcvitals_path")
    path = conf.get("BASIC", "tcvitals_path")
    args = dict(fdate=curr_fdate, stormId=conf.get("BASIC", "storm_id"))
    path = path.format(**args)
    trk = trkutils.get_track_data(path)
    return trk.fcst_date_dict
    
    
def _create_map(conf, plotSettings, dataset, log=None, curr_fdate=None):
    """
    Create and return a Basemap object set according to options specified
    in ConfigParser object `conf'. Use extents from Nio dataset `dataset' 
    as a backup.
    :conf: ConfigParser object containing map settings
    :plotSettings: Dict containing plot-specific config settings (field, tcvitals, etc)
    :curr_fdate: Current forecast date (init_date + fhr). This is only needed 
                 if doing storm-centric plots (to get the domain center from tc vitals
    :type curr_fdate: datetime.datetime
    """
    global _trk_fdates
    if log is None: log=_default_log()
    confmap = lambda param: conf.get("map_settings", param)
    extents_setting = plotSettings["map_extents_setting"]
    if extents_setting == "storm_centric":
        if _trk_fdates is None:
            _trk_fdates = _get_track_fdates_dict(conf, curr_fdate)
        try:
            width,length = [float(plotSettings[x]) for x in ("map_width","map_length")]
        except KeyError:
            raise Exception("If giving tcvitals_path to follow, you must specify"
                            " map_width and map_height as well".format())
        #import pdb ; pdb.set_trace()
        if curr_fdate in _trk_fdates:
            nearest = curr_fdate
        elif "find_nearest_fdate" in plotSettings and plotSettings["find_nearest_fdate"]:
            nearest = min(_trk_fdates.keys(), key=lambda d: abs(d - curr_fdate) )
            log.info("Current forecast date {0} is not in tc vitals. Using "
                     " nearest date: {1}".format(curr_fdate, nearest))
        else:
            raise Exception("Current forecast date {0} is not in tc vitals"
                            .format(curr_fdate))
        cenlat = _trk_fdates[nearest].lat
        cenlon = _trk_fdates[nearest].lon
        # pm_relative for basemap
        if cenlon > 180. or cenlon < -180.:
            cenlon += 180. ; cenlon %= 360. ; cenlon -= 180.
        # TODO : what if it's out of bounds (e.g. > 90N)
        south_lat = cenlat - length/2. ; north_lat = cenlat + length/2.
        west_lon = cenlon - width/2. ; east_lon = cenlon + width/2.
        assert abs(south_lat) < 90 ; assert abs(north_lat) < 90
        assert abs(west_lon) <= 180. ; assert abs(east_lon) <= 180.
    elif extents_setting == "config":
        extents = confmap("map_extents")
        extents = extents.replace("'", "").replace('"', '')
        try:
            extents = confmap('map_extents').replace("'","").replace('"','')
            toks = extents.strip().split()
            south_lat = float(toks[0])
            west_lon = float(toks[1])
            north_lat = float(toks[2])
            east_lon = float(toks[3])
        except:
            log.info("Error parsing 'map_extents' section, will use dataset")
    elif extents_setting == "data":
        (south_lat, west_lon, north_lat, east_lon) = \
                    _get_extents_from_dataset(dataset, log) 
    else:
        raise Exception("Bad map_extents_setting: {0}".format(extents_setting))

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
_trk_fdates = None
def _parse_contour_levels(contourStr):
    """
    Given a string with a range, convert it to a list of values. If 
    The range can be "min:max" or "min:max:step". In the former case, np.linspace
    will be used to calculate the range. In the latter case, np.arange will be
    used.
    """
    if contourStr == IGNORED_OPTION: return None
    r = [float(r) for r in contourStr.split(":")]
    if len(r) == 2:
        return np.linspace(r[0],r[1])
    elif len(r) == 3:
        return np.arange(r[0],r[1],r[2])
    else:
        raise Exception("Invalid range given: {0}".format(contourStr))

    
def _parse_level_idc(csvString):
    """
    Given a string containing comma-separated list of integers or ranges, return a 
    list of integers corresponding to the indices specified.
    Ranges use the standard format understood by Python's range()
    """
    el_list = []
    items = csvString.split(",")
    for item in items:
        try:
            item = int(item)
            el_list.append(item)
        except ValueError:
            toks = [int(x) for x in item.split(":")]
            if len(toks) == 2:
                interval = 1
            elif len(toks) == 3:
                interval = toks[2]
            else:
                raise Exception("Unsupported configuration value: {0}"
                                .format(item))
            el_list.extend(range(toks[0],toks[1], interval))
    return el_list



def _add_colorbars_to_plot(plot_settings, basemap, cs, plotsPerFigure, ax=None, fig=None):
    """
    Add colorbars to designated areas of figure/plot, depending on 
    configuration settings specified in <plot_settings>
    """
    #global subplot_counter
    if ax is None: ax  = plt.gca()
    if fig is None: fig = plt.gcf()
    
    if plot_settings is not None:
        #fig_colorbar_bottom = strtobool(plot_settings["fig_colorbar_bottom"])
        #fig_colorbar_right = strtobool(plot_settings["fig_colorbar_right"])
        plot_colorbar_bottom_freq = int(plot_settings["plot_colorbar_bottom_freq"])
        plot_colorbar_right_freq = int(plot_settings["plot_colorbar_right_freq"])
        cbar_on_top = strtobool(plot_settings["cbar_on_top"])
    else:
        #fig_colorbar_bottom = False
        #fig_colorbar_right = False
        plot_colorbar_bottom_freq = -9999
        plot_colorbar_right_freq = -9999
        cbar_on_top = False
    #import pdb ; pdb.set_trace()
    if subplot_counter % plot_colorbar_right_freq == 0:
        cb = basemap.colorbar(cs, "right", pad=0.20)
        if "right_colorbar_label" in plot_settings:
            s = plot_settings["right_colorbar_label"].format(**conversion_args)
            cb.set_label(s, rotation=270, labelpad=15, color=CBAR_LABEL_COLOR)
    if subplot_counter % plot_colorbar_bottom_freq == 0:
        cb = basemap.colorbar(cs, "bottom", pad=0.25)
        if "bottom_colorbar_label" in plot_settings:
            s = plot_settings["bottom_colorbar_label"].format(**conversion_args)
            cb.set_label(s, labelpad=1, color=CBAR_LABEL_COLOR)
    # import pdb ; pdb.set_trace()
    if subplot_counter == 1:
        # Calculate amount of padding to put on top
        #top = (1 - plotsPerFigure/8.)
        #top = 0.75
        top = 1
        inches_in_plot = plotsPerFigure * INCHES_PER_FIGURE
        inches_per_point = 0.01388889 # const
        padding = 0.032 # account for figure padding TODO: Retrieve this: currently guessing
        if plots_per_figure < 3:
            padding = 0.094    # works well with 1 plot/figure and x-axis text
        top -= padding
        if True: # Assume we are plotting a title
            title_height = FIGURE_TITLE_SIZE * inches_per_point
            title_height *= 2  # ASSUME : 2 lines title for now
            top -=  (title_height/inches_in_plot)
        if cbar_on_top:
            cbar_hgt = 0.5 # inch
            #padding = inches_in_plot * 0.01 # 0.022 # TODO: Retrieve this
            top -= (cbar_hgt / inches_in_plot) 
            # - (.004*inches_in_plot)
            #top = 6. / (plotsPerFigure*4.)
            log.debug("top = '{0}', title_height = '{1}'".format(top, title_height))
        #else:
        #    top -= (0.2 / inches_in_plot)
        fig.subplots_adjust(top=top) # leave .25 for title and cbar
        bot = top

        # TODO : finetune this offset/"margin"
        if plots_per_figure > 3:
            offset = 0.0025 # works well with 4 plots per page
        else:
            offset = 0.0325 # works well with 1 plot per page, no x-axis text
            if len(plot_settings["x_axis_text"]) > 0:
                #offset += 0.015
                offset += 0.01 # less => moves cbar down

        bot = top + (offset*inches_in_plot) 

        if cbar_on_top:
            #                     left  bot   wid   hgt %
            cb_ax = fig.add_axes([0.1, bot, 0.8, 0.001])
            cb_ax.set_axis_bgcolor("red")
            padding = 0.0 # This takes away from all corners!
            cb = basemap.colorbar(cs, "top", pad=padding, ax=cb_ax, size=0.30)
            if "top_colorbar_label" in plot_settings:
                cb.set_label(plot_settings["top_colorbar_label"], 
                             color=CBAR_LABEL_COLOR)

            cb_ax.get_xaxis().set_visible(False)
            cb_ax.get_yaxis().set_visible(False)

        plt.sca(ax) # further drawings on regular axes
    '''
    l, b, w, h = plt.gca().get_position().bounds
    ll, bb, ww, hh = cb.ax.get_position().bounds
    cb.ax.set_position([ll, b + 0.4*h, ww, h*1.2])
    '''

def _create_ffmpeg_animation(pattern, outfile, fps=5):
    """ Create mp4 movie out of the image files in the list passed in """
    cmd = [g_paths["ffmpeg"], "-r", str(fps), "-i", pattern, "-vcodec", "mpeg4", 
           "-y", outfile]
    subprocess.check_call(cmd)
    #ffmpeg  -i _tmp_mslp_filled.png%07d.png  -vcodec mpeg4 -y jza.mp4

def save_fig(specdata, filename, fhr, settings, animatedGif_file, mp4_file, fig=None):
    """
    Save a file to disk and close the plot. If filename extension is "PDF", 
    use PdfPages to save the figure, consequently saving multiple plots to the 
    same file. For other file types, if a file of the same name already exists,
    rename to <filename base>_i.<filename ext>. Where 'i' is the next available
    file name.
    This function also adds the title to the figure
    @param filename The name of the output file
    @param settings Dictionary containing settings specified in the "BASIC"
                    section of the main configuration object as well as the 
                    following dynamic parameters:
                        -"lev" - pressure-level
                        -"fhr" - forecast hour (float)
                        -"sd" - datetime of start/initialization date
    @param fig the Figure object. Default: use pyplot::gcf()
    """
    #import pdb ; pdb.set_trace()
    global __already_saved, __pdf, __animation_img, _animations, _mp4_frame_ctr
    if fig is None: fig = plt.gcf()
    if not "__already_saved" in globals():
        __already_saved = []
    if not "__pdf" in globals():
        __pdf = {}
    if not "__animation_img" in globals():
        __animation_img = {} # map `filename' (which does not have index suffix) to list of images
    settings.update(dict(specdata.inputspec.items("BASIC")))
    plot_title = settings["title"].format(**settings).decode("string-escape")
    fig.suptitle(plot_title, size=FIGURE_TITLE_SIZE, color=settings["fig_title_color"])
    base = os.path.splitext(filename)[0]
    ext = os.path.splitext(filename)[1]
   
    #plt.subplots_adjust(left=0.005, right=0.995)
    #plt.tight_layout()

    # TODO Adjust margins - need to account for titles and color bars
    # Note that this does not address the issue of thin/square domains having wide
    # L/R margins due to short height vs. wide title/cbar
    #plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    
    # Since we'll be purging old files, ensure they are image files
    assert ext.lower() in (".png", ".pdf", ".jpg", ".jpeg", ".gif", ".tiff", ".bmp")
    if not base in __already_saved:
        if os.path.exists(filename):
            log.info("Will delete existing images")
            os.unlink(filename)
        if filename.endswith("pdf"):
            log.debug("Creating PdfPages for '{fil}'".format(fil=filename))
            __pdf[filename] = PdfPages(filename)
            __already_saved.append(base)

    if filename.endswith("pdf"):
        log.debug("Will append to existing PDF")
        __pdf[filename].savefig(fig, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close()
    else:
        # not PDF
        if not base in __already_saved:
            __already_saved.append(base)
            #import pdb ; pdb.set_trace()
            files = glob.glob(base + "*" + ext)
            if len(files) > 0:
                log.info("Will delete existing output files: {ar}"
                         .format(ar=files))
                for fil in files: os.unlink(fil)
        if append_suffix_to_images:
            i = 1 # first is simply named <filename>
            newNamePattern = "{base}_{num:02d}{ext}"
            while True:
                newName = newNamePattern.format(base=base, ext=ext, num=i)
                if not os.path.exists(newName):
                    break
                assert i < 300 
                i += 1
        else:
            newName = filename
        plt.savefig(newName, facecolor=fig.get_facecolor(), edgecolor='none', 
                    bbox_inches='tight')

    if animatedGif_file :
        # For MPL.animation method:
        #if not filename in g_animation_img
            #_animations[filename] = ImageMagickFileWriter()
            #_animations[filename].setup(fig, filename+"_noctxmgr.gif", 80, 
            #              frame_prefix="_tmp_"+filename)
        #_animations[filename].fig = fig
        #_animations[filename].grab_frame()
        #plt.close()

        # If using a format supported by `convert' for gif,
        # no need to create a new figure. Otherwise, create PNG
        b,x = os.path.splitext(animatedGif_file)
        if append_suffix_to_images:
            anim_frame_img = b + "{0:05d}".format(i) 
        else:
            anim_frame_img = b + "{0:05d}".format(fhr)
        # NOTE: gif not supported for mp4 frames
        #import pdb ;pdb.set_trace()
        if ext[1:] in ['png', 'jpeg', 'ppm', 'tiff', 'sgi', 'bmp', 'gif']:
            anim_frame_img += ext
            try:
                os.symlink(newName, anim_frame_img)
            except OSError:
                log.debug("Removing existing link: " + anim_frame_img)
                os.unlink(anim_frame_img)
                os.symlink(newName, anim_frame_img)
        else:
            anim_frame_img = newName.replace(ext, ".png")
            plt.savefig(anim_frame_img, bbox_inches='tight')
        if anim_frame_img != newName:
            g_temp_img_files.append(anim_frame_img)
        if not animatedGif_file in g_animation_img:
            g_animation_img[animatedGif_file] = []
        g_animation_img[animatedGif_file].append(anim_frame_img)
    if mp4_file:
        b,x = os.path.splitext(mp4_file)
        anim_frame_img = b + "{0:05d}".format(_mp4_frame_ctr)
        _mp4_frame_ctr += 1
        # NOTE: duplicate code from here ...
        if ext[1:] in ['png', 'jpeg', 'ppm', 'tiff', 'sgi', 'bmp']:
            anim_frame_img += ext
            ext_frame = ext
            try:
                os.symlink(newName, anim_frame_img)
            except OSError:
                log.debug("Removing existing link: " + anim_frame_img)
                os.unlink(anim_frame_img)
                os.symlink(newName, anim_frame_img)
        else:
            anim_frame_img = newName.replace(ext, ".png")
            ext_frame = ".png"
            plt.savefig(anim_frame_img, bbox_inches='tight')
        if anim_frame_img != newName:
            g_temp_img_files.append(anim_frame_img)
            # ... to here 
        g_animation_img[mp4_file] = b + "%05d" + ext_frame # ffmpeg takes in pattern

    plt.close()

    #plt.tight_layout()
    #plt.savefig(filename, facecolor='blue', edgecolor='red')
#__already_saved = []
#__pdf = {}
# ffmpeg requires input files to be named sequentially, which ain't the case
# if using forecast offset as suffix
_mp4_frame_ctr = 0 
def _add_text_to_plot(plot_settings, basemap, fhr, fdate, fig=None, ax=None):
    """
    Add text to axes according to settings in <plot_settings>. The parameters
    in plot_settings will translate the following key words:
        -fhr = ``fhr'' (passed-in forecast hour)
        -fdate = The forecast date as a datetime
    TODO :
     - position of "left_column_text" if there is a colorbar
     - top-positioned y label will overlap with bottom-oriented colorbar label
    """
    if ax is None: ax  = plt.gca()
    if fig is None: fig = plt.gcf()
    log.debug("Adding text to plot")
    log.debug("plot_settings dump: {}".format(plot_settings))
    key_val = lambda d,k,default: d[k] if k in d.keys() else default
    setting_val = lambda k,default: key_val(plot_settings, k, default) 
    x_axis_text = setting_val("x_axis_text", None)
    y_axis_text = setting_val("y_axis_text", None)
    x_axis_position = setting_val("x_axis_text_location", "top")
    y_axis_position = setting_val("y_axis_text_location", "right")
    left_column_text = setting_val("left_column_text", None)
    #left_column_text  = "0 dark"
    left_column_text_offset = int(setting_val("left_column_text_offset", 40))
    left_column_text_size = int(setting_val("left_column_text_size", 22))  
    left_column_text_color = setting_val("left_column_text_color", "blue")
    # interpolate text
    parms = { "fhr":fhr, "fdate":fdate}
    if y_axis_text: y_axis_text = y_axis_text.format(**parms)
    if x_axis_text: x_axis_text = x_axis_text.format(**parms)
    if left_column_text: left_column_text = left_column_text.format(**parms)
            
    # Insert left-column text
    (loLat, hiLat) = ax.get_ylim()
    (leftLon, rightLon) = ax.get_xlim()
    if left_column_text:
        ax.text(leftLon-left_column_text_offset, loLat+(hiLat-loLat)/2, 
                left_column_text, fontsize=left_column_text_size, 
                color=left_column_text_color)
    # Insert axis labels
    if x_axis_text: 
        ax.set_xlabel(x_axis_text, color=left_column_text_color) # TODO : Use xlabel color form config
        ax.xaxis.set_label_position(x_axis_position)
    if y_axis_text:
        ax.set_ylabel(y_axis_text)
        ax.yaxis.set_label_position(y_axis_position)
      
def _get_levels_to_plot(plotSettings):
    """
    Determine either indices corresponding to the vertical level to
    plot or the pressure level to plot, depending on the values 
    set in `plotSettings'. Either "level_idc_to_plot" to plot OR 
    "pressure_levels_to_plot_hPa" must be set. If the former, those
    vertical indices will be ultimately be plotted. If the latter,
    those pressure levels, or those nearest to them, will ultimately be
    plotted, depnding on the "level_lookup_setting".
    "pressure_levels_to_plot_hPa should be CSV of pressure levels to plot, 
    in millibars/hPa.
    @return a 2-tuple: (level_idc_to_plot, plevs_to_plot); one of
             these will be None, depending on which option is set. 
    """
    level_idc_to_plot = None
    plevs_to_plot = None
    #import pdb ; pdb.set_trace()
    if plotSettings.has_key("level_idc_to_plot"):
        idc = plotSettings["level_idc_to_plot"]
        #assert isinstance(idc, list)
        level_idc_to_plot = _parse_level_idc(idc)
    if plotSettings.has_key("pressure_levels_to_plot_hPa"):
        if level_idc_to_plot is not None:
            raise Exception("Config error: either `level_idc_to_plot _or_ "
                            "pressure_levels_to_plot_hPa must be set".format())
        pLevs = plotSettings["pressure_levels_to_plot_hPa"].strip().split(",")
        plevs_to_plot = np.array([float(lev) for lev in pLevs])
    #assert (level_idc_to_plot,plevs_to_plot) != (None,None)
    # -> OK if plotting 2-d data only
    return (level_idc_to_plot, plevs_to_plot)


def _get_units(dataset, fieldName):
    """
    Get units of variable named "fieldName" from Nio dataset ``dataset''.
    First, try to get from the variable attributes. If that fails, use
    the modelspec
    """
    try:
        return getattr(dataset.variables[fieldName], "units")
    except AttributeError:
        #import pdb ; pdb.set_trace()
        unit_mappings = dict(__input_src_conf.items("units"))
        #field_mappings = dict(__input_src_conf.items("field_mappings"))
        if fieldName in unit_mappings.keys():
            return unit_mappings[fieldName]
        # see if there's a mapping to the standard name
        for k,v in _model2standard_fields.iteritems():
            if v == fieldName:
                if k in unit_mappings.keys():
                    return unit_mappings[k]
        raise Exception("Variable '{v}' does not have  a 'units' attribute "
                            " or an entry under [units] in modelspec file {fil}"
                            .format(v=v, fil=_get_modelspec_file()))

def _plot_wind_vectors(u2d, v2d, lats, lons, resolution, plotFunc, basemap,
                       plotSettings, ax=None, fig=None, log=None, dimOrder=None,
                       fieldUnits=None):
    """
    Plot wind barbs or arrows for a slice of flow field data
    :param u2d: Array of U-wind values
    :param v2d: Array of V-wind values
    :param resolution: Resolution in degrees of dataset    
    @param dataset Nio dataset containing wind fields
    @param basemap Basemap object to use for obtaining projection
    @param latVarName Name of the latitude variable in the dataset.
           May be used to determine the resolution of the data - 
           see get_resolution()
    @param plotType Type of plot to generate - either "arrows" or "barbs"
    @param level Vertical level to plot. Should be None to plot a 
            height-above-surface value or the level value to plot
            3-D winds at the given level.
    @param ax Axes to draw on. If None, use gca
    @param fig Figure object to use. If None, use gcf
    @param fieldUnits Units of given wind field values
    @param log Logger to use for logging. If None, use _default_log
    """
    if fig is None: fig = plt.gcf()
    if ax is None: ax = plt.gca()
    if log is None: log = log
    kwargs = {}
    (x, y) = basemap(lons, lats)
    
    # Do the quiver. NMM-B does not produce initial u10, so check for that
    if not np.all(u2d == 0):
        # Get data, according to arrow_data_algorithm
        alg = plotSettings["arrow_data_algorithm"].strip() 
        degreesPerArrow=float(plotSettings["degrees_per_arrow"])
        degs_lat_per_arrow = degs_lon_per_arrow = degreesPerArrow
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if "arrows_per_inch_lat" in plotSettings:
            log.debug("arrows_per_inch_lat will override degrees_per_arrow")
            arrows_per_inch_lat = float(plotSettings["arrows_per_inch_lat"])
            figsize_lat_degs = abs(basemap.llcrnrlat - basemap.urcrnrlat)
            figsize_lat_inches = bbox.height
            lat_degs_per_inch = figsize_lat_degs / figsize_lat_inches 
            degs_lat_per_arrow = lat_degs_per_inch / arrows_per_inch_lat 
            # TODO : Same for lon
        if "arrows_per_inch_lon" in plotSettings:
            log.debug("arrows_per_inch_lon will override degrees_per_arrow")
            arrows_per_inch_lon = float(plotSettings["arrows_per_inch_lon"])
            figsize_lon_degs  = abs(basemap.llcrnrlon - basemap.urcrnrlat)
            figsize_lon_inches = bbox.width
            lon_degs_per_inch = figsize_lon_degs / figsize_lon_inches 
            degs_lon_per_arrow = lon_degs_per_inch / arrows_per_inch_lon
        
        if alg == "skip":
            #import pdb ; pdb.set_trace()
            #skip = degreesPerArrow / resolution
            skipY = degs_lat_per_arrow / resolution 
            skipX = degs_lon_per_arrow / resolution 
            if dimOrder == "ij":
                u2d = u2d[::skipX,::skipY]
                v2d = v2d[::skipX,::skipY]
            elif dimOrder == "ji":
                u2d = u2d[::skipY,::skipX]
                v2d = v2d[::skipY,::skipX]
            else:
                raise Exception("unexpected dimension order for u2d/v2d")
            if x.ndim == 1:
                (x,y) = (x[::skipX], y[::skipY])
            elif x.ndim == 2:
                if dimOrder == "ij":
                    x = x[::skipX,::skipY]
                    y = y[::skipX,::skipY]
                else:
                    x = x[::skipY,::skipX]
                    y = y[::skipY,::skipX]
            else:
                raise Exception("Unexpected dimensionality for x/y")
        else:
            raise Exception("Unknown arrow_data_algorithm: '{0}'".format(alg))

        # Set color/cmap
        do_cmap = False
        if "arrow_color_map" in plotSettings:
            do_cmap = True
            kwargs["cmap"] = getattr(plt.cm, "jet")
            #kwargs["cmap"] = plt.cm.jet
            #kwargs["cmap"] = "jet"
        else:
            if plotFunc.__name__ == "quiver":
                key = "color"
            else:
                key = "barbcolor"
            try:
                kwargs[key] = plotSettings["arrow_color"]
            except KeyError:
                kwargs[key] = "black"

        # Plot
        if plotFunc.__name__ == "quiver":
            if do_cmap:
                args = (x, y, u2d, v2d, np.sqrt(u2d**2 + v2d**2))
            else:
                args = (x, y, u2d, v2d)
            kwargs.update({"scale":400, "linewidth":1.0})
            #basemap.quiver(x, y, u10, v10)
        elif plotFunc.__name__ == "barbs":
            #basemap.barbs(x, y, u10, v10, scale=400)
            #args = (x, y, u2d, v2d, length=5, barbcolor='black', flagcolor='r',
            #        linewidth=0.5)
            kwargs.update({"length":5, "flagcolor":'r', "linewidth":0.5})
            if do_cmap:
                # pass in speed
                args = (x, y, u2d, v2d, np.sqrt(u2d**2 + v2d**2))
                #kwargs["barb_increments"] = dict(half=10, full=20, flag=60)
            else:
                args = (x, y, u2d, v2d)
        else:
            raise Exception("Unknown wind field plot type")
        #import pdb ; pdb.set_trace()
        p = plotFunc(*args, **kwargs)
        # Create quiverkey or colorbar, if applicable
        if plotFunc.__name__ == "quiver":
            if "quiverkey_color" in plotSettings:
                color = plotSettings["quiverkey_color"]
                if color and not do_cmap:
                    #qk = plt.quiverkey(p, 0.8, 1.05, 10, r'$10 \frac{m}{s}$', 
                    #                     L->R  B->T
                    #qk = plt.quiverkey(p, 1.0,  1.05, 10, r'$10 \frac{m}{s}$', 
                    s = "10 " + fieldUnits
                    qk = plt.quiverkey(p, 1.0,  1.05, 10, s,
                                       labelpos='E', color=color, 
                                       labelcolor=color)
                elif do_cmap:
                    cb = basemap.colorbar(p, cmap=kwargs["cmap"])
                    s = "Wind speed " + fieldUnits
                    cb.set_label(s, rotation=270, labelpad=15, color=CBAR_LABEL_COLOR)
                    
        # Colorbar for barbs
        if plotFunc.__name__ == "barbs":
            if do_cmap:
                # NOTE: basemap.barbs returns a 2-tuple: MPL::barbs for SH and
                #       another for NH
                cb = basemap.colorbar(p[1], cmap=kwargs["cmap"])
                s = "Wind speed " + fieldUnits
                cb.set_label(s, rotation=270, labelpad=15, color=CBAR_LABEL_COLOR)


def _contour_plot_helper(lons, lats, var2d, basemap, plot_func, plotsPerFigure,
                         contour_levels=None, colors=None, fontsize=10, ax=None,
                         plot_settings=None, fig=None, map_settings=None,
                         basic_settings=None):
    """
    Do the actual plotting of the data, given array of lat, lons, and a 2-d slab
    of data. Also, draw and decorate the map according to ``map_settings'' and
    add text according to ``plot_settings''
    :param lons: longitude array
    :param lats: latitude array
    :param var2d: Either an array (if doing contour plot) or a list of 2-arrays 
                  if plotting barbs/arrows (element one is the u-wind array and 
                  eleent two is the v-wind array)
    :param basemap: Basemap object (ie canvas)
    :param plot_func: Function pointer to plotting function being used
    """
    if ax is None: ax  = plt.gca()
    if fig is None: fig = plt.gcf()
    kwargs = {}
    (x, y) = basemap(lons, lats)
    if "max" in contour_levels: 
        contour_levels = contour_levels.replace("max", str(var2d.max()))
    if "min" in contour_levels: 
        contour_levels = contour_levels.replace("min", str(var2d.min()))
    # Set colormap if specified. First see if it is a known cmap (string).
    # If not, try to import it
    if "colormap" in plot_settings:
        cmap = plot_settings["colormap"]
        if cmap in cmap_d.keys():
            kwargs["cmap"] = cmap
        else:
            # import function from lib as string: path.to.lib.field_name
            try:
                lib_and_func = cmap
                lib_name = lib_and_func[0:lib_and_func.rindex(".")]
                func_name = lib_and_func[lib_and_func.rindex(".")+1:]
                lib = importlib.import_module(lib_name)
                kwargs["cmap"] = getattr(lib, func_name)
            except Exception as e:
                log.error("Unable to import colormap '{0}'. If specifying a "
                          "module and function, it should be in the format "
                          "'path.to.lib.field_name' (e.g. matplotlib.cm.bone)"
                          .format(lib_and_func))
                raise e

            
    contour_levels = _parse_contour_levels(contour_levels) 
    cs = plot_func(x, y, var2d, contour_levels, colors=colors, **kwargs)
    
    if plot_func in (plt.contourf, basemap.contourf):
        _add_colorbars_to_plot(plot_settings, basemap, cs, plotsPerFigure, ax, fig)
    elif plot_func in (plt.contour, basemap.contour):
        plt.clabel(cs, inline=1, fmt="%3d", inline_spacing=1, fontsize=fontsize)

    # Final processing
    #ax.set_axis_bgcolor((1, 0, 0))

def _set_field_transforms(field, fieldCtr, plotSettings):
    """
	Set the field transforms for the given field based on plotSettings.
    plotSettings is a dictionary containing various settings; each setting
    may have settings for multiple fields, so only read the one at 
    index `fieldCtr'.
    """
    kv = [ ("convert_units", "units", str), ("multiply_by", "multiply_by", float),
           ("equation", "equation", str) ]
    for transName,settingName,typeCast in kv:
        try:    
            settingVal = plotSettings[settingName].split(",")[fieldCtr]
            # TODO? Have set_transform() for field since _specdata is private
            #field.set_transform(transName, typeCast(settingVal))
            field._specdata.set_field_transform(transName, typeCast(settingVal))
        except KeyError:
            log.debug("'{setting}' not specified in plotSettings. Not set."
                      .format(setting=settingName))

def _set_specFile_transforms(specFile, plotSettings):
    """
    Set default transforms that should be performed on all fields
    """
    # Set the coord transforms
    # For Basemap, grid data must go from -180->180 (GFS goes 0->360)             
    specFile.set_coord_transform("pm_relative_longitude", True) 
    # Need lat and lon to be 2D for contour plots to work in Basemap
    specFile.set_coord_transform("make_2d", True)
    # Assume we won't be using time dim
    #specFile.set_transform("squeeze_out_time", True)
        
def __get_plot_settings(plotSettings, fieldIdx):
    """ Helper function to read settings """
    # each of these may be a csv, one for each product within the figure
    plot_type = plotSettings["type"].split(",")[fieldIdx].strip()
    clevs_list = plotSettings["contour_levels"].split(",")[fieldIdx].strip()
    
    colors = None
    if plot_type == "line":
        colors = plotSettings["colors"].split(",")[fieldIdx]
        plot_func = basemap.contour
    elif plot_type == 'filled':
        plot_func = basemap.contourf
    elif plot_type == "arrows":
        plot_func = basemap.quiver
    elif plot_type == "barbs":
        plot_func = basemap.barbs
    else:
        raise Exception("Invalid plot_type: `{pt}'. Must be one of:"
                        "'line', 'filled', 'arrows', 'quiver'"
                        .format(pt=plot_type))
    return (plot_type, clevs_list, plot_func, colors)
    
def _get_var_data_sans_time(var, noTime=True):
    """ 
    Extract data from Nio/netCDF variables var. Remove the time dimension
    if present.
    """
    varVals = var[:].squeeze()
    if varVals.ndim != var[:].ndim:
        if var.dimensions[0].lower() != "time":
            raise Exception("'Time' dimension must be first")
    return varVals

def _get_plot_sets(plotSettings):
    """
    Given a dictionary with plot settings, return a "plot set" that describes 
    fields to be plotted.
    :param plotSettings: Dictionary with the following keys:
      - field_name - name of the field 
      - units - units to convert field's data to, if necessary
      - level_type - type of level ('isobaric', 'sfcDelta', None)
      - level_value(s) - value of the level. Only isobaric levels may have
                         multiple values as a CSV 
      - level_idc - Data indices. Either values or indices may be given. Not both.
    :type plotSettings: dict 
    :returns: list -- a 2-d array of 3-tuples. There is an array for each isobaric
    level to be plotted; the array consists of an array of ``field sets," which 
    is a 5-tuple (field_name, units, levelType, levelValue, levelIdx). 
    Either levelValue or levelIdx must be None. If both are None, it should be 
    a 2-d field, in which case levelValue should also be None
     e.g. If creating two figures, one with T@800 and PHIS and another
          with T@500 and PHIS:
         [ [(T, 'K', 'isobaric', 800, None), (PHIS, 'm', None, None, None)], 
           [(T, 'K', 'isobaric', 500, None), (PHIS, None, None, None)] ]     
    """
    field_names = [f.strip() for f in plotSettings["fields"].split(",")] 
    lev_types = [lt.strip() for lt in plotSettings["level_types"].split(",")]
    fld_units = [u.strip() for u in plotSettings["units"].split(",")]
    assert len(field_names) == len(lev_types) == len(fld_units)
    isobaric_levs = None
    if "isobaric_levels_hPa" in plotSettings.keys():
        isobaric_levs = plotSettings["isobaric_levels_hPa"].split(",")
    level_idc = None
    if "level_idc_to_plot" in plotSettings.keys():
        level_idc = plotSettings["level_idc_to_plot"].split(",")
    if level_idc is not None and isobaric_levs is not None:
        raise Exception("You can specify either 'level_idc_to_plot' OR " 
                        "'isobaric_levels_hPa'. Not both".format())
    multi_level_fieldsets = [] # isobaric or index
    single_level_fieldsets = [] # 2d, sfcDelta, etc.
    # Create fieldsets; each fieldset is a list of (5-element) lists. The fieldset 
    # encapsulates all fields to be plotted in a single image. The inner 
    # 5-element-list encapsulates the attributes of one of the fields.
    # placeholders for "isobaric" and "index" levels, which is why the inner
    # list is not a tuple
    for fieldName,levType,units in zip(field_names, lev_types, fld_units):
        if levType == "isobaric":
            assert level_idc is None
            if isobaric_levs is None:
                raise Exception("Field '{f}' specified as isobaric, but "
                                "'isobaric_levels_hPa' not specified"
                                .format(f=fieldName))
            multi_level_fieldsets.append([fieldName, units, levType, "placeholder", None])
        elif levType == "sfcDelta":
            delta = float(plotSettings["surface_delta_m"])
            single_level_fieldsets.append((fieldName, units, levType, delta, None))
        elif levType in (None,"-"):  #2d
            single_level_fieldsets.append((fieldName, units, None, None, None))
        elif levType == "index":
            assert isobaric_levs is None
            assert level_idc is not None
            multi_level_fieldsets.append([fieldName, units, levType, "placeholder", "placeholder"])
        else:
            raise ValueError("Unknown level type '{0}'".format(levType))
    # Now create the array of arrays - one for each level to be plotted. Set
    # the 'placeholder' values for isobaric/index levels
    if isobaric_levs is not None:
        ret = []
        #import pdb ; pdb.set_trace()
        for levCtr,levValue in enumerate(isobaric_levs):
            # level value is the third item in the fieldset tuple 
            isobaric_fieldsets = copy.deepcopy(multi_level_fieldsets)
            for fieldset in isobaric_fieldsets:
                assert fieldset[2] == "isobaric" # all multilevel should be isobaric at this point
                fieldset[3] = float(levValue)
                fieldset = tuple(fieldset) # see comment about having to make it a list for placeholder
            isobaric_fieldsets.extend(single_level_fieldsets)
            ret.append(isobaric_fieldsets)
        return ret
    elif level_idc is not None:
        ret = []
        for lev in level_idc:
            # level value,index are the 3rd,4th items in the fieldset tuple 
            threeD_fieldsets = copy.deepcopy(multi_level_fieldsets)
            for fieldset in threeD_fieldsets:
                fieldset[3] = float(lev)
                fieldset[4] = int(lev)
                fieldset = tuple(fieldset) # see comment about having to make it a list for placeholder
            threeD_fieldsets.extend(single_level_fieldsets)
            ret.append(threeD_fieldsets)
        return ret 
    else: # only single-level fieldsets
        return [single_level_fieldsets]
    
def make_plot(specdata, plotSettings, basemap, plotsPerFigure, log=None, 
              create_figure=False, basic_settings=None, fhr=None, fdate=None):
    """
    Create a plot/subplot for a particular field (or set of fields) at a 
    particular level. Actually, it can be multiple levels, but each level
    will be done in a separate figure.
    Technically speaking, we will only produce the request of one plot config
    section in each call to this function. So while there may be multiple 
    figures generated, they all pertain to the same set of fields.
    @param specdata SpecifiedForecastData object from which to retrive field
    @param basemap Basemap object to use as a canvas. 
    @param plotsPerFigure Numb of subplots to create on each Figure. The 
           size of the Figure will be a function of this number. One Figure
           is put on each output file.
    @param create_figure If True, create a new Figure to put the subplot in
    @param log If passed in, use it as the logger. Otherwise, use _default_log
    @param basic_settings Dictionary containing settings from the "[BASIC]" 
           section of the config file as well as some dynamic settings. See
           the documentation for save_fig() for more details.
    @param fhr Forecast hour (int)
    @param fdate Forecast date (datetime.datetime)
    @param plotSettings Dictionary containing various plot-related settings;
           most of these come from the plot configuration section of the 
           main config file(s). See the sample configuration file (test.conf)
           for details. Here is an abbreviated synopsis:
      Expected keys in `plotSettings':
        fields - List of field(s) to plot. Given name(s)
        type - the plot type, "line" for line contour; "filled" for filled contour
        levels - The contour level range, using the same format expected by 
                 numpy.arange 
        colors - For line contours, this is the color to use. Ignored for filled
                 contours. Only a single color is supported.
        colorbar_label - Label to put by the color bar. -9999=no label
        file_name - The name to used for the produced plot, Extension 
                    determines type.
        pressure_lookup_setting - The setting to use when looking for a pressure 
                                  level (when given pressure_levels_to_plot)
         Possible values:
           # 0 = No lookup, read the data file to determine (SLOW)
           # 1 = just trust value specified in the model's config file ;
           # 2 = Look for it in the first value of the output's PRESSSURE variable
           # 3 = In addition to (2), also ensure it's isobaric
           # This options overrides the setting in the model's configuratino

        GENERAL INFORMATION
            -The values for fields, standard_fields, type, and colors may be a comma-separated
             list (i.e. to plot multiple fields on the same canvas). The length
             of each list _must_ be the same. For fields for which an option does
             not apply (e.g. 'color' for a filled contour), the "ignore" value
             (-9999) should be used.
            -The "standard_name" is a name that applies to any model. There must
             be a mapping from it to the model's native field name
        
       [default_plot_settings]
       plot_colorbar_bottom_freq = 4
       plot_colorbar_right_freq = 4

       [mslpAndTemp]
       # Generate contour plot of MSLP
       fields = mslp, T
       contour_levels= 980:1044:4, 268:310:3
       colors = black, -9999
       file_name = 'mslp_and_temp.png'
       type = contour, contourf
       labels = yes, -9999
       colorbar_label = -9999, "Temperature"
       level_idc_to_plot = 0, 2:9:4
       pressure_levels_to_plot = 250, 500, 850
       plot_colorbar_bottom_freq = 8

    """
    global _figures, conversion_args
    if log is None: log = _default_log()
    # ensure uniform number of entries
    # Determine fields to plot (returns ConformedField elements)
    #fields = _get_fields_from_settings(specNio, plotSettings) 
    # TODO derived:  Need some additional abstraction here to deal with derived fields
    #  Maybe have a Field class
    out_filename = plotSettings["file_name"]
    slice_retrieval_setting = plotSettings["slab_retrieval_setting"]
    log.debug("Generating file: '{fil}'".format(fil=out_filename))
    
    if create_figure:
        _figures = [] # reset array for new set of figures
    
    #import pdb ; pdb.set_trace()
    plot_sets = _get_plot_sets(plotSettings) # , field_names, level_types)
    # e.g. [ [(T, 'K', 'isobaric', 800, None), (PHIS, 'm', None, None, None)],
    #        [(T, 'K', 'isobaric', 500, None), (PHIS, 'm', None, None, None)] ]

    # Iterate through plotsets
    for levCtr,plotSet in enumerate(plot_sets):
        # NOTE: Each plotset consists of a set of fields to be plotted on the
        # same figure. Ipso facto, just one level, but an arbitrary number
        # of fields are plotted.

        ### Create figure if necessary. The _figures array will be of length
        ### ``number-of-levels'' and indexed by i. Figure name will be <i>
        ##log.info("Plotting all fields at {levType} ~= {lev}"
                 ##.format(levType="pLev" if pressure_vals_given else "levIdx",
                         ##lev = lev_to_plot))
        if create_figure:
            log.debug("Creating Figure for plotset {0}".format(plotSet))
            _figures.append(plt.figure(levCtr))
            assert len(_figures) == levCtr + 1
            fig_length = INCHES_PER_FIGURE * plotsPerFigure
            _figures[levCtr].set_size_inches(FIGURE_WIDTH, fig_length)
            _figures[levCtr].patch.set_facecolor(plotSettings["fig_background_color"])
        else:
            plt.figure(levCtr)
        
        ax = plt.gcf().add_subplot(plotsPerFigure, 1, subplot_counter)

        # Iterate through fields in the plotset
        lev_units = { "isobaric":"hPa", "index":None, "sfcDelta":"m" } # TODO : others
        curr_filename = "None" # hack (see below)
        for fieldCtr,(fieldName, units, levType, levVal, levIdx) in enumerate(plotSet):
            field_settings = dict(levType=levType, levValue=levVal, units=units,
                                  levIdx=levIdx, fieldName=fieldName)
            log.info("Plotting field {fld}".format(fld=fieldName))
            (plot_type, clevs_list, plot_func, colors) = \
                     __get_plot_settings(plotSettings, fieldCtr)
            #import pdb ; pdb.set_trace()
            ## new (TODO? Can this go in a separate func?)
            if levType is None:
                levUnits=None
            else:
                levUnits = lev_units[levType]
            #import pdb ; pdb.set_trace()
            conversion_args = dict(levVal=levVal, fhr=fhr, fieldUnits=units)
            conversion_args.update(field_settings)
            if levVal is not None and int(levVal) == levVal:  
                # we can remove decimal from file name
                conversion_args["levVal"] = int(levVal)
                conversion_args["levValue"] = int(levVal)
            
            # Hack: If plotting multiple fields/levels on the same figure, the 
            # fieldName, levValue, units may change. Stick with the first 
            # filename unless there is a 'None' in it
            if  "None" in curr_filename:
                curr_filename = plot_settings["file_name"].format(**conversion_args)

            # Set varVals2d and (`fields'(vector plots) or `field' (contour plots))
            if plot_func.__name__ in ("barbs", "quiver"):
                fields = (specdata.get_field("x_wind", units, levType=levType),
                          specdata.get_field("y_wind", units, levType=levType))
                varVals2d = [None,None] ; coords = [None,None] ; lev = [None,None]
                #import pdb ; pdb.set_trace()
                for i,field in enumerate(fields):
                    _set_field_transforms(field, fieldCtr, plotSettings)
                    cslice = field.get_horz_slice(levVal, levUnits=levUnits, 
                                                  levIdx=levIdx, outDimOrder="ji", 
                                                   #lookupSetting=lev_lookup_setting,
                                                   sliceRetrievalSetting=slice_retrieval_setting,
                                                   fcstOffset=None)
                    varVals2d[i] = cslice.data
                    coords[i] = cslice.coords
                    lev[i] = cslice.speclev
                # ensure coords and levs are the same for U and V; then flatten
                assert len(varVals2d) == 2 and coords[0] == coords[1] and lev[0] == lev[1]
                coords = coords[0]
                lev = lev[0]
                (lats,lons) = coords.latlons
                # !! Field data may be transformed; only access thru varVals2d
                (dx,dy) = fields[0].resolution
                _plot_wind_vectors(varVals2d[0], varVals2d[1], lats, lons,
                                   dx, plot_func, basemap, plotSettings,
                                   ax=ax, fig=_figures[levCtr], log=log, 
                                   dimOrder=cslice.dim_order, fieldUnits=units)
                # TODO : The setting of padding/margins is currently done in _add_colorbars_to_plot,
                #        since they are co-dependant. That's why we call it here for other plot types
                _add_colorbars_to_plot(plot_settings, basemap, None, plotsPerFigure)

            elif plot_func.__name__ in ("contour", "contourf"):
                field = specdata.get_field(fieldName, units, levType=levType)
                _set_field_transforms(field, fieldCtr, plotSettings)
                cslice =  field.get_horz_slice(levVal, levUnits=levUnits, levIdx=levIdx,
                                         outDimOrder="ji", 
                                         #lookupSetting=lev_lookup_setting,
                                         sliceRetrievalSetting=slice_retrieval_setting, 
                                         fcstOffset=None)
                varVals2d = cslice.data
                (lats,lons) = cslice.coords.latlons
                # !! Field data may be transformed; only access thru varVals2d 
                _contour_plot_helper(
                        lons, lats, varVals2d, basemap,
                        plot_func, plotsPerFigure, 
                        contour_levels=clevs_list, colors=colors, 
                        fontsize=10, ax=ax, plot_settings=plot_settings, 
                        map_settings=dict(conf.items("map_settings")),
                        basic_settings=basic_settings, fig=_figures[levCtr])
            else:
                raise Exception("Unknown plot function")
            
            ### TODO derived:
            #args = {}
            ###if derived:
            ###    func = libstring_to_func(field)
            ###    varDims = ...
            ###    args = ...
            ###    varVals = # func(args)
            ### else:
            ##'''
            
            # TODO
                ##factor = plot_settings["multiply_by"].split(",")[fieldIdx].strip()
                ##if factor not in ("-", "1", "1.0", "", " "):
                    ##log.debug("Multiplying values by factor of {0}".format(factor))
                    ##varVals2d *= float(factor)

            #import pdb ; pdb.set_trace()
            fig = plt.gcf()
            # TODO ? Seems if a fill color for the continents is given, it 
            # does not paint the contour over it. Possible solutions:
            # decorate_map() before calling plot_func OR use z-levs
            map_settings=dict(conf.items("map_settings"))
            map_settings.update({"ax":ax})
            #import pdb ; pdb.set_trace()
            map_bling.decorate_map(basemap, **map_settings)
            
            _add_text_to_plot(plot_settings, basemap, fhr, fdate, fig, ax)
 

        # See if we want animated images
        if not "animated_gif_file" in plotSettings:
            animated_gif_file = None
        else:
            animated_gif_file = plotSettings["animated_gif_file"].format(**conversion_args)
        if not "mp4_file" in plotSettings:
            mp4_file = None
        else:
            mp4_file = plotSettings["mp4_file"].format(**conversion_args)

        # Save figure if it's time
        if subplot_counter == plots_per_figure \
           or fhr == final_forecast_hour \
           or ("current_file" in globals() and current_file == final_file_to_process)\
           or currDate == all_dates[-1]:
                settings = copy.copy(basic_settings)
                settings.update(conversion_args)
                settings.update(plot_settings)
                save_fig(specdata, curr_filename, fhr, settings, 
                         animated_gif_file, mp4_file)

# Global list of figures. There will be one figure for each 
# level and for each plot config (i.e. each section in the config file). Since
# we iteratre through forecast hours outside of this function, this list is
# kept as a global variable.
_figures = None

#   
# MAIN
#
if __name__ == "__main__":
    global subplot_counter, __pdf, final_forecast_hour
    global append_suffix_to_images, current_file 
    global g_animation_img, g_temp_img_files
    
    g_animation_img = {}
    g_temp_img_files = [] # will purge at end

    confbasic = lambda param: conf.get("BASIC", param)
    (infiles, config_files, log_level) = __parse_args()
    #import pdb ; pdb.set_trace()
    log = _default_log(log2stdout=log_level, name="omappy")
    conf = ConfigParser()
    conf.optionxform = str # retain capitalization
    conf.read(config_files)
    g_paths = dict(conf.items("paths"))

    # Get start date, duration, and frequency
    start_date = datestr_to_datetime(confbasic("start_date"))
    first_fhr = int(confbasic("first_fhr"))
    first_fhr_secs = first_fhr * 3600
    duration = tdelta(hours = float(confbasic("duration")))
    frequency = tdelta(hours = float(confbasic("interval")))
    #final_forecast_hour = duration.total_seconds() / 3600.
    periods = range(first_fhr_secs, int(duration.total_seconds()+1), 
                       int(frequency.total_seconds()))
    all_fhrs = [ x/3600. for x in periods]
    final_forecast_hour = all_fhrs[-1]
    #dates = [start_date + tdelta(seconds=curr) for curr in periods]
   
    # set initial conversion args
    conversion_args = {"init_date":start_date, 
                       "fdate":start_date + tdelta(hours=first_fhr), 
                       "fhr":first_fhr}

    # IF files were passed in, use them, otherwise determine from config
    if infiles is not None:
        log.debug("Input files to process: {0}".format(infiles))
        log.warn("Files may not be processed in chronological order, but "
                 " this program assumes that they are".format())
        # TODO ? Let user pass in inspec
        inspecs = specdata_utils.filepath2inpsec(infiles[0], INSPEC_PATH)
        # TODO : sort files by fhr
        #infiles, final_forecast_hour = rearrange_infiles_by_date(infiles) TODO
        final_file_to_process = infiles[-1]
        topdir = os.path.dirname(infiles[0])
    else: # files not passed in cmdline
        # TODO - clean this up and document which config parameters are needed in which cases
        #import pdb ; pdb.set_trace()
        log.debug("Input files not passed in, will use config setting")
        #input_file_source = globals()[confbasic("model")]
        topdir = confbasic("path").format(**conversion_args) 
        #import pdb ; pdb.set_trace()
        try:
            inspecs = [confbasic("inspec")]
        except ConfigParser.NoOptionError:
            log.info("No 'inspec' option specified in config. Will try guessing"
                     " which could lead to problems down the line")
            first_fil = confbasic("first_file").format(**conversion_args)
            first_fil_path = os.path.join(topdir, first_fil)
            # TODO ? fix this function to account for new directory scheme for inspecs
            inspecs = specdata_utils.filepath2inputspec(first_fil_path, INSPEC_PATH)
        #infiles = _get_files_in_path(inspecs, topdir, start_date, duration, frequency)
    
    
    desired_plots = confbasic("plots").split(", ")
    date_range = range(first_fhr_secs, int(duration.total_seconds())+1, 
                       int(frequency.total_seconds()))
    #import pdb ; pdb.set_trace()
    all_dates = [start_date + tdelta(seconds=curr) for curr in date_range]
    for plotConfSection in desired_plots:
        log.debug("Processing plot '{sec}'".format(sec=plotConfSection))
        basic_settings = {"sd":start_date} # need a datetime
        basic_settings.update(dict(conf.items("BASIC")))
        plot_settings = dict(conf.items("default_plot_settings"))
        plot_settings.update(dict(conf.items(plotConfSection)))
        plots_per_figure = int(plot_settings["plots_per_figure"])
        if "domain" in plot_settings.keys():
            domain = int(plot_settings["domain"])
        else:
            domain = 1
        basic_settings["title"] = plot_settings["title"]

        reset_counter = True

        fhr = tdelta(hours=first_fhr)
        for currDate in all_dates:    
            
            fhrInt = int(fhr.total_seconds() / 3600)
            # If there will be more than one file and the outfile does
            # not have forecast offset, append numerical suffix to names
            fname = plot_settings["file_name"]
            if "fhr" in fname or "fmin" in fname:
                append_suffix_to_images = False
            elif plots_per_figure < len(all_fhrs):
                append_suffix_to_images = True
            else:
                append_suffix_to_images = False

            if reset_counter:
                #plt.clf()
                #fig = plt.figure()
                #fig.set_size_inches(8.0, INCHES_PER_FIGURE*plots_per_figure)
                create_figure = True
                subplot_counter = 0
                reset_counter = False
            subplot_counter += 1

            # hack: inspecs is modified in SpecifiedForecastDatast (to include topdir)
            inspecs_copy = copy.copy(inspecs)
            specdata = SpecData(inspecs_copy, topdir, start_date, 
                                fcst_offset=fhr, domain=domain,
                                inspecsTopdir=INSPEC_PATH, #TODO get from config
                                log=log,
                                #fieldTransforms=
                                #coordTransforms=
                               )


            _set_specFile_transforms(specdata, plot_settings)
            log.info("Generating plot for {0}".format(currDate))
            #import pdb ; pdb.set_trace()
            # TODO / BUG? : conversion_args is getting overwritten after the first file
            #  as a result, it no longer has fdate. Is this expected?
            fdate = start_date + fhr 
            #basemap = _create_map(conf, plot_settings, dataset=None, curr_fdate=conversion_args["fdate"]) # TODO ? Can we create this outside?
            basemap = _create_map(conf, plot_settings, dataset=None, curr_fdate=fdate) # TODO ? Can we create this outside?
            # TODO : pass in dataset so that it uses the dataset extents if the extents are not specified in the config. This will require writing the function _get_extents_from_dataset()
            # Actually, we should leverage the Specdata for this since we need access to the HorizontalCoordinate for the specific Field - but it should be after performing transforms
            # -> gotta consider fact that domain can move

            #import pdb ; pdb.set_trace()
            if subplot_counter == plots_per_figure:
                reset_counter = True
            log.debug("subplot_ctr = {}".format(subplot_counter))
            # Note that the value of subplot counter stays the same, even though
            # the call to make_plot() will plot multiple levels for 3-d fields.
            # This is fine since there are separate Figures for each level
            make_plot(specdata, plot_settings, basemap, plots_per_figure, 
                      log=log, create_figure=create_figure, 
                      basic_settings=basic_settings, fhr=fhrInt, fdate=currDate)
            create_figure = False
            fhr += frequency

        # If there were any PDFs, close the PdfPages 
    # BUG : This doesn't work if I pass in the file(s) to process TODO
    #  -> actually, the underlying problem is that save_fig is not called
    # because last_fhr is not set
    # Also, we shouldn't even be going here nor creating __pdf if we aren't 
    # working with PDFs
    # UPDATE : Added the final_file_to_process workaround for now
    for k in __pdf.keys(): __pdf[k].close()
    #import pdb ; pdb.set_trace()
    for animFileName,v in g_animation_img.iteritems():
        if animFileName.endswith("gif"):
           imgList = v
           delay=str(int(plot_settings["gif_delay"]))
           cmd = ["convert", "-delay", delay, "-loop", "0"]
           cmd.extend(imgList)
           cmd.append(animFileName)  #TODO read from config (so key will have to be used in save_fig as well)
           subprocess.check_call(cmd)
        elif animFileName.endswith("mp4"):
            pattern = v
            _create_ffmpeg_animation(pattern, animFileName, 
                                     fps=int(plot_settings["mp4_fps"]))
        else:
            raise Exception("Unknown extension for animation: " + animFileName)

    # Purge any temp. files/links for animations
    for fil in g_temp_img_files:
        os.unlink(fil)

    # This is the MovieWriter version. Not going with this one since it 
    # may create images unnecessarily (since individual PDFs are usually created anyway)
    # * May make sense if we have optional filename AND/or different plots_per_figure for animations and "postage stamps"
    #if "_animations" in globals():
    #   for k in _animations.keys():
    #        _animations[k].finish()

