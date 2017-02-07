#!/usr/bin/env python
"""
Tile images vertically or horizontally. 
Output is controlled by the parameters under "Set stuff"
Images are expected to have
the same names (<img_type>+<img_suffix>)  and be placed under 
different directories (<dataset_dirs>).
If the make_pdf parameter is True, combine individual files
into a single PDF after processing all files.
This script should work on most systems and does not depend on any omapy stuff.
However, for movies to work, you need ffmpeg installed with the necessary codecs.
REQUIREMENTS
 - The 'convert' command should be in user PATH
 - For make_pdf, pdftk should be installed
 - FOr make_movie: ffmpeg (set path below) with needed codecs

"""

import os
import subprocess
import glob

#
# CONSTANTS
#

# Set path to pdftk and corresponding GCJ path. 
# These are the paths on Theia.
# On Storm, pdftk is installed natively.
# These are only needed if make_pdf is True
PDFTK = "/home/Javier.Delgado/bin/pdftk"
GCJ_PATH = "/home/Javier.Delgado/libs_for_apps/pdftk/lib64"
ffmpeg_path = "/home/Javier.Delgado/scratch/apps_tmp/ffmpeg-git-20161211-64bit-static/ffmpeg"

##
# SETTINGS
##
append="+append" # use -append for vertical tiling
# Extention of input data to be tiled
inputs_extension = "png"
# Directories containing images
dataset_dirs = ["bsnr", "g5nr"]
# 3-d levels to plot (will be used to interpolate `img_types'
standard_levels = [250,500, 850]

# Pattern of image file names
'''
img_types = ["spechumd_{lev}".format(lev=lev) for lev in standard_levels]
img_types = ["mslp"] # _filled"]
img_types = ["air_temp_{lev}".format(lev=lev) for lev in standard_levels]
'''
img_types = ["prate"]

# Numbers on file names of images to process 
img_numbers = range(0, 240, 3)
#img_numbers = range(0, 24, 3)
# Make a single PDF file with one page per figure
make_pdf = False
# Make movie(s) of tiled images
make_movie = True
# Make animated gif of tiled images
make_gif = True
# Delay paramter to use when generating the animated gif
# (Tip: 20 for hourly, 40 for 3 hourly, xx for 6 hourly)
gif_delay = 40
# Suffix of image file name
#img_suffix = ".{ext}" # If just one image in set
#img_suffix = "_{idx:04d}.{ext}" # no 'fhr' in filename
img_suffix = "_f{idx:04d}.{ext}" # yes 'fhr' in filename
# List of movies to create,specifically what codec to use and what extension to use for each
# the file name (before extention) will be `img_type`
movie_sets = [ ("mpeg2video", "mpg"), ("mpeg4", "mp4") ]

#
# MAIN
#
if __name__ == "__main__":
   
    # Tile images
    for i in img_numbers:
        #import pdb ; pdb.set_trace()
        for img_type in img_types:
            file_name = img_type+img_suffix.format(idx=i, ext=inputs_extension) 
            in_paths = [os.path.join(ds, file_name) for ds in dataset_dirs]
            tiled_outfile = img_type + img_suffix.format(idx=i, ext=inputs_extension)
            
            # Tile using convert
            cmd_args = ["convert"]  
            cmd_args.extend(in_paths) 
            cmd_args.extend([append, tiled_outfile])
            subprocess.check_call(cmd_args)

            if make_pdf:
                if inputs_extension == "pdf": continue
                args = ['convert', tiled_outfile, tiled_outfile+".pdf"]
                subprocess.check_call(args)
    # Make PDF "presentation"
    if make_pdf:
        os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":" + GCJ_PATH
        for img_type in img_types:
            outfile = "comparo_"+img_type+".pdf"
            cmd_args = [PDFTK]
            cmd_args.extend(sorted(glob.glob(img_type+"*.pdf")))
            cmd_args.extend(["cat", "output", outfile])
            subprocess.check_call(cmd_args)
    # Make movie
    if make_movie or make_gif:
        for img_type in img_types:
            imgList = []
            fps = 2 # for 6 hourly
            fps = 5 # for hourly
            fps = 3 # 3-hourly
            #b + "%05d" + ext_frame 
            # ffmpeg requires files to be sequential
            if img_numbers[1] - img_numbers[0] != 1:
                for i,imgNumber in enumerate(img_numbers):
                    tile = img_type + img_suffix.format(idx=imgNumber, ext=inputs_extension)
                    link = "tmp_ffmpeg_{0:04d}.png".format(i)
                    os.symlink(tile, link)
                    imgList.append(link)
                pattern = "tmp_ffmpeg_%04d.png"
            else:
                #pattern = img_type +  "_%04d." + inputs_extension # should correspond to `tiled_outfile'
                pattern = img_type +  "_f%04d." + inputs_extension # should correspond to `tiled_outfile'
                for i,imgNumber in enumerate(img_numbers):
                    imgList.append(img_type + img_suffix.format(idx=imgNumber, ext=inputs_extension))
            if make_movie:
                for (codec,ext) in movie_sets:
                    outfile = img_type + "." + ext
                    cmd = [ffmpeg_path, "-r", str(fps), "-i", pattern, "-vcodec", codec,
                           "-y", outfile]
                    subprocess.check_call(cmd)
            if make_gif:
                #import pdb ; pdb.set_trace()
                cmd = ["convert", "-delay", str(gif_delay), "-loop", "0"]
                cmd.extend(imgList)
                cmd.append(img_type + ".gif") 
                subprocess.check_call(cmd)
            if img_numbers[1] - img_numbers[0] != 1:
                for ii in range(i+1):
                    os.unlink("tmp_ffmpeg_{0:04d}.png".format(ii))
