import os
from datetime import datetime as dtime
from datetime import timedelta as tdelta
import numpy as np

#from PyNIO import Nio
#import objects as specdata_objects
from objects import SpecifiedForecastDataset
#import utils as specdata_utils
import ESMF

# enable logging
esmpy = ESMF.Manager(debug=True)

INSPEC_TOPDIR = "./inspec"
if __name__ == "__main__":

    field = "isobaric.500.hPa.air_temperature.K"
    src_inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    src_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_now/"

    #dest_inspec = ["grib/upp/upp_nmb_grib2_multifile.conf"]
    #dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_5000x2500/postprd_now"
    dest_inspec = ["grib/upp/upp_nmb_grib1.conf"]
    #dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/2j_800x800/postprd_grib1"
    dest_topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/_new_beta_pl/800x800/postprd.orig"
    
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

    #import pdb ; pdb.set_trace()
    regrid = ESMF.Regrid(src_esmf_field, dest_esmf_field, 
                         src_mask_values=src_field.missing_value, 
                         dst_mask_values=dest_field.missing_value,
                         unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE)
    #jzanote : need the unmapped_action setting or it will fail due to their 
    # being dest. points that do not map to source points
    dest_esmf_field = regrid(src_esmf_field, dest_esmf_field)
    
    #diff = src_esmf_field.data[:] - dest_esmf_field.data[:]
    #print diff.min(), diff.max(), diff.mean()

    src_data = np.ma.masked_equal(src_esmf_field.data[:], src_field.missing_value)
    dest_data = np.ma.masked_equal(dest_esmf_field.data[:], dest_field.missing_value)
    diff = src_data - dest_data
    print diff.shape, diff.min(), diff.max(), diff.mean()

    #src_data.inputspec

"""

** think we gotta mask missing points before subtracting

(Pdb) diff = src_esmf_field.data[:] - dest_esmf_field.data[:]
(Pdb) diff.min()
-1.0000000200408773e+20
(Pdb) diff.max()
1.0000000200408773e+20
(Pdb) diff.mean()
-8.9218751788022029e+17
(Pdb) diff2 = diff[diff < 1e19]
(Pdb) diff2 = diff2[diff2 > -1e19]
(Pdb) diff2.shape
(602506,)
# was the same after this
(Pdb) diff2 = diff2[diff2 > 0.0001]
(Pdb) diff2.shape
(197799,)
# much smaller after his
(Pdb) diff.shape
(800, 800)
(Pdb) 800*800
640000
(Pdb) diff2.min()
0.009063720703125
(Pdb) diff2.max()
1.0000000200408773e+20
"""
