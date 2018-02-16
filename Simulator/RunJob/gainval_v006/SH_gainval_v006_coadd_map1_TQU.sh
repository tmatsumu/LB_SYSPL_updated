#!/bin/sh
python /home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator//src/fits_sep2oneTQU.py /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006/coadd_map/coadd_map1/mapT.fits /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006/coadd_map/coadd_map1/mapQ.fits /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006/coadd_map/coadd_map1/mapU.fits /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006/coadd_map/coadd_map1/mapTQU.fits
rm -f /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006/coadd_map/coadd_map1/*.npy
wait
