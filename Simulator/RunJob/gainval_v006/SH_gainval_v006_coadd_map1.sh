#!/bin/sh
python /home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator//src/run_coaddmaps.py /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/gainval_v006 coadd_map1 256 "select * from LBSimPtg where id < 366" /home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator//xml/xml_par_gainval_v006.xml & 
wait
