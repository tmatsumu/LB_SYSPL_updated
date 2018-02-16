#!/bin/sh
python /home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator//src/run_coaddmaps.py /group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3//SimedMaps/RunLog/02365_140GHz_alpha65beta30_prec93_p1rpm_10Hz_alpha0 coadd_map1 512 "select * from LBSimPtg where id < 366" /home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator//xml/xml_par_02365_140GHz_alpha65beta30_prec93_p1rpm_10Hz_alpha0.xml & 
wait
