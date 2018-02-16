#!/bin/sh

dir_bin=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/

# LiteBIRD sidelobe source 

# (option of shell input)
# (title)
# (sampleing [Hz])
# (sec in a day)
# (precession angle [degs])
# (precession period [min])
# (spin axis angle [degs])
# (spin axis period [rpm])
# (sim run period [days])
# (nside for projection)

date=20161223_debug100
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30.  90. 65. 0.1 365 128 &

date=20161223_debug101
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30.  92. 65. 0.1 365 128 &

date=20161223_debug102
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30.  94. 65. 0.1 365 128 &

date=20161223_debug103
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30.  96. 65. 0.1 365 128 &

date=20161223_debug104
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30.  98. 65. 0.1 365 128 &

date=20161223_debug105
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 30. 100. 65. 0.1 365 128 &

exit
#++++++++
