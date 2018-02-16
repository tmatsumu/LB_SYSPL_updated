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



#date=20161220_debug45
#python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 5580. 30. 96. 65. -0.1 1 256 &

#date=20161220_debug46
#python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 11160. 30. 96. 65. -0.1 1 256 &

#date=20161220_debug47
#python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 16740. 30. 96. 65. -0.1 1 256 &

#date=20161220_debug48
#python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 22320. 30. 96. 65. -0.1 1 256 &

date=20161220_debug49
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 1 256 &

date=20161220_debug50
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 2 256 &

date=20161220_debug51
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 7 256 &

date=20161220_debug52
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 30 256 &

date=20161220_debug53
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 60 256 &

date=20161220_debug54
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 180 256 &

date=20161220_debug55
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 365 256 &

exit
#++++++++
