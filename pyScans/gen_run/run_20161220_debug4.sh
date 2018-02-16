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



date=20161220_debug56
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 5580. 30. 5760. 65. 0.1 1 128 &

date=20161220_debug57
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 11160. 30. 5760. 65. 0.1 1 128 &

date=20161220_debug58
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 16740. 30. 5760. 65. 0.1 1 128 &

date=20161220_debug59
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 22320. 30. 5760. 65. 0.1 1 128 &

date=20161220_debug60
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 1 128 &

date=20161220_debug61
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 2 128 &

date=20161220_debug62
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 7 128 &

date=20161220_debug63
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 30 128 &

date=20161220_debug64
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 60 128 &

date=20161220_debug65
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 180 128 &

date=20161220_debug66
python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 5760. 65. 0.1 365 128 &

exit
#++++++++
