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


date=20161220_debug34
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 5580. 30. 96. 65. 0.1 1 256 &

date=20161220_debug35
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 11160. 30. 96. 65. 0.1 1 256 &

date=20161220_debug36
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 16740. 30. 96. 65. 0.1 1 256 &

date=20161220_debug37
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 22320. 30. 96. 65. 0.1 1 256 &

date=20161220_debug38
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 1 256 &

date=20161220_debug39
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 2 256 &

date=20161220_debug40
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 7 256 &

date=20161220_debug41
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 30 256 &

date=20161220_debug42
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 60 256 &

date=20161220_debug43
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 180 256 &

date=20161220_debug44
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. 0.1 365 256 &

#++++++++

date=20161220_debug45
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 5580. 30. 96. 65. -0.1 1 256 &

date=20161220_debug46
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 11160. 30. 96. 65. -0.1 1 256 &

date=20161220_debug47
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 16740. 30. 96. 65. -0.1 1 256 &

date=20161220_debug48
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 22320. 30. 96. 65. -0.1 1 256 &

date=20161220_debug49
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 1 256 &

date=20161220_debug50
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 2 256 &

date=20161220_debug51
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 7 256 &

date=20161220_debug52
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 30 256 &

date=20161220_debug53
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 60 256 &

date=20161220_debug54
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 180 256 &

date=20161220_debug55
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 19 86400. 30. 96. 65. -0.1 365 256 &

exit
#++++++++

#date=20161220_debug01
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 5580. 30. 93. 65. 0.1 1 256 &

#date=20161220_debug02
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 11160. 30. 93. 65. 0.1 1 256 &

#date=20161220_debug03
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 16740. 30. 93. 65. 0.1 1 256 &

#date=20161220_debug04
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 22320. 30. 93. 65. 0.1 1 256 &

#date=20161220_debug05
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 1 256 &

#date=20161220_debug06
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 2 256 &

#date=20161220_debug07
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 7 256 &

#date=20161220_debug08
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 30 256 &

#date=20161220_debug09
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 60 256 &

#date=20161220_debug10
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 180 256 &

#date=20161220_debug11
#bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 93. 65. 0.1 365 256 &

#+++++++++++++++++++++

date=20161220_debug12
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 5580. 30. 96. 65. 0.1 1 256 &

date=20161220_debug13
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 11160. 30. 96. 65. 0.1 1 256 &

date=20161220_debug14
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 16740. 30. 96. 65. 0.1 1 256 &

date=20161220_debug15
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 22320. 30. 96. 65. 0.1 1 256 &

date=20161220_debug16
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 1 256 &

date=20161220_debug17
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 2 256 &

date=20161220_debug18
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 7 256 &

date=20161220_debug19
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 30 256 &

date=20161220_debug20
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 60 256 &

date=20161220_debug21
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 180 256 &

date=20161220_debug22
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 30. 96. 65. 0.1 365 256 &


#++++

date=20161220_debug23
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 5580. 30. 96. 65. 0.1 1 256 &

date=20161220_debug24
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 11160. 30. 96. 65. 0.1 1 256 &

date=20161220_debug25
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 16740. 30. 96. 65. 0.1 1 256 &

date=20161220_debug26
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 22320. 30. 96. 65. 0.1 1 256 &

date=20161220_debug27
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 1 256 &

date=20161220_debug28
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 2 256 &

date=20161220_debug29
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 7 256 &

date=20161220_debug30
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 30 256 &

date=20161220_debug31
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 60 256 &

date=20161220_debug32
bsub -q s -o ${date}_out.o python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 180 256 &

date=20161220_debug33
bsub -q s -o ${date}_out.o python  ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 11 86400. 30. 96. 65. 0.1 365 256 &
