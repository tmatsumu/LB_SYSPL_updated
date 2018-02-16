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

dir_in=/group/cmb/litebird/simdata/Scans/dataout/

file[0]=LB_L2_20161223_debug100_samplerate_SCANSPEC_1440min_30.0degs_90.0min_65.0degs_0.1rpm_365day_nside128_20Hz
file[1]=LB_L2_20161223_debug101_samplerate_SCANSPEC_1440min_30.0degs_92.0min_65.0degs_0.1rpm_365day_nside128_20Hz
file[2]=LB_L2_20161223_debug102_samplerate_SCANSPEC_1440min_30.0degs_94.0min_65.0degs_0.1rpm_365day_nside128_20Hz
file[3]=LB_L2_20161223_debug103_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_365day_nside128_20Hz
file[4]=LB_L2_20161223_debug104_samplerate_SCANSPEC_1440min_30.0degs_98.0min_65.0degs_0.1rpm_365day_nside128_20Hz
file[5]=LB_L2_20161223_debug105_samplerate_SCANSPEC_1440min_30.0degs_100.0min_65.0degs_0.1rpm_365day_nside128_20Hz

for ((i=1;i<2;i++)) ; do
    python ${dir_bin}/src/main_plotMollviewMap_from_fits.py $dir_in/${file[i]}/ps/ ${file[i]}_nhits 0 0
done

