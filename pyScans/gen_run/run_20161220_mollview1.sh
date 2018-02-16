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

file[0]=LB_L2_20161220_debug34_samplerate_SCANSPEC_93min_30.0degs_96.0min_65.0degs_0.1rpm_1day_nside256_19Hz
file[1]=LB_L2_20161220_debug35_samplerate_SCANSPEC_186min_30.0degs_96.0min_65.0degs_0.1rpm_1day_nside256_19Hz
file[2]=LB_L2_20161220_debug36_samplerate_SCANSPEC_279min_30.0degs_96.0min_65.0degs_0.1rpm_1day_nside256_19Hz
file[3]=LB_L2_20161220_debug37_samplerate_SCANSPEC_372min_30.0degs_96.0min_65.0degs_0.1rpm_1day_nside256_19Hz
file[4]=LB_L2_20161220_debug38_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_1day_nside256_19Hz
file[5]=LB_L2_20161220_debug39_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_2day_nside256_19Hz
file[6]=LB_L2_20161220_debug40_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_7day_nside256_19Hz
file[7]=LB_L2_20161220_debug41_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_30day_nside256_19Hz
file[8]=LB_L2_20161220_debug42_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_60day_nside256_19Hz
file[9]=LB_L2_20161220_debug43_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_180day_nside256_19Hz
file[10]=LB_L2_20161220_debug44_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_0.1rpm_365day_nside256_19Hz
file[11]=LB_L2_20161220_debug45_samplerate_SCANSPEC_93min_30.0degs_96.0min_65.0degs_-0.1rpm_1day_nside256_19Hz
file[12]=LB_L2_20161220_debug46_samplerate_SCANSPEC_186min_30.0degs_96.0min_65.0degs_-0.1rpm_1day_nside256_19Hz
file[13]=LB_L2_20161220_debug47_samplerate_SCANSPEC_279min_30.0degs_96.0min_65.0degs_-0.1rpm_1day_nside256_19Hz
file[14]=LB_L2_20161220_debug48_samplerate_SCANSPEC_372min_30.0degs_96.0min_65.0degs_-0.1rpm_1day_nside256_19Hz
file[15]=LB_L2_20161220_debug49_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_1day_nside256_19Hz
file[16]=LB_L2_20161220_debug50_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_2day_nside256_19Hz
file[17]=LB_L2_20161220_debug51_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_7day_nside256_19Hz
file[18]=LB_L2_20161220_debug52_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_30day_nside256_19Hz
file[19]=LB_L2_20161220_debug53_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_60day_nside256_19Hz
file[20]=LB_L2_20161220_debug54_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_180day_nside256_19Hz
file[21]=LB_L2_20161220_debug55_samplerate_SCANSPEC_1440min_30.0degs_96.0min_65.0degs_-0.1rpm_365day_nside256_19Hz

for ((i=1;i<2;i++)) ; do
    python ${dir_bin}/src/main_plotMollviewMap_from_fits.py $dir_in/${file[i]}/ps/ ${file[i]}_nhits 0 0
done

