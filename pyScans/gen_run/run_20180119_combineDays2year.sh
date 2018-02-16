#!/bin/sh

echo 
echo '========================='
echo 'BEGINNING OF SCRIPT'
echo ''

dir_pyScans=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.4/pyScans/

#sub_dir=LB_L2_20180118_P2G_case1_samplerate_SCANSPEC_1440min_45.0degs_90.6min_40.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case2_samplerate_SCANSPEC_1440min_30.0degs_90.6min_55.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case3_samplerate_SCANSPEC_1440min_5.0degs_90.6min_90.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case4_samplerate_SCANSPEC_1440min_7.0degs_90.6min_95.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case5_samplerate_SCANSPEC_1440min_45.0degs_90.6min_40.0degs_0.1rpm_365day_nside256_20Hz
#sub_dir=LB_L2_20180118_P2G_case6_samplerate_SCANSPEC_1440min_30.0degs_90.6min_65.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case7_samplerate_SCANSPEC_1440min_45.0degs_90.6min_50.0degs_0.1rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20180118_P2G_case8_samplerate_SCANSPEC_1440min_45.0degs_90.6min_50.0degs_0.1rpm_365day_nside256_20Hz/
sub_dir=LB_L2_20180118_P2G_case9_samplerate_SCANSPEC_1440min_65.0degs_90.6min_30.0degs_0.1rpm_365day_nside256_10Hz/

filename_db=$dir_pyScans/dataout/$sub_dir/ptg/LBPTG_latlon_eclip

filename_out=$dir_pyScans/dataout/$sub_dir/ptg/sample_boresight_365days
bsub -q s -o run_main_combined_day2all.o python $dir_pyScans/src/main_combined_day2all.py $filename_db 'select * from LBSimPtg where id < 366' $filename_out &

ls ${filename_out}*

echo 
echo 'END OF SCRIPT'
echo '========================='
echo ''

exit