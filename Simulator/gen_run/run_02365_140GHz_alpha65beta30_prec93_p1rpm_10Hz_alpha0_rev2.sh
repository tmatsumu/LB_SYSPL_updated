#!/bin/sh

waittime=1s

# 
run_name=02365_140GHz_alpha65beta30_prec93_p1rpm_10Hz_alpha0_rev2
#sleep $waittime | ./run_LB_SYSPL_v4.2.sh $1 $run_name 'select * from LBSimPtg where id > 180 and id < 366' 'select * from detector where detid < 370' $2 $3
sleep $waittime | ./run_LB_SYSPL_v4.3.sh $1 $run_name 'select * from LBSimPtg where id < 366' 'select * from detector where detid < 2' $2 $3

exit
