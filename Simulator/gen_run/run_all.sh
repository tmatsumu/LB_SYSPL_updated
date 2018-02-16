#!/bin/sh

#command=init
#command=mm
#command=coadd
command=clean_npy
#command=clean_npz
#command=clean_tod
#command=anafast
#command=plot_anafast

./run_gainval_v001.sh $command 0 &
./run_gainval_v002.sh $command 0 &
./run_gainval_v003.sh $command 0 &
./run_gainval_v004.sh $command 0 &
./run_gainval_v005.sh $command 0 &
./run_gainval_v006.sh $command 0 &

./run_gainval_v012.sh $command 0 &
./run_gainval_v013.sh $command 0 &
./run_gainval_v014.sh $command 0 &
./run_gainval_v015.sh $command 0 &
./run_gainval_v016.sh $command 0 &

./run_gainval_v022.sh $command 0 &
./run_gainval_v023.sh $command 0 &
./run_gainval_v024.sh $command 0 &
./run_gainval_v025.sh $command 0 &
./run_gainval_v026.sh $command 0 &
