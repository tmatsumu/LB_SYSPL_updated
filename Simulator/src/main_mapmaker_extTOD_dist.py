import numpy as np
import healpy as h
import math as m
from numpy import fft
from numpy import random
import Gen_Mapmaking_inputs as gmi
import ctypes
import fileinput
import sys
import os
import re
import time
import glob
import ReadMapMakeXml as rxml
import string
from multiprocessing import Pool
from multiprocessing import Process
import copy
import lib_mapmaker as lmm
import time
import pickle
import gen_qsub_PsuedoCl as gqp
import lib_qsub as lqsub
import sqlite3 as sq
import lib_gain as lib_g

pi = np.pi
radeg=(180./pi)

#######################################################################################################
#######################################################################################################
print ""
print ""
print ""
print "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#"
print "#### LiteBIRD SIMULATOR                                                        ####"
print "####  written by T. Matsumura (2013-12-29)                                     ####"
print "####    revision: 2013-12-29, copied from main_mapmaker_dist.py                ####"
print "####                          take the external TOD and make maps              ####"
print "#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#"
print ""
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#------------------------------------------------------------
runtime_init = time.time()

xml_filename = sys.argv[1]
print "[main_mapmaker_dist.py]: START", xml_filename
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
#------------------------------------------------------------
print xml_input['file_input_maps2']

#------------------------------------------------------------
tod_label = sys.argv[3]
#------------------------------------------------------------
readDB = lmm.readDB()
readDB.xml=xml_input
readDB.sq_command = sys.argv[2]
CES = readDB.read_lb_observation()
readDB.display_all(CES)
nb_scanset = len(CES['id'])

readDB.sq_command = 'select * from detector'
readDB.filename = xml_input['file_fpdb_mmin']
fpdb_mmin = readDB.read_LBFPDB()
#------------------------------------------------------------


#------------------------------------------------------------
shuffle = lmm.shuffle_pixlist()
shuffle.boloid = fpdb_simin
pix_list = shuffle.gen_pixlist()
np.save(xml_input["dir_simedmap"]+'/processedPix',pix_list)

fpdb_input4mm = shuffle.gen_fpdb4mm(fpdb_mmin)
#------------------------------------------------------------


#------------------------------------------------------------
# CONSTRUCT THE SIMULATION INPUTS
inputs = {}
option_debug = True
muellermatrix = lmm.read_MuellerMatrix(xml_input['file_muellermatrix'],option_debug=option_debug)
if option_debug==True:
    num = len(fpdb_simin['detid'])
    muellermatrix['detid'] = fpdb_simin['detid']
    muellermatrix['detname'] = fpdb_simin['detname']
    muellermatrix['T1'] = np.ones(num)*muellermatrix['T1']
    muellermatrix['T2'] = np.ones(num)*muellermatrix['T2']
    muellermatrix['trans'] = np.ones(num)*muellermatrix['trans']
    muellermatrix['rho'] = np.ones(num)*muellermatrix['rho']
    muellermatrix['cosdel'] = np.ones(num)*muellermatrix['cosdel']
inputs["muellermatrix"] = shuffle.gen_muellermatrix(muellermatrix)

inputs["TQU"] = xml_input['TQU']
inputs["fpdb_list_mmin"] = fpdb_input4mm

inputs["NoiseInput"] = lmm.read_simnoise(xml_input['file_input_noise'],option_debug=True)
inputs["pix_list"] = pix_list 
if 'poly' in xml_input['choice']:
    inputs["poly"] = xml_input['poly']
inputs["nside"] = xml_input['nside']
inputs["runtime_init"] = runtime_init
inputs["silent"] = xml_input['silent']
inputs["run_type"] = xml_input['run_type']
inputs["pixelmapio"] = xml_input['pixelmapio']
#------------------------------------------------------------


#------------------------------------------------------------
# SET UP THE OUTPUT DIRECTORY
print ""
print "dir_simedmap:", xml_input["dir_simedmap"]
print "=============================================================="
print "the # of scansets: ", nb_scanset
print "the # of pixels: ", len(pix_list)
print "the total # of scansets (scanset x pix):", nb_scanset*len(pix_list)
print "=============================================================="
print ""

#------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#------------------------------------------------------------
# DISTRIBUTE THE SIMULTION INPUTS AS A BATCH JOB
# AND
# GENERATE THE QSUB SCRIPT FOR EACH CES/PIX
lmm.info('########### Parent Process ID ############', runtime_init)
subscan_total = 0
for i_scanset in range(0,nb_scanset):
    if 'E2G' in xml_input['coord']:
        if "sim" == xml_input['run_type']:
            print CES['dir_ptg'][i_scanset]
            ptg_idx = lmm.read_LBptg_idx(CES['dir_ptg'][i_scanset]+'.db')
            subscan_total += len(ptg_idx['s_idx'])
print ''
print xml_input['gain_type']

gain_cl = lib_g.gen_gain4mm()
gain_cl.sqlite_command = 'select * from GainParams;'
gain_cl.filename = xml_input['db_gain']
gain_cl.gain_type = xml_input['gain_type']
gain_cl.pix_list = pix_list
gain_in = gain_cl.prep_relgain4mm(subscan_total)

file_inputnpy_arr = []
i_subscan = 0
f_subscan = 0
for i_scanset in range(0,nb_scanset):
    print "i_scanset", i_scanset+1, '/', nb_scanset
    if 'E2G' in xml_input['coord']:
        if "sim" == xml_input['run_type']:
            print CES['dir_ptg'][i_scanset]
            ptg_package = lmm.read_LBptg(CES['dir_ptg'][i_scanset]+'.npz') 
            ptg_idx = lmm.read_LBptg_idx(CES['dir_ptg'][i_scanset]+'.db')
            f_subscan += len(ptg_idx['s_idx'])
    
    inputs["relgain"] = gain_in[:,i_subscan:f_subscan]
    inputs["gain_type"] = xml_input['gain_type']
    inputs["out_dir"] = xml_input["dir_simedmap"]+"/day"+str(i_scanset)+"/"
    os.popen("mkdir -p "+inputs["out_dir"])
    inputs["out_dir_ptg"] = xml_input["dir_simedmap"]
    inputs["ptg_package"] = ptg_package
    inputs["ptg_idx"] = ptg_idx
    inputs["id"] = CES["id"][i_scanset]

    file_inputnpy = xml_input["dir_simedmap"]+"/input_npy/inputs_"+str(i_scanset)
    np.save(file_inputnpy, inputs)
    file_inputnpy_arr.append(file_inputnpy)
    i_subscan += len(ptg_idx['s_idx'])

gen_qsub_mm = lqsub.gen_qsub()
gen_qsub_mm.file_inputnpy = file_inputnpy_arr
gen_qsub_mm.dir_simulator = xml_input['dir_simulator']
gen_qsub_mm.runID = xml_input['runID']
gen_qsub_mm.dir_out = xml_input['dir_simedmap']
gen_qsub_mm.mode = xml_input['debug']
gen_qsub_mm.machine = xml_input['machine']
out_qsub = gen_qsub_mm.gen_qsub_mm()

print '[main_mapmaker_dist.py]: Run mode', xml_input['debug']
for qsub_file in out_qsub:
    print qsub_file
    if xml_input['debug']!='skip_qsub':
        print ""
        print ""
        print "========================================================================"
        print "= Submit qsub "
        print qsub_file
        os.popen(qsub_file)
        print "========================================================================"
        print ""
        print ""

print "[MAIN] End for loop in main_mapmaker.py", '/', nb_scanset
sys.exit()

