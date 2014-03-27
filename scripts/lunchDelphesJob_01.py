#!/usr/bin/python
#**********************************************************************************************
#
#	this script lunch a job for each file contained in the directory specified with -d
#	
#**********************************************************************************************

import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import string


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def submitJob (jobID, fileName, queue, mqqCut, lheDir):
    jobname = 'jDel_'+queue+'_'+jobID+'.sh'
    fileName = fileName[:len(fileName)-1]
    rootName = fileName[:string.find(fileName,'.lhe')]+'.root'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('cmsStage -f /store/lhe/'+lheDir+'/'+fileName+ ' . \n\n')
    f.write ('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
    f.write ('cd /afs/cern.ch/user/s/spigazzi/work/EXOVBF/CMSSW_6_2_0/src/ \n')
    f.write ('eval `scramv1 runtime -sh` \n')
    f.write ('source /afs/cern.ch/user/s/spigazzi/work/EXOVBF/setup_slc6.sh \n')
    f.write ('cd - \n\n')
    f.write ('/afs/cern.ch/user/s/spigazzi/work/EXOVBF/Delphes-Simulation/DelphesPythia8 /afs/cern.ch/user/s/spigazzi/work/EXOVBF/Delphes-Simulation/cards/delphes_card_CMS_PileUp.tcl '+ fileName + ' ' + rootName + ' ' + mqqCut + '\n\n')
    f.write ('cmsStage -f ' + rootName + ' /store/user/govoni/Delphes/'+lheDir+'/')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    getstatusoutput ('bsub -q ' + queue + ' ' + '-u simone.pigazzini@cern.ch ' + jobname)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'boh')
    parser.add_argument ('-d', '--lheDir' , default = 'newDir' , help='sub directory of the lhe file (newDir)')
    parser.add_argument ('-q', '--queue' , default = '8nh', help='batch queue (1nd)')
    parser.add_argument ('-c', '--mqqCut' , default = '150', help='cut in mqq (150)')
    parser.add_argument ('-s', '--start' , default = 0, type=int, help='start file')
    parser.add_argument ('-e', '--end' , default = -1, type=int, help='end file')

    args = parser.parse_args ()

    # better crate the directory before running the jobs!
    #getstatusoutput ('cmsMkdir /store/user/govoni/Delphes/'+args.lheDir)
    getstatusoutput ('ls ../eosdir/cms/store/lhe/'+args.lheDir+'/ > lhe_file_'+args.queue+'.list')

    listf = open ('lhe_file_'+args.queue+'.list','r') 
    i = 0
    print 'submitting jobs to queue', args.queue
    for line in listf:
        if (i > args.end & args.end != -1) : 
            break        
        if (i >= args.start) :
            fileName = line
            submitJob (str (i), fileName, args.queue, args.mqqCut, args.lheDir) 
            i = i+1
        else :
            i = i+1

    listf.close()

