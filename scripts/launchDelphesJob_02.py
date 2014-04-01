#!/usr/bin/python


import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def submitJob (jobID, queue, events, mqqCut, fileName):
    startEvent = (int(jobID)*int(events))
    jobname = 'jDel_'+queue+'_'+jobID+'.sh'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('cmsStage -f /store/user/rgerosa/DelphesAnalysis/kk_graviton_corrected/'+fileName+'.lhe . \n\n')
    f.write ('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
    f.write ('cd /afs/cern.ch/user/s/spigazzi/work/EXOVBF/CMSSW_6_2_0/src/ \n')
    f.write ('eval `scramv1 runtime -sh` \n')
    f.write ('source /afs/cern.ch/user/s/spigazzi/work/EXOVBF/setup_slc6.sh \n')
    f.write ('cd - \n\n')
    f.write ('/afs/cern.ch/user/s/spigazzi/work/EXOVBF/Delphes-Simulation/DelphesPythia8 /afs/cern.ch/user/s/spigazzi/work/EXOVBF/Delphes-Simulation/cards/delphes_card_CMS_PileUp_signal.tcl '+ fileName+'.lhe' + ' ' + fileName+'_'+jobID+'.root' + ' ' + mqqCut + ' ' + str(startEvent) + ' ' + events + '\n\n')
    f.write ('cmsStage -f ' + fileName+'_'+jobID+'.root ' + '/store/user/govoni/Delphes/signal_lvj/')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    getstatusoutput ('bsub -q ' + queue + ' ' + '-u simone.pigazzini@cern.ch ' + jobname)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'boh')
    parser.add_argument('-j', '--jobsNum' , default = 1, help='number of sub job')
    parser.add_argument('-f', '--fileName' , default = 'nofile', help='lhe file to be processed without .lhe')
    parser.add_argument('-q', '--queue' , default = '8nh', help='batch queue (1nd)')
    parser.add_argument('-c', '--mqqCut' , default = '150', help='cut in mqq (150)')
    parser.add_argument('-n', '--events' , default = '10000', help='total number of events per job(10000)')
    parser.add_argument('-o', '--offset' , default = 0, type=int, help='job numbering offset')
    
    args = parser.parse_args ()

    print 'submitting', args.jobsNum, 'jobs to queue', args.queue
    for i in range (0, int(args.jobsNum)):
        submitJob (str (i+args.offset), args.queue, args.events, args.mqqCut, args.fileName)
        
