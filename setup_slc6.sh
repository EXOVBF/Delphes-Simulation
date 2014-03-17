#!/bin/sh

source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh 
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh

export SCRAM_ARCH=slc6_amd64_gcc472
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt66/pythia8/183/x86_64-slc6-gcc46-opt/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc46-opt/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/cern.ch/user/l/lbrianza/work/SIMONE/delphes_code/:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt66/pythia8/183/x86_64-slc6-gcc46-opt/bin/:$PATH
export PYTHIA8=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt66/pythia8/183/x86_64-slc6-gcc46-opt
export PYTHIA8DATA=/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt66/pythia8/183/x86_64-slc6-gcc46-opt/xmldoc
