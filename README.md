Delphes-Simulation
==================

* mods for boosted WW analysis:

```
- modules:  -TreeWriter
            -Isolation
            -EnergySmearing
            -MomentumSmearing
            -FastJetFinder
- classes:  -DelphesClasses
- external: -Fastejet ---> added fastjet/contrib
            -ExRootAnalysis ---> modified ExRootBranch & ExRootEntry 
```

* Compile on slc6

```
- source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh 
- cmsrel CMSSW_6_2_0 (or other) + cmsenv
- source setup_slc6.sh
- ./configure
- make -j 16
```
* Execute on slc6
 
- lhe input with plain root tree output (std::vector<float>):

```
- ./DelphesPythia8 cards/delphes_card_CMS_PileUp.tcl input_file.lhe output_file.root startEvent endEvent
```

- HepMC input with heavy delphes tree output:

```
- ./DelphesHepMC cards/delphes_card_CMS_PileUp.tcl output_file.root input_file.hepmc 
```
