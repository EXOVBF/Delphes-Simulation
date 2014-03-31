DELPHES-SIMULATION:
====================

  * this repo contains modified version of Delphes.
  * The main difference is the output format: plain root tree of std::vector<float> instead of a root tree with the Delphes' classes, this means that the output is much lighter (~1kb/event). 
  * It also implement the calculation of nsubjettines for CA8 jets and pruned mass for all jets.


  * Installation steps:
```
  - cmsrel CMSSW_6_2_0
  - cd CMSSW_6_2_0
  - cmsenv
  - cd -
  - source setup_slc6.sh
  - make -j 16 -B (-B to rebuild everythings from the scratch)
```

  * Run Delphes with hepmc input and plain root output:
```
  - DelphesHepMC cards/delphes_card_CMS_PileUp.tcl input.hepmc output.root
```
  - you can also specify: 
        - start event number as 5° argument,
        - number of event to be processed 6°.
  
  * Run Delphes with lhe file input and root plain tree output:
```
  - DelphesPythia8 cards/delphes_card_CMS_PileUp.tcl input.hepmc output.root
```
  - you can also specify: 
        - Mqq_cut at lhe event level (VBF cut) as 4° argument,
        - start event number as 5°,
        - number of event to be processed 6°.
  - this program do all the stuff, shower+hadronization with Pythia8 and detector simulation with Delphes
  
  * the other readers are the standard Delphes codes.
  * to get back the standard heavy root file output with Delphes classes download Delphes from its site and replace modules/TreeWriter.cc and (.h) with the original one.
      
