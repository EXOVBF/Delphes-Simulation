DELPHES-SIMULATION:
====================

  * this repo contains modified version of Delphes.
  * The main difference is the output format: plain root tree of std::vector<float> instead of a root tree with the Delphes' classes, this means that the output is much lighter (~1kb/event). 
  * It also implement the calculation of nsubjettines for CA8 jets and pruned mass for all jets.


  * Installation steps:
```
  - source setup_slc6.sh
  - make -j -B (-B to rebuild everythings from the scratch)

```

  * Run Delphes with lhe file input and root plain tree output:
```
  - DelphesPythia8 cards/delphes_card_CMS_PileUp.tcl input.hepmc output.root

  - you can also specify: 
        - Mqq_cut at lhe event level (VBF cut) as 4° argument,
        - start event number as 5°,
        - number of event to be processed 6°.
  - this program do all the stuff, shower+hadronization with Pythia8 and detector simulation with Delphes
```
  
  * the other readers are the standard Delphes codes (they may not work with the mods).
  * to get back the standard heavy root file output with Delphes classes download Delphes from its site and replace modules/TreeWriter.cc and (.h) with the original one.
      
