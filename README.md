DELPHES-SIMULATION:
====================

 * this repo contains modified version of Delphes.
 * The main difference is the output format: plain root tree of std::vector<float> instead of a root tree with the Delphes' classes, this means that the output is much lighter (~1kb/event). 
 * It also implement the calculation of nsubjettines for CA8 jets and pruned mass for all jets.
