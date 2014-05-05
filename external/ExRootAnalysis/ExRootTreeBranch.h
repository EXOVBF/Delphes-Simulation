#ifndef ExRootTreeBranch_h
#define ExRootTreeBranch_h

/** \class ExRootTreeBranch
 *
 *  Class handling object creation.
 *  It is also used for output ROOT tree branches
 *
 *  $Date: 2008-06-04 13:57:27 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "Rtypes.h"

#include <vector>

class TTree;
class TClonesArray;

class ExRootTreeBranch
{
public:

  ExRootTreeBranch(const char *name, TClass *cl, TTree *tree = 0);
  ExRootTreeBranch(const char *name, TTree *tree = 0); // mod  
  ~ExRootTreeBranch();

  TObject *NewEntry();
  std::vector<float>* NewFloatEntry(); // mod  
  void Clear();

private:

  Int_t fSize, fCapacity; //!
  std::vector<float>* fDataFloat; // mod  
  TClonesArray *fData; //!
};

#endif /* ExRootTreeBranch */

