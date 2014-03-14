#ifndef TreeWriter_h
#define TreeWriter_h

/** \class TreeWriter
 *
 *  Fills ROOT tree branches.
 *
 *  $Date: 2013-11-08 17:09:58 +0100 (Fri, 08 Nov 2013) $
 *  $Revision: 1323 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <vector>

class TClass;
class TObjArray;
class TRefArray;

class Candidate;
class ExRootTreeBranch;

class TreeWriter: public DelphesModule
{
public:

  TreeWriter();
  ~TreeWriter();

  void Init();
  void Process();
  void Finish();

private:

  void ProcessLeptons(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);
  void ProcessParticles(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);
  void ProcessJets(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);
  void ProcessMissingET(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);
  void ProcessRho(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);
  void ProcessVertices(std::vector<ExRootTreeBranch*> branchVector, std::vector<TObjArray*> arrayVector);

#ifndef __CINT__
  typedef void (TreeWriter::*TProcessMethod)(std::vector<ExRootTreeBranch*>, std::vector<TObjArray*>); //!

  typedef std::map< std::vector<ExRootTreeBranch*>, std::pair< TProcessMethod, std::vector<TObjArray*> > > TBranchMap; //!

  TBranchMap fBranchMap; //!

  std::map< TClass *, TProcessMethod > fClassMap; //!
#endif

  ClassDef(TreeWriter, 1);
  int filep,fnlep;
  int fjet_type;
};

#endif
