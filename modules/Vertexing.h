#ifndef Vertexing_h
#define Vertexing_h

/** \class Merger
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  $Date: 2013-02-09 18:32:10 +0100 (Sat, 09 Feb 2013) $
 *  $Revision: 894 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>

class TIterator;
class TObjArray;

class Vertexing: public DelphesModule
{
public:

  Vertexing();
  ~Vertexing();

  void Init();
  void Process();
  void Finish();

private:

  TObjArray *fTrackInputArray; //!
  TObjArray *fVertexInputArray; //!
  TObjArray *fVertexOutputArray; //!
  
  float fpt_min; 
  
  ClassDef(Vertexing, 1)
};

#endif
