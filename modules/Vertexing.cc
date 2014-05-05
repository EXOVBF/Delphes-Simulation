
/** \class Vertexing
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  $Date: 2013-04-26 12:39:14 +0200 (Fri, 26 Apr 2013) $
 *  $Revision: 1099 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Vertexing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------

Vertexing::Vertexing():
    fpt_min(0)
{
}

//------------------------------------------------------------------------------

Vertexing::~Vertexing()
{
}

//------------------------------------------------------------------------------

void Vertexing::Init()
{
  // import arrays with output from other modules

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "tracks"));
  
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "vertices"));
  
  fpt_min = GetDouble("pt_min",0.1);
}

//------------------------------------------------------------------------------

void Vertexing::Finish()
{
  if(fTrackInputArray) delete fTrackInputArray;
  if(fVertexInputArray) delete fVertexInputArray;
}

//------------------------------------------------------------------------------

void Vertexing::Process()
{
  Candidate *candidate_track=0, *candidate_vertex=0;
  TIter track_iter(fTrackInputArray);
  TIter vertex_iter(fVertexInputArray);
  vector<double> sumPt;
  int vertex_id_tmp=0;
  
  for(int i=0; i<fVertexInputArray->GetEntries(); i++)
  {
    sumPt.push_back(0);
  }

  // loop over all input tracks arrays
  track_iter.Reset();
  vertex_iter.Reset();
  while((candidate_track = static_cast<Candidate*>(track_iter.Next())))
  {
    vertex_id_tmp = candidate_track->VertexID_gen;
    Double_t pt_tmp = (candidate_track->Momentum).Pt();
    if( abs(pt_tmp) > fpt_min )
    {
      sumPt.at(vertex_id_tmp) = sumPt.at(vertex_id_tmp) + pt_tmp*pt_tmp;
    }
  }
  while((candidate_vertex = static_cast<Candidate*>(vertex_iter.Next())))
  {
    vertex_id_tmp = candidate_vertex->VertexID_gen;
    candidate_vertex->sumPtSquare = sumPt.at(vertex_id_tmp);
  }
}

//------------------------------------------------------------------------------
