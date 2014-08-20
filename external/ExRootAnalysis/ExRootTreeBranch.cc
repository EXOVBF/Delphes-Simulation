
/** \class ExRootTreeBranch
 *
 *  Class handling object creation
 *  It is also used for output ROOT tree branches
 *
 *  $Date: 2008-06-04 13:57:56 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

ExRootTreeBranch::ExRootTreeBranch(const char *name, TClass *cl, TTree *tree) :
    fSize(0), fCapacity(1),fDataFloat(0), fData(0)
{
    stringstream message;
//  cl->IgnoreTObjectStreamer();
//  cl->BypassStreamer();
    fData = new TClonesArray(cl, fCapacity);

    if(fData)
    {
	fData->SetName(name);
	fData->ExpandCreateFast(fCapacity);
	fData->Clear();
	if(tree)
	{
	    tree->Branch(name, &fData, 64000);
	    tree->Branch(TString(name) + "_size", &fSize, TString(name) + "_size/I");
	}
    }
    else
    {
	message << "can't create TClonesArray for branch '" << name << "'";
	throw runtime_error(message.str());
    }
}

//------------------------------------------------------------------------------
// construct a plain branch

ExRootTreeBranch::ExRootTreeBranch(const char *name, TTree *tree) :
    fDataFloat(0), fData(0)
{
    if(tree)
    {
	tree->Branch(name, &fDataFloat);
    }
}

//------------------------------------------------------------------------------

ExRootTreeBranch::~ExRootTreeBranch()
{
    if(fData) delete fData;
//    if(fDataFloat) delete fDataFloat;
}

//------------------------------------------------------------------------------

TObject *ExRootTreeBranch::NewEntry()
{
    if(!fData) return 0;

    if(fSize >= fCapacity)
    {
	if(fCapacity < 10) fCapacity = 10;
	else if(fCapacity < 30) fCapacity = 30;
	else if(fCapacity < 100) fCapacity = 100;
	else if(fCapacity < 250) fCapacity = 250;
	else fCapacity *= 2;

	fData->ExpandCreateFast(fCapacity);

	fData->Clear();
	fData->ExpandCreateFast(fSize);
    }
  
    return fData->AddrAt(fSize++);
}

//------------------------------------------------------------------------------

vector<float>* ExRootTreeBranch::NewFloatEntry()
{
    fDataFloat->clear();
    return fDataFloat;
}

//------------------------------------------------------------------------------

void ExRootTreeBranch::Clear()
{
    fSize = 0;
    if(fData) fData->Clear();
}

//------------------------------------------------------------------------------

