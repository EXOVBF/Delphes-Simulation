
/** \class TreeWriter
 *
 *  Fills ROOT tree branches.
 *
 *  $Date: 2014-01-03 16:41:28 +0100 (Fri, 03 Jan 2014) $
 *  $Revision: 1348 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \and you
 */

#include "modules/TreeWriter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TROOT.h"
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

TreeWriter::TreeWriter()
{
}

//------------------------------------------------------------------------------

TreeWriter::~TreeWriter()
{
}

//------------------------------------------------------------------------------

void TreeWriter::Init()
{ 
//  fClassMap[Photon::Class()] = &TreeWriter::ProcessPhotons;
  fClassMap[Electron::Class()] = &TreeWriter::ProcessLeptons;
  fClassMap[Muon::Class()] = &TreeWriter::ProcessLeptons;
  fClassMap[GenParticle::Class()] = &TreeWriter::ProcessParticles;
  fClassMap[Jet::Class()] = &TreeWriter::ProcessJets;
  fClassMap[MissingET::Class()] = &TreeWriter::ProcessMissingET;
  fClassMap[Rho::Class()] = &TreeWriter::ProcessRho;
  fClassMap[Vertex::Class()] = &TreeWriter::ProcessVertices;  

  TBranchMap::iterator itBranchMap;
  map< TClass *, TProcessMethod >::iterator itClassMap;

  // read branch configuration and
  // import array with output from filter/classifier/jetfinder modules

  ExRootConfParam param = GetParam("Branch");
  Long_t i, size;
  TString branchName, branchClassName, branchInputArray;
  TClass *branchClass;
  TObjArray *array;
  size = param.GetSize();
  vector<ExRootTreeBranch*> branchVector[size];
  vector<TObjArray*> arrayVector[size];
  TProcessMethod methodVector[size];
  int nFill=0;
  fnlep=0;
  filep=0;
  
  size = param.GetSize();
  for(i = 0; i < size/3; ++i)
  {
    branchInputArray = param[i*3].GetString();
    branchName = param[i*3 + 1].GetString();
    branchClassName = param[i*3 + 2].GetString();

    branchClass = gROOT->GetClass(branchClassName);

    if(!branchClass)
    {
      cout << "** ERROR: cannot find class '" << branchClassName << "'" << endl;
      continue;
    }
    
    itClassMap = fClassMap.find(branchClass);
    if(itClassMap == fClassMap.end())
    {
      cout << "** ERROR: cannot create branch for class '" << branchClassName << "'" << endl;
      continue;
    }
    
    array = ImportArray(branchInputArray);
    
    //--- lvj preselection
    // MAYBE
    //--- Leptons
    if(strcmp(branchName.Data(),"Electron")==0 || strcmp(branchName.Data(),"Muon")==0)
    {
      if(fnlep==0)
      {
        ExRootTreeBranch* branch_lep_0 = NewFloatBranch("lep_pt");
        ExRootTreeBranch* branch_lep_1 = NewFloatBranch("lep_eta");
        ExRootTreeBranch* branch_lep_2 = NewFloatBranch("lep_phi");
        ExRootTreeBranch* branch_lep_3 = NewFloatBranch("lep_flv");
        ExRootTreeBranch* branch_lep_4 = NewFloatBranch("lep_isolation");
        ExRootTreeBranch* branch_lep_5 = NewFloatBranch("lep_E_no_smearing");
        ExRootTreeBranch* branch_lep_6 = NewFloatBranch("lep_E_with_smearing");
        ExRootTreeBranch* branch_lep_7 = NewFloatBranch("lep_Hcal_over_Ecal");
        ExRootTreeBranch* branch_lep_8 = NewFloatBranch("lep_number");
        ExRootTreeBranch* branch_lep_9 = NewFloatBranch("gen_lep_pt");
        ExRootTreeBranch* branch_lep_10 = NewFloatBranch("gen_lep_eta");
        ExRootTreeBranch* branch_lep_11 = NewFloatBranch("gen_lep_phi");
        branchVector[nFill].push_back(branch_lep_0);
        branchVector[nFill].push_back(branch_lep_1);
        branchVector[nFill].push_back(branch_lep_2);
        branchVector[nFill].push_back(branch_lep_3);
        branchVector[nFill].push_back(branch_lep_4);
        branchVector[nFill].push_back(branch_lep_5);
        branchVector[nFill].push_back(branch_lep_6);
        branchVector[nFill].push_back(branch_lep_7);
        branchVector[nFill].push_back(branch_lep_8);
        branchVector[nFill].push_back(branch_lep_9);
        branchVector[nFill].push_back(branch_lep_10);
        branchVector[nFill].push_back(branch_lep_11);
        arrayVector[nFill].push_back(array);
        methodVector[nFill] = itClassMap->second;
        filep=nFill;
        nFill++;
      }
      else
      {
        arrayVector[filep].push_back(array);
      }      
      fnlep++;
    }
    //--- GenParticles -> just to tag leptons from tau decay
    if(strcmp(branchName.Data(),"Particle")==0)
    {
      ExRootTreeBranch* branchTau = NewFloatBranch("gen_was_tau");
      branchVector[nFill].push_back(branchTau);
      //--- gen MET
      ExRootTreeBranch* branchMET = NewFloatBranch("gen_MET");
      ExRootTreeBranch* branchEta = NewFloatBranch("gen_MET_eta");
      ExRootTreeBranch* branchPhi = NewFloatBranch("gen_MET_phi");
      branchVector[nFill].push_back(branchMET);
      branchVector[nFill].push_back(branchEta);
      branchVector[nFill].push_back(branchPhi);
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    } 
    //--- MET
    if(strcmp(branchName.Data(),"MissingET")==0)
    {
      ExRootTreeBranch* branchMET = NewFloatBranch("MET");
      ExRootTreeBranch* branchPhi = NewFloatBranch("MET_phi");
      branchVector[nFill].push_back(branchMET);
      branchVector[nFill].push_back(branchPhi);
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- jets
    if(branchName.Contains("jet"))
    {
      ExRootTreeBranch* branch_jet_0 = NewFloatBranch("number_"+branchName);
      ExRootTreeBranch* branch_jet_1 = NewFloatBranch(branchName+"_pt");
      ExRootTreeBranch* branch_jet_2 = NewFloatBranch(branchName+"_eta");
      ExRootTreeBranch* branch_jet_3 = NewFloatBranch(branchName+"_phi");
      ExRootTreeBranch* branch_jet_4 = NewFloatBranch(branchName+"_mass");
      ExRootTreeBranch* branch_jet_5 = NewFloatBranch(branchName+"_mass_pruned");
      ExRootTreeBranch* branch_jet_6 = NewFloatBranch(branchName+"_btag");
      ExRootTreeBranch* branch_jet_7 = NewFloatBranch(branchName+"_Hcal_over_Ecal");
      branchVector[nFill].push_back(branch_jet_0);
      branchVector[nFill].push_back(branch_jet_1);
      branchVector[nFill].push_back(branch_jet_2);
      branchVector[nFill].push_back(branch_jet_3);
      branchVector[nFill].push_back(branch_jet_4);
      branchVector[nFill].push_back(branch_jet_5);
      branchVector[nFill].push_back(branch_jet_6);
      branchVector[nFill].push_back(branch_jet_7);
      if(branchName.Contains("CA8"))
      {
        ExRootTreeBranch* branch_jet_8 = NewFloatBranch(branchName+"_tau1");
        ExRootTreeBranch* branch_jet_9 = NewFloatBranch(branchName+"_tau2");
        ExRootTreeBranch* branch_jet_10 = NewFloatBranch(branchName+"_tau3");
        branchVector[nFill].push_back(branch_jet_8);
        branchVector[nFill].push_back(branch_jet_9);
        branchVector[nFill].push_back(branch_jet_10);
      }
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- Rho
    if(strcmp(branchName.Data(),"Rho")==0)
    {
      ExRootTreeBranch* branchRho = NewFloatBranch("Rho_PU");
      ExRootTreeBranch* branchRhoFW = NewFloatBranch("Rho_PU_forward");
      branchVector[nFill].push_back(branchRho);
      branchVector[nFill].push_back(branchRhoFW);
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- Vertices
    if(strcmp(branchName.Data(),"Vertex")==0)
    {
      ExRootTreeBranch* branchReco = NewFloatBranch("nPV");
      ExRootTreeBranch* branchGen = NewFloatBranch("gen_nPV");
      branchVector[nFill].push_back(branchReco);
      branchVector[nFill].push_back(branchGen);
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
  }
  for(int iFill=0; iFill<nFill; iFill++)
  {
    fBranchMap.insert(make_pair(branchVector[iFill], make_pair(methodVector[iFill], arrayVector[iFill])));
  }
}

//------------------------------------------------------------------------------

void TreeWriter::Finish()
{
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessParticles(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  TLorentzVector met4vect, momentum_tmp;
  Candidate *candidate = 0;
  int tau_count=0;
  vector<float> *was_tau;
  vector<float> *met_gen,*met_phi_gen,*met_eta_gen;
  
  was_tau = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  met_gen = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  met_eta_gen = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  met_phi_gen = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
    
  // loop over all particles -> identifies taus and get gen_MET ==> v(s)
  iterator.Reset();
  int n=0;
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    if(TMath::Abs(candidate->PID) == 16 || TMath::Abs(candidate->PID) == 14 || TMath::Abs(candidate->PID) == 12)
    {
      momentum_tmp = candidate->Momentum;
      met4vect += momentum_tmp;
      n++;
    }
    if(TMath::Abs(candidate->PID) == 16)
    {
      tau_count++;
    }
  }
  was_tau->push_back(tau_count);
  Double_t signPz = (momentum_tmp.Pz() >= 0.0) ? 1.0 : -1.0;
  Double_t cosTheta = TMath::Abs(momentum_tmp.CosTheta());
  met_gen->push_back(momentum_tmp.Et());
  met_eta_gen->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum_tmp.Eta()));
  //met_eta_gen->push_back(momentum_tmp.Eta());
  met_phi_gen->push_back(momentum_tmp.Phi());
}
//------------------------------------------------------------------------------

void TreeWriter::ProcessVertices(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  Candidate *candidate=0;
  vector<float> *gen_nPV, *nPV, Z;
  int vertexR=0,vertexG=0;

  nPV = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  gen_nPV = (vector<float>*)((branchVector.at(1))->NewFloatEntry());

  // loop over all vertices
  // -> vertex closer then 100 mircon in the z direction are considered unresolved
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    TLorentzVector position = candidate->Position;
    Z.push_back(position.Z());
    vertexG++;
  }
  sort(Z.begin(), Z.end());
  for(int i=1; i<Z.size(); i++)
  {
    if(TMath::Abs(Z.at(i) - Z.at(i-1)) > 0.1)
    {
      vertexR++;
    }
  }
  gen_nPV->push_back(vertexG);
  nPV->push_back(vertexR);
}

//----------------------------------------------------------------------------------------
/*
void TreeWriter::ProcessPhotons(vector<ExRootTreeBranch*> branchVector, TObjArray *array)
{
  TIter iterator(array);
  Candidate *candidate = 0;
  Photon *entry = 0;
  Double_t pt, signPz, cosTheta, eta, rapidity;
  const Double_t c_light = 2.99792458E8;
  
  array->Sort();

  // loop over all photons
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    TIter it1(candidate->GetCandidates());
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;
    

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    entry = static_cast<Photon*>(branch->NewEntry());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;
    entry->E = momentum.E();
    
    entry->T = position.T()*1.0E-3/c_light;
   
    entry->EhadOverEem = candidate->Eem > 0.0 ? candidate->Ehad/candidate->Eem : 999.9;

    FillParticles(candidate, &entry->Particles);
  }
}
*/
//------------------------------------------------------------------------------
// Leptons

void TreeWriter::ProcessLeptons(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  Candidate *candidate = 0;
  Candidate *gen_candidate = 0;
  int number_lep=0;
  unsigned int iArray=0;
  vector<float> *pt, *eta, *phi, *flv, *iso, *Ein, *Eout, *HcalEcal, *nLep;
  vector<float> *pt_gen, *eta_gen, *phi_gen;

  pt = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  eta = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  phi = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  flv = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  iso = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  Ein = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  Eout = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  HcalEcal = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  nLep = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
  pt_gen = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
  eta_gen = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
  phi_gen = (vector<float>*)((branchVector.at(11))->NewFloatEntry());
  
  while(iArray < arrayVector.size())
  {
    TIter iterator = arrayVector.at(iArray);
    (arrayVector.at(iArray))->Sort();
    iterator.Reset();
    while((candidate = static_cast<Candidate*>(iterator.Next())))
    {
      //--- lep
      TLorentzVector momentum = candidate->Momentum;
      Double_t signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
      Double_t cosTheta = TMath::Abs(momentum.CosTheta());
      pt->push_back(momentum.Pt());
      eta->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum.Eta()));
      phi->push_back(momentum.Phi());
      iso->push_back(candidate->Isolation);
      Ein->push_back(candidate->P_in);
      Eout->push_back(candidate->P_out);
      HcalEcal->push_back((candidate->Ehad)/(candidate->Eem));
      flv->push_back(TMath::Abs(candidate->PID));
      //--- gen lep
      gen_candidate = static_cast<Candidate*>(candidate->GetCandidates()->At(0));
      momentum = gen_candidate->Momentum;
      signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
      cosTheta = TMath::Abs(momentum.CosTheta());
      pt_gen->push_back(momentum.Pt());
      eta_gen->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum.Eta()));
      phi_gen->push_back(momentum.Phi());
      number_lep++;
    }
    iArray++;
  }
  nLep->push_back(number_lep);   
}

//------------------------------------------------------------------------------
// Jets

void TreeWriter::ProcessJets(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  Candidate *candidate = 0, *constituent=0;
  int jets_count=0;
  float hcalEnergy=0,ecalEnergy=0;
  vector<float> *pt, *eta, *phi, *mass, *mass_pruned, *btag, *HcalEcal, *tau1, *tau2, *tau3, *n_jets;
  
  n_jets = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  pt = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  eta = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  phi = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  mass = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  mass_pruned = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  btag = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  HcalEcal = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  if(branchVector.size()==11)
  {
    tau1 = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
    tau2 = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
    tau3 = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
  }
    
  arrayVector.at(0)->Sort();

  // loop over all jets
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    const TLorentzVector &momentum = candidate->Momentum;
    Double_t signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    Double_t cosTheta = TMath::Abs(momentum.CosTheta());
    pt->push_back(momentum.Pt());
    eta->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum.Eta()));
    phi->push_back(momentum.Phi());
    mass->push_back(momentum.M());
    mass_pruned->push_back(candidate->PrunedMass);
    btag->push_back(candidate->BTag);
    //--- loop on the costituens to get Hcal/Ecal
    TIter itConstituents(candidate->GetCandidates());    
    while((constituent = static_cast<Candidate*>(itConstituents.Next())))
    {
      ecalEnergy += constituent->Eem;
      hcalEnergy += constituent->Ehad;
    }
    HcalEcal->push_back(ecalEnergy > 0.0 ? hcalEnergy/ecalEnergy : 999.9);
    //--- only for CA8 jets
    if(branchVector.size()==11)
    {
      tau1->push_back(candidate->tau1);
      tau2->push_back(candidate->tau2);
      tau3->push_back(candidate->tau3);
    }
    jets_count++;
  }
  n_jets->push_back(jets_count);
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMissingET(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  Candidate *candidate = 0;
  vector<float> *met,*met_phi;

  met = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  met_phi = (vector<float>*)((branchVector.at(1))->NewFloatEntry());

  // get the first entry
  if((candidate = static_cast<Candidate*>((arrayVector.at(0))->At(0))))
  {
    //--- MET
    TLorentzVector momentum = candidate->Momentum;
    met->push_back(momentum.Pt());
    met_phi->push_back((-momentum).Phi());
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessRho(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  vector<float> *rho_tmp,*rho_fw_tmp;

  // loop over all rho
  iterator.Reset();
  rho_tmp = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  rho_fw_tmp = (vector<float>*)((branchVector.at(1))->NewFloatEntry());

  rho_tmp->push_back((static_cast<Candidate*>(iterator.Next()))->Momentum.E());
  rho_fw_tmp->push_back((static_cast<Candidate*>(iterator.Next()))->Momentum.E());
}

//------------------------------------------------------------------------------

void TreeWriter::Process()
{
  TBranchMap::iterator itBranchMap;
  vector<ExRootTreeBranch*> branchVector;
  TProcessMethod method;
  vector<TObjArray*> arrayVector;

  for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
  {
    branchVector = itBranchMap->first;
    method = itBranchMap->second.first;
    arrayVector = itBranchMap->second.second;

    (this->*method)(branchVector, arrayVector);
  }
}

//------------------------------------------------------------------------------
