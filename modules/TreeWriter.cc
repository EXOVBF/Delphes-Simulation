
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
  fsignal = GetBool("isSignal",false);  
  
  fClassMap[Photon::Class()] = &TreeWriter::ProcessPhotons;
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
        ExRootTreeBranch* branch_lep_7 = NewFloatBranch("lep_n_particle_cone");
        ExRootTreeBranch* branch_lep_8 = NewFloatBranch("lep_number");
        ExRootTreeBranch* branch_lep_9 = NewFloatBranch("lep_X_vertex");
        ExRootTreeBranch* branch_lep_10 = NewFloatBranch("lep_Y_vertex");
        ExRootTreeBranch* branch_lep_11 = NewFloatBranch("lep_Z_vertex");
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
      //--- number of generated leptons
      ExRootTreeBranch* branchNlep = NewFloatBranch("gen_lep_number");
      branchVector[nFill].push_back(branchNlep);
      //--- identify leptons from tau
      ExRootTreeBranch* branchTau = NewFloatBranch("gen_was_tau");
      branchVector[nFill].push_back(branchTau);
      //--- gen MET
      ExRootTreeBranch* branchMET = NewFloatBranch("gen_MET");
      ExRootTreeBranch* branchEta = NewFloatBranch("gen_MET_eta");
      ExRootTreeBranch* branchPhi = NewFloatBranch("gen_MET_phi");
      branchVector[nFill].push_back(branchMET);
      branchVector[nFill].push_back(branchEta);
      branchVector[nFill].push_back(branchPhi);      
      //--- genXmass
      if(fsignal)
      {
        ExRootTreeBranch* branchX_pt = NewFloatBranch("gen_X_pt");
        ExRootTreeBranch* branchX_eta = NewFloatBranch("gen_X_eta");
        ExRootTreeBranch* branchX_phi = NewFloatBranch("gen_X_phi");
        ExRootTreeBranch* branchX_m = NewFloatBranch("gen_X_mass");
        branchVector[nFill].push_back(branchX_pt);
        branchVector[nFill].push_back(branchX_eta);
        branchVector[nFill].push_back(branchX_phi);
        branchVector[nFill].push_back(branchX_m);
      }
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
      ExRootTreeBranch* branch_jet_6 = NewFloatBranch(branchName+"_jetcharge");
      ExRootTreeBranch* branch_jet_7 = NewFloatBranch(branchName+"_btag");
      ExRootTreeBranch* branch_jet_8 = NewFloatBranch(branchName+"_Hcal_over_Ecal");
      branchVector[nFill].push_back(branch_jet_0);
      branchVector[nFill].push_back(branch_jet_1);
      branchVector[nFill].push_back(branch_jet_2);
      branchVector[nFill].push_back(branch_jet_3);
      branchVector[nFill].push_back(branch_jet_4);
      branchVector[nFill].push_back(branch_jet_5);
      branchVector[nFill].push_back(branch_jet_6);
      branchVector[nFill].push_back(branch_jet_7);
      branchVector[nFill].push_back(branch_jet_8);
      if(branchName.Contains("CA8"))
      {
        ExRootTreeBranch* branch_jet_9 = NewFloatBranch(branchName+"_tau1");
        ExRootTreeBranch* branch_jet_10 = NewFloatBranch(branchName+"_tau2");
        ExRootTreeBranch* branch_jet_11 = NewFloatBranch(branchName+"_tau3");
        branchVector[nFill].push_back(branch_jet_9);
        branchVector[nFill].push_back(branch_jet_10);
        branchVector[nFill].push_back(branch_jet_11);
      }
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- Rho
    if(strcmp(branchName.Data(),"Rho")==0)
    {
      ExRootTreeBranch* branchRho_0 = NewFloatBranch("Rho_PU_barrel");
      ExRootTreeBranch* branchRho_1 = NewFloatBranch("Rho_PU_endcap");
      branchVector[nFill].push_back(branchRho_0);
      branchVector[nFill].push_back(branchRho_1);
      if(array->GetSize())
      {
        ExRootTreeBranch* branchRho_2 = NewFloatBranch("Rho_PU_foward");
        branchVector[nFill].push_back(branchRho_2);
      }
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- Vertices
    if(strcmp(branchName.Data(),"Vertex")==0)
    {
      ExRootTreeBranch* branchRecoN = NewFloatBranch("nPV");
      ExRootTreeBranch* branchRecoX = NewFloatBranch("vertex_X");
      ExRootTreeBranch* branchRecoY = NewFloatBranch("vertex_Y");
      ExRootTreeBranch* branchRecoZ = NewFloatBranch("vertex_Z");
      ExRootTreeBranch* branchRecoS = NewFloatBranch("sum_pt_square");
      ExRootTreeBranch* branchGenN = NewFloatBranch("gen_nPV");
      ExRootTreeBranch* branchGenPU = NewFloatBranch("gen_is_PU");
      branchVector[nFill].push_back(branchRecoN);
      branchVector[nFill].push_back(branchRecoX);
      branchVector[nFill].push_back(branchRecoY);
      branchVector[nFill].push_back(branchRecoZ);
      branchVector[nFill].push_back(branchRecoS);
      branchVector[nFill].push_back(branchGenN);
      branchVector[nFill].push_back(branchGenPU);
      arrayVector[nFill].push_back(array);
      methodVector[nFill] = itClassMap->second;
      nFill++;
    }
    //--- Photons 
    if(strcmp(branchName.Data(),"Photon")==0)
    {
      ExRootTreeBranch* branch_pho_0 = NewFloatBranch("photon_pt");
      ExRootTreeBranch* branch_pho_1 = NewFloatBranch("photon_eta");
      ExRootTreeBranch* branch_pho_2 = NewFloatBranch("photon_phi");
      ExRootTreeBranch* branch_pho_3 = NewFloatBranch("photon_E");
      ExRootTreeBranch* branch_pho_4 = NewFloatBranch("photon_isolation");
      ExRootTreeBranch* branch_pho_5 = NewFloatBranch("photon_n_particle_cone");
      ExRootTreeBranch* branch_pho_6 = NewFloatBranch("photon_number");
      ExRootTreeBranch* branch_pho_7 = NewFloatBranch("photon_Hcal_over_Ecal");
      branchVector[nFill].push_back(branch_pho_0);
      branchVector[nFill].push_back(branch_pho_1);
      branchVector[nFill].push_back(branch_pho_2);
      branchVector[nFill].push_back(branch_pho_3);
      branchVector[nFill].push_back(branch_pho_4);
      branchVector[nFill].push_back(branch_pho_5);
      branchVector[nFill].push_back(branch_pho_6);
      branchVector[nFill].push_back(branch_pho_7);
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
// GenParticles
// all the trick in here is done to manage the mess in the Pythia8 event listing.
// due tu showering the same particle may appears multiple times, in sometimes the good
// one is the first, sometimes the last...

void TreeWriter::ProcessParticles(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  TLorentzVector met4vect, momentum_tmp;
  Candidate *candidate = 0;
  int lep_count=0, tau_count=0;
  vector<int> W1_code, W2_code;
  vector<float> *lep_gen,*was_tau;
  vector<float> *met_gen,*met_phi_gen,*met_eta_gen;
  vector<float> *x_pt, *x_phi, *x_eta, *x_mass;

  lep_gen = (vector<float>*)((branchVector.at(0))->NewFloatEntry());  
  was_tau = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  met_gen = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  met_eta_gen = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  met_phi_gen = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  if(fsignal)
  {
    x_pt = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
    x_eta = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
    x_phi = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
    x_mass = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
    // initialization
    x_pt->push_back(0);
    x_eta->push_back(0);
    x_phi->push_back(0);
    x_mass->push_back(0);
  }
    
  // loop over all particles -> identifies taus and get gen_MET ==> v(s)
  iterator.Reset();
  int i=0;
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    int PID_tmp = candidate->PID;
    if(PID_tmp == 24)
    {
      W1_code.push_back(i);
    }
    if(PID_tmp == -24)
    {
      W2_code.push_back(i);
    }
    i++;
  }
  if( W1_code.size() == 0 )
  {
    W1_code.push_back(-3);
  }
  if( W2_code.size() == 0 )
  {
    W2_code.push_back(-3);
  }
  iterator.Reset();  
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    int PID_tmp = TMath::Abs(candidate->PID);
    int M1 = candidate->M1 + 1;
    int D1 = candidate->D1 + 1;
    int D2 = candidate->D2 + 1;
    if((PID_tmp == 16 || PID_tmp == 14 || PID_tmp == 12) && candidate->Status > 0 && candidate->IsPU == 0)
    {
      momentum_tmp = candidate->Momentum;
      met4vect += momentum_tmp;
    }
    if((PID_tmp == 11 || PID_tmp == 13 || PID_tmp == 15) && (M1 == W1_code.back() || M1 == W2_code.back()))
    {
      lep_count++;       
    }
    if(PID_tmp == 15 && (M1 == W1_code.back() || M1 == W2_code.back()))
    {
      tau_count++;
    }
    //--- graviton gen infos, find the right one linking it with thw 2 W
    if(PID_tmp == 39 && ((D1 == W1_code.front() && D2 == W2_code.front()) || (D2 == W1_code.front() && D1 == W2_code.front()) ))
    {
      x_pt->at(0) = (candidate->Momentum).Pt();
      x_eta->at(0) = (candidate->Momentum).Eta();
      x_phi->at(0) = (candidate->Momentum).Phi();
      x_mass->at(0) = (candidate->Momentum).M();
    }
  }
  lep_gen->push_back(lep_count);
  was_tau->push_back(tau_count);
  Double_t signPz = (momentum_tmp.Pz() >= 0.0) ? 1.0 : -1.0;
  Double_t cosTheta = TMath::Abs(momentum_tmp.CosTheta());
  met_gen->push_back(momentum_tmp.Et());
  met_eta_gen->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum_tmp.Eta()));
  met_phi_gen->push_back(momentum_tmp.Phi());
}
//------------------------------------------------------------------------------
// Vertices

void TreeWriter::ProcessVertices(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator1(arrayVector.at(0));
  TIter iterator2(arrayVector.at(0));
  Candidate *candidate1=0;
  Candidate *candidate2=0;
  vector<float> *nPV, *X, *Y, *Z, *S, *gen_nPV, *gen_PU;
  vector<int> skip;
  int vertexR=0,i=0;
  float X_tmp, Y_tmp, Z_tmp, S_tmp, PU_tmp;

  nPV = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  X = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  Y = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  Z = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  S = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  gen_nPV = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  gen_PU = (vector<float>*)((branchVector.at(6))->NewFloatEntry());  

  // loop over all vertices
  // -> vertex closer then 100 mircon in the z direction are considered unresolved
  iterator1.Reset();
  while((candidate1 = static_cast<Candidate*>(iterator1.Next())))
  {
    i++;
    PU_tmp = 1;
    X_tmp = 0;
    Y_tmp = 0;    
    Z_tmp = 0;
    S_tmp = 0;
    int j=0;
    int n=1;
    if( find(skip.begin(), skip.end(), i) != skip.end() )
    {
      continue;
    }
    TLorentzVector position1 = candidate1->Position;
    iterator2.Reset();
    while((candidate2 = static_cast<Candidate*>(iterator2.Next())))
    {
      j++;
      TLorentzVector position2 = candidate2->Position;
      if( TMath::Abs(position1.Z() - position2.Z()) < 0.1 && candidate2->VertexID_gen != candidate1->VertexID_gen)
      {
        if(candidate2->IsPU == 0)
        {
          PU_tmp = 0;
        }
        X_tmp = X_tmp + position2.X();
        Y_tmp = Y_tmp + position2.Y();   
        Z_tmp = Z_tmp + position2.Z();
        S_tmp = S_tmp + candidate2->sumPtSquare;
        skip.push_back(j);
        n++;
      }
    }
    //---reco
    if(candidate1->IsPU == 0)
    {
      PU_tmp = 0;
    }
    X_tmp = X_tmp + position1.X();
    X_tmp = X_tmp / n;
    Y_tmp = Y_tmp + position1.Y();
    Y_tmp = Y_tmp / n;
    Z_tmp = Z_tmp + position1.Z();
    Z_tmp = Z_tmp / n;
    S_tmp = S_tmp + candidate1->sumPtSquare;
    X->push_back(X_tmp);
    Y->push_back(Y_tmp);
    Z->push_back(Z_tmp);
    S->push_back(S_tmp);
    gen_PU->push_back(PU_tmp);
    vertexR++;
  }
  gen_nPV->push_back((arrayVector.at(0))->GetEntries());
  nPV->push_back(vertexR);
}

//----------------------------------------------------------------------------------------
// Photons

void TreeWriter::ProcessPhotons(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  Candidate *candidate = 0;
  int number=0;
  vector<float> *pt, *eta, *phi, *E, *iso, *nPartCone, *nPho, *HcalEcal;

  pt = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  eta = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  phi = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  E = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  iso = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  nPartCone = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  nPho = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  HcalEcal = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  
  arrayVector.at(0)->Sort();

  // loop over all photons
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
      TLorentzVector momentum = candidate->Momentum;
      Double_t signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
      Double_t cosTheta = TMath::Abs(momentum.CosTheta());
      float EhadOverEem = candidate->Eem > 0.0 ? candidate->Ehad/candidate->Eem : 999.9;
      
      pt->push_back(momentum.Pt());
      eta->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum.Eta()));
      phi->push_back(momentum.Phi());
      E->push_back(momentum.E());
      iso->push_back(candidate->Isolation);
      nPartCone->push_back((candidate->ParticleInCone));
      HcalEcal->push_back(EhadOverEem);
      number++;
  }
  nPho->push_back(number);
}

//------------------------------------------------------------------------------
// Leptons

void TreeWriter::ProcessLeptons(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  Candidate *candidate = 0;
  Candidate *gen_candidate = 0;
  int number_lep=0;
  unsigned int iArray=0;
  vector<float> *pt, *eta, *phi, *flv, *iso, *Ein, *Eout, *nPartCone, *nLep; 
  vector<float> *X, *Y, *Z;

  pt = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  eta = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  phi = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  flv = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  iso = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  Ein = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  Eout = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  nPartCone = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  nLep = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
  X = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
  Y = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
  Z = (vector<float>*)((branchVector.at(11))->NewFloatEntry());
  
  while(iArray < arrayVector.size())
  {
    TIter iterator = arrayVector.at(iArray);
    (arrayVector.at(iArray))->Sort();
    iterator.Reset();
    while((candidate = static_cast<Candidate*>(iterator.Next())))
    {
      TLorentzVector momentum = candidate->Momentum;
      Double_t signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
      Double_t cosTheta = TMath::Abs(momentum.CosTheta());

      pt->push_back(momentum.Pt());
      eta->push_back((cosTheta == 1.0 ? signPz*999.9 : momentum.Eta()));
      phi->push_back(momentum.Phi());
      iso->push_back(candidate->Isolation);
      Ein->push_back(candidate->P_in);
      Eout->push_back(candidate->P_out);
      nPartCone->push_back((candidate->ParticleInCone));
      flv->push_back(candidate->PID);
      gen_candidate = static_cast<Candidate*>(candidate->GetCandidates()->At(0));
      X->push_back((gen_candidate->Position).X());
      Y->push_back((gen_candidate->Position).Y());
      Z->push_back((gen_candidate->Position).Z());
      number_lep++;
    }
    iArray++;
  }
  nLep->push_back(number_lep);
//  cout << "leptons number:  " << number_lep << endl;   
}

//------------------------------------------------------------------------------
// Jets

void TreeWriter::ProcessJets(vector<ExRootTreeBranch*> branchVector, vector<TObjArray*> arrayVector)
{
  TIter iterator(arrayVector.at(0));
  Candidate *candidate = 0, *constituent=0;
  int jets_count=0;
  float hcalEnergy=0,ecalEnergy=0;
  vector<float> *pt, *eta, *phi, *mass, *mass_pruned, *jetcharge, *btag, *HcalEcal, *tau1, *tau2, *tau3, *n_jets;
  
  n_jets = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  pt = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  eta = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  phi = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  mass = (vector<float>*)((branchVector.at(4))->NewFloatEntry());
  mass_pruned = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  jetcharge   = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  btag = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  HcalEcal = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
  if(branchVector.size()==11)
  {
    tau1 = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
    tau2 = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
    tau3 = (vector<float>*)((branchVector.at(11))->NewFloatEntry());
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
    jetcharge->push_back(candidate->jetCharge);
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
      tau1->push_back(candidate->Tau[0]);
      tau2->push_back(candidate->Tau[1]);
      tau3->push_back(candidate->Tau[2]);
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
  Candidate* candidate=0;
  vector<float> *rho_tmp[3];
  int i=0;

  // loop over all rho
  iterator.Reset();
  rho_tmp[0] = (vector<float>*)((branchVector.at(0))->NewFloatEntry());
  rho_tmp[1] = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  if(branchVector.size()==3)
  {
    rho_tmp[2] = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  }
  while((candidate = static_cast<Candidate*>(iterator.Next())))
  {
    rho_tmp[i]->push_back(candidate->Momentum.E());
    i++;
  }
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
