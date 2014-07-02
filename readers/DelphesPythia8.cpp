#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

#include <signal.h>

#include "Pythia.h"

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//**************************************************************************************************

bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut, vector<ExRootTreeBranch*> branchVector)
{
    vector<TLorentzVector> outPartons;
    int leptons=0;
    int Mjj_check=0;
    
    //---search for final state partons in the event---
    for(int iPart = 0; iPart < LHE_event->size(); iPart++)
    {
        vector<float> particle = LHE_event->at(iPart);
        TLorentzVector tmp4vect;
        tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));
        if(particle.at(1) == 1 && ((abs(particle.at(0)) > 0 && abs(particle.at(0)) < 7) || abs(particle.at(0)) == 21))
        {
            outPartons.push_back(tmp4vect);
        }
        if(particle.at(1) == 1 && (abs(particle.at(0)) == 11 || abs(particle.at(0)) == 13 || abs(particle.at(0)) == 15))
        {
            leptons++;
        }
    }
    //---reject fully hadronic events
    if(leptons < 1)
    {
        return false;
    }
    //---apply the VBF preselection cut on Mjj---
    for(int iPartons = 0; iPartons < outPartons.size(); iPartons++)
    {
        TLorentzVector tmp4vect = outPartons.at(iPartons);
        for(int jQuark = iPartons+1; jQuark < outPartons.size(); jQuark++)
        {
            tmp4vect += outPartons.at(jQuark);
            if(tmp4vect.M() > Mjj_cut)
            {
                Mjj_check = 1;
            }
        }
    }
    if( outPartons.size() > 2 && Mjj_check == 0 )
    {
        return false;
    }
    //---if preselection is passed save gen_lhe infos---
    vector<int> W_codes;
    vector<float> *n_lep_gen;
    vector<float> *lep_pt, *lep_eta, *lep_phi, *lep_flv;
    vector<float> *nu_pt, *nu_eta, *nu_phi, *nu_flv;
    vector<float> *p_pt, *p_eta, *p_phi, *p_flv, *p_fW;
    vector<float> *W_pt, *W_phi, *W_eta, *W_m, *W_pid;
    vector<float> *x_pt, *x_phi, *x_eta, *x_m;

    n_lep_gen = (vector<float>*)((branchVector.at(0))->NewFloatEntry());  
    lep_pt = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
    lep_eta = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
    lep_phi = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
    lep_flv = (vector<float>*)((branchVector.at(4))->NewFloatEntry());    
    nu_pt = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
    nu_eta = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
    nu_phi = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
    nu_flv = (vector<float>*)((branchVector.at(8))->NewFloatEntry());
    p_pt = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
    p_eta = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
    p_phi = (vector<float>*)((branchVector.at(11))->NewFloatEntry());
    p_flv = (vector<float>*)((branchVector.at(12))->NewFloatEntry());
    p_fW = (vector<float>*)((branchVector.at(13))->NewFloatEntry());
    W_pt = (vector<float>*)((branchVector.at(14))->NewFloatEntry());
    W_eta = (vector<float>*)((branchVector.at(15))->NewFloatEntry());
    W_phi = (vector<float>*)((branchVector.at(16))->NewFloatEntry());
    W_m = (vector<float>*)((branchVector.at(17))->NewFloatEntry());
    W_pid = (vector<float>*)((branchVector.at(18))->NewFloatEntry());
    if( branchVector.size() > 19 )
    {
        x_pt = (vector<float>*)((branchVector.at(19))->NewFloatEntry());
        x_eta = (vector<float>*)((branchVector.at(20))->NewFloatEntry());
        x_phi = (vector<float>*)((branchVector.at(21))->NewFloatEntry());
        x_m = (vector<float>*)((branchVector.at(22))->NewFloatEntry());    
    }
    //---loop on lhe events particle searching for W's---
    for(int iPart = 0; iPart < LHE_event->size(); iPart++)
    { 
        vector<float> particle = LHE_event->at(iPart);
        TLorentzVector tmp4vect;
        tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));   
        if(particle.at(1) == 2 && abs(particle.at(0)) == 24)
        { 
            W_pt->push_back(tmp4vect.Pt());
            W_eta->push_back(tmp4vect.Eta());
            W_phi->push_back(tmp4vect.Phi());
            W_m->push_back(tmp4vect.M());
            W_pid->push_back(particle.at(0));
            W_codes.push_back(iPart+1);
        }
    }
    //---loop on lhe events particle---
    for(int iPart = 0; iPart < LHE_event->size(); iPart++)
    { 
        vector<float> particle = LHE_event->at(iPart);
        TLorentzVector tmp4vect;
        tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));
        //---Graviton
        if(branchVector.size() > 19 && particle.at(1) == 2 && abs(particle.at(0)) == 39)
        {
            x_pt->push_back(tmp4vect.Pt());
            x_eta->push_back(tmp4vect.Eta());
            x_phi->push_back(tmp4vect.Phi());
            x_m->push_back(tmp4vect.M());
        }
        //---Partons
        if(particle.at(1) == 1 && ((abs(particle.at(0)) > 0 && abs(particle.at(0)) < 7) || abs(particle.at(0)) == 21))
        {
            p_pt->push_back(tmp4vect.Pt());
            p_eta->push_back(tmp4vect.Eta());
            p_phi->push_back(tmp4vect.Phi());
            p_flv->push_back(particle.at(0));
            if( particle.at(2) == particle.at(3) && (particle.at(2) == W_codes.front() || particle.at(2) == W_codes.back()) )
            {
                p_fW->push_back(1);
            }
            else
            {
                p_fW->push_back(0);
            }
        }
        //---charged leptons
        if(particle.at(1) == 1 && (abs(particle.at(0)) == 11 || abs(particle.at(0)) == 13 || abs(particle.at(0)) == 15))
        {
            lep_pt->push_back(tmp4vect.Pt());
            lep_eta->push_back(tmp4vect.Eta());
            lep_phi->push_back(tmp4vect.Phi());
            lep_flv->push_back(particle.at(0));
        }
        //---neutrinos
        if(particle.at(1) == 1 && (abs(particle.at(0)) == 12 || abs(particle.at(0)) == 14 || abs(particle.at(0)) == 16))
        {
            nu_pt->push_back(tmp4vect.Pt());
            nu_eta->push_back(tmp4vect.Eta());
            nu_phi->push_back(tmp4vect.Phi());
            nu_flv->push_back(particle.at(0));
        }
    }
    n_lep_gen->push_back(leptons);
    
    return true;
}

//**************************************************************************************************

void ConvertInput(Long64_t eventCounter, Pythia8::Pythia *pythia,
		  ExRootTreeBranch *branch, DelphesFactory *factory,
		  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray,
		  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
    int i;

    HepMCEvent *element=0;
    Candidate *candidate=0;
    TDatabasePDG *pdg=0;
    TParticlePDG *pdgParticle=0;
    Int_t pdgCode;

    Int_t pid, status;
    Double_t px, py, pz, e, mass;
    Double_t x, y, z, t;

    // event information
    element = static_cast<HepMCEvent *>(branch->NewEntry());

    element->Number = eventCounter;

    element->ProcessID = pythia->info.code();
    element->MPI = 1;
    element->Weight = pythia->info.weight();
    element->Scale = pythia->info.QRen();
    element->AlphaQED = pythia->info.alphaEM();
    element->AlphaQCD = pythia->info.alphaS();

    element->ID1 = pythia->info.id1();
    element->ID2 = pythia->info.id2();
    element->X1 = pythia->info.x1();
    element->X2 = pythia->info.x2();
    element->ScalePDF = pythia->info.QFac();
    element->PDF1 = pythia->info.pdf1();
    element->PDF2 = pythia->info.pdf2();

    element->ReadTime = readStopWatch->RealTime();
    element->ProcTime = procStopWatch->RealTime();

    pdg = TDatabasePDG::Instance();

    for(i = 0; i < pythia->event.size(); ++i)
    {
        Pythia8::Particle &particle = pythia->event[i];
        
        pid = particle.id();
        status = pythia->event.statusHepMC(i);

        px = particle.px(); py = particle.py(); pz = particle.pz(); e = particle.e(); mass = particle.m();
        x = particle.xProd(); y = particle.yProd(); z = particle.zProd(); t = particle.tProd();

        candidate = factory->NewCandidate();

        candidate->PID = pid;
        pdgCode = TMath::Abs(candidate->PID);

        candidate->Status = status;

        candidate->M1 = particle.mother1() - 1;
        candidate->M2 = particle.mother2() - 1;

        candidate->D1 = particle.daughter1() - 1;
        candidate->D2 = particle.daughter2() - 1;

        pdgParticle = pdg->GetParticle(pid);
        candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
        candidate->Mass = mass;

        candidate->Momentum.SetPxPyPzE(px, py, pz, e);

        candidate->Position.SetXYZT(x, y, z, t);

        allParticleOutputArray->Add(candidate);

        if(status == 1)
        {
            stableParticleOutputArray->Add(candidate);
        }
        else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
        {
            partonOutputArray->Add(candidate);
        }
    }
}

//****************************************************************************************

static bool interrupted = false;

void SignalHandler(int sig)
{
    interrupted = true;
}

//****************************************************************************************
// main

int main(int argc, char *argv[])
{
//----------------------------------------------------------------------------------------
// definitions

    char appName[] = "DelphesPythia8";
    stringstream message;
    string inputFile;
    TFile *outputFile = 0;
    TStopwatch readStopWatch, procStopWatch;
    ExRootTreeWriter *treeWriter = 0;
    ExRootTreeWriter *treeHepMC = 0;    
    ExRootTreeBranch *branchEvent = 0;
    vector<ExRootTreeBranch*> branchGen;
    ExRootConfReader *confReader = 0;
    Delphes *modularDelphes = 0;
    DelphesFactory *factory = 0;
    TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
    Long64_t eventCounter, errorCounter;

    Pythia8::Pythia *pythia = 0;

//----------------------------------------------------------------------------------------

    if(argc < 4)
    {
        cout << "------------------------------------Manual----------------------------------------" << endl;
        cout << "- Usage: " << appName << " config_file" << " lhe_file" << " output_file" << " Mjj_cut" << " start" << " number" << " signal" << endl;
        cout << "- config_file          ->  configuration file in Tcl format" << endl;
        cout << "- input_file           ->  lhe file for Pythia8" << endl;
        cout << "- output_file          ->  output file in ROOT format" << endl;
        cout << "- Mjj_cut  (optional)  ->  cut on Mjj in GeV -- default = 150 GeV" << endl;
        cout << "- start    (optional)  ->  number of starting event" << endl;
        cout << "- number   (optional)  ->  number of total events to be processed" << endl;
        cout << "- signal   (optional)  ->  1 if the sample is a graviton signal sample" << endl;
        cout << "----------------------------------------------------------------------------------" << endl << endl;  
        return 1;
    }

    signal(SIGINT, SignalHandler);

    gROOT->SetBatch();

    int appargc = 1;
    char *appargv[] = {appName};
    TApplication app(appName, &appargc, appargv);

//--------------------------------------------------------------------------------------------------
// do the simulation

    try
    {
        inputFile = argv[2];
        outputFile = TFile::Open(argv[3], "RECREATE");

        if(outputFile == NULL)
        {
            message << "can't create output file " << argv[3];
            throw runtime_error(message.str());
        }
        int fsignal=0;        
        if( argc >= 8)
        {
            fsignal = atoi(argv[7]);
        }
        //--- create output tree ---
        treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
        //--- create gen (lhe level) branch ---
        //--- lep number
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_number"));      
        //--- leptons
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pt"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_eta"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_phi"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_flv"));
        //--- neutrinos
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pt"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_eta"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_phi"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_flv"));
        //--- gen partons infos
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_p_pt"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_p_eta"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_p_phi"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_p_flv"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_p_from_W"));
        //--- gen W infos
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_pt"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_eta"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_phi"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_mass"));
        branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_pid"));
        //--- gen Graviton infos
        if(fsignal == 1)
        {
            branchGen.push_back(treeWriter->NewFloatBranch("lhe_X_pt"));
            branchGen.push_back(treeWriter->NewFloatBranch("lhe_X_eta"));
            branchGen.push_back(treeWriter->NewFloatBranch("lhe_X_phi"));
            branchGen.push_back(treeWriter->NewFloatBranch("lhe_X_mass"));
        }
        
        //--- deals with the HepMc output of Pythia8 ---> no need to store it
        treeHepMC = new ExRootTreeWriter();
        branchEvent = treeHepMC->NewBranch("Event", HepMCEvent::Class());

        //----- Delphes init -----        
        confReader = new ExRootConfReader;
        confReader->ReadFile(argv[1]);

        modularDelphes = new Delphes("Delphes");
        modularDelphes->SetConfReader(confReader);
        modularDelphes->SetTreeWriter(treeWriter);

        factory = modularDelphes->GetFactory();
        allParticleOutputArray = modularDelphes->ExportArray("allParticles");
        stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
        partonOutputArray = modularDelphes->ExportArray("partons");

        modularDelphes->InitTask();
        //------------------------------------------------------------------------------------------
        //----- Initialize pythia -----
        pythia = new Pythia8::Pythia;        
        if(pythia == NULL)
        {
            throw runtime_error("can't create Pythia instance");
        }
        std::string sSeed = "0";
	float Mjj_cut = 150;
	if (argc >= 5)
	{
	    Mjj_cut = atoi(argv[4]);
	}
	int startEvent = 0, nEvent = -1;
        if (argc >= 6) 
        {
            startEvent = atoi(argv[5]);
            //sSeed = argv[5];
        }
        if (argc >= 7) 
        {
            nEvent = atoi(argv[6]);                                                     
        }
        //--- Initialize Les Houches Event File run. List initialization information.
        std::string sfile = "Beams:LHEF ="+inputFile;
        std::string sRandomSeed = "Random:seed = "+sSeed;
        //--- random seed from start event number
        pythia->readString("Random:setSeed = on");
	    pythia->readString("HadronLevel:Hadronize = on");
	    pythia->readString(sRandomSeed.c_str());        
	    pythia->readString("Beams:frameType = 4");
	    pythia->readString(sfile.c_str());

	    pythia->init();
	    //------------------------------------------------------------------------------------------
	    //----- initialize fast lhe file reader for preselection -----
	    ifstream inputLHE (inputFile.c_str(), ios::in);
	    int count=-1;
	    int skippedCounter = 0;
	    char buffer[256];
	    vector< vector<float> > LHE_event;
	    //--- start from specified event
	    while(count < startEvent)
	    {
		inputLHE.getline(buffer,256);
		while(strcmp(buffer, "<event>") != 0)
		{
		    inputLHE.getline(buffer,258);
		}
		count++;
	    }
	    if(pythia->LHAeventSkip(startEvent))
	    {
		cout << "### skipped first " << startEvent << " events" << endl;
	    }
	    //------------------------------------------------------------------------------------------
	    ExRootProgressBar progressBar(-1);
	    // Loop over all events
	    errorCounter = 0;
	    eventCounter = 0;
	    modularDelphes->Clear();
	    readStopWatch.Start();
	    while(!interrupted)
	    {
		if(eventCounter >= nEvent && nEvent != -1)
		{
		    break;
		}
		LHE_event.clear();  
		//--- end of file check
		if(strcmp(buffer, "</LesHouchesEvents>") == 0)        
		{
		    cout << "### end of LHE file reached! ### " << endl;
		    break;
		}//--- reads and store the event
		else if(strcmp(buffer, "<event>") == 0)
		{   
		    int nPart;
		    float info_tmp;
		    inputLHE >> nPart;
		    for(int i=0; i<5; i++)
		    {
			inputLHE >> info_tmp;
		    }
		    for(int i=0; i<nPart; i++)
		    {
			vector<float> LHE_particle;
			LHE_particle.clear();                
			for(int j=0; j<13; j++)
			{
			    inputLHE >> info_tmp;
			    LHE_particle.push_back(info_tmp);
			}
			LHE_event.push_back(LHE_particle);
		    }
		    //--- process only selected events
		    if(lhe_event_preselection(&LHE_event, Mjj_cut, branchGen))
		    {
			if(!pythia->next())
			{
			    //--- If failure because reached end of file then exit event loop
			    if (pythia->info.atEndOfFile())
			    {
				cerr << "Aborted since reached end of Les Houches Event File" << endl;
				break;
			    }
			    //--- keep trace of faulty events
			    errorCounter++;
			}
			readStopWatch.Stop();
			//--- delphes simulation fase
			procStopWatch.Start();
			ConvertInput(eventCounter, pythia, branchEvent, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray, &readStopWatch, &procStopWatch);
			modularDelphes->ProcessTask();
			procStopWatch.Stop();
			//--- filling the output tree
			treeWriter->Fill();
    
			//--- logistic 
			//treeWriter->Clear();
			modularDelphes->Clear();
			readStopWatch.Start();
		    }
		    else
		    {
			if(pythia->LHAeventSkip(1))
			{
			    skippedCounter++;
			}
			else
			{
			    cout << "### ERROR: couldn't skip event" << endl;
			}
		    }
		    eventCounter++;
		    progressBar.Update(eventCounter, eventCounter);
        	}
        	inputLHE.getline(buffer,256);
	    }	
	    progressBar.Update(eventCounter, eventCounter, kTRUE);
	    progressBar.Finish();

	    cout << "--------------------Statistics---------------------" << endl;
	    cout << "-#######  read events:        " << eventCounter << endl; 
	    cout << "-#######  failed events:      " << errorCounter << endl;
	    cout << "-#######  skipped events:     " << skippedCounter << endl;
	    cout << "---------------------------------------------------" << endl;
        
	    modularDelphes->FinishTask();
	    treeWriter->Write();
    
	    cout << endl <<  "** Exiting..." << endl;

	    delete pythia;
	    delete confReader;
	    delete outputFile;

	    return 0;
    }
    catch(runtime_error &e)
    {
        if(treeWriter) delete treeWriter;
        if(outputFile) delete outputFile;
        cerr << "** ERROR: " << e.what() << endl;
        return 1;
    }
}
