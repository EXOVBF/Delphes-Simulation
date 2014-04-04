#include <stdexcept>
#include <iostream>
#include <sstream>

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

bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut)
{
    vector<TLorentzVector> outQuark;
    int leptons=0;
    
    //--- search for final state quark in the event
    for(int iPart = 0; iPart < LHE_event->size(); iPart++)
    {
        vector<float> particle = LHE_event->at(iPart);
        if(particle.at(1) == 1 && abs(particle.at(0)) > 0 && abs(particle.at(0)) < 7)
        {
            TLorentzVector tmp4vect;
            tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));
            outQuark.push_back(tmp4vect);
        }
        if(particle.at(1) == 1 && (abs(particle.at(0)) == 11 || abs(particle.at(0)) == 13 || abs(particle.at(0)) == 15))
        {
            leptons++;
        }
    }
    // reject fully hadronic events
    if(leptons < 1)
    {
        return false;
    }
    // keep events with 1 or 2 jets
    if(outQuark.size() < 3)
    {
        return true;
    }
    //--- apply the VBF preselection cut on Mjj
    for(int iQuark = 0; iQuark < outQuark.size(); iQuark++)
    {
        TLorentzVector tmp4vect = outQuark.at(iQuark);
        for(int jQuark = iQuark+1; jQuark < outQuark.size(); jQuark++)
        {
            tmp4vect += outQuark.at(jQuark);
            if(tmp4vect.M() > Mjj_cut)
            {
                return true;
            }
        }
    }
    return false;
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
        cout << "- Usage: " << appName << " config_file" << " lhe_file" << " output_file" << " Mjj_cut" << " start" << " number" << endl;
        cout << "- config_file          ->  configuration file in Tcl format" << endl;
        cout << "- input_file           ->  lhe file for Pythia8" << endl;
        cout << "- output_file          ->  output file in ROOT format" << endl;
        cout << "- Mjj_cut  (optional)  ->  cut on Mjj in GeV -- default = 150 GeV" << endl;
        cout << "- start    (optional)  ->  number of starting event" << endl;
        cout << "- number   (optional)  ->  number of total events to be processed" << endl;
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

        treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
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
                if(lhe_event_preselection(&LHE_event, Mjj_cut))
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
        //delete modularDelphes;
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
