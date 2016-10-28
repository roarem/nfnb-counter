#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include "Timer.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#define bcorr 0
#define NBNF 1 
    ///Dependent on NBNF///
    #define NBNFRegular 1
    #define NBNFSingle 1
    #define NBNFDouble 1
    #define NPOM 1  
        //Dependent on NPOM//
        #define NPOMptcut 1
        #define NPOMnsd 0

class Count 
{
    public:
        Count (const char* nbnfout, std::string datapath,double numberofevents);
        void ReadAndCount();
        void Progress(int eventnr);

        std::string data_loc;
        std::string finalpr_loc;
        std::string B_MULT_loc;
        std::string NPOM_loc;
        double number_of_events; 

        Timer timer;

        #if NBNF
        void InitializeNBNF();

        float NBins = 240;
        float start = -0.5;
        float stop  = NBins + start;
        TFile *output;
        const char* NBNFFilename;

        #if NBNFRegular
        std::vector<const char*> HistNamesReg;
        std::vector<TH1F*> ALLREG;
        std::vector<TH1F*> DIVREG;
        std::vector<TH1F*> NFREG;

        TTree *ALLREGTree; 
        TTree *DIVREGTree; 
        TTree *NFREGTree; 
        #endif

        #if NBNFSingle
        std::vector<const char*> HistNamesSin;
        std::vector<TH1F*> ALLSIN;
        std::vector<TH1F*> DIVSIN;
        std::vector<TH1F*> NFSIN;
        std::vector<TH1F*> NBSIN;

        TTree *ALLSINTree; 
        TTree *DIVSINTree; 
        TTree *NFSINTree; 
        TTree *NBSINTree; 
        #endif

        #if NBNFDouble
        std::vector<const char*> HistNamesDou;
        std::vector<TH1F*> ALLDOU;
        std::vector<TH1F*> DIVDOU;
        std::vector<TH1F*> NFDOU;
        std::vector<TH1F*> NBDOU;

        TTree *ALLDOUTree; 
        TTree *DIVDOUTree; 
        TTree *NFDOUTree; 
        TTree *NBDOUTree; 
        #endif

        #if NPOM
        std::vector<std::vector<TH1F*>> NPOMSH;
        std::vector<std::vector<TH1F*>> NPOMSH_nf;

        TTree *NPOMTree;
        #endif
        #endif

        #if bcorr
        void BcorrCheck(int EVENTNR, double eta);
        void Bcorrgap();

        int bcorr_count     = 0;
        double bcorr_Nevents   = 0;

        std::vector<std::vector<double>> eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(4,0));
        std::vector<std::vector<double>> temp_eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(2,0));
        std::vector<std::vector<double>> bne =
            std::vector<std::vector<double>>(9,std::vector<double>(2,0));
        #endif
};
