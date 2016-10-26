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

#define NBNF 1 
#define NPOM 0  
#define NPOMptcut 0
#define NPOMnsd 0
#define bcorr 0

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
        float NBins = 240;
        float start = -0.5;
        float stop  = NBins + start;
        void InitializeNBNF();

        const char* NBNFFilename;
        std::vector<const char*> HistNames;
        std::vector<TH1F*> ALL;
        std::vector<TH1F*> DIV;
        std::vector<TH1F*> NF;

        TTree *ALLTree; 
        TTree *DIVTree; 
        TTree *NFTree; 

        TFile *output;

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


