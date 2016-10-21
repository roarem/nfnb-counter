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

#define NBNF 0 
#define NPOM 0  
#define ptcut 0
#define nsd 0
#define bcorr 1

class Count 
{
    public:
        Count (std::string datapath,double numberofevents);
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

        //const char* NBNFFilename = "7TeV_4M.root";
        const char* NBNFFilename = "7TeV_4M_nsd.root";
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
        void BcorrCheck(double eta);

        int bcorr_count     = 0;
        int bcorr_Nevents   = 0;

        std::vector<std::vector<int>> bcorr_nfnb_event =
            std::vector<std::vector<int>>(9,std::vector<int>(2,0));
        std::vector<std::vector<int>> bcorr_nfnb = 
            std::vector<std::vector<int>>(8,std::vector<int>(4,0));
        #endif
};


