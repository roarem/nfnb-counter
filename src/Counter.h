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

        const char* NBNFFilename = "7TeV_4M.root";
        std::vector<const char*> HistNames;
        std::vector<TH1F*> ALL;
        std::vector<TH1F*> DIV;
        std::vector<TH1F*> NF;

        TTree *ALLTree; 
        TTree *DIVTree; 
        TTree *NFTree; 

        TFile *output;

        #if NPOM
        std::vector<std::vector<TH1F*>> NPOMS;
        std::vector<std::vector<TH1F*>> NPOMS_nf;
        //std::vector<std::vector<TH1F*>> NPOMH;
        //std::vector<std::vector<TH1F*>> NPOMH_nf;

        TTree *NPOMTree;
        #endif
        #endif
        #if bcorr
        void BcorrCheck(int eta);

        int bcorr_count     = 0;
        int bcorr_Nevents   = 0;

        std::vector<std::vector<int>> bcorr_nfnb = std::vector<std::vector<int>>(8,std::vector<int>(2,0));
        //std::vector<std::vector<int>> bcorr1 = std::vector<std::vector<int>>(1,std::vector<int>(2,0));
        //std::vector<std::vector<int>> bcorr2 = std::vector<std::vector<int>>(2,std::vector<int>(2,0));
        //std::vector<std::vector<int>> bcorr3 = std::vector<std::vector<int>>(3,std::vector<int>(2,0));
        //std::vector<std::vector<int>> bcorr7 = std::vector<std::vector<int>>(7,std::vector<int>(2,0));
        #endif
};


