#pragma once

//#include "Eigen/Core"
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

#define NBNF true
#define NPOM true 
#define bcorr false

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

        const char* NBNFFilename = "test.root";
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
};


