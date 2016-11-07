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

        float NBins = 600;
        float start = -0.5;
        float stop  = NBins + start;
        TFile *output;
        const char* NBNFFilename;
        std::vector<const char*> folders = {"ptcut","all","nsd","sin","dou","NBNFSin","NBNFDou"};
        std::vector<std::string> prefix = {"ptcut","all","nsd","sin","dou"};
        std::vector<int> count_this  = std::vector<int>(prefix.size(),0);

        #if NBNFRegular
        //folders.push_back("NBNFRegular");
        //const char* REGFolder = "NBNFRegular";
        std::vector<TH1F*> NBNFREG;
        std::vector<TH1F*> NFREG;
        std::vector<TH1F*> NBREG;
        #endif

        #if NBNFSingle
        //folders.push_back("NBNFSingle");
        //const char* SINFolder = "NBNFSingle";
        std::vector<TH1F*> NBNFSIN;
        std::vector<TH1F*> NFSIN;
        std::vector<TH1F*> NBSIN;
        #endif

        #if NBNFDouble
        //folders.push_back("NBNFDouble");
        //const char* DOUFolder = "NBNFDouble";
        std::vector<TH1F*> NBNFDOU;
        std::vector<TH1F*> NFDOU;
        std::vector<TH1F*> NBDOU;
        #endif

        #if NPOM
        //std::vector<const char*> NPOMFolders = 
        //    std::vector<const char*>{"NPOMptcut","NPOMall","NPOMnsd","NPOMdou","NPOMsin"};
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH;
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH_nf;
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH_nb;
        #endif
        #endif

        #if bcorr
        void BcorrCheck(int EVENTNR, double eta);
        void Bcorrgap();

        int bcorr_count         = 0;
        double bcorr_Nevents    = 0;

        std::vector<std::vector<double>> eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(4,0));
        std::vector<std::vector<double>> temp_eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(2,0));
        std::vector<std::vector<double>> bne =
            std::vector<std::vector<double>>(9,std::vector<double>(2,0));
        #endif
};
