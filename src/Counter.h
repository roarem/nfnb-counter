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
  #define bcorr 1

class Count 
{
    public:
        Count (const char* nbnfout, const char* bcorrout, std::string datapath,double numberofevents);
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
        void Sin_Dou(int nbnf_index,float psrap_abs,int IDIAG,int ICHJ);
        void Non_sin_diff(int nbnf_index,float psrap_abs,int IDIAG,int ICHJ);
        void eta_pt_cut(int nbnf_index,float psrap_abs,float p_T, int ICHJ);
        void Filler(int npoms, int npomh);
        void Writer();

        float NBins = 600;
        float start = -0.5;
        float stop  = NBins + start;

        int nch = 0;
        int nf_nb [10] = {0,0,0,0,0,0,0,0,0,0};
        int nf_nb_sin [8] = {0,0,0,0,0,0,0,0};
        int counted_sin [4] = {0,0,0,0};
        int nf_nb_dou [8] = {0,0,0,0,0,0,0,0};
        int counted_dou [4] = {0,0,0,0};

        TFile *output;
        const char* NBNFFilename;
        std::vector<const char*> folders = {"ptcut","all","nsd","sin","dou","NBNFSin",
					    "NBNFDou","multi"};
        std::vector<std::string> prefix = {"ptcut","all","nsd","sin","dou"};
        std::vector<int> count_this  = std::vector<int>(prefix.size(),0);

        TH1F* N_CH;

        std::vector<TH1F*> NBNFREG;
        std::vector<TH1F*> NFREG;
        std::vector<TH1F*> NBREG;

        std::vector<TH1F*> NBNFSIN;
        std::vector<TH1F*> NFSIN;
        std::vector<TH1F*> NBSIN;

        std::vector<TH1F*> NBNFDOU;
        std::vector<TH1F*> NFDOU;
        std::vector<TH1F*> NBDOU;

        std::vector<std::vector<TH1F*>> NPOM_NCH;
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH;
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH_nf;
        std::vector<std::vector<std::vector<TH1F*>>> NPOMSH_nb;
        #endif//NBNF

        #if bcorr
        void bcorr_initialize();
        void BcorrCheck(int EVENTNR, double eta);
        void Bcorrgap();

        int bcorr_count         = 0;
        double bcorr_Nevents    = 0;
        const char* bcorrFilename;

        std::vector<std::vector<double>> eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(4,0));
        std::vector<std::vector<double>> temp_eta_gaps = 
            std::vector<std::vector<double>>(13,std::vector<double>(2,0));
        std::vector<std::vector<double>> bne =
            std::vector<std::vector<double>>(9,std::vector<double>(2,0));
        #endif//bcorr
};
