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
#define PI 3.14159265358979323846264338327950288

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
        void Sin_Dou(int nbnf_index,float psrap_abs,int IDIAG);
        void Non_sin_diff(int nbnf_index,float psrap_abs,int IDIAG);
        void eta_pt_cut(int nbnf_index,float psrap_abs,float p_T);
        void Filler(int npoms, int npomh);
        void Writer();

        float NBins = 600;
        float start = -0.5;
        float stop  = NBins + start;


        TFile *output;
        const char* NBNFFilename;
        std::vector<const char*> folders = {"ptcut","all","nsd","sin","dou","NBNFSin",
					    "NBNFDou","multi"};
        std::vector<std::string> prefix = {"ptcut","all","nsd","sin","dou"};
        std::vector<int> count_this  = std::vector<int>(prefix.size(),0);
	
        int nch = 0;
	std::vector<int> nf_nb		= std::vector<int>(10,0);
	std::vector<int> nf_nb_sin	= std::vector<int>(16,0); 
        std::vector<int> counted_sin	= std::vector<int>(8,0); 
        std::vector<int> nf_nb_dou  	= std::vector<int>(16,0); 
        std::vector<int> counted_dou	= std::vector<int>(8,0); 

	int folders_size     = (int)folders.size();
	int prefix_size      = (int)prefix.size();
	int count_this_size  = (int)count_this.size();
	int counted_sin_size = (int)counted_sin.size();
	int counted_dou_size = (int)counted_dou.size();
	int NBNFSIN_size     = 0;
	int NBNFDOU_size     = 0;
	int NPOM_NCH_size    = 0;
	int NPOM_NCHi_size   = 0;
	int NPOMSH_size	     = 0;
	int NPOMSHk_size     = 0;
	int NPOMSHki_size    = 0;


        TH1F* N_CH;

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
	void PhiCheck(double pxi, double pyi, double eta_abs, int FB);
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

	int eta_gap_index_size = 6;
	int phi_count_size     = 8;

	std::vector<std::vector<std::vector<double>>> phieta_temp =
	    std::vector<std::vector<std::vector<double>>>(2,
		    std::vector<std::vector<double>>(eta_gap_index_size,
			std::vector<double>(phi_count_size,0)));

	std::vector<std::vector<double>> phi_eta_bcorr = 
	    std::vector<std::vector<double>>(8,std::vector<double>(6,0));

        #endif//bcorr
};
