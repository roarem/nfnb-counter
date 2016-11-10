#include "Counter.h"

Count::Count(const char* nbnfout, std::string datapath,double numberofevents)
#if NBNF
    : NBNFFilename(nbnfout)
#endif//NBNF
{
    data_loc            = datapath;
    finalpr_loc         = datapath + "finalpr.data";
    B_MULT_loc          = datapath + "B_MULT";
    NPOM_loc            = datapath + "NPOM.dat";
    number_of_events    = numberofevents;

    #if NBNF
    output = new TFile(NBNFFilename,"recreate");
    for (int i=0 ; i<folders.size() ; i++)
        output->mkdir(folders[i]);
    InitializeNBNF();
    #endif//NBNF
    #if bcorr
    output->mkdir(bcorr_folder);
    bcorr_initialize();
    #endif//bcorr
}

#if NBNF
void Count::InitializeNBNF()
{
    std::vector<const char*> HistNamesSin;
    HistNamesSin.push_back("ALL Single");
    HistNamesSin.push_back("|\\eta|<1 Single");
    HistNamesSin.push_back("0.3<p_{T}<1.5 Single");

    NBNFSIN.push_back(new TH1F("allSinFull",HistNamesSin[0],NBins,start,stop));
    NFSIN. push_back(new TH1F("nfFullSin",HistNamesSin[0],NBins,start,stop));
    NBSIN. push_back(new TH1F("nbFullSin",HistNamesSin[0],NBins,start,stop));
    NBNFSIN.push_back(new TH1F("allEta2Sin",HistNamesSin[0],NBins,start,stop));
    NFSIN. push_back(new TH1F("nfEta2Sin",HistNamesSin[0],NBins,start,stop));
    NBSIN. push_back(new TH1F("nbEta2Sin",HistNamesSin[0],NBins,start,stop));
    NBNFSIN.push_back(new TH1F("allEta1Sin",HistNamesSin[0],NBins,start,stop));
    NFSIN. push_back(new TH1F("nfEta1Sin",HistNamesSin[0],NBins,start,stop));
    NBSIN. push_back(new TH1F("nbEta1Sin",HistNamesSin[0],NBins,start,stop));
    NBNFSIN.push_back(new TH1F("allEta05Sin",HistNamesSin[0],NBins,start,stop));
    NFSIN. push_back(new TH1F("nfEta05Sin",HistNamesSin[0],NBins,start,stop));
    NBSIN. push_back(new TH1F("nbEta05Sin",HistNamesSin[0],NBins,start,stop));
    
    for(int i=0 ; i<4 ; i++)
    {
        // Setting errorbar calculations
        NBNFSIN[i]->Sumw2(true);
        NFSIN[i]->Sumw2(true);
        NBSIN[i]->Sumw2(true);
    }

    std::vector<const char*> HistNamesDou;
    HistNamesDou.push_back("ALL Double");
    HistNamesDou.push_back("|\\eta|<1 Double");
    HistNamesDou.push_back("0.3<p_{T}<1.5 Double");
    
    NBNFDOU.push_back(new TH1F("allDouFull",HistNamesDou[0],NBins,start,stop));
    NFDOU. push_back(new TH1F("nfFullDou",HistNamesDou[0],NBins,start,stop));
    NBDOU. push_back(new TH1F("nbFullDou",HistNamesDou[0],NBins,start,stop));
    NBNFDOU.push_back(new TH1F("allEta2Dou",HistNamesDou[0],NBins,start,stop));
    NFDOU. push_back(new TH1F("nfEta2Dou",HistNamesDou[0],NBins,start,stop));
    NBDOU. push_back(new TH1F("nbEta2Dou",HistNamesDou[0],NBins,start,stop));
    NBNFDOU.push_back(new TH1F("allEta1Dou",HistNamesDou[0],NBins,start,stop));
    NFDOU. push_back(new TH1F("nfEta1Dou",HistNamesDou[0],NBins,start,stop));
    NBDOU. push_back(new TH1F("nbEta1Dou",HistNamesDou[0],NBins,start,stop));
    NBNFDOU.push_back(new TH1F("allEta05Dou",HistNamesDou[0],NBins,start,stop));
    NFDOU. push_back(new TH1F("nfEta05Dou",HistNamesDou[0],NBins,start,stop));
    NBDOU. push_back(new TH1F("nbEta05Dou",HistNamesDou[0],NBins,start,stop));
    
    for(int i=0 ; i<4 ; i++)
    {
        // Setting errorbar calculations
        NBNFDOU[i]->Sumw2(true);
        NFDOU[i]->Sumw2(true);
        NBDOU[i]->Sumw2(true);
    }

    for (int i=0 ; i<prefix.size() ; i++)
    {
        NPOMSH.push_back(std::vector<std::vector<TH1F*>>());
        NPOMSH_nf.push_back(std::vector<std::vector<TH1F*>>());
        NPOMSH_nb.push_back(std::vector<std::vector<TH1F*>>());
        for (int j=0 ; j<25 ; j++)
        {
            NPOMSH[i].push_back(std::vector<TH1F*>());
            NPOMSH_nf[i].push_back(std::vector<TH1F*>());
            NPOMSH_nb[i].push_back(std::vector<TH1F*>());

            for (int k=0 ; k<25 ; k++)
            {
                char number_string [7];
                sprintf(number_string,"_%02d_%02d",j,k);
                number_string[6] = '\0';
                std::string temp  = prefix[i]+"_NPOM"+std::string(number_string);
                std::string temp1 = prefix[i]+"_npom"+std::string(number_string);
                std::string temp2 = prefix[i]+"_NPOM_NF"+std::string(number_string);
                std::string temp3 = prefix[i]+"_npom_nf"+std::string(number_string);
                std::string temp4 = prefix[i]+"_NPOM_NB"+std::string(number_string);
                std::string temp5 = prefix[i]+"_npom_nb"+std::string(number_string);

                NPOMSH[i][j].push_back(new TH1F(temp.c_str(),temp1.c_str(),NBins,start,stop));
                NPOMSH_nf[i][j].push_back(new TH1F(temp2.c_str(),temp3.c_str(),NBins,start,stop));
                NPOMSH_nb[i][j].push_back(new TH1F(temp4.c_str(),temp5.c_str(),NBins,start,stop));

                // Setting errorbar calculations
                NPOMSH[i][j][k]->Sumw2(true);
                NPOMSH_nf[i][j][k]->Sumw2(true);
                NPOMSH_nb[i][j][k]->Sumw2(true);
            }
        }
    }
}
#endif//NBNF
#if bcorr
void Count::bcorr_initialize()
{
    for (int i=0 ; i<13 ; i++)
    {
        bcorr_hists.push_back(std::vector<TH1F*>());
        for (int j=0 ; j<4 ; j++)
        {
            char number_string [6];
            sprintf(number_string,"_%d_%d",i,j);
            number_string[5] = '\0';
            std::string temp = "BCORR_"+std::string(number_string);
            std::string temp2 = "bcorr_"+std::string(number_string);
            bcorr_hists[i].push_back(new TH1F(temp.c_str(),temp2.c_str(),NBins,start,stop));
            bcorr_hists[i][j]->Sumw2(true);
        }
    }
}
#endif//bcorr

void Count::ReadAndCount()
{
    timer.startTimer();

    #if NBNF
    int npoms, npomh;
    std::ifstream NPOMFile(NPOM_loc.c_str());
    #endif//NBNF

    //B_MULT columns
    int   EVENTNR;  // Event number
    int   PARTNR;   // Particles in event
    float temp;     // unknown

    //finalpr columns
    float FREEZJ;   // Time event is completed
    float XXJ;      // position x at freeze time
    float YYJ;      // position y at freeze time
    float ZZJ;      // position z at freeze time
    float EPAT;     // energy at freeze time
    float PXJ;      // momentum p_x at freeze time
    float PYJ;      // momentum p_y at freeze time
    float PZJ;      // momentum p_z at freeze time
    float AMJ;      // mass
    float IDENT;    // identity of particle
    float IDIAG;    // Type of event
    float IBJ;      // baryon number
    float ISJ;      // strangeness
    float ICHJ;     // charge
    float TFORMJ;   // formation time
    float XXJI;     // formatiom position x
    float YYJI;     // formation position y
    float ZZJI;     // formation position z
    float IORIGJ;   // origin of particle
    float TFORMRJUK;// unknown
    
    std::string finalprLine, B_MULTLine, NPOMLine;
    std::ifstream finalprFile(finalpr_loc.c_str());
    std::ifstream B_MULTFile(B_MULT_loc.c_str());
    while (std::getline(B_MULTFile,B_MULTLine))
    {
        std::istringstream aa(B_MULTLine);
        aa >> EVENTNR >> temp >> PARTNR;

        #if NBNF
        std::getline(NPOMFile,NPOMLine);
        std::istringstream cc(NPOMLine);
        cc >> npoms >> npomh;
        #endif//NBNF
        
        if(EVENTNR%100==0)
            Progress(EVENTNR);
        //if (EVENTNR%100!=0) 
        //    continue;
        
        for (int i=0 ; i<PARTNR ; i++)
        {
            std::getline(finalprFile,finalprLine);
            std::istringstream bb(finalprLine);
            bb >> FREEZJ >> XXJ >> YYJ >> ZZJ >> EPAT >> PXJ >> PYJ >> PZJ >> AMJ 
               >> IDENT >> IDIAG >> IBJ >> ISJ >> ICHJ >> TFORMJ >> XXJI >> YYJI 
               >> ZZJI >> IORIGJ >> TFORMRJUK;

            // Checks if there is anything happening
            if (i==0 and IDIAG==4)
            {
                std::getline(finalprFile,finalprLine);
                break;
            }
            else
            {
                #if NBNF
                const double p_abs      = std::sqrt(PXJ*PXJ + PYJ*PYJ + PZJ*PZJ);
                const double p_T        = std::sqrt(PXJ*PXJ + PYJ*PYJ);
                const double rap        = 0.5*std::log((EPAT+PZJ)/(EPAT-PZJ));
                const double psrap      = 0.5*std::log((p_abs+PZJ)/(p_abs-PZJ));
                const double psrap_abs  = std::abs(psrap);
                const int nbnf_index = (psrap<0);

                Sin_Dou(nbnf_index,psrap_abs,IDIAG,ICHJ);
                Non_sin_diff(nbnf_index,psrap_abs,IDIAG,ICHJ);
                eta_pt_cut(nbnf_index,psrap_abs,p_T,ICHJ);

                nf_nb[nbnf_index+2] += 1;
                if (ICHJ != 0)
                    count_this[1] = 1;

                #endif//NBNF

                #if bcorr
                if (psrap_abs < 1 and p_T > 0.05) 
                {
                    if (ICHJ != 0)
                        bcorr_count = 1;
                    if (p_T > 0.3 and p_T < 1.5)
                        BcorrCheck(EVENTNR,psrap);
                }
                #endif//bcorr
            }
        }

        #if NBNF
        Filler(npoms,npomh);
        // Resets nf and nb counters
        for (int i=0 ; i<count_this.size() ; i++)
            nf_nb[2*i] = nf_nb[2*i+1] = count_this[i] = 0;
        for(int i=0 ; i<4 ; i++)
            nf_nb_sin[2*i] = nf_nb_sin[2*i+1] = counted_sin[i] = 0;
        for(int i=0 ; i<4 ; i++)
            nf_nb_dou[2*i] = nf_nb_dou[2*i+1] = counted_dou[i] = 0;

        #endif//NBNF

        #if bcorr
        if (bcorr_count)
        {
            Bcorrgap();
            for(int i=0 ; i<9 ; i++)
                std::fill(bne[i].begin(),bne[i].end(),0);
            bcorr_count = 0;
            bcorr_Nevents += 1;
        }
        #endif//bcorr
    }
    #if NBNF
    Writer();
    #endif//NBNF
    std::cout << std::endl;
    #if bcorr
    output->cd(bcorr_folder);
    for (int i=0 ; i<13 ; i++)
    {
        for (int j=0 ; j<4 ; j++)
        {
            bcorr_hists[i][j]->Write();
        }
    }
	std::ofstream bcorr_file("/home/roar/master/qgsm_analysis_tool/ana/out/7000_4M_bcorr.csv");
    bcorr_file << bcorr_Nevents << std::endl;
    double b_corr = 0;
    for (int i=0 ; i<13 ; i++)
    {
        b_corr = (eta_gaps[i][3] - eta_gaps[i][0]*eta_gaps[i][1]/(double)bcorr_Nevents)/
                 (eta_gaps[i][2] - eta_gaps[i][0]*eta_gaps[i][0]/(double)bcorr_Nevents);
        bcorr_file << b_corr << std::endl;
    }
	bcorr_file.close();
    #endif//bcorr
    output->Close();
}

#if bcorr
void Count::Bcorrgap()
{

    for(int i=0 ; i<7 ; i++)
    {
        temp_eta_gaps[i][0] = bne[i][0]+bne[i+1][0];
        temp_eta_gaps[i][1] = bne[i][1]+bne[i+1][1];
    }
    int j = 0;
    for(int i=7 ; i<10 ; i++)
    {
        temp_eta_gaps[i][0]  = bne[j][0]+bne[j+1][0]+bne[j+2][0]+bne[j+3][0];
        temp_eta_gaps[i][1]  = bne[j][1]+bne[j+1][1]+bne[j+2][1]+bne[j+3][1];
        j+=2;
    }
    j = 0;
    for(int i=10 ; i<12 ; i++)
    {
        temp_eta_gaps[i][0] = bne[j][0]+bne[j+1][0]+bne[j+2][0]+bne[j+3][0]+
                              bne[j+4][0]+bne[j+5][0];
        temp_eta_gaps[i][1] = bne[j][1]+bne[j+1][1]+bne[j+2][1]+bne[j+3][1]+
                              bne[j+4][1]+bne[j+5][1];
        j+=2;

    }

    temp_eta_gaps[12][0] = bne[0][0]+bne[1][0]+bne[2][0]+bne[3][0]+
                           bne[4][0]+bne[5][0]+bne[6][0]+bne[7][0];
    temp_eta_gaps[12][1] = bne[0][1]+bne[1][1]+bne[2][1]+bne[3][1]+
                           bne[4][1]+bne[5][1]+bne[6][1]+bne[7][1];

    for (int k=0 ; k<13 ; k++)
    {
        bcorr_hists[k][0]->Fill(temp_eta_gaps[k][0]);
        bcorr_hists[k][1]->Fill(temp_eta_gaps[k][1]);
        bcorr_hists[k][2]->Fill(temp_eta_gaps[k][0]*temp_eta_gaps[k][0]);
        bcorr_hists[k][3]->Fill(temp_eta_gaps[k][0]*temp_eta_gaps[k][1]);
        eta_gaps[k][0] += temp_eta_gaps[k][0];
        eta_gaps[k][1] += temp_eta_gaps[k][1];
        eta_gaps[k][2] += temp_eta_gaps[k][0]*temp_eta_gaps[k][0];
        eta_gaps[k][3] += temp_eta_gaps[k][0]*temp_eta_gaps[k][1];
    }
} 

void Count::BcorrCheck(int EVENTNR,double eta)
{
    int nfnbi = (eta<0);
    double eta10 = std::abs(eta)*10;
    for(int i=8 ; i>-1 ; i--)
    {
        if(eta10>=i)
        {
            bne[i][nfnbi] += 1;
            break;
        }
    }
}
#endif//bcorr
#if NBNF
void Count::eta_pt_cut(int nbnf_index, float psrap_abs, float p_T, int ICHJ)
{
    if (psrap_abs < 1 and ICHJ != 0)
    {
        count_this[0] = 1;
        if (psrap_abs > 0.2 and psrap_abs < 0.8)
        {
            if (p_T > 0.3 and p_T < 1.5)
                nf_nb[nbnf_index] += 1;
        }
    }
}
void Count::Non_sin_diff(int nbnf_index,float psrap_abs,int IDIAG,int ICHJ)
{
    if (IDIAG !=1 and IDIAG !=6 and IDIAG != 10)
    {
        nf_nb[nbnf_index+4] += 1;
        if (ICHJ != 0)
        {
            count_this[2] = 1;
        }
    }
}

void Count::Sin_Dou(int nbnf_index,float psrap_abs,int IDIAG,int ICHJ)
{
    int i_start = 0;
    if(psrap_abs<0.5)
        i_start = 0;
    else if(psrap_abs<1)
        i_start = 1;
    else if(psrap_abs<2)
        i_start = 2;
    else
        i_start = 3;
    
    if (IDIAG==1 or IDIAG==6 or IDIAG==10)
    {
        for(int i=i_start ; i<4 ; i++)
        {
            nf_nb_sin[2*i+nbnf_index] += 1;
            counted_sin[i] = 1;
        }
        if (ICHJ != 0)
            count_this[3] = 1;
    }
    
    if(IDIAG==11)
    {
        for(int i=i_start ; i<4 ; i++)
        {
            nf_nb_dou[2*i+nbnf_index] += 1;
            counted_dou[i] = 1;
        }
        if (ICHJ!=0)
            count_this[4] = 1;
    }
}

void Count::Filler(int npoms,int npomh)
{
    for(int i=0 ; i<4 ; i++)
    {
        if (counted_sin[i])
        {
            NBNFSIN[i]->Fill(nf_nb_sin[2*i],nf_nb_sin[2*i+1]);
            NFSIN[i]->Fill(nf_nb_sin[2*i]);
            NBSIN[i]->Fill(nf_nb_sin[2*i+1]);
        }
    }
    
    for(int i=0 ; i<4 ; i++)
    {
        if (counted_dou[i])
        {
            NBNFDOU[i]->Fill(nf_nb_dou[2*i],nf_nb_dou[2*i+1]);
            NFDOU[i]->Fill(nf_nb_dou[2*i]);
            NBDOU[i]->Fill(nf_nb_dou[2*i+1]);
        }
    }
    
    for (int i=0 ; i<count_this.size() ; i++)
    {
        if (count_this[i])
        {
            Nevents[i] += 1;
            NPOMSH[i][npoms][npomh]->Fill(nf_nb[2*i],nf_nb[2*i+1]);
            NPOMSH_nf[i][npoms][npomh]->Fill(nf_nb[2*i]);
            NPOMSH_nb[i][npoms][npomh]->Fill(nf_nb[2*i+1]);
        }
    }
}

void Count::Writer()
{
    output->cd(folders[5]);
    for(int i=0 ; i<4 ; i++)
    {
        NBNFSIN[i]->Write();
        NFSIN [i]->Write();
        NBSIN [i]->Write();
    }

    output->cd(folders[6]);
    for(int i=0 ; i<4 ; i++)
    {
        NBNFDOU[i]->Write();
        NFDOU [i]->Write();
        NBDOU [i]->Write();
    }

    for(int k=0 ; k<prefix.size() ; k++)
    {
        output->cd(folders[k]);
        for(int i=0 ; i<25 ; i++)
        {
            for(int j=0 ; j<25 ; j++)
            {
                NPOMSH[k][i][j]->Write();
                NPOMSH_nf[k][i][j]->Write();
                NPOMSH_nb[k][i][j]->Write();
            }
        }
    }
}

#endif//NBNF

void Count::Progress(int eventnr)
{
    timer.stopTimer();
    int *returnTime = new int[3];
    returnTime = timer.elapsedTimeClock();
    printf("\r %3d%% %02dh %02dm %02ds  ",(int)(eventnr/number_of_events*100),
                                       returnTime[0],returnTime[1],returnTime[2]);
}
