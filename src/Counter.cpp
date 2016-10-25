#include "Counter.h"

Count::Count(std::string datapath,double numberofevents)
{
    data_loc    = datapath;
    finalpr_loc = datapath + "finalpr.data";
    B_MULT_loc  = datapath + "B_MULT";
    NPOM_loc    = datapath + "NPOM.dat";
    number_of_events = numberofevents;
    #if NBNF
    InitializeNBNF();
    #endif
}

#if NBNF
void Count::InitializeNBNF()
{

    output = new TFile(NBNFFilename,"recreate");

    // Creates NFNB trees and histograms
    std::vector<const char*> branchNamesALL = {"all","nsd","etalim","ptcut"}; 
    std::vector<const char*> branchNamesDIV = {"all_div","nsd_div","etalim_div","ptcut_div"}; 
    std::vector<const char*> branchNamesNF  = {"all_nf","nsd_nf","etalim_nf","ptcut_nf"}; 
    HistNames.push_back("ALL");
    HistNames.push_back("Non-single diffraction");
    HistNames.push_back("|\\eta|<1");
    HistNames.push_back("0.3<p_{T}<1.5");

    ALL.push_back(new TH1F("all",HistNames[0],NBins,start,stop));
    ALL.push_back(new TH1F("nsd",HistNames[1],NBins,start,stop));
    ALL.push_back(new TH1F("etalim",HistNames[2],NBins,start,stop));
    ALL.push_back(new TH1F("ptcut",HistNames[3],NBins,start,stop));
    DIV.push_back(new TH1F("all_div",HistNames[0],NBins,start,stop));
    DIV.push_back(new TH1F("nsd_div",HistNames[1],NBins,start,stop));
    DIV.push_back(new TH1F("etalim_div",HistNames[2],NBins,start,stop));
    DIV.push_back(new TH1F("ptcut_div",HistNames[3],NBins,start,stop));
    NF. push_back(new TH1F("nf",HistNames[0],NBins,start,stop));
    NF. push_back(new TH1F("nsd_nf",HistNames[1],NBins,start,stop));
    NF. push_back(new TH1F("etalim_nf",HistNames[2],NBins,start,stop));
    NF. push_back(new TH1F("ptcut_nf",HistNames[3],NBins,start,stop));

    ALLTree = new TTree("H_all","all");
    DIVTree = new TTree("H_nsd","nsd");
    NFTree  = new TTree("H_nf","nf");
    
    for(int i=0 ; i<4 ; i++)
    {
        // Setting errorbar calculations
        ALL[i]->Sumw2(true);
        DIV[i]->Sumw2(true);
        NF [i]->Sumw2(true);

        ALLTree->Branch(branchNamesALL[i],ALL[i]);
        DIVTree->Branch(branchNamesDIV[i],DIV[i]);
        NFTree-> Branch(branchNamesNF [i],NF [i]);
    }
#if NPOM
    // Creates NPOM tree and histograms
    NPOMTree = new TTree("NPOM","npom");
    for (int i=0 ; i<25 ; i++)
    {
        NPOMSH.push_back(std::vector<TH1F*>());
        NPOMSH_nf.push_back(std::vector<TH1F*>());

        for (int j=0 ; j<25 ; j++)
        {
            char number_string [7];
            sprintf(number_string,"_%02d_%02d",i,j);
            number_string[6] = '\0';
            std::string temp  = "NPOM"+std::string(number_string);
            std::string temp1 = "npom"+std::string(number_string);
            std::string temp2 = "NPOM_NF"+std::string(number_string);
            std::string temp3 = "npom_nf"+std::string(number_string);

            NPOMSH[i].push_back(new TH1F(temp.c_str(),temp1.c_str(),NBins,start,stop));
            NPOMSH_nf[i].push_back(new TH1F(temp2.c_str(),temp3.c_str(),NBins,start,stop));

            // Setting errorbar calculations
            NPOMSH[i][j]->Sumw2(true);
            NPOMSH_nf[i][j]->Sumw2(true);

            NPOMTree->Branch(temp1.c_str(),NPOMSH[i][j]);
            NPOMTree->Branch(temp3.c_str(),NPOMSH_nf[i][j]);
        }
    }
#endif
}
#endif


void Count::ReadAndCount()
{
    timer.startTimer();

    #if NBNF
    int nf_nb [8] = {0,0,0,0,0,0,0,0};
    int counted [4] = {0,0,0,0};
    #endif
    #if NPOM
    int npoms, npomh;
    std::ifstream NPOMFile(NPOM_loc.c_str());
    #endif
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
    //testfil.open("900_bcorr.out");
    
    while (std::getline(B_MULTFile,B_MULTLine))
    {
        std::istringstream aa(B_MULTLine);
        aa >> EVENTNR >> temp >> PARTNR;
        #if NPOM
        std::getline(NPOMFile,NPOMLine);
        std::istringstream cc(NPOMLine);
        cc >> npoms >> npomh;
        #endif
        
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

            // Checks if there is nothing happening
            if (i==0 and IDIAG==4)
            {
                std::getline(finalprFile,finalprLine);
                break;
            }
            else
            {
                if (ICHJ != 0)
                {
                    const double p_abs      = std::sqrt(PXJ*PXJ + PYJ*PYJ + PZJ*PZJ);
                    const double p_T        = std::sqrt(PXJ*PXJ + PYJ*PYJ);
                    const double rap        = 0.5*std::log((EPAT+PZJ)/(EPAT-PZJ));
                    const double psrap      = 0.5*std::log((p_abs+PZJ)/(p_abs-PZJ));
                    const double psrap_abs  = std::abs(psrap);

                    #if NBNF
                    int nbnf_index = (psrap < 0);
                    nf_nb[nbnf_index] += 1;
                    counted[0] = 1;

                    // Non-single diffraction
                    if (IDIAG !=1 and IDIAG !=6 and IDIAG != 10)
                    {
                        nf_nb[nbnf_index+2] += 1;
                        counted[1] = 1;
                    }
                    if (psrap_abs < 1)
                    {
                        nf_nb[nbnf_index+4] += 1;
                        counted[2] = counted[3] = 1;
                        if (psrap_abs > 0.2 and psrap_abs < 0.8)
                        {
                            if (p_T > 0.3 and p_T < 1.5)
                                nf_nb[nbnf_index+6] += 1;
                        }
                    }
                    #endif
                #if bcorr
                if (psrap_abs < 1 and p_T > 0.05) 
                {
                    bcorr_count = 1;
                    if (p_T > 0.3 and p_T < 1.5)
                        BcorrCheck(EVENTNR,psrap);
                }

                #endif
                }
            }
        }

        #if NBNF
        // Counts for different conditions
        for(int i=0 ; i<4 ; i++)
        {
            if (counted[i])
            {
                ALL[i]->Fill(nf_nb[2*i],nf_nb[2*i+1]);
                DIV[i]->Fill(nf_nb[2*i],nf_nb[2*i+1]);
                NF [i]->Fill(nf_nb[2*i]);
            }
        }
        #if NPOM
        // Counts if abs(eta) < 1
        #if ptcut
        if (counted[3])
        {
            NPOMSH[npoms][npomh]->Fill(nf_nb[6],nf_nb[7]);
            NPOMSH_nf[npoms][npomh]->Fill(nf_nb[6]);
        }
        #elif nsd
        if (counted[1])
        {
            NPOMSH[npoms][npomh]->Fill(nf_nb[2],nf_nb[3]);
            NPOMSH_nf[npoms][npomh]->Fill(nf_nb[2]);
        }
        #endif
        
        #endif
        // Resets nf and nb counters
        for(int i=0 ; i<4 ; i++)
            nf_nb[2*i] = nf_nb[2*i+1] = counted[i] = 0;
        #endif
        #if bcorr
        if (bcorr_count)
        {
            Bcorrgap();
            for(int i=0 ; i<9 ; i++)
                std::fill(bne[i].begin(),bne[i].end(),0);
            bcorr_count = 0;
            bcorr_Nevents += 1;
        }
        #endif
    }
    #if NBNF
    ALLTree->Fill();
    DIVTree->Fill();
    NFTree ->Fill();
    NPOMTree->Fill();

    for(int i=0 ; i<4 ; i++)
        DIV[i]->Divide(NF[i]);

    output->Write();
    output->Close();
    #endif
    std::cout << std::endl;
    #if bcorr
	std::ofstream bcorr_file("/home/roar/master/qgsm_analysis_tool/ana/out/bcorr.csv");
    bcorr_file << bcorr_Nevents << std::endl;
    double b_corr = 0;
    //std::cout << bcorr_Nevents << std::endl;
    for (int i=0 ; i<13 ; i++)
    {
        b_corr = (eta_gaps[i][3] - eta_gaps[i][0]*eta_gaps[i][1]/(double)bcorr_Nevents)/
                 (eta_gaps[i][2] - eta_gaps[i][0]*eta_gaps[i][0]/(double)bcorr_Nevents);
        bcorr_file << b_corr << std::endl;
        //std::cout << "("<<eta_gaps[i][3] << " - " <<eta_gaps[i][0]*eta_gaps[i][1]/(double)bcorr_Nevents <<
        //    ")/("<<eta_gaps[i][2] << " - " <<eta_gaps[i][0]*eta_gaps[i][0]/(double)bcorr_Nevents<<")"<<std::endl;
        //std::cout << eta_gaps[i][0] << "  " << eta_gaps[i][1] << "  " << eta_gaps[i][2] << 
        //             "  " << eta_gaps[i][3] << std::endl;
        //std::cout << "="<<b_corr << std::endl<<std::endl;
    }
	bcorr_file.close();
    #endif
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
        eta_gaps[k][0] += temp_eta_gaps[k][0];
        eta_gaps[k][1] += temp_eta_gaps[k][1];
        eta_gaps[k][2] += temp_eta_gaps[k][0]*temp_eta_gaps[k][0];
        eta_gaps[k][3] += temp_eta_gaps[k][0]*temp_eta_gaps[k][1];
    }
} 
void Count::BcorrCheck(int EVENTNR,double eta)
{
    //testfil << EVENTNR << "," << eta<< std::endl;
    
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
#endif

void Count::Progress(int eventnr)
{
    timer.stopTimer();
    int *returnTime = new int[3];
    returnTime = timer.elapsedTimeClock();
    printf("\r %3d%% %02dh %02dm %02ds  ",(int)(eventnr/number_of_events*100),
                                       returnTime[0],returnTime[1],returnTime[2]);
}
