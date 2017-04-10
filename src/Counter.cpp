#include "Counter.h"

Count::Count(const char* nbnfout,const char* bcorrout, std::string datapath,double numberofevents)
#if NBNF
    : NBNFFilename(nbnfout)
#endif//NBNF
#if bcorr
    ,bcorrFilename(bcorrout)
#endif//bcorr
{
    data_loc            = datapath;
    finalpr_loc         = datapath + "finalpr.data";
    B_MULT_loc          = datapath + "B_MULT";
    NPOM_loc            = datapath + "NPOM.dat";
    number_of_events    = numberofevents;

    output = new TFile(NBNFFilename,"recreate");
    #if NBNF
    for (int i=0 ; i<folders_size ; i++)
        output->mkdir(folders[i]);
    InitializeNBNF();
    #endif//NBNF
}

#if NBNF
void Count::InitializeNBNF()
{
    
    std::string HistTitleSin = "single diffraction";
    std::string HistTitleDou = "double diffraction";
    for (int i=1 ; i<4 ; i++)
    {
        char number_string [4];
    	sprintf(number_string,"%02d",i);
    	number_string[3] = '\0';
    	std::string temp0  = "SIN_NBNF_"+std::string(number_string);
    	std::string temp1  = "SIN_NF_"+std::string(number_string);
    	std::string temp2  = "SIN_NB_"+std::string(number_string);
        NBNFSIN.push_back(new TH1F(temp0.c_str(),HistTitleSin.c_str(),NBins,start,stop));
        NFSIN.push_back(new TH1F(temp1.c_str(),HistTitleSin.c_str(),NBins,start,stop));
        NBSIN.push_back(new TH1F(temp2.c_str(),HistTitleSin.c_str(),NBins,start,stop));

        // Setting errorbar calculations
        NBNFSIN[i-1]->Sumw2(true);
        NFSIN[i-1]->Sumw2(true);
        NBSIN[i-1]->Sumw2(true);

    }

    std::string temp0   = "DOU_NBNF";
    std::string temp1   = "DOU_NF_";
    std::string temp2   = "DOU_NB_";
    NBNFDOU = new TH1F(temp0.c_str(),HistTitleDou.c_str(),NBins,start,stop);
    NFDOU   = new TH1F(temp1.c_str(),HistTitleDou.c_str(),NBins,start,stop);
    NBDOU   = new TH1F(temp2.c_str(),HistTitleDou.c_str(),NBins,start,stop);
    
    NBNFDOU->Sumw2(true);
    NFDOU->Sumw2(true);
    NBDOU->Sumw2(true);

    for (int i=0 ; i<prefix_size ; i++)
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

    N_CH = new TH1F("nch","NCH",NBins,start,stop);
    for (int i=0 ; i<25 ; i++)
    {
        NPOM_NCH.push_back(std::vector<TH1F*>());
        for (int j=0 ; j<25 ; j++)
        {
            char number_string [7];
            sprintf(number_string,"_%02d_%02d",i,j);
            number_string[6] = '\0';
            std::string temp  = "multi_NPOM"+std::string(number_string);
            std::string temp1 = "multi_npom"+std::string(number_string);
            NPOM_NCH[i].push_back(new TH1F(temp.c_str(),temp1.c_str(),NBins,start,stop));
            NPOM_NCH[i][j]->Sumw2(true);
        }
    }
    
   NBNFSIN_size	    = (int)NBNFSIN.size(); 
   //NBNFDOU_size     = (int)NBNFDOU.size(); 
   NPOM_NCH_size    = (int)NPOM_NCH.size();
   NPOM_NCHi_size   = (int)NPOM_NCH[0].size();
   NPOMSH_size	    = (int)NPOMSH.size();
   NPOMSHk_size	    = (int)NPOMSH[0].size();
   NPOMSHki_size    = (int)NPOMSH[0][0].size();

}
#endif//NBNF

void Count::ReadAndCount()
{
    std::cout << "using data from " << data_loc << " and storing in " << NBNFFilename << std::endl;

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
    float FREEZJ   = 0;// Time event is completed
    float XXJ      = 0;// position x at freeze time
    float YYJ      = 0;// position y at freeze time
    float ZZJ      = 0;// position z at freeze time
    float EPAT     = 0;// energy at freeze time
    float PXJ      = 0;// momentum p_x at freeze time
    float PYJ      = 0;// momentum p_y at freeze time
    float PZJ      = 0;// momentum p_z at freeze time
    float AMJ      = 0;// mass
    float IDENT    = 0;// identity of particle
    float IDIAG    = 0;// Type of event
    float IBJ      = 0;// baryon number
    float ISJ      = 0;// strangeness
    float ICHJ     = 0;// charge
    float TFORMJ   = 0;// formation time
    float XXJI     = 0;// formatiom position x
    float YYJI     = 0;// formation position y
    float ZZJI     = 0;// formation position z
    float IORIGJ   = 0;// origin of particle
    float TFORMRJUK= 0;// unknown
    
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
                const double p_abs      = std::sqrt(PXJ*PXJ + PYJ*PYJ + PZJ*PZJ);
                const double p_T        = std::sqrt(PXJ*PXJ + PYJ*PYJ);
                const double rap        = 0.5*std::log((EPAT+PZJ)/(EPAT-PZJ));
                const double psrap      = 0.5*std::log((p_abs+PZJ)/(p_abs-PZJ));
                const double psrap_abs  = std::abs(psrap);
                const int nbnf_index    = (psrap<0);

                #if NBNF
		if(i==0)
		{
		    if (IDIAG != 1 and IDIAG != 6 and IDIAG != 10 and IDIAG != 4)
		    {
			count_this[2] = 1; //nsd
		    }
		}

		if(ICHJ!=0)
		{
		    Sin_Dou(nbnf_index,psrap,IDIAG);
            Non_sin_diff(nbnf_index,psrap_abs,IDIAG);
            eta_pt_cut(nbnf_index,psrap_abs,p_T);
		    PhiCheck(PXJ,PYJ,psrap_abs,nbnf_index);
		}

        if (ICHJ != 0)
		{
		    nf_nb[nbnf_index+2] += 1;
                    count_this[1] = 1;   //all charged
		}
                #endif//NBNF

                #if bcorr
                if (psrap_abs < 1 and p_T > 0.05) 
                {
                    if (ICHJ != 0)
		            {
                        bcorr_count = 1;
			            if (p_T > 0.3 and p_T < 1.5)
                            BcorrCheck(EVENTNR,psrap);
		            }
                }
                #endif//bcorr
            }
        }

        #if NBNF
        Filler(npoms,npomh);
        // Resets nf and nb counters
        nch = 0;
        for (int i=0 ; i<count_this_size ; i++)
            nf_nb[2*i] = nf_nb[2*i+1] = count_this[i] = 0;

        counted_sin1 = counted_sin6 = counted_sin10 = counted_dou = 0;
        for(int i=0 ; i<counted_sin_size1 ; i++)
            nf_nb_sin1[i] = nf_nb_sin6[i] = nf_nb_sin10[i] = 0;
        for(int i=0 ; i<counted_dou_size ; i++)
            nf_nb_dou[i] = 0;
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
    std::cout << std::endl;
    #if NBNF
    std::cout << "writing nbnf to root file" << std::endl;
    Writer();
    #endif//NBNF
    #if bcorr
    std::cout << "writing bcorr to file" << std::endl;
    std::ofstream bcorr_file(bcorrFilename);
    bcorr_file << bcorr_Nevents << std::endl;
    double b_corr = 0;
    for (int i=0 ; i<13 ; i++)
    {
        b_corr = (eta_gaps[i][3] - eta_gaps[i][0]*eta_gaps[i][1]/(double)bcorr_Nevents)/
                 (eta_gaps[i][2] - eta_gaps[i][0]*eta_gaps[i][0]/(double)bcorr_Nevents);
        bcorr_file << b_corr << std::endl;
	//std::cout << "gap: " << i << "  nfnb: " <<eta_gaps[i][3] << " nf**2: " << eta_gaps[i][2] << 
	//    " nf: " << eta_gaps[i][0] << " nb: " << eta_gaps[i][1] << std::endl;
    }
    bcorr_file.close();
    
    #endif//bcorr
    std::cout << "closing root file" << std::endl;
    output->Close();
    //std::cout << "destroying all new" << std::endl;
    //Destroy();
    std::cout << "done and done, bye!" << std::endl;
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
void Count::PhiCheck(double pxi, double pyi, double eta_abs, int FB)
{
    if (eta_abs<1.2)
    {
	int eta_gap_index;
	double pi_8 = PI/8.0;
	double eta10 = eta_abs*10;

	for (int i=eta_gap_index_size-1 ; i>-1 ; i--)
	{
	    if (eta10>i)
	    {
		eta_gap_index=i;
		break;
	    }
	}

	double phi_abs	= std::atan(std::abs(pyi)/std::abs(pxi));
	double temp_phi = 0;

	if (pxi>0 and pyi>0)
	    temp_phi = phi_abs;
	else if (pxi<0 and pyi>0)
	    temp_phi = PI-phi_abs;
	else if (pxi<0 and pyi<0)
	    temp_phi = PI+phi_abs;
	else if (pxi>0 and pyi<0 and phi_abs>pi_8)
	    temp_phi = 2*PI-phi_abs;
	else if (pxi>0 and pyi<0 and phi_abs<pi_8)
	    temp_phi = -phi_abs;
	
	double phi_counter = pi_8;
	for (int i=0 ; i<phi_count_size ; i++)
	{
	    if (temp_phi<phi_counter) 
	    {
		phieta_temp[FB][eta_gap_index][i] += 1;
		break;
	    }
	    phi_counter += 2*pi_8;
	}
    }
}
#endif//bcorr
#if NBNF
void Count::eta_pt_cut(int nbnf_index, float psrap_abs, float p_T)
{
    if (psrap_abs < 1)
    {
        count_this[0] = 1; //ptcut
        if (psrap_abs > 0.2 and psrap_abs < 0.8)
        {
            if (p_T > 0.3 and p_T < 1.5)
                nf_nb[nbnf_index] += 1;
        }
    }
}
void Count::Non_sin_diff(int nbnf_index,float psrap_abs,int IDIAG)
{
    nf_nb[nbnf_index+4] += 1;
    nch += 1;
}

void Count::Sin_Dou(int nbnf_index,float psrap,int IDIAG)
{
    int psrap_limit     = 10;
    float psrap_shift   = psrap*10;
    float psrap_count   = -50;

    if (psrap > -psrap_limit and psrap < psrap_limit)
    {
        for (int i = 0 ; i<counted_sin_size1 ; i++)
        {
            if (psrap_shift > psrap_count and psrap_shift < psrap_count+1)
            {
                if (IDIAG==1)
                {
                    counted_sin1 = 1;
                    nf_nb_sin1[i] += 1;
                }
                else if (IDIAG==6)
                {
                    counted_sin6 = 1;
                    nf_nb_sin6[i] += 1;
                }
                else if (IDIAG==10)
                {
                    counted_sin10 = 1;
                    nf_nb_sin10[i] += 1;
                }
                else if (IDIAG==11)
                {
                    counted_dou = 1;
                    nf_nb_dou[i] += 1;
                }
                break;
            }
            psrap_count += 1;
        }
    }
}


void Count::Filler(int npoms,int npomh)
{
    if (count_this[2])
    {
        NPOM_NCH[npoms][npomh]->Fill(nf_nb[4]+nf_nb[5]);
        N_CH->Fill(nch);
    }
    
    int diff_size = counted_sin_size1/2;
    for (int i=0 ; i<diff_size ; i++)
    {
        if (counted_sin1)
        {
            if(nf_nb_sin1[i+diff_size]){std::cout << nf_nb_sin1[i+diff_size] << std::endl;}
            NBNFSIN[0]->Fill(nf_nb_sin1[i],nf_nb_sin1[i+diff_size]);
            NFSIN[0]->Fill(nf_nb_sin1[i]);
            NBSIN[0]->Fill(nf_nb_sin1[i+diff_size]);
        }
        if (counted_sin6)
        {
            NBNFSIN[1]->Fill(nf_nb_sin6[i],nf_nb_sin6[i+diff_size]);
            NFSIN[1]->Fill(nf_nb_sin6[i]);
            NBSIN[1]->Fill(nf_nb_sin6[i+diff_size]);
        }
        if (counted_sin10)
        {
            NBNFSIN[2]->Fill(nf_nb_sin10[i],nf_nb_sin10[i+diff_size]);
            NFSIN[2]->Fill(nf_nb_sin10[i]);
            NBSIN[2]->Fill(nf_nb_sin10[i+diff_size]);
        }
        
        
        for(int i=0 ; i<diff_size ; i++)
        {
            if (counted_dou)
            {
                NBNFDOU->Fill(nf_nb_dou[i],nf_nb_dou[i+diff_size]);
                NFDOU->Fill(nf_nb_dou[i]);
                NBDOU->Fill(nf_nb_dou[i+diff_size]);
            }
        }
    }
    
    for (int i=0 ; i<count_this_size ; i++)
    {
        if (count_this[i])
        {
            NPOMSH[i][npoms][npomh]->Fill(nf_nb[2*i],nf_nb[2*i+1]);
            NPOMSH_nf[i][npoms][npomh]->Fill(nf_nb[2*i]);
            NPOMSH_nb[i][npoms][npomh]->Fill(nf_nb[2*i+1]);
        }
    }
}

void Count::Writer()
{
    output->cd(folders[7]);
    N_CH->Write();
    for(int i=0 ; i<NPOM_NCH_size ; i++)
    {
        for(int j=0 ; j<NPOM_NCHi_size ; j++)
            NPOM_NCH[i][j]->Write();
    }

    output->cd(folders[5]);
    for(int i=0 ; i<NBNFSIN_size ; i++)
    {
        NBNFSIN[i]->Write();
        NFSIN  [i]->Write();
        NBSIN  [i]->Write();
    }

    output->cd(folders[6]);
    //for(int i=0 ; i<NBNFDOU_size ; i++)
    //{
        NBNFDOU->Write();
        NFDOU->Write();
        NBDOU->Write();
    //}

    for(int k=0 ; k<prefix_size ; k++)
    {
        output->cd(folders[k]);
        for(int i=0 ; i<NPOMSHk_size ; i++)
        {
            for(int j=0 ; j<NPOMSHki_size ; j++)
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
    //delete returnTime;
}
