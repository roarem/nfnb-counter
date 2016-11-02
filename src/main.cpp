#include "Counter.h"
#include <iostream>


int main() 
{
    //std::string datapath = "/home/roar/DISKS/BIG1/19-31_oct/4mln/no_DD/7TeV/build/data/";
    std::string datapath = "/home/roar/DISKS/BIG1/19-31_oct/4mln/code_recieved_2810/7000/build/data/"; 
    //std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/7000/data/";
    const char* nbnfout = "7TeV_4M.root";
    //const char* nbnfout = "7TeV_4M_nsd.root";
    double number_of_events = 4000000.0;

    //std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/900/data/";
    //const char* nbnfout = "900GeV_1M.root"; 
    //const char* nbnfout = "900GeV_1M_nsd.root"; 
    //double number_of_events = 1000000.0;
    
    //std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/13000/data/";
    //const char* nbnfout = "13TeV_1M.root";
    //const char* nbnfout = "13TeV_1M_nsd.root";
    //double number_of_events = 1000000.0;

    Count NBNFandNPOM(nbnfout=nbnfout,datapath=datapath,number_of_events=number_of_events);

    NBNFandNPOM.ReadAndCount();

    return 0;
}
    
