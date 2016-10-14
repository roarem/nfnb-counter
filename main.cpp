#include "Counter.h"
#include <iostream>

int main() 
{

    std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/7000/data/";
    const char* outfile  = "7TeV_4M.root";
    double number_of_events = 4000000.0;
    Count NBNFandNPOM(datapath=datapath,number_of_events=number_of_events);
    float NBin=4; float start=0; float stop=3;
    NBNFandNPOM.InitializeNBNF();
    NBNFandNPOM.ReadAndCount();
    return 0;
}
    
