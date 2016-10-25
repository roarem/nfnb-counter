#include "Counter.h"
#include <iostream>


int main() 
{

    //std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/7000/data/";
    //double number_of_events = 4000000.0;

    std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/900/data/";
    double number_of_events = 1000000.0;
    
    //std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/13000/data/";
    //double number_of_events = 1000000.0;

    Count NBNFandNPOM(datapath=datapath,number_of_events=number_of_events);
    
    NBNFandNPOM.ReadAndCount();

    return 0;
}
    
