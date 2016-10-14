#include "Counter.h"
#include <iostream>
//#include "TRandom1.h"

int main() 
{

    std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/7000/data/";
    const char* outfile  = "7TeV_4M.root";
    double number_of_events = 4000000.0;
    Count NBNFandNPOM(datapath=datapath,number_of_events=number_of_events);
    float NBin=4; float start=0; float stop=3;
    NBNFandNPOM.InitializeNBNF();
    NBNFandNPOM.ReadAndCount();
    //NBNFandNPOM.ReadAndCount(outfile);

    //std::cout << NBNFandNPOM.B_MULT_loc << std::endl;
    /*
    TRandom1* myrand = new TRandom1();
    for (int i=0 ; i<100000 ; i++)
        NBNFandNPOM.ALL->Fill(myrand->Gaus(5,1));
    std::cout << NBNFandNPOM.ALL->GetBinContent(0)<<std::endl;
    std::cout << NBNFandNPOM.ALL->GetBinContent(1)<<std::endl;
    std::cout << NBNFandNPOM.ALL->GetBinContent(2)<<std::endl;
    std::cout << NBNFandNPOM.ALL->GetBinContent(3)<<std::endl;
    */
    return 0;
}
    
