#include "Counter.h"
#include <iostream>

int main(int argc, char* argv[]) 
{
    std::string input = std::string(argv[1]);
    if (argc<2)
    {
	std::cout << "Need input, yo" << std::endl;
	exit(1);
    } 
    else if (input==std::string("900"))
    {
	//std::string datapath = "/home/roar/master/qgsm_lhc_analyze/data/rawData/900/data/";
    	std::string datapath = "/home/roar/DISKS/1/19-31_oct/4mln/code_recieved_2810/900/build/data/"; 
    	const char* nbnfout = "900_4M.root"; 
    	const char* bcorrout= "900_4M_bcorr.csv";
    	double number_of_events = 4000000.0;
	Count NBNFandNPOM(nbnfout=nbnfout,bcorrout=bcorrout,datapath=datapath,number_of_events=number_of_events);
    	NBNFandNPOM.ReadAndCount();
    }
    else if (input==std::string("2760"))
    {
	std::string datapath = "/home/roar/DISKS/1/19-31_oct/4mln/code_recieved_2810/2760/build/data/"; 
    	const char* nbnfout = "2760_4M.root";
    	const char* bcorrout= "2760_4M_bcorr.csv";
    	double number_of_events = 4000000.0;
	Count NBNFandNPOM(nbnfout=nbnfout,bcorrout=bcorrout,datapath=datapath,number_of_events=number_of_events);
    	NBNFandNPOM.ReadAndCount();
    }
    else if (input==std::string("7000"))
    {
    	std::string datapath = "/home/roar/DISKS/1/19-31_oct/4mln/code_recieved_2810/7000/build/data/"; 
    	//std::string datapath = "/home/roar/DISKS/1/30_aug/7000/data/";
    	const char* nbnfout = "7000_4M.root";
    	//const char* nbnfout = "7000_4M.nsd_count.root";
    	//const char* nbnfout = "7000_4M.test.root";
    	const char* bcorrout= "7000_4M_bcorr.csv";
    	double number_of_events = 4000000.0;
	Count NBNFandNPOM(nbnfout=nbnfout,bcorrout=bcorrout,datapath=datapath,number_of_events=number_of_events);
    	NBNFandNPOM.ReadAndCount();
    }
    else if (input==std::string("13000"))
    {
	std::string datapath = "/home/roar/DISKS/1/13000_attempts/";
    	const char* nbnfout = "13000_4M.root";
    	const char* bcorrout= "13000_4M_bcorr.csv";
    	double number_of_events = 4580316.0;
	Count NBNFandNPOM(nbnfout=nbnfout,bcorrout=bcorrout,datapath=datapath,number_of_events=number_of_events);
    	NBNFandNPOM.ReadAndCount();
    }
    else if (input==std::string("alt"))
    {
	std::string datapath = "/home/roar//DISKS/1/13000_attempts/900/200k/";
    	const char* nbnfout = "900_200k.root";
    	const char* bcorrout= "900_200k_bcorr.csv";
    	double number_of_events = 200000.0;
	Count NBNFandNPOM(nbnfout=nbnfout,bcorrout=bcorrout,datapath=datapath,number_of_events=number_of_events);
    	NBNFandNPOM.ReadAndCount();
    }
    else
    {
	std::cout << "Wrong input, yo" << std::endl;
	exit(1);
    }

    return 0;
}
    
