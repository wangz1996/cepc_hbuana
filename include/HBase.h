#ifndef HBASE_HH
#define HBASE_HH

#include <TH2D.h>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class HBase{
public:
	//Constructor, destructor and instance of base class
	HBase();
	virtual ~HBase();

protected:
	//Protected member functions
	virtual void ReadTree(const TString &fname,const TString &tname); //Read TTree from ROOT files
	virtual void ReadList(const string &_list); // Read the file list and save to the protected vector
	virtual void CreateFile(const TString &_outname); // Create output file
	virtual void Init(const TString &_outname);// Initialize derived members
	
	// Protected member variables
	vector<string>	list;
	TFile *fin;	//TFile pointer to open files
	TTree *tin;	//Read Tree from fin
	TFile *fout;	//Create output files
	TTree *tout;	//Create output trees
	//Double_t        cycleID;
	//Double_t        triggerID;
	Int_t        cycleID;
	Int_t        triggerID;
	vector<int>     *cellIDs;
	vector<int>     *BCIDs;
	vector<int>     *hitTags;
	vector<int>     *gainTags;
	vector<double>  *charges;
	vector<double>  *times;

};

#endif
