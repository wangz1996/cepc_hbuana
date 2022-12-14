#include "HBase.h"

using namespace std;

HBase::HBase() : fin(0),fout(0),tin(0),tout(0)
{
	list.clear();
	cout<<"HBase class instance initialized."<<endl;
}

HBase::~HBase()
{
	cout<<"Base destructor called"<<endl;
	fout->Close();
	//fin->Close();
}

void HBase::Init(const TString &_outname)
{
	CreateFile(_outname);
}

void HBase::CreateFile(const TString &_outname)
{
	fout = new TFile(TString(_outname),"RECREATE");
}


void HBase::ReadList(const string &_list)
{
	ifstream data(_list);
	while(!data.eof())
	{
		string temp;
		data>>temp;
		if(temp=="")continue;
		list.push_back(temp);
	}
}

void HBase::ReadTree(const TString &fname,const TString &tname)
{
	cout<<"Reading tree "<<fname<<endl;
	fin = TFile::Open(TString(fname),"READ");
	tin = (TTree*)fin->Get(TString(tname));
	tin->SetBranchAddress("cycleID", &cycleID);
	tin->SetBranchAddress("triggerID", &triggerID);
	tin->SetBranchAddress("cellIDs", &cellIDs);
	tin->SetBranchAddress("BCIDs", &BCIDs);
	tin->SetBranchAddress("hitTags", &hitTags);
	tin->SetBranchAddress("gainTags", &gainTags);
	tin->SetBranchAddress("charges", &charges);
	tin->SetBranchAddress("times", &times);
}
