#include "yaml-cpp/yaml.h"
#include "JetClusterManager.h"
#include "TFile.h"
#include <ctime>

using namespace std;

int main(int argc, char* argv[])
{
	time_t time1,time2;
	clock_t startTime,endTime;
	float diff_time;
	time(&time1);
	startTime = clock();
	//...
	string file,tree;
	for(int i=1;i<argc;i++)
	{
		if(string(argv[i])=="-f")
		{
			file = string(argv[i+1]);
		}
		if(string(argv[i])=="-t")
		{
			tree = string(argv[i+1]);
		}
	}
	ROOT::DisableImplicitMT();
	ROOT::RDataFrame df(tree,file);
	typedef struct Jet
	{
		double e=0.;
		double xwidth=0.;
		double ywidth=0.;
		double zwidth=0.;
	}Jet;
	string outname=file;
	outname = outname.substr(outname.find_last_of('/')+1);
	outname = "jet_"+outname;
	JetClusterManager *jm = new JetClusterManager(100.,700.);
	vector<string> columns={"Hit_X","Hit_Y","Hit_Z","Hit_Energy"};
	TFile *f=TFile::Open(TString(file),"READ");
	TFile *fout=new TFile(TString(outname),"RECREATE");
	if(!f){cout<<"cannot open file"<<endl;throw file;}
	TTree *t=(TTree*)f->Get("T");
	vector<double> *Hit_X;Hit_X=0;
	vector<double> *Hit_Y;Hit_Y=0;
	vector<double> *Hit_Z;Hit_Z=0;
	vector<double> *Hit_Energy;Hit_Energy=0;
	t->SetBranchAddress("Hit_X",&Hit_X);
	t->SetBranchAddress("Hit_Y",&Hit_Y);
	t->SetBranchAddress("Hit_Z",&Hit_Z);
	t->SetBranchAddress("Hit_Energy",&Hit_Energy);
	TTree *tout=(TTree*)t->CloneTree(0);

	int njet;
	vector<double> *xwidth=0;
	vector<double> *ywidth=0;
	vector<double> *zwidth=0;
	vector<double> *j_e=0;
	tout->Branch("njet",&njet);
	tout->Branch("xwidth",&xwidth);
	tout->Branch("ywidth",&ywidth);
	tout->Branch("zwidth",&zwidth);
	tout->Branch("j_e",&j_e);
	for(int i=0;i<t->GetEntries();i++)
	{
		if(i%1000==0)cout<<i<<endl;
		t->GetEntry(i);
		jm->Clear();
		xwidth->clear();
		ywidth->clear();
		zwidth->clear();
		j_e->clear();
		for(int j=0;j<Hit_X->size();j++)
		{
			jm->AddPar(Hit_X->at(j),Hit_Y->at(j),Hit_Z->at(j),Hit_Energy->at(j));
		}
		jm->DoClustering();
		njet = jm->GetNJet();
		//cout<<i<<" "<<Hit_X->size()<<" "<<njet<<" "<<endl;
		for(int j=0;j<njet;j++)
		{
			xwidth->emplace_back(jm->GetXWidth(j));
			ywidth->emplace_back(jm->GetYWidth(j));
			zwidth->emplace_back(jm->GetZWidth(j));
			j_e->emplace_back(jm->GetEnergy(j));
		}
		tout->Fill();
	}
	//f->Close();
	fout->cd();
	tout->Write();
	fout->Close();
	
	////////////////////////
	endTime = clock();
	time(&time2);
	diff_time = difftime(time2,time1);
	cout<<"Running(CPU) time: "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s."<<endl;
	cout<<"Actual time: "<<diff_time<<" s."<<endl;
	return 1;
}
