#ifndef PIDTOOL_HH
#define PIDTOOL_HH
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <cstdlib>
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "HTTool.h"
//#include "Fit3D.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/TMVAMultiClassGui.h"
using namespace std;
using namespace ROOT;

class PIDTool
{
public:
	PIDTool();
	~PIDTool();

	Int_t GenNtuple(const string& file, const string& tree);

	void AddTrainSignal(const string& file, const string& tree, const string& particle_name)
	{
		train_signal.insert(pair<pair<TString, TString>, TString>(pair<TString, TString>(TString(file), TString(tree)), TString(particle_name)));
	}

    void AddTrainBkg(const string& file, const string& tree, const string& particle_name)
    {
        train_bkg.insert(pair<pair<TString, TString>, TString>(pair<TString, TString>(TString(file), TString(tree)), TString(particle_name)));
    }

	void AddTestSignal(const string& file, const string& tree, const string& particle_name)
	{
		test_signal.insert(pair<pair<TString, TString>, TString>(pair<TString, TString>(TString(file), TString(tree)), TString(particle_name)));
	}

	void AddTestBkg(const string& file, const string& tree, const string& particle_name)
	{
		test_bkg.insert(pair<pair<TString, TString>, TString>(pair<TString, TString>(TString(file), TString(tree)), TString(particle_name)));
	}

	void AddVar(const string& v, const Char_t& type)
	{
		var.insert(pair<TString, Char_t>(TString(v), type));
	}

	Int_t TrainBDT();
	Int_t BDTNtuple(const string& fname, const string& tname);

	void Clear()
	{
		var.clear();
		train_signal.clear();
        train_bkg.clear();
        test_signal.clear();
		test_bkg.clear();
	}


private:
	map<TString, Char_t> var;
	map<pair<TString, TString>, TString> train_signal;
	map<pair<TString, TString>, TString> train_bkg;
	map<pair<TString, TString>, TString> test_signal;
	map<pair<TString, TString>, TString> test_bkg;
};

#endif
