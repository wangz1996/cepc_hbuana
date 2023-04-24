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
	int GenNtuple(const string &file,const string &tree);
	void AddSignal(const string &file,const string &tree,const string &particle_name)
	{
		signal.insert(pair<pair<TString,TString>,TString>(pair<TString,TString>(TString(file),TString(tree)),TString(particle_name)));
	}
	void AddVar(const string &v,const char &type)
	{
		var.insert(pair<TString,char>(TString(v),type));
	}
	int TrainBDT();
	int BDTNtuple(const string &fname,const string &tname);
	void Clear()
	{
		var.clear();
		signal.clear();
	}
private:
	map<TString,char> var;
	map<pair<TString,TString>,TString> signal;
	
};
#endif
