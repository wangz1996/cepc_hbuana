#ifndef PIDTOOL_HH
#define PIDTOOL_HH
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <ROOT/RDataFrame.hxx>
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
using namespace std;
class PIDTool
{
public:
	PIDTool();
	~PIDTool();
	int GenNtuple(const string &file,const string &tree);
	void AddSig(const string &file,const string &tree)
	{
		sig.insert(pair<TString,TString>(TString(file),TString(tree)));
	}
	void AddBkg(const string &file,const string &tree)
	{
		bkg.insert(pair<TString,TString>(TString(file),TString(tree)));
	}
	void AddVar(const string &v,const char &type)
	{
		var.insert(pair<TString,char>(TString(v),type));
	}
	int TrainBDT(const string &bdtname);
	void Clear()
	{
		var.clear();
		sig.clear();
		bkg.clear();
	}
private:
	map<TString,char> var;
	map<TString,TString> sig;
	map<TString,TString> bkg;
	
};
#endif
