#ifndef HTTOOL_HH
#define HTTOOL_HH
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
#include <unordered_map>
using namespace std;

class HTTool
{
public:
	HTTool(const vector<int> &id,const vector<double> &x,const vector<double> &y,const vector<double> &z,const vector<double> &e);
	~HTTool();
	int GetNtrack(){return ntrack;};
	vector<double> GetHclX(){return hcx;};
	vector<double> GetHclY(){return hcy;};
	vector<double> GetHclZ(){return hcz;};
	vector<double> GetHclE(){return hce;};
	
	
private:
	int ntrack=0;
	typedef struct hit
	{
		int id=0;
		double x=0.;
		double y=0.;
		double z=0.;
		double e=0.;
	}hit;
	typedef struct hcl
	{
		double x=0.;
		double y=0.;
		double z=0.;
		double e=0.;
	}hcl;
	TH2D *hht;
	unordered_map<int,vector<hit*>> layer_hits;
	vector<hcl*> Hough_Cluster;
	vector<vector<hit*>> MergeAdjacent(vector<hit*> hits);
	vector<double> hcx;
	vector<double> hcy;
	vector<double> hcz;
	vector<double> hce;
};

#endif
