#ifndef JETCLUSTERMANAGER_HH
#define JETCLUSTERMANAGER_HH

#include <TH2D.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include <fstream>
#include <string>
#include <atomic>
#include <mutex>
#include <future>
#include <thread>
#include <utility>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <ROOT/RDataFrame.hxx>
#include "TROOT.h"
#include <algorithm>

using namespace std;

class JetClusterManager{

private:
	typedef struct JetPar
	{
		double x=0.;
		double y=0.;
		double z=0.;
		double e=0.;
	}JetPar;
public:
	JetClusterManager(const double &_r,const double &_a);
	~JetClusterManager(){};
	virtual void AddPar(const double &_x,
						const double &_y,
						const double &_z,
						const double &_e)
	{
		JetPar tmp = {_x,_y,_z,_e};
		_ParCollection.emplace_back(tmp);
	}
	virtual void DoClustering();
	virtual void Print() const;
	virtual void Clear()
	{
		NJet=0;
		vector<int>().swap(_npar);
		vector<double>().swap(_energy);
		vector<JetPar>().swap(_ParCollection);
		vector<vector<JetPar>>().swap(_par);
		vector<double>().swap(_x_width);
		vector<double>().swap(_y_width);
		vector<double>().swap(_z_width);
	}
	virtual double GetEnergy(const int &ijet=0) const
	{
		return _energy.at(ijet);
	}
	virtual int GetNPar(const int &ijet=0) const
	{
		return _npar.at(ijet);
	}
	virtual int GetNJet(){return NJet;}
	virtual double GetXWidth(const int &ijet=0){return _x_width.at(ijet);}
	virtual double GetYWidth(const int &ijet=0){return _y_width.at(ijet);}
	virtual double GetZWidth(const int &ijet=0){return _z_width.at(ijet);}
	double GetParX(const int &ijet,const int &ipar){return _par.at(ijet).at(ipar).x;}
	double GetParY(const int &ijet,const int &ipar){return _par.at(ijet).at(ipar).y;}
	double GetParZ(const int &ijet,const int &ipar){return _par.at(ijet).at(ipar).z;}

private:
	double D(const JetPar &a,const JetPar &b) const;
	double D(const JetPar &a,const vector<JetPar> &b) const;
	double XWidth(const vector<JetPar> &a) const;
	double YWidth(const vector<JetPar> &a) const;
	double ZWidth(const vector<JetPar> &a) const;
	double EnergyLeft() const;
	void SetR(const double &_r){_R = _r;}
	void SetA(const double &_a){_A = _a;}
	int NParLeft() const
	{
		return _ParCollection.size();
	}
	JetPar FindCore();
	vector<int> _npar={};
	vector<double> _energy={};
	vector<vector<JetPar>> _par={};
	vector<double> _x_width={};
	vector<double> _y_width={};
	vector<double> _z_width={};
	int NJet=0;
	double _R=0.;
	double _A=0.;
	vector<JetPar> _ParCollection;

public:
	//Delete Copy constructor
	JetClusterManager(const JetClusterManager &) = delete;
	JetClusterManager &operator=(JetClusterManager const &) = delete;

	
private:
	
};

#endif
