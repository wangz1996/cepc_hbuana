#ifndef PEDESTALMANAGER_HH
#define PEDESTALMANAGER_HH

#include "HBase.h"
#include <TH2D.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <fstream>
#include <string>
#include <atomic>
#include <mutex>
#include <future>
#include <thread>
#include <utility>
#include <chrono>

using namespace std;

class PedestalManager : public HBase{

public:
	static PedestalManager* CreateInstance();
	static void DeleteInstance();
	~PedestalManager();

private:
	PedestalManager();

public:
	//Delete Copy constructor
	PedestalManager(const PedestalManager &) = delete;
	PedestalManager &operator=(PedestalManager const &) = delete;

	void Init(const TString &_outname);
	int AnaPedestal(const std::string &list,const int &sel_hittag);
	void Setmt(bool mt){usemt = mt;};
	
private:
	//using HBase::HBase;
	mutex mtx;
	bool usemt;
	vector<int> vec_cellid;
	std::unique_ptr<TH2D> htimepeak;
	std::unique_ptr<TH2D> htimerms;
	std::unique_ptr<TH2D> hchargepeak;
	std::unique_ptr<TH2D> hchargerms;
	unordered_map<int,TH2D*> map_layer_timepeak;
	unordered_map<int,TH2D*> map_layer_timerms;
	unordered_map<int,TH2D*> map_layer_chargepeak;
	unordered_map<int,TH2D*> map_layer_chargerms;
	unordered_map<int,TH1D*> map_cellid_htime;
	unordered_map<int,TH1D*> map_cellid_hcharge;
	double time_min=1000.,time_max=0.;
	double charge_min=1000.,charge_max=0.;
	//Branch Name
	double time_peak=0.,time_rms=0.,charge_peak=0.,charge_rms=0.;
	int _cellid;
	
	void SaveCanvas(TH2D* h,const TString &name);
};

extern PedestalManager *_instance;

#endif
