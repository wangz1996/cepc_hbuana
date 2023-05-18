#include "PedestalManager.h"
#include "TStyle.h"
#include <TH2.h>
#include <memory>
#include <iostream>
#include <TCanvas.h>
#include <sstream>
#include <algorithm>

using namespace std;

PedestalManager *_instance = nullptr;
//Get Instance Class
PedestalManager *PedestalManager::CreateInstance()
{
	if(_instance == nullptr)
	{
		_instance = new PedestalManager();
	}
	return _instance;
}

void PedestalManager::DeleteInstance()
{
	if(_instance == nullptr)
	{

	}
	else
	{
		delete _instance;
	}
}

PedestalManager::PedestalManager()
{
	list.clear();
	cout<<"PedestalManager class instance initialized."<<endl;
}

void PedestalManager::Init(const TString &_outname)
{
	CreateFile(_outname);
	tout = new TTree("pedestal","Pedestal");
	cellIDs = 0;
	BCIDs = 0;
	hitTags = 0;
	gainTags = 0;
	charges = 0;
	times = 0;
	tout->Branch("cellid",&_cellid);
	tout->Branch("time_peak",&time_peak);
	tout->Branch("time_rms",&time_rms);
	tout->Branch("charge_peak",&charge_peak);
	tout->Branch("charge_rms",&charge_rms);
	htimepeak=std::make_unique<TH2D>("htimepeak","HighGain Peak",360,0,360,36,0,36);
	htimerms=std::make_unique<TH2D>("htimerms","HighGain RMS",360,0,360,36,0,36);
	hchargepeak=std::make_unique<TH2D>("hchargepeak","LowGain Peak",360,0,360,36,0,36);
	hchargerms=std::make_unique<TH2D>("hchargerms","LowGain RMS",360,0,360,36,0,36);
	int ini_cellid=0;
	for(int i_layer=0;i_layer<40;i_layer++)
	{
		for(int i_chip=0;i_chip<9;i_chip++)
		{
			for(int i_chn=0;i_chn<36;i_chn++)
			{
				ini_cellid=i_layer*100000+i_chip*10000+i_chn;
				vec_cellid.push_back(ini_cellid);
				TString timename="htime_"+TString(to_string(ini_cellid).c_str());
				map_cellid_htime[ini_cellid] = new TH1D(timename,timename,1500,0,1500);
				TString chargename="hcharge_"+TString(to_string(ini_cellid).c_str());
				map_cellid_hcharge[ini_cellid] = new TH1D(chargename,chargename,1600,0,1600);
			}
		}
	}
	for(int i=0;i<40;i++)
	{
		TString name_timepeak="htimepeak_"+TString(to_string(i).c_str());
		map_layer_timepeak[i] = new TH2D(name_timepeak,name_timepeak,9,0,9,36,0,36);
		TString name_timerms="htimerms_"+TString(to_string(i).c_str());
		map_layer_timerms[i] = new TH2D(name_timerms,name_timerms,9,0,9,36,0,36);
		TString name_chargepeak="hchargepeak_"+TString(to_string(i).c_str());
		map_layer_chargepeak[i] = new TH2D(name_chargepeak,name_chargepeak,9,0,9,36,0,36);
		TString name_chargerms="hchargerms_"+TString(to_string(i).c_str());
		map_layer_chargerms[i] = new TH2D(name_chargerms,name_chargerms,9,0,9,36,0,36);
	}
	cout<<"Initialization done"<<endl;
}

int PedestalManager::AnaPedestal(const std::string &_list,const int &sel_hittag)
{
	cout<<"Starting Ana"<<endl;
	cout<<"Ana preparation done"<<endl;
	ReadList(_list); // read file list _list to list
	cout<<"read list done"<<endl;
	if(usemt){
		ROOT::EnableImplicitMT();
		ROOT::EnableThreadSafety();
		const int nthreads = 10;
		thread t[nthreads];
		vector<vector<string>> listth;
		int nfile_pert = list.size()/nthreads +1;
		if(list.size()%nthreads == 0)nfile_pert--;
		auto f = [this,sel_hittag](vector<string> list_tmp)
		{
			for(auto tmp:list_tmp)
			{
				string skipchannel = tmp;
				int dac_chn=-1;// Which channel should not be used here for pedestal analysis
				if(sel_hittag == 1)
				{
				skipchannel = skipchannel.substr(skipchannel.find_last_of('/')+1);
				skipchannel = skipchannel.substr(skipchannel.find("chn")+3);
				skipchannel = skipchannel.substr(0,skipchannel.find_last_of('_'));
				dac_chn = stoi(skipchannel);
				}
				//this->ReadTree(TString(tmp.c_str()),"Cosmic_Event");
				TFile *_fin = TFile::Open(TString(tmp.c_str()),"READ");
				TTree *_tin = (TTree*)_fin->Get("Raw_Hit");
				Double_t        _cycleID;
				Double_t        _triggerID; 
				vector<int>     *_cellIDs;
				vector<int>     *_BCIDs;
				vector<int>     *_hitTags;
				vector<int>     *_gainTags;
				vector<double>  *_charges;
				vector<double>  *_times;
				_cellIDs=0;_BCIDs=0;_hitTags=0;_gainTags=0;_charges=0;_times=0;
				_tin->SetBranchAddress("cycleID", &_cycleID);
				_tin->SetBranchAddress("triggerID", &_triggerID);
				_tin->SetBranchAddress("cellIDs", &_cellIDs);
				_tin->SetBranchAddress("BCIDs", &_BCIDs);
				_tin->SetBranchAddress("hitTags", &_hitTags);
				_tin->SetBranchAddress("gainTags", &_gainTags);
				_tin->SetBranchAddress("charges", &_charges);
				_tin->SetBranchAddress("times", &_times);
				int Nentry = _tin->GetEntries();
				for(int ientry=0;ientry<Nentry;ientry++)
				{
					_tin->GetEntry(ientry);
					for(int i=0;i<_hitTags->size();i++)
					{
						if(_hitTags->at(i)!=sel_hittag)continue;
						int cellid = _cellIDs->at(i);
						int channel = cellid%100;
						int memo = (cellid%10000)/100;
						if(memo !=0 )continue;
						if(dac_chn==channel)continue;
						int layer = cellid/1e5;
						if(_times->at(i)<time_min)time_min=_times->at(i);
						if(_times->at(i)>time_max)time_max=_times->at(i);
						if(_charges->at(i)<charge_min)charge_min=_charges->at(i);
						if(_charges->at(i)>charge_max)charge_max=_charges->at(i);
						mtx.lock();
						map_cellid_htime[cellid]->Fill(_times->at(i));
						map_cellid_hcharge[cellid]->Fill(_charges->at(i));
						mtx.unlock();
					}
				}
			//cout<<tmp<<endl;
			}
		};
		for(int ith=0;ith<nthreads;ith++)
		{
			int begin = ith*nfile_pert;
			int end = begin + nfile_pert - 1;
			end = end>list.size()-1? list.size()-1 : end;
			vector<string> tmp_list;
			for(int i=begin;i<=end;i++)tmp_list.push_back(list.at(i));
			listth.push_back(tmp_list);
		}
		for(int ith=0;ith<nthreads;ith++)
		{
			async(launch::async,f,listth.at(ith));
			//t[ith] = thread(f,listth.at(ith));
		}
		//for(int i=0;i<nthreads;i++)
		//{
		//	if(t[i].joinable())t[i].join();
		//}
	}
	else
	{
		for_each(list.begin(),list.end(),[this,sel_hittag](string tmp)
				{
				string skipchannel = tmp;
				int dac_chn=-1;// Which channel should not be used here for pedestal analysis
				if(sel_hittag == 1)
				{
				skipchannel = skipchannel.substr(skipchannel.find_last_of('/')+1);
				skipchannel = skipchannel.substr(skipchannel.find("chn")+3);
				skipchannel = skipchannel.substr(0,skipchannel.find_last_of('_'));
				dac_chn = stoi(skipchannel);
				}
				this->ReadTree(TString(tmp.c_str()),"Raw_Hit");
				int Nentry = tin->GetEntries();
				for(int ientry=0;ientry<Nentry;ientry++)
				{
					tin->GetEntry(ientry);
					for(int i=0;i<hitTags->size();i++)
					{
						if(hitTags->at(i)!=sel_hittag)continue;
						int cellid = cellIDs->at(i);
						int channel = cellid%100;
						int memo = (cellid%10000)/100;
						if(memo !=0 )continue;
						if(dac_chn==channel)continue;
						int layer = cellid/1e5;
						if(times->at(i)<time_min)time_min=times->at(i);
						if(times->at(i)>time_max)time_max=times->at(i);
						if(charges->at(i)<charge_min)charge_min=charges->at(i);
						if(charges->at(i)>charge_max)charge_max=charges->at(i);
						map_cellid_htime[cellid]->Fill(times->at(i));
						map_cellid_hcharge[cellid]->Fill(charges->at(i));
					}
				}
				cout<<tmp<<endl;
				});
	}
	// Analysis done
	//
	// Fill the output tree
	for_each(vec_cellid.begin(),vec_cellid.end(),[this](int i)->void
			{
			_cellid = i ;
			time_peak=map_cellid_htime[i]->GetBinCenter(map_cellid_htime[i]->GetMaximumBin());
			time_rms=map_cellid_htime[i]->GetRMS();
			charge_peak=map_cellid_hcharge[i]->GetBinCenter(map_cellid_hcharge[i]->GetMaximumBin());
			charge_rms=map_cellid_hcharge[i]->GetRMS();
			tout->Fill();
			});
	cout<<"Out Tree Filled"<<endl;
	fout->cd();
	tout->Write();
	cout<<"time min: "<<time_min<<" max: "<<time_max<<endl;
	cout<<"charge min: "<<charge_min<<" max: "<<charge_max<<endl;

	// Save hists into the output file:
	// mode_name = times or charges
	// tmp_map is the map from cellid to TH1D to loop over
	// tmp_layer_timepeak and tmp_layer_timerms is the map from layer to peak and rms
	// hpeak and hrms are the general TH2D for peak and rms
	// alias is the alternative name. high or low
	auto f_save = [this](TString mode_name,unordered_map<int,TH1D*> tmp_map,unordered_map<int,TH2D*> tmp_layer_timepeak,unordered_map<int,TH2D*> tmp_layer_timerms,std::unique_ptr<TH2D> &hpeak,std::unique_ptr<TH2D> &hrms,TString alias)
	{
		fout->mkdir(TString(mode_name));
		fout->cd(TString(mode_name));
		for(int i=0;i<40;i++)gDirectory->mkdir(TString("layer_")+TString(to_string(i).c_str()));
		for(auto i:tmp_map)
		{
			int cellid=i.first;
			int layer = cellid/1e5;
			int channel = cellid%100;
			int chip = (cellid%100000)/10000;
			tmp_layer_timepeak[layer]->Fill(chip,channel,i.second->GetBinCenter(i.second->GetMaximumBin()));
			tmp_layer_timerms[layer]->Fill(chip,channel,i.second->GetRMS());      
			hpeak->Fill(layer*9+chip,channel,i.second->GetBinCenter(i.second->GetMaximumBin()));
			hrms->Fill(layer*9+chip,channel,i.second->GetRMS());
			TString dir_name = TString(mode_name+"/layer_") + TString(to_string(layer).c_str());
			fout->cd(dir_name);
			i.second->Write();
		}
		for(int i=0;i<40;i++)
		{
			TString dir_name = TString(mode_name+"/layer_") + TString(to_string(i).c_str());
			fout->cd(dir_name);
			tmp_layer_timepeak[i]->Write();
			tmp_layer_timerms[i]->Write();
			//this->SaveCanvas(tmp_layer_timepeak[i],alias+TString("gain_peak_")+to_string(i).c_str());
			//this->SaveCanvas(tmp_layer_timerms[i],alias+TString("gain_rms_")+to_string(i).c_str());
		}
		cout<<mode_name<<" done"<<endl;
	};
	f_save("times",map_cellid_htime,map_layer_timepeak,map_layer_timerms,htimepeak,htimerms,"high");
	f_save("charges",map_cellid_hcharge,map_layer_chargepeak,map_layer_chargerms,hchargepeak,hchargerms,"low");
	fout->cd("");
	htimepeak->Write();
	htimerms->Write();
	hchargepeak->Write();
	hchargerms->Write();
	//fout->Close();
	cout<<"2D hists written"<<endl;
	return 0;
}

void PedestalManager::SaveCanvas(TH2D* h,const TString &name)
{
	gStyle->SetPaintTextFormat("4.1f");
	std::unique_ptr<TCanvas> c1=std::make_unique<TCanvas>(TString(name),TString(name),1024,768);
	c1->cd();
	gStyle->SetOptStat("");
	h->Draw("colztext");
	TString outname=name+".png";
	c1->SaveAs(TString(outname));
}

PedestalManager::~PedestalManager()
{
	cout<<"Pedestal destructor called"<<endl;
}
