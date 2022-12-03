#include "config.h"
#include "DatManager.h"
#include "DacManager.h"
#include "PedestalManager.h"
#include <fstream>

using namespace std;

void Config::Parse(const string config_file)
{
	conf = YAML::LoadFile(config_file);
	//Printing version number
	cout<<"HBUANA Version: "<<conf["hbuana"]["version"].as<std::string>()<<endl;
	cout<<"HBUANA Github Repository: "<<conf["hbuana"]["github"].as<std::string>()<<endl;
}

int Config::Run()
{
	if(conf["DAT-ROOT"]["on-off"].as<bool>())
	{
		cout<<"DAT mode: ON"<<endl;
		if(conf["DAT-ROOT"]["file-list"].as<std::string>()=="" || conf["DAT-ROOT"]["output-dir"].as<std::string>()=="")
		{
			cout<<"ERROR: Please specify file list or output-dir for dat files"<<endl;
		}
		else
		{
			ifstream dat_list(conf["DAT-ROOT"]["file-list"].as<std::string>());
			DatManager dm;
			while(!dat_list.eof())
			{
				string dat_temp;
				dat_list >> dat_temp;
				if(dat_temp=="")continue;
				dm.Decode(dat_temp,conf["DAT-ROOT"]["output-dir"].as<std::string>());
			}
		}
	}
	if(conf["Pedestal"]["on-off"].as<bool>())
	{
		PedestalManager::CreateInstance();
		cout<<"Pedestal mode: ON"<<endl;
		if(conf["Pedestal"]["Cosmic"]["on-off"].as<bool>())
		{
			cout<<"Pedestal mode for cosmic events: ON"<<endl;
			_instance->Init(conf["Pedestal"]["Cosmic"]["output-file"].as<string>().c_str());
			_instance->Setmt(conf["Pedestal"]["Cosmic"]["usemt"].as<bool>());
			_instance->AnaPedestal(conf["Pedestal"]["Cosmic"]["file-list"].as<std::string>(),0);
			PedestalManager::DeleteInstance();
		}
		if(conf["Pedestal"]["DAC"]["on-off"].as<bool>())
		{
			cout<<"Pedestal mode for DAC events: ON"<<endl;
			_instance->Init(conf["Pedestal"]["DAC"]["output-file"].as<string>().c_str());
			_instance->AnaPedestal(conf["Pedestal"]["DAC"]["file-list"].as<std::string>(),1);
			PedestalManager::DeleteInstance();
		}
	}
	if(conf["Calibration"]["on-off"].as<bool>())
	{
		if(conf["Calibration"]["Cosmic"]["on-off"].as<bool>())
		{
			cout<<"Cosmic calibration mode:ON"<<endl;
			DacManager dacmanager("cosmic_calib.root");
			dacmanager.SetPedestal(conf["Calibration"]["Cosmic"]["ped-file"].as<string>().c_str());
			dacmanager.AnaDac(conf["Calibration"]["Cosmic"]["file-list"].as<std::string>(),"cosmic");
		}
		if(conf["Calibration"]["DAC"]["on-off"].as<bool>())
		{
			cout<<"DAC Calibration mode:ON"<<endl;
			DacManager dacmanager("dac_calib.root");
			dacmanager.SetPedestal(conf["Calibration"]["DAC"]["ped-file"].as<string>().c_str());
			dacmanager.AnaDac(conf["Calibration"]["DAC"]["file-list"].as<std::string>().c_str(),"dac");
		}
	}
	return 1;
}

void Config::Print()
{
	YAML::Node example_config = YAML::LoadFile("/mnt2/SJTU/wangzhen/hbuana/config/config.yaml");
	ofstream fout("./config.yaml");
	fout << example_config;
	fout.close();
}


Config::Config() : conf(0)
{
}

Config::~Config()
{
}
