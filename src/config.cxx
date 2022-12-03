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
	ofstream fout("./config.yaml");
	fout << "#Printing version information"<<endl;
	fout << "hbuana:"<<endl;
	fout << "\t version: 1.0"<<endl;
	fout << "\t github: git@github.com:wangz1996/cepc_hbuana.git"<<endl;
	fout << "\n"<<endl;
	fout << "#Dat file to ROOT Decoder"<<endl;
	fout << "DAT-ROOT:"<<endl;
	fout << "\t on-off: False"<<endl;
	fout << "\t file-list: list_dat.txt"<<endl;
	fout << "\t output-dir: ./"<<endl;
	fout << "\n"<<endl;
	fout << "#Pedestal analyse manager"<<endl;
	fout << "Pedestal: "<<endl;
	fout << "\t on-off: False"<<endl;
	fout << "\t #If work in cosmic mode (hittag==0)"<<endl;
	fout << "\t Cosmic:"<<endl;
	fout << "\t \t on-off: False"<<endl;
	fout << "\t \t file-list: list.txt"<<endl;
	fout << "\t \t output-file: cosmic_pedestal.root"<<endl;
	fout << "\t \t usemt: False"<<endl;
	fout << "\t #If work in DAC mode (hittag==1 and skip the calibration channel)"<<endl;
	fout << "\t DAC:"<<endl;
	fout << "\t \t on-off: False"<<endl;
	fout << "\t \t file-list: list.txt"<<endl;
	fout << "\t \t output-file: dac_pedestal.root"<<endl;
	fout << "\n"<<endl;
	fout << "#DAC Calibration Manager"<<endl;
	fout << "Calibration:"<<endl;
	fout << "\t on-off: False"<<endl;
	fout << "\t #If work in cosmic mode"<<endl;
	fout << "\t Cosmic:"<<endl;
	fout << "\t \t on-off: False"<<endl;
	fout << "\t \t file-list: list.txt"<<endl;
	fout << "\t \t ped-file: pedestal_cosmic.root"<<endl;
	fout << "\t #If work in DAC mode"<<endl;
	fout << "\t DAC:"<<endl;
	fout << "\t \t on-off: False"<<endl;
	fout << "\t \t file-list: list.txt"<<endl;
	fout << "\t \t ped-file: pedestal_dac.root"<<endl;
	fout.close();
}


Config::Config() : conf(0)
{
}

Config::~Config()
{
}
