#ifndef CONFIG_HH
#define CONFIG_HH
#include <iostream>
#include <string>
#include <vector>
#include "yaml-cpp/yaml.h"

using namespace std;

class Config
{
public:

	YAML::Node conf;

	Config();
	~Config();

	virtual void Print();
	virtual void Parse(const string config_file);
	virtual int Run();
};



#endif
