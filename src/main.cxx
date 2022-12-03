#include "yaml-cpp/yaml.h"
#include "config.h"
#include "TFile.h"
#include <ctime>

using namespace std;

int main(int argc, char* argv[])
{
	time_t time1,time2;
	clock_t startTime,endTime;
	float diff_time;
	time(&time1);
	startTime = clock();
	Config config;
	//Initialize a config parser
	//...
	for(int i=1;i<argc;i++)
	{
		if(string(argv[i])=="-c")
		{
			string config_file="";
			config_file=string(argv[i+1]);
			config.Parse(config_file);
			config.Run();
		}
		else if(string(argv[i])=="-x")
		{
			config.Print();
		}
	}
	////////////////////////
	endTime = clock();
	time(&time2);
	diff_time = difftime(time2,time1);
	cout<<"Running(CPU) time: "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s."<<endl;
	cout<<"Actual time: "<<diff_time<<" s."<<endl;
	return 1;
}
