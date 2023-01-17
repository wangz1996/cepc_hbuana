#include "PIDTool.h"

using namespace std;

int main(int argc,char* argv[])
{
	string file="",tree="";
	int train=0;
	for(int i=1;i<argc;i++)
	{
		if(string(argv[i])==string("-f"))
		{
			file = string(argv[i+1]);
		}
		else if(string(argv[i])==string("-t"))
		{
			tree = string(argv[i+1]);
		}
		else if(string(argv[i])==string("-train"))
		{
			train=1;
		}
	}
	cout<<file<<" "<<tree<<endl;
	PIDTool *pt = new PIDTool();
	if(file!="" && tree!="")
	{
		pt->GenNtuple(file,tree);
	}
	if(train==1)
	{
		pt->AddVar("xwidth",'D');
		pt->AddVar("ywidth",'D');
		pt->AddVar("zwidth",'D');
		pt->AddVar("Edep",'D');
		pt->AddVar("shower_start",'I');
		pt->AddVar("shower_layer_ratio",'D');
		pt->AddSig("pid_pion.root","T");
		pt->AddBkg("pid_e.root","T");
		pt->AddBkg("pid_muon.root","T");
		pt->TrainBDT("pion");
	}
	delete pt;
	return 0;
}
