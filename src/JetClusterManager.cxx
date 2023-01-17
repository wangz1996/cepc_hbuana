#include "JetClusterManager.h"

using namespace std;

JetClusterManager::JetClusterManager(const double &_r,const double &_a)
{
	this->SetR(_r);
	this->SetA(_a);
}

void JetClusterManager::DoClustering()
{
	while(EnergyLeft()>1. || NParLeft()>=1)
	{
		JetPar ParCore = FindCore(); // Find the core
		_par.emplace_back(vector<JetPar>()); // Init a empty collection
		_par.back().emplace_back(ParCore); // Add the core
		vector<int> h2erase={};
		while(NParLeft())
		{
			for(int j=0;j<_ParCollection.size();j++)
			{
				if(D(_ParCollection.at(j),_par.back())<_A)
				{
					_par.back().emplace_back(_ParCollection.at(j));
					h2erase.emplace_back(j);
				}
			}
			if(h2erase.size()==0)break;
			else
			{
				for(int i=h2erase.size()-1;i>=0;i--)
				{
					_ParCollection.erase(_ParCollection.begin()+h2erase.at(i));
				}
				h2erase.clear();
			}
		}
		NJet++;
		_npar.emplace_back(_par.back().size());
		_energy.emplace_back(0.);
		for_each(_par.back().begin(),_par.back().end(),[this](const JetPar &i)
		{
			this->_energy.back()+=i.e;
		});
		_x_width.emplace_back(XWidth(_par.back()));
		_y_width.emplace_back(YWidth(_par.back()));
		_z_width.emplace_back(ZWidth(_par.back()));
	}
}

double JetClusterManager::D(const JetPar &a,const JetPar &b=JetPar()) const
{
	double Rij_2 = (a.x-b.x)*(a.x-b.x)
				  +(a.y-b.y)*(a.y-b.y)
				  +(a.z-b.z)*(a.z-b.z);
	return Rij_2/_R;
}

double JetClusterManager::XWidth(const vector<JetPar> &a) const
{
	TH1D *h1=new TH1D("h1","test",100,-400,400);
	for(auto i:a)
	{
		h1->Fill(i.x,i.e);
	}
	double tmp = h1->GetRMS();
	delete h1;
	return tmp;
}

double JetClusterManager::YWidth(const vector<JetPar> &a) const
{
	TH1D *h1=new TH1D("h1","test",100,-400,400);
	for(auto i:a)
	{
		h1->Fill(i.y,i.e);
	}
	double tmp = h1->GetRMS();
	delete h1;
	return tmp;
}

double JetClusterManager::ZWidth(const vector<JetPar> &a) const
{
	TH1D *h1=new TH1D("h1","test",100,-400,400);
	for(auto i:a)
	{
		h1->Fill(i.z,i.e);
	}
	double tmp = h1->GetRMS();
	delete h1;
	return tmp;
}

double JetClusterManager::D(const JetPar &a,const vector<JetPar> &b) const
{
	double min=999.;
	for(auto i:b)
	{
		double tmp=D(a,i);
		min = (tmp<min)?tmp:min;
	}
	return min;
}

JetClusterManager::JetPar JetClusterManager::FindCore() 
{
	JetPar Parbeam;
	double min=99999.;
	int ier=0;
	for(int i=0;i<_ParCollection.size();i++)
	{
		if(D(_ParCollection.at(i))<min)
		{
			Parbeam = _ParCollection.at(i);
			ier = i;
		}
	}
	auto iter = _ParCollection.erase(_ParCollection.begin()+ier);
	return Parbeam;
}

double JetClusterManager::EnergyLeft() const
{
	double sum=0.;
	for_each(_ParCollection.begin(),_ParCollection.end(),[&sum](const JetPar &i)->void
	{
		sum+=i.e;
	});
	return sum;
}

void JetClusterManager::Print() const
{
	TH1D *h1=new TH1D("h1","test",100,0,1000);
	for(auto i=_ParCollection.begin()+1;i!=_ParCollection.end();i++)
	{
		cout<<D(*i,*(i-1))<<endl;
		h1->Fill(D(*i,*(i-1)));
	}
	TCanvas *c1=new TCanvas("c1","test",1024,768);
	c1->cd();
	gStyle->SetOptStat("");
	h1->Draw("");
	c1->SaveAs("hit0.png");
}
