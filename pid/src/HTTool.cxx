#include "HTTool.h"

using namespace std;

HTTool::HTTool(const vector<int> &id,const vector<double> &x,const vector<double> &y,const vector<double> &z,const vector<double> &e)
{
	//cout<<"Init"<<endl;
	for(int i=0;i<40;i++)
	{
		layer_hits[i] = vector<HTTool::hit*>();
	}
	for(int j=0;j<id.size();j++)
	{
		int layer=id.at(j)/10000;
		HTTool::hit* tmp_hit = new HTTool::hit();
		tmp_hit->id=id.at(j);
		tmp_hit->x=x.at(j);
		tmp_hit->y=y.at(j);
		tmp_hit->z=z.at(j);
		tmp_hit->e=e.at(j);
		layer_hits[layer].emplace_back(tmp_hit);
	}
	for(int i=0;i<40;i++)
	{
		vector<vector<HTTool::hit*>> tmpcol = this->MergeAdjacent(layer_hits[i]);
		for(auto t:tmpcol)
		{
			if(t.size()>=10)continue;
			else
			{
				HTTool::hcl *tmp_hcl = new HTTool::hcl();
				double tmp_x=0.;
				double tmp_y=0.;
				double tmp_z=0.;
				double tmp_e=0.;
				for(auto z:t)
				{
					tmp_x+=z->x*z->e;
					tmp_y+=z->y*z->e;
					tmp_z+=z->z*z->e;
					tmp_e+=z->e;
				}
				tmp_x/=tmp_e;
				tmp_y/=tmp_e;
				tmp_z/=tmp_e;
				tmp_hcl->x=tmp_x;
				tmp_hcl->y=tmp_y;
				tmp_hcl->z=tmp_z;
				tmp_hcl->e=tmp_e;
				hcx.emplace_back(tmp_x);
				hcy.emplace_back(tmp_y);
				hcz.emplace_back(tmp_z);
				hce.emplace_back(tmp_e);
				if(tmp_e>0.05)Hough_Cluster.emplace_back(tmp_hcl);
			}
		}
	}
	int nbineta=100;
	hht=new TH2D("hht","test",nbineta,0.,3.14,160,-10,150);
	for(auto hh:Hough_Cluster)
	{
		for(int itheta=0;itheta<nbineta;itheta++)
        {
            double theta = itheta*3.1415/double(nbineta);
            double rho = hh->z/10.*TMath::Cos(theta) +
                        hh->x/10.*TMath::Sin(theta);
            hht->Fill(theta,rho);
        }
	}
	map<std::pair<int,int>,int> map_track;
	auto ftrack = [&map_track,this](pair<int,int> pos)->bool
	{
		for(auto i:map_track)
		{
			pair<int,int> posi = i.first;
			double d = sqrt(pow(posi.first-pos.first,2)+pow(posi.second-pos.second,2));
			if(d<10)
			{
				if(this->hht->GetBinContent(posi.first+1,posi.second+1)<this->hht->GetBinContent(pos.first+1,pos.second+1))
				{
					map_track.erase(map_track.find(posi));
					map_track.insert(pair<pair<int,int>,int>(pair<int,int>(pos.first,pos.second),1));
				}
				return true;
			}
		}
		return false;
	};
	for(int i=0;i<hht->GetNbinsX();i++)
	{
		for(int j=0;j<hht->GetNbinsY();j++)
		{
			if(hht->GetBinContent(i+1,j+1)>4)
			{
				if(!ftrack(pair<int,int>(i,j)))
				{
					map_track.insert(pair<pair<int,int>,int>(pair<int,int>(i,j),1));
				}
			}
		}
	}
	ntrack = map_track.size();
	//TCanvas *c1=new TCanvas("c1","test",1024,768);
	//c1->cd();
	//gStyle->SetOptStat("");
	//hht->Draw("colztext");
	//c1->SaveAs("test.png");
}
HTTool::~HTTool()
{
	vector<hcl*>().swap(Hough_Cluster);
	layer_hits.clear();
	delete hht;
}


vector<vector<HTTool::hit*>> HTTool::MergeAdjacent(vector<HTTool::hit*> hits)
{
	vector<vector<HTTool::hit*>> vechit;
	auto fadj=[](const HTTool::hit* _hit,const vector<HTTool::hit*> _hits)->bool
	{
		int m = ((_hit->id%10000)/100);
		int n = ((_hit->id%10000)%100);
		for(auto i:_hits)
		{
			int im = (i->id%10000)/100;
			int in = (i->id%10000)%100;
			if((m==im && abs(n-in)==1) || 
				(n==in && abs(m-im)==1) ||
				(abs(m-im)==1 && abs(n-in)==1 && false))
			{
				return true;
			}
		}
		return false;
	};
	while(hits.size()>0)
	{
		vechit.emplace_back(vector<HTTool::hit*>({hits.at(0)}));
		hits.erase(hits.begin());
		vector<int> hit2erase;
		while(1)
		{
			vector<int>().swap(hit2erase);
			for(int i=0;i<hits.size();i++)
			{
				if(fadj(hits.at(i),vechit.back()))
				{
					vechit.back().emplace_back(hits.at(i));
					hit2erase.emplace_back(i);
				}
			}
			if(hit2erase.size()==0)break;
			else
			{
				for(int i=(hit2erase.size()-1);i>=0;i--)
				{
					hits.erase(hits.begin()+i);
				}
			}
		}
	}
	return vechit;
}
