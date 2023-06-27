#include "HTTool.h"

using namespace std;

HTTool::HTTool(const vector<Double_t>& x, const vector<Double_t>& y, const vector<Double_t>& z, const vector<Double_t>& e)
{
    const Int_t nlayer = 40;
    const Int_t thick = 30;

	//cout << "Initialising..." << endl;

	for (Int_t i = 0; i < nlayer; i++)
		layer_hits[i] = vector<HTTool::hit*>();

    for (Int_t j = 0; j < x.size(); j++)
	{
        Int_t layer = z.at(j) / thick;

		HTTool::hit* tmp_hit = new HTTool::hit();

		tmp_hit->x = x.at(j);
		tmp_hit->y = y.at(j);
		tmp_hit->z = z.at(j);
		tmp_hit->e = e.at(j);

		layer_hits[layer].emplace_back(tmp_hit);
	}

	for (Int_t i = 0; i < nlayer; i++)
	{
		vector<vector<HTTool::hit*>> tmpcol = this->MergeAdjacent(layer_hits[i]);

		for (auto t : tmpcol)
		{
			if (t.size() >= 5)
                continue;
			else
			{
				HTTool::hcl* tmp_hcl = new HTTool::hcl();

				Double_t tmp_x = 0.0;
				Double_t tmp_y = 0.0;
				Double_t tmp_z = 0.0;
				Double_t tmp_e = 0.0;

				for (auto z : t)
				{
					tmp_x += z->x * z->e;
					tmp_y += z->y * z->e;
					tmp_z += z->z * z->e;
					tmp_e += z->e;
				}

				tmp_x /= tmp_e;
				tmp_y /= tmp_e;
				tmp_z /= tmp_e;

				tmp_hcl->x = tmp_x;
				tmp_hcl->y = tmp_y;
				tmp_hcl->z = tmp_z;
				tmp_hcl->e = tmp_e;

				hcx.emplace_back(tmp_x);
				hcy.emplace_back(tmp_y);
				hcz.emplace_back(tmp_z);
				hce.emplace_back(tmp_e);

				if (tmp_e > 0.05)
                    Hough_Cluster.emplace_back(tmp_hcl);
			}
		}
	}

	Int_t nbineta = 100;
	hht = new TH2D("hht", "test", nbineta, 0.0, 3.14, 160, -10, 150);

	for (auto hh : Hough_Cluster)
	{
		for (Int_t itheta = 0; itheta < nbineta; itheta++)
        {
            Double_t theta = TMath::Pi() * itheta / Double_t(nbineta);
            Double_t rho = 0.1 * hh->z * TMath::Cos(theta) + 0.1 * hh->x * TMath::Sin(theta);

            hht->Fill(theta, rho);
        }
	}

	map<std::pair<Int_t, Int_t>, Int_t> map_track;

	auto ftrack = [&map_track, this] (pair<Int_t, Int_t> pos)->Bool_t
	{
		for (auto i : map_track)
		{
			pair<Int_t, Int_t> posi = i.first;
			Double_t d = sqrt(pow(posi.first - pos.first, 2) + pow(posi.second - pos.second, 2));
			if (d < 10)
			{
				if (this->hht->GetBinContent(posi.first + 1, posi.second + 1) < this->hht->GetBinContent(pos.first + 1, pos.second + 1))
				{
					map_track.erase(map_track.find(posi));
					map_track.insert(pair<pair<Int_t, Int_t>, Int_t>(pair<Int_t, Int_t>(pos.first, pos.second), 1));
				}
				return true;
			}
		}
		return false;
	};

	for (Int_t i = 0; i < hht->GetNbinsX(); i++)
		for (Int_t j = 0; j < hht->GetNbinsY(); j++)
			if (hht->GetBinContent(i + 1, j + 1) > 10)
				if (!ftrack(pair<Int_t, Int_t>(i, j)))
					map_track.insert(pair<pair<Int_t, Int_t>, Int_t>(pair<Int_t, Int_t>(i, j), 1));

	ntrack = map_track.size();

	//TCanvas* c1 = new TCanvas("c1", "test", 1024, 768);
	//c1->cd();
	//gStyle->SetOptStat("");
	//hht->Draw("colz");
	//c1->SaveAs("ht.png");
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

	auto fadj = [] (const HTTool::hit* _hit, const vector<HTTool::hit*> _hits)->Bool_t
	{
        const Double_t bias = 342.55;
        const Double_t width = 40.3;

        Int_t m = (_hit->y + bias) / width;
        Int_t n = (_hit->x + bias) / width;

		for (auto i : _hits)
		{
            Int_t im = (i->y + bias) / width;
            Int_t in = (i->x + bias) / width;

			if ((m == im && abs(n - in) == 1) || (n == in && abs(m - im) == 1) || (abs(m - im) == 1 && abs(n - in) == 1 && true))
				return true;
		}
		return false;
	};

	while (hits.size() > 0)
	{
		vechit.emplace_back(vector<HTTool::hit*>({hits.at(0)}));
		hits.erase(hits.begin());
		vector<Int_t> hit2erase;

		while(1)
		{
			vector<Int_t>().swap(hit2erase);

			for (Int_t i = 0; i < hits.size(); i++)
				if (fadj(hits.at(i), vechit.back()))
				{
					vechit.back().emplace_back(hits.at(i));
					hit2erase.emplace_back(i);
				}

			if (hit2erase.size() == 0)
                break;
			else
				for (Int_t i = (hit2erase.size() - 1); i >= 0; i--)
					hits.erase(hits.begin() + i);
		}
	}
	return vechit;
}
