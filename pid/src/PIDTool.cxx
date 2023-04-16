#include "PIDTool.h"

PIDTool::PIDTool()
{
}

PIDTool::~PIDTool()
{
}

int PIDTool::GenNtuple(const string &file,const string &tree)
{
    const Int_t nlayer = 40;
    const Int_t thick = 25;
	//ROOT::EnableImplicitMT();
	ROOT::DisableImplicitMT();
	ROOT::RDataFrame *dm=new ROOT::RDataFrame(tree,file);
	string outname = file;
    outname = outname.substr(outname.find_last_of('/')+1);
	outname = "pid_"+outname;
	auto fout = dm->Define("xwidth", [] (vector<Double_t> Hit_X)
    {
        TH1D* h1 = new TH1D("h1", "", 100, -400, 400);
        for (Double_t i : Hit_X)
            h1->Fill(i);
        Double_t xwidth = h1->GetRMS();
        delete h1;
        return xwidth;
    }, {"Hit_X"})
    .Define("ywidth", [] (vector<Double_t> Hit_Y)
    {
        TH1D* h1 = new TH1D("h1", "", 100, -400, 400);
        for (Double_t i : Hit_Y)
            h1->Fill(i);
        Double_t ywidth = h1->GetRMS();
        delete h1;
        return ywidth;
    }, {"Hit_Y"})
    .Define("zwidth", [] (vector<Double_t> Hit_Z)
    {
        TH1D* h1 = new TH1D("h1", "", 100, 0, 1200);
        for (Double_t i : Hit_Z)
            h1->Fill(i);
        Double_t zwidth = h1->GetRMS();
        delete h1;
        return zwidth;
    }, {"Hit_Z"})
    .Define("Edep", [] (vector<Double_t> Hit_Energy)
    {
        Double_t sum = 0.0;
        for (Double_t i : Hit_Energy)
            sum += i;
        return sum;
    }, {"Hit_Energy"})
	.Define("layer_hitcell", [] (vector<Int_t> CellID, vector<Double_t> Hit_Energy)
	{
		vector<Int_t> layer_HitCell(nlayer);
		for(Int_t i = 0; i < CellID.size(); i++)
		{
			Int_t layer = CellID.at(i) / 100000;
			if (Hit_Energy.at(i) < 0.1)
                continue;
			layer_HitCell.at(layer)++;
		}
		return layer_HitCell;
	}, {"CellID", "Hit_Energy"})
	.Define("layer_rms", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> CellID, vector<Double_t> Hit_Energy)->vector<Double_t>
	{
		vector<Double_t> layer_rms(nlayer);
		vector<TH2D*> hvec;
		for (Int_t i = 0; i < nlayer; i++)
            hvec.emplace_back(new TH2D("h" + TString(to_string(i)) + "_rms", "Layer RMS", 100, -400, 400, 100, -400, 400));
		for (Int_t i = 0; i < Hit_X.size(); i++)
		{
            Int_t layer = CellID.at(i) / 100000;
            hvec.at(layer)->Fill(Hit_X.at(i), Hit_Y.at(i), Hit_Energy.at(i));
		}
		for (Int_t i = 0; i < hvec.size(); i++)
		{
			if (hvec.at(i)->GetEntries() < 4)
				layer_rms.at(i) = 0.;
			else
				layer_rms.at(i) = hvec.at(i)->GetRMS();
			delete hvec.at(i);
		}
		vector<TH2D*>().swap(hvec);
		return layer_rms;
	}, {"Hit_X", "Hit_Y", "CellID", "Hit_Energy"})
	.Define("shower_start", [] (vector<Int_t> layer_hitcell, vector<Double_t> layer_rms)
	{
		Int_t shower_start = nlayer + 2;
		for (Int_t i = 0; i < layer_hitcell.size() - 3; i++)
		{
			if (layer_hitcell.at(i) >= 4 && layer_rms.at(i) < 50.0 && layer_hitcell.at(i + 1) >= 4 && layer_hitcell.at(i + 2) >= 4 && layer_hitcell.at(i + 3) >= 4)
			{
				shower_start = i;
				break;
			}
		}
		return shower_start;
	}, {"layer_hitcell", "layer_rms"})
	.Define("shower_end", [] (vector<Int_t> layer_hitcell, vector<Double_t> layer_rms, Int_t shower_start)
	{
		Int_t shower_end = nlayer + 2;
		if (shower_start == nlayer + 2)
            return shower_end;
		for (Int_t i = shower_start; i < layer_hitcell.size() - 3; i++)
		{
			if (layer_hitcell.at(i) < 4 && layer_hitcell.at(i + 1) < 4)
			{
				shower_end = i;
				break;
			}
		}
        if (shower_end == nlayer + 2)
            shower_end = nlayer;
		return shower_end;
	}, {"layer_hitcell", "layer_rms", "shower_start"})
    .Define("radius_shower", [] (vector<Double_t> hit_x, vector<Double_t> hit_y, vector<Double_t> hit_z, Int_t beginning, Int_t ending)
    {
        LongDouble_t d2 = 0;
        Int_t hits = 0;
        const Int_t n = hit_x.size();
        for (Int_t i = 0; i < n; i++)
        {
            if ((hit_z.at(i) - 1.5) / thick >= beginning && (hit_z.at(i) - 1.5) / thick < ending)
            {
                hits++;
                d2 += TMath::Power(hit_x.at(i), 2) + TMath::Power(hit_y.at(i), 2);
            }
        }
        if (hits == 0)
            return 0.0;
        else
        {
            Double_t radius = TMath::Sqrt(d2 / hits);
            return radius;
        }
    }, {"Hit_X", "Hit_Y", "Hit_Z", "shower_start", "shower_end"})
	.Define("layer_xwidth", [] (vector<Double_t> Hit_X, vector<Int_t> Hit_PSDID)
	{
		vector<Double_t> layer_xwidth(nlayer);
		vector<TH1D*> h;
		for (Int_t l = 0; l < nlayer; l++)
			h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -400, 400));
		for (Int_t i = 0; i < Hit_PSDID.size(); i++)
		{
			Int_t layer = Hit_PSDID.at(i) / 100000;
			h.at(layer)->Fill(Hit_X.at(i));
		}
		for (Int_t i = 0; i < h.size(); i++)
		{
			layer_xwidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_xwidth;
	}, {"Hit_X", "CellID"})
	.Define("layer_ywidth", [] (vector<Double_t> Hit_Y, vector<Int_t> Hit_PSDID)
	{
		vector<Double_t> layer_ywidth(nlayer);
		vector<TH1D*> h;
		for (Int_t l = 0; l < nlayer; l++)
			h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -400, 400));
		for (Int_t i = 0; i < Hit_PSDID.size(); i++)
		{
			Int_t layer = Hit_PSDID.at(i) / 100000;
			h.at(layer)->Fill(Hit_Y.at(i));
		}
		for (Int_t i = 0; i < h.size(); i++)
		{
			layer_ywidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_ywidth;
	}, {"Hit_Y", "CellID"})
	.Define("shower_layer", [] (vector<Int_t> Hit_PSDID, vector<Double_t> layer_xwidth, vector<Double_t> layer_ywidth)
	{
		Double_t shower_layer = 0;
		for (Int_t i = 0; i < nlayer; i++)
		{
			if (layer_xwidth.at(i) > 60 && layer_ywidth.at(i) > 60)
				shower_layer++;
		}
		return shower_layer;
	}, {"CellID", "layer_xwidth", "layer_ywidth"})
	.Define("hit_layer", [] (vector<Int_t> Hit_PSDID)
	{
		Double_t hit_layer = 0;
		unordered_map<Int_t, Int_t> map_layer_hit;
		for (Int_t i : Hit_PSDID)
		{
			Int_t layer = i / 100000;
			map_layer_hit[layer]++;
		}
		for (Int_t i = 0; i < nlayer; i++)
		{
			if (map_layer_hit.count(i) > 0)
                hit_layer++;
		}
		return hit_layer;
	}, {"CellID"})
	.Define("shower_layer_ratio", "shower_layer / hit_layer")
	.Define("shower_density", [] (vector<Int_t> CellID, vector<Double_t> Hit_Energy)
	{
		Double_t shower_density = 0.0;
		unordered_map<Int_t, Int_t> map_CellID;
		for (Int_t i : CellID)
			map_CellID[i] = 1;
		for(Int_t i = 0; i < CellID.size(); i++)
		{
			if (Hit_Energy.at(i) < 0.1)
                continue;
			Int_t layer = CellID.at(i) / 100000;
			Int_t x = (CellID.at(i) % 10000) / 100;
			Int_t y = CellID.at(i) % 100;
			for (Int_t _l = layer - 1; _l <= layer+1; _l++)
			{
				for (Int_t _x = x-1; _x <= x + 1; _x++)
				{
					for (Int_t _y = y-1; _y <= y+1; _y++)
					{
						Int_t _id = _l * 10000 + _x * 100 + _y;
						shower_density += map_CellID[_id];
					}
				}
			}
		}
		shower_density /= CellID.size();
		return shower_density;
	}, {"CellID", "Hit_Energy"})
	.Define("shower_length", [] (vector<Double_t> layer_rms, Int_t shower_start)
	{
		Double_t shower_length = 0.0;
		Double_t startz = 300.0 + 1.5 + shower_start * thick;
		Int_t max_layer = 0;
		Double_t max_rms = 0.0;
		for (Int_t i = 0; i < layer_rms.size(); i++)
		{
			if (layer_rms.at(i) > max_rms)
			{
				max_layer = i;
				max_rms = layer_rms.at(i);
			}
		}
		//auto maxPosition = max_element(layer_rms.begin() + shower_start, layer_rms.end());
		//Int_t max_layer = maxPosition - layer_rms.begin();
		shower_length = (max_layer - shower_start) * thick;
		if (shower_start == nlayer + 2)
            shower_length = 0.0;
		return shower_length;
	}, {"layer_rms", "shower_start"})
	.Define("hclx", [] (vector<Int_t> CellID, vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
	{
		vector<Double_t> hclx;
		HTTool *httool = new HTTool(CellID, Hit_X, Hit_Y, Hit_Z, Hit_Energy);
		hclx = httool->GetHclX();
		delete httool;
		return hclx;
	}, {"CellID", "Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
	.Define("hcly", [] (vector<Int_t> CellID, vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
	{
		vector<Double_t> hcly;
		HTTool* httool = new HTTool(CellID, Hit_X, Hit_Y, Hit_Z, Hit_Energy);
		hcly = httool->GetHclY();
		delete httool;
		return hcly;
	}, {"CellID", "Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
	.Define("hclz", [] (vector<Int_t> CellID, vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
	{
		vector<Double_t> hclz;
		HTTool* httool = new HTTool(CellID, Hit_X, Hit_Y, Hit_Z, Hit_Energy);
		hclz = httool->GetHclZ();
		delete httool;
		return hclz;
	}, {"CellID", "Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
	.Define("hcle", [](vector<Int_t> CellID, vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
	{
		vector<Double_t> hcle;
		HTTool* httool = new HTTool(CellID, Hit_X, Hit_Y, Hit_Z, Hit_Energy);
		hcle = httool->GetHclE();
		delete httool;
		return hcle;
	}, {"CellID", "Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
	.Define("ntrack", [] (vector<Int_t> CellID, vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Hit_Energy)
	{
		Int_t ntrack = 0;
		HTTool* httool = new HTTool(CellID, Hit_X, Hit_Y, Hit_Z, Hit_Energy);
		ntrack = httool->GetNtrack();
		delete httool;
		return ntrack;
	}, {"CellID", "Hit_X", "Hit_Y", "Hit_Z", "Hit_Energy"})
    /*
    .Define("time_mean_hit", [] (vector<Double_t> hit_time)
    {
        Double_t tot = 0;
        const Int_t hits = hit_time.size();
        for (Double_t& i : hit_time)
            tot += i;
        return tot / hits;
    }, {"Hit_Time"})
    .Define("time_rms_hit", [] (vector<Double_t> hit_time)
    {
        Double_t tot2 = 0;
        const Int_t hits = hit_time.size();
        for (Double_t& i : hit_time)
            tot2 += i * i;
        return TMath::Sqrt(tot2 / hits);
    }, {"Hit_Time"})
    .Define("time_mean_shower", [] (vector<Double_t> hit_time, vector<Double_t> hit_z, Int_t beginning, Int_t ending)
    {
        Double_t tot = 0;
        Int_t hits = 0;
        const Int_t n = hit_time.size();
        for (Int_t i = 0; i < n; i++)
        {
            if (hit_z.at(i) / thick >= beginning && hit_z.at(i) / thick < ending)
            {
                hits++;
                tot += hit_time.at(i), 2;
            }
        }
        if (hits == 0)
            return 0.0;
        else
            return tot / hits;
    }, {"Hit_Time", "shower_start", "shower_end"})
    .Define("time_rms_shower", [] (vector<Double_t> hit_time, vector<Double_t> hit_z, Int_t beginning, Int_t ending)
    {
        Double_t tot2 = 0;
        Int_t hits = 0;
        const Int_t n = hit_time.size();
        for (Int_t i = 0; i < n; i++)
        {
            if (hit_z.at(i) / thick >= beginning && hit_z.at(i) / thick < ending)
            {
                hits++;
                tot2 += TMath::Power(hit_time.at(i), 2);
            }
        }
        if (hits == 0)
            return 0.0;
        else
            return tot2 / hits;
    }, {"Hit_Time", "shower_start", "shower_end"})
    */
	//.Range(1)
    .Snapshot(tree,outname);
	delete dm;
	return 1;
}

int PIDTool::TrainBDT()
{
   TMVA::Tools::Instance();

	unordered_map<TTree*,TString> tsignal;
	for(auto i:signal)
	{
		TFile *f=TFile::Open(i.first.first,"READ");
		TTree *t=(TTree*)f->Get(i.first.second);
		tsignal.insert(pair<TTree*,TString>(t,i.second));
	}
   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVAMulticlass.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
	for(auto i:var)
	{
		dataloader->AddVariable(i.first,i.second);
	}
   // You can add an arbitrary number of signal or background trees
   for(auto i:tsignal)dataloader->AddTree(i.first,i.second);
   dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );

	factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   //if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
	
}

int PIDTool::BDTNtuple(const string &fname,const string &tname)
{
	ROOT::EnableImplicitMT();
	string outname = fname;
	outname = outname.substr(outname.find_last_of('/')+1);
	outname = "bdt_"+outname;
	// This loads the library
	TMVA::Tools::Instance();

	// Default MVA methods to be trained + tested
	std::map<std::string,int> Use;

	// Cut optimisation
    Use["BDTG"]            = 1;
	std::cout << "==> Start TMVAMulticlassApplication" << std::endl;
	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
	Float_t        bdt_xwidth;
	Float_t        bdt_ywidth;
	Float_t        bdt_zwidth;
	Float_t        bdt_Edep;
	Float_t			bdt_shower_start;
	Float_t        bdt_shower_layer_ratio;
	Float_t        bdt_shower_density;
	Float_t        bdt_shower_length;
	Float_t		bdt_ntrack;
	reader->AddVariable("Edep",&bdt_Edep);
	reader->AddVariable("ntrack",&bdt_ntrack);
	reader->AddVariable("shower_density",&bdt_shower_density);
	reader->AddVariable("shower_layer_ratio",&bdt_shower_layer_ratio);
	reader->AddVariable("shower_length",&bdt_shower_length);
	reader->AddVariable("shower_start",&bdt_shower_start);
	reader->AddVariable("xwidth",&bdt_xwidth);
	reader->AddVariable("ywidth",&bdt_ywidth);
	reader->AddVariable("zwidth",&bdt_zwidth);

	reader->BookMVA("BDTG method",TString("dataset/weights/TMVAMulticlass_BDTG.weights.xml"));
	cout<<"Booked"<<endl;
	vector<string> rdf_input={"Edep","ntrack","shower_density","shower_layer_ratio","shower_length","shower_start","xwidth","ywidth","zwidth"};
	ROOT::RDataFrame df(tname,fname);
	auto bdtout = df
	.Define("BDT_pi_plus",[&](double e,int n,double d,double lr,double l,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[1];
	},rdf_input)
	.Define("BDT_e_plus",[&](double e,int n,double d,double lr,double l,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[2];
	},rdf_input)
	.Define("BDT_mu_plus",[&](double e,int n,double d,double lr,double l,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[3];
	},rdf_input)
	.Define("bdt_proton",[&](double e,int n,double d,double lr,double l,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[0];
	},rdf_input)
	.Snapshot(tname,outname);
	return 1;
}
