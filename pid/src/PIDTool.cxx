#include "PIDTool.h"

Int_t NHScaleV2(RVec<Double_t> const& pos_x, RVec<Double_t> const& pos_y, RVec<Double_t> const& pos_z, Int_t const& RatioX, Int_t const& RatioY, Int_t const& RatioZ)
{
    Int_t ReScaledNH = 0;
    Int_t tmpI = 0;
    Int_t tmpJ = 0;
    Int_t tmpK = 0;
    Double_t tmpEn = 0;
    Int_t NewCellID0 = 0;
    const Int_t NumHit = pos_x.size();

    std::map<Double_t, Double_t> testIDtoEnergy;

    for (Int_t i = 0; i < NumHit; i++)
    {
        Double_t x = pos_x.at(i);
        Double_t y = pos_y.at(i);
        Double_t z = pos_z.at(i);

        tmpI = (Int_t(x / 40.2996) + Int_t(fabs(x) / x)) / RatioX;
        tmpJ = (Int_t(y / 39.9874) + Int_t(fabs(y) / y)) / RatioY;
        tmpK = Int_t(z / 26) / RatioZ;
        tmpEn = 1;

        NewCellID0 = (tmpK << 24) + (tmpJ << 12) + tmpI;

        if (testIDtoEnergy.find(NewCellID0) == testIDtoEnergy.end())
            testIDtoEnergy[NewCellID0] = tmpEn;
        else
            testIDtoEnergy[NewCellID0] += tmpEn;
    }

    ReScaledNH = testIDtoEnergy.size();
    return ReScaledNH;
}

PIDTool::PIDTool()
{
}

PIDTool::~PIDTool()
{
}

int PIDTool::GenNtuple(const string &file,const string &tree)
{
	//ROOT::EnableImplicitMT();
	ROOT::DisableImplicitMT();
	ROOT::RDataFrame *dm=new ROOT::RDataFrame(tree,file);
	string outname = file;
    outname = outname.substr(outname.find_last_of('/')+1);
	outname = "pid_"+outname;
	auto fout = dm
	->Define("xwidth",[](vector<double> Hit_X)
    {
        TH1D *h1=new TH1D("h1","test",100,-400,400);
        for(auto i:Hit_X)h1->Fill(i);
        double xwidth = h1->GetRMS();
        delete h1;
        return xwidth;
    },{"Hit_X"})
    .Define("ywidth",[](vector<double> Hit_Y){
        TH1D *h1=new TH1D("h1","test",100,-400,400);
        for(auto i:Hit_Y)h1->Fill(i);
        double ywidth = h1->GetRMS();
        delete h1;
        return ywidth;
    },{"Hit_Y"})
    .Define("zwidth",[](vector<double> Hit_Z){
        TH1D *h1=new TH1D("h1","test",100,-0,1200);
        for(auto i:Hit_Z)h1->Fill(i);
        double zwidth = h1->GetRMS();
        delete h1;
        return zwidth;
    },{"Hit_Z"})
    .Define("Edep",[](vector<double> Hit_Energy){
        double sum=0.;
        for(auto i:Hit_Energy)sum+=i;
        return sum;
    },{"Hit_Energy"})
	.Define("layer_hitcell",[](vector<int> CellID,vector<double> Hit_Energy)
	{
		vector<int> layer_HitCell(40);
		for(int i=0;i<CellID.size();i++)
		{
			int layer = CellID.at(i)/100000;
			if(Hit_Energy.at(i)<0.1)continue;
			layer_HitCell.at(layer)++;
		}
		return layer_HitCell;
	},{"CellID","Hit_Energy"})
	.Define("layer_rms",[](vector<double> Hit_X,vector<double> Hit_Y,vector<int> CellID,vector<double> Hit_Energy)->vector<double>
	{
		vector<double> layer_rms(40);
		vector<TH2D*> hvec;
		for(int i=0;i<40;i++)hvec.emplace_back(new TH2D("h"+TString(to_string(i))+"_rms","Layer RMS",100,-400,400,100,-400,400));
		for(int i=0;i<Hit_X.size();i++)
		{
			int layer = CellID.at(i)/100000;
			hvec.at(layer)->Fill(Hit_X.at(i),Hit_Y.at(i),Hit_Energy.at(i));
		}
		for(int i=0;i<hvec.size();i++)
		{
			if(hvec.at(i)->GetEntries()<4)
			{
				layer_rms.at(i) = 0.;
			}
			else
			{
				layer_rms.at(i) = hvec.at(i)->GetRMS();
			}
			delete hvec.at(i);
		}
		vector<TH2D*>().swap(hvec);
		return layer_rms;
	},{"Hit_X","Hit_Y","CellID","Hit_Energy"})
	.Define("shower_start",[](vector<int> layer_hitcell,vector<double> layer_rms)
	{
		int shower_start=42;
		for(int i=0;i<layer_hitcell.size()-3;i++)
		{
			if(layer_hitcell.at(i)>=4 && layer_rms.at(i)<50. && layer_hitcell.at(i+1)>=4 && layer_hitcell.at(i+2)>=4 && layer_hitcell.at(i+3)>=4)
			{
				shower_start = i;
				break;
			}
		}
		return shower_start;
	},{"layer_hitcell","layer_rms"})
	.Define("shower_end",[](vector<int> layer_hitcell,vector<double> layer_rms,int shower_start)
	{
		int shower_end=42;
		if(shower_start==42)return shower_end;
		for(int i=shower_start;i<layer_hitcell.size()-3;i++)
		{
			if(layer_hitcell.at(i)<4 && layer_hitcell.at(i+1)<4)
			{
				shower_end = i;
				break;
			}
		}
        if (shower_end == 42)
            shower_end = 40;
		return shower_end;
	},{"layer_hitcell","layer_rms","shower_start"})
	.Define("layer_xwidth",[](vector<double> Hit_X,vector<int> Hit_PSDID)
	{
		vector<double> layer_xwidth(40);
		vector<TH1D*> h;
		for(int l=0;l<40;l++)
		{
			h.emplace_back(new TH1D(TString("h")+TString(to_string(l)),"test",100,-400,400));
		}
		for(int i=0;i<Hit_PSDID.size();i++)
		{
			int layer = Hit_PSDID.at(i)/100000;
			h.at(layer)->Fill(Hit_X.at(i));
		}
		for(int i=0;i<h.size();i++)
		{
			layer_xwidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_xwidth;
	},{"Hit_X","CellID"})
	.Define("layer_ywidth",[](vector<double> Hit_Y,vector<int> Hit_PSDID)
	{
		vector<double> layer_ywidth(40);
		vector<TH1D*> h;
		for(int l=0;l<40;l++)
		{
			h.emplace_back(new TH1D(TString("h")+TString(to_string(l)),"test",100,-400,400));
		}
		for(int i=0;i<Hit_PSDID.size();i++)
		{
			int layer = Hit_PSDID.at(i)/100000;
			h.at(layer)->Fill(Hit_Y.at(i));
		}
		for(int i=0;i<h.size();i++)
		{
			layer_ywidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_ywidth;
	},{"Hit_Y","CellID"})
	.Define("shower_layer",[](vector<int> Hit_PSDID,vector<double> layer_xwidth,vector<double> layer_ywidth)
	{
		double shower_layer=0;
		for(int i=0;i<40;i++)
		{
			if(layer_xwidth.at(i)>60 && layer_ywidth.at(i)>60)
			{
				shower_layer++;
			}
		}
		return shower_layer;
	},{"CellID","layer_xwidth","layer_ywidth"})
	.Define("hit_layer",[](vector<int> Hit_PSDID)
	{
		double hit_layer=0;
		unordered_map<int,int> map_layer_hit;
		for(auto i:Hit_PSDID)
		{
			int layer = i/100000;
			map_layer_hit[layer]++;
		}
		for(int i=0;i<40;i++)
		{
			if(map_layer_hit.count(i)>0)hit_layer++;
		}
		return hit_layer;
	},{"CellID"})
	.Define("shower_layer_ratio","shower_layer/hit_layer")
	.Define("shower_density",[](vector<int> CellID,vector<double> Hit_Energy,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z)
	{
		double shower_density=0.;
		unordered_map<int,int> map_CellID;
		double xend = -342.551;
		double xgap = 40.2996;
		double yend = -340.557;
		double ygap = 39.9874;
		for(int i=0;i<Hit_X.size();i++)
		{
			int ix = round((Hit_X.at(i)-xend)/xgap);
        	int iy = round((Hit_Y.at(i)-yend)/ygap);
			int layer = CellID.at(i)/100000;
			int cellid = layer*10000+ix*100 + iy;
			int memo = (CellID.at(i)%10000)/100;
			if(memo!=0)continue;
			map_CellID[cellid]=1;
		}
		int nhit = map_CellID.size();
		for(auto i:map_CellID)
		{
			int layer = i.first/10000;
			int x = (i.first%10000)/100;
			int y = i.first%100;
			for(int il=layer-1;il<=layer+1;il++)
			{
				if(il<0 || il>39)continue;
				for(int ix=x-1;ix<=x+1;ix++)
				{
					if(ix<0 || ix>17)continue;
					for(int iy=y-1;iy<=y+1;iy++)
					{
						if(iy<0 || iy>17)continue;
						int tmp = il*10000+ix*100+iy;
						if(map_CellID.count(tmp))shower_density++;
						//if(1)shower_density++;
					}
				}
			}
		}
		shower_density/=nhit;
		//cout<<"DEBUG: "<<shower_density<<endl;
		return shower_density;
	},{"CellID","Hit_Energy","Hit_X","Hit_Y","Hit_Z"})
	.Define("shower_length",[](vector<double> layer_rms,int shower_start,int shower_end)
	{
		double shower_length=0.;
		double startz = 300. + 1.5+shower_start*25.;
		int max_layer =0;
		double max_rms=0.;
		//for(int i=0;i<layer_rms.size();i++)
		//{
		//	if(layer_rms.at(i)>max_rms)
		//	{
		//		max_layer=i;
		//		max_rms=layer_rms.at(i);
		//	}
		//}
		//auto maxPosition = max_element(layer_rms.begin()+shower_start, layer_rms.end());
		//int max_layer= maxPosition - layer_rms.begin();
		shower_length = (shower_end-shower_start)*25.;
		if(shower_start==42)shower_length=0.;
		return shower_length;
	},{"layer_rms","shower_start","shower_end"})
    .Define("shower_radius", [] (vector<Double_t> hit_x, vector<Double_t> hit_y, vector<Double_t> hit_z, Int_t beginning, Int_t ending)
    {
        Double_t d2 = 0;
        Int_t hits = 0;
        const Int_t n = hit_x.size();
        for (Int_t i = 0; i < n; i++)
        {
            if (hit_z.at(i) / 25.0 >= beginning && hit_z.at(i) / 25.0 < ending)
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
    .Define("FD_2D", [] (RVec<Double_t> const& pos_x, RVec<Double_t> const& pos_y, RVec<Double_t> const& pos_z)
    {
        Double_t fd = 0;
        const Int_t num = 10;
        const Int_t nhit = pos_x.size();
        Int_t NResizeHit[num] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        Int_t scale[num] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 20 };
        for (Int_t i = 0; i < num; i++)
        {
            NResizeHit[i] = NHScaleV2(pos_x, pos_y, pos_z, scale[i], scale[i], 1);
            fd += 0.1 * TMath::Log((Double_t)nhit / NResizeHit[i]) / TMath::Log((Double_t)scale[i]);
        }
        if (pos_x.size() == 0)
            fd = -1;
        return fd;
    }, {"Hit_X", "Hit_Y", "Hit_Z"})
	.Define("hclx",[](vector<int> CellID,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z,vector<double> Hit_Energy)
	{
		vector<double> hclx;
		HTTool *httool = new HTTool(CellID,Hit_X,Hit_Y,Hit_Z,Hit_Energy);
		hclx = httool->GetHclX();
		delete httool;
		return hclx;
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy"})
	.Define("hcly",[](vector<int> CellID,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z,vector<double> Hit_Energy)
	{
		vector<double> hcly;
		HTTool *httool = new HTTool(CellID,Hit_X,Hit_Y,Hit_Z,Hit_Energy);
		hcly = httool->GetHclY();
		delete httool;
		return hcly;
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy"})
	.Define("hclz",[](vector<int> CellID,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z,vector<double> Hit_Energy)
	{
		vector<double> hclz;
		HTTool *httool = new HTTool(CellID,Hit_X,Hit_Y,Hit_Z,Hit_Energy);
		hclz = httool->GetHclZ();
		delete httool;
		return hclz;
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy"})
	.Define("hcle",[](vector<int> CellID,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z,vector<double> Hit_Energy)
	{
		vector<double> hcle;
		HTTool *httool = new HTTool(CellID,Hit_X,Hit_Y,Hit_Z,Hit_Energy);
		hcle = httool->GetHclE();
		delete httool;
		return hcle;
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy"})
	.Define("ntrack",[](vector<int> CellID,vector<double> Hit_X,vector<double> Hit_Y,vector<double> Hit_Z,vector<double> Hit_Energy)
	{
		int ntrack=0;
		HTTool *httool = new HTTool(CellID,Hit_X,Hit_Y,Hit_Z,Hit_Energy);
		ntrack = httool->GetNtrack();
		delete httool;
		return ntrack;
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy"})
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
	Float_t        bdt_shower_radius;
	Float_t        bdt_shower_length;
	Float_t		bdt_ntrack;
	reader->AddVariable("Edep",&bdt_Edep);
	reader->AddVariable("ntrack",&bdt_ntrack);
	reader->AddVariable("shower_density",&bdt_shower_density);
	reader->AddVariable("shower_layer_ratio",&bdt_shower_layer_ratio);
	reader->AddVariable("shower_length",&bdt_shower_length);
	reader->AddVariable("shower_radius",&bdt_shower_radius);
	reader->AddVariable("shower_start",&bdt_shower_start);
	reader->AddVariable("xwidth",&bdt_xwidth);
	reader->AddVariable("ywidth",&bdt_ywidth);
	reader->AddVariable("zwidth",&bdt_zwidth);

	reader->BookMVA("BDTG method",TString("/lustre/collider/wangzhen/cepc/ahcal/cepc_hbuana/run/pid/dataset/weights/TMVAMulticlass_BDTG.weights.xml"));
	cout<<"Booked"<<endl;
	vector<string> rdf_input={"Edep","ntrack","shower_density","shower_layer_ratio","shower_length","shower_radius","shower_start","xwidth","ywidth","zwidth"};
	ROOT::RDataFrame df(tname,fname);
	auto bdtout = df
	.Define("BDT_pi_plus",[&](double e,int n,double d,double lr,double l,double r,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_radius = r;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[0];
	},rdf_input)
	.Define("BDT_e_plus",[&](double e,int n,double d,double lr,double l,double r,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_radius = r;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[1];
	},rdf_input)
	.Define("BDT_mu_plus",[&](double e,int n,double d,double lr,double l,double r,int s,double x,double y,double z)
	{
		bdt_xwidth = x;
		bdt_ywidth = y;
		bdt_zwidth = z;
		bdt_Edep   = e;
		bdt_shower_start = s;
		bdt_shower_layer_ratio = lr;
		bdt_shower_density = d;
		bdt_shower_radius = r;
		bdt_shower_length = l;
		bdt_ntrack = n;
		return (reader->EvaluateMulticlass( "BDTG method" ))[2];
	},rdf_input)
	//.Define("bdt_proton",[&](double e,int n,double d,double lr,double l,int s,double x,double y,double z)
	//{
	//	bdt_xwidth = x;
	//	bdt_ywidth = y;
	//	bdt_zwidth = z;
	//	bdt_Edep   = e;
	//	bdt_shower_start = s;
	//	bdt_shower_layer_ratio = lr;
	//	bdt_shower_density = d;
	//	bdt_shower_length = l;
	//	bdt_ntrack = n;
	//	return (reader->EvaluateMulticlass( "BDTG method" ))[0];
	//},rdf_input)
	.Snapshot(tname,outname);
	return 1;
}
