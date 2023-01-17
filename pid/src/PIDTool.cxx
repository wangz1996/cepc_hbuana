#include "PIDTool.h"

PIDTool::PIDTool()
{
}

PIDTool::~PIDTool()
{
}

int PIDTool::GenNtuple(const string &file,const string &tree)
{
	ROOT::RDataFrame *dm=new ROOT::RDataFrame(tree,file);
	string outname = file;
    outname = outname.substr(outname.find_last_of('/')+1);
	outname = "pid_"+outname;
	auto fout = dm->Define("xwidth",[](vector<double> Hit_X)
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
	.Define("layer_hitcell",[](vector<int> Hit_PSDID)
	{
		vector<int> layer_HitCell(40);
		unordered_map<int,int> map_cell_hit; // PSDID, Count
		for(auto i:Hit_PSDID)
		{
			int layer = i/10000;
			if(map_cell_hit.count(i)==0)
			{
				layer_HitCell.at(layer)++;
			}
			map_cell_hit[i]++;
		}
		return layer_HitCell;
	},{"Hit_PSDID"})
	.Define("shower_start",[](vector<int> layer_hitcell)
	{
		int shower_start=38;
		for(int i=0;i<layer_hitcell.size()-3;i++)
		{
			if(layer_hitcell.at(i)>=4 && layer_hitcell.at(i+1)>=4 && layer_hitcell.at(i+2)>=4 && layer_hitcell.at(i+3)>=4)
			{
				shower_start = i;
				break;
			}
		}
		return shower_start;
	},{"layer_hitcell"})
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
			int layer = Hit_PSDID.at(i)/10000;
			h.at(layer)->Fill(Hit_X.at(i));
		}
		for(int i=0;i<h.size();i++)
		{
			layer_xwidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_xwidth;
	},{"Hit_X","Hit_PSDID"})
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
			int layer = Hit_PSDID.at(i)/10000;
			h.at(layer)->Fill(Hit_Y.at(i));
		}
		for(int i=0;i<h.size();i++)
		{
			layer_ywidth.at(i) = h.at(i)->GetRMS();
			delete h.at(i);
		}
		vector<TH1D*>().swap(h);
		return layer_ywidth;
	},{"Hit_Y","Hit_PSDID"})
	.Define("shower_layer",[](vector<int> Hit_PSDID,vector<double> layer_xwidth,vector<double> layer_ywidth)
	{
		double shower_layer=0;
		for(int i=0;i<40;i++)
		{
			if(layer_xwidth.at(i)>50 && layer_ywidth.at(i)>50)
			{
				shower_layer++;
			}
		}
		return shower_layer;
	},{"Hit_PSDID","layer_xwidth","layer_ywidth"})
	.Define("hit_layer",[](vector<int> Hit_PSDID)
	{
		double hit_layer=0;
		unordered_map<int,int> map_layer_hit;
		for(auto i:Hit_PSDID)
		{
			int layer = i/10000;
			map_layer_hit[layer]++;
		}
		for(int i=0;i<40;i++)
		{
			if(map_layer_hit.count(i)>0)hit_layer++;
		}
		return hit_layer;
	},{"Hit_PSDID"})
	.Define("shower_layer_ratio","shower_layer/hit_layer")
    .Snapshot(tree,outname);
	delete dm;
	return 1;
}

int PIDTool::TrainBDT(const string &bdtname)
{
   TMVA::Tools::Instance();

   std::map<std::string,int> Use;

   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost

	vector<TTree*> tsig;
	vector<TTree*> tbkg;
	for(auto i:sig)
	{
		TFile *f=TFile::Open(i.first,"READ");
		TTree *t=(TTree*)f->Get(i.second);
		tsig.emplace_back(t);
	}
	for(auto i:bkg)
	{
		TFile *f=TFile::Open(i.first,"READ");
		TTree *t=(TTree*)f->Get(i.second);
		tbkg.emplace_back(t);
	}
   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA_"+bdtname+".root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_"+TString(bdtname));
	for(auto i:var)
	{
		dataloader->AddVariable(i.first,i.second);
	}
   //dataloader->AddVariable( "xwidth", 'D' );
   //dataloader->AddVariable( "ywidth", 'D' );
   //dataloader->AddVariable( "zwidth", 'D' );
   //dataloader->AddVariable( "Edep ", 'D' );

   // You can add an arbitrary number of signal or background trees
   for(auto i:tsig)dataloader->AddSignalTree(i,1.0);
	for(auto i:tbkg)dataloader->AddBackgroundTree(i,1.0);

   dataloader->SetBackgroundWeightExpression( "1" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

	int nsig=0;for(auto i:tsig)nsig+=i->GetEntries();
	int nbkg=0;for(auto i:tbkg)nbkg+=i->GetEntries();
	int nTrain_Signal = 0.7*nsig;
	int nTrain_Bkg = 0.7*nbkg;
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal="+TString(to_string(nTrain_Signal))+":nTrain_Background="+TString(to_string(nTrain_Bkg))+":SplitMode=Random:NormMode=NumEvents:!V" );

   // Cut optimisation
   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   // --------------------------------------------------------------------------------------------------

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
