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

        tmpI = (Int_t ((x + 342.55) / 40.3) + Int_t(TMath::Abs(x) / x)) / RatioX;
        tmpJ = (Int_t ((y + 342.55) / 40.3) + Int_t(TMath::Abs(y) / y)) / RatioY;
        tmpK = (Int_t)(z / 30) / RatioZ;
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
    const Int_t nlayer = 40;
    const Int_t thick = 30;
    //ROOT::EnableImplicitMT();
    ROOT::DisableImplicitMT();
    ROOT::RDataFrame *dm = new ROOT::RDataFrame(tree, file);
    string outname = file;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "pid_" + outname;
    auto fout = dm->Define("nhits", "(Int_t) Hit_X.size()")
    .Define("layer", [] (vector<Double_t> Hit_Z, Int_t nhits)
    {
        vector<Int_t> layer;
        for (Int_t i = 0; i < nhits; i++)
            layer.emplace_back((Int_t) Hit_Z.at(i) / thick);
        return layer;
    }, {"Hit_Z", "nhits"})
    .Define("xwidth", [] (vector<Double_t> Hit_X)
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
    .Define("Edep", [] (vector<Double_t> Digi_Hit_Energy)
    {
        Double_t sum = 0;
        for (Double_t i : Digi_Hit_Energy)
            sum += i;
        return sum;
    }, {"Digi_Hit_Energy"})
    .Define("Emean", "Edep / nhits")
    .Define("layer_hitcell", [] (vector<Int_t> layer, vector<Double_t> Digi_Hit_Energy, Int_t nhits)
    {
        vector<Int_t> layer_HitCell(nlayer);
        for (Int_t i = 0; i < nhits; i++)
        {
            Int_t ilayer = layer.at(i);
            if (Digi_Hit_Energy.at(i) < 0.1)
                continue;
            layer_HitCell.at(ilayer)++;
        }
        return layer_HitCell;
    }, {"layer", "Digi_Hit_Energy", "nhits"})
    .Define("layer_rms", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Digi_Hit_Energy, Int_t nhits)->vector<Double_t>
    {
        vector<Double_t> layer_rms(nlayer);
        vector<TH2D*> hvec;
        for (Int_t i = 0; i < nlayer; i++)
            hvec.emplace_back(new TH2D("h" + TString(to_string(i)) + "_rms", "Layer RMS", 100, -400, 400, 100, -400, 400));
        for (Int_t i = 0; i < nhits; i++)
        {
            Int_t ilayer = layer.at(i);
            hvec.at(ilayer)->Fill(Hit_X.at(i), Hit_Y.at(i), Digi_Hit_Energy.at(i));
        }
        for (Int_t i = 0; i < hvec.size(); i++)
        {
            if (hvec.at(i)->GetEntries() < 4)
                layer_rms.at(i) = 0.0;
            else
                layer_rms.at(i) = hvec.at(i)->GetRMS();
            delete hvec.at(i);
        }
        vector<TH2D*>().swap(hvec);
        return layer_rms;
    }, {"Hit_X", "Hit_Y", "layer", "Digi_Hit_Energy", "nhits"})
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
    .Define("shower_radius", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, Int_t beginning, Int_t ending, Int_t nhits)
    {
        LongDouble_t d2 = 0;
        Int_t hits = 0;
        for (Int_t i = 0; i < nhits; i++)
        {
            if (layer.at(i) >= beginning && layer.at(i) < ending)
            {
                hits++;
                d2 += TMath::Power(Hit_X.at(i), 2) + TMath::Power(Hit_Y.at(i), 2);
            }
        }
        if (hits == 0)
            return 0.0;
        Double_t radius = TMath::Sqrt(d2 / hits);
        return radius;
    }, {"Hit_X", "Hit_Y", "layer", "shower_start", "shower_end", "nhits"})
    .Define("layer_xwidth", [] (vector<Double_t> Hit_X, vector<Int_t> layer, Int_t nhits)
    {
        vector<Double_t> layer_xwidth(nlayer);
        vector<TH1D*> h;
        for (Int_t l = 0; l < nlayer; l++)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -400, 400));
        for (Int_t i = 0; i < nhits; i++)
        {
            Int_t ilayer = layer.at(i);
            h.at(ilayer)->Fill(Hit_X.at(i));
        }
        for (Int_t i = 0; i < h.size(); i++)
        {
            layer_xwidth.at(i) = h.at(i)->GetRMS();
            delete h.at(i);
        }
        vector<TH1D*>().swap(h);
        return layer_xwidth;
    }, {"Hit_X", "layer", "nhits"})
    .Define("layer_ywidth", [] (vector<Double_t> Hit_Y, vector<Int_t> layer, Int_t nhits)
    {
        vector<Double_t> layer_ywidth(nlayer);
        vector<TH1D*> h;
        for (Int_t l = 0; l < nlayer; l++)
            h.emplace_back(new TH1D(TString("h") + TString(to_string(l)), "test", 100, -400, 400));
        for (Int_t i = 0; i < nhits; i++)
        {
            Int_t ilayer = layer.at(i);
            h.at(ilayer)->Fill(Hit_Y.at(i));
        }
        for (Int_t i = 0; i < h.size(); i++)
        {
            layer_ywidth.at(i) = h.at(i)->GetRMS();
            delete h.at(i);
        }
        vector<TH1D*>().swap(h);
        return layer_ywidth;
    }, {"Hit_Y", "layer", "nhits"})
    .Define("shower_layer", [] (vector<Double_t> layer_xwidth, vector<Double_t> layer_ywidth)
    {
        Double_t shower_layer = 0;
        for (Int_t i = 0; i < nlayer; i++)
        {
            if (layer_xwidth.at(i) > 60 && layer_ywidth.at(i) > 60)
                shower_layer++;
        }
        return shower_layer;
    }, {"layer_xwidth", "layer_ywidth"})
    .Define("hit_layer", [] (vector<Int_t> layer)
    {
        Double_t hit_layer = 0;
        unordered_map<Int_t, Int_t> map_layer_hit;
        for (Double_t i : layer)
            map_layer_hit[i]++;
        for (Int_t i = 0; i < nlayer; i++)
        {
            if (map_layer_hit.count(i) > 0)
                hit_layer++;
        }
        return hit_layer;
    }, {"layer"})
    .Define("shower_layer_ratio", "shower_layer / hit_layer")
    .Define("shower_density", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Digi_Hit_Energy, Int_t nhits)
    {
        const Double_t bias = 342.55;
        const Double_t width = 40.3;
        Double_t shower_density = 0.0;
        unordered_map<Int_t, Int_t> map_CellID;
        for (Int_t j = 0; j < nhits; j++)
        {
            Int_t x = round((Hit_X.at(j) + bias) / width);
            Int_t y = round((Hit_Y.at(j) + bias) / width);
            Int_t z = layer.at(j);
            Int_t index = z * 100000 + x * 100 + y;
            map_CellID[index] += 1;
        }
        for (Int_t i = 0; i < nhits; i++)
        {
            if (Digi_Hit_Energy.at(i) == 0.0)
                continue;
            Int_t x = round((Hit_X.at(i) + bias) / width);
            Int_t y = round((Hit_Y.at(i) + bias) / width);
            Int_t z = layer.at(i);
            for (Int_t ix = x - 1; ix <= x + 1; ix++)
            {
                if (ix < 0 || ix > 17)
                    continue;
                for (Int_t iy = y - 1; iy <= y + 1; iy++)
                {
                    if (iy < 0 || iy > 17)
                        continue;
                    Int_t tmp = z * 100000 + ix * 100 + iy;
                    shower_density += map_CellID[tmp];
                }
            }
        }
        shower_density /= nhits;
        return shower_density;
    }, {"Hit_X", "Hit_Y", "layer", "Digi_Hit_Energy", "nhits"})
    .Define("clusterE1E9", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Int_t> layer, vector<Double_t> Digi_Hit_Energy, Int_t nhits)
    {
        if (nhits == 0)
            return 0.0;
        const Double_t bias = 342.55;
        const Double_t width = 40.3;
        Double_t clusterE1E9 = 0.0;
        unordered_map<Int_t, Double_t> map_cellid;
        for (Int_t i = 0; i < nhits; i++)
        {
            Int_t x = round((Hit_X.at(i) + bias) / width);
            Int_t y = round((Hit_Y.at(i) + bias) / width);
            Int_t z = layer.at(i);
            Int_t index = z * 100000 + x * 100 + y;
            map_cellid[index] += Digi_Hit_Energy.at(i);
        }
        for (Int_t j = 0; j < nhits; j++)
        {
            if (Digi_Hit_Energy.at(j) == 0.0)
                continue;
            Int_t x = round((Hit_X.at(j) + bias) / width);
            Int_t y = round((Hit_Y.at(j) + bias) / width);
            Int_t z = layer.at(j);
            Int_t index = z * 100000 + x * 100 + y;
            Double_t tempE1 = map_cellid[index];
            Double_t tempE9 = map_cellid[index];
            for (Int_t ix = x - 1; ix <= x + 1; ix++)
            {
                if (ix < 0 || ix > 17)
                    continue;
                for (Int_t iy = y - 1; iy <= y + 1; iy++)
                {
                    if (iy < 0 || iy > 17)
                        continue;
                    Int_t tmp = z * 100000 + ix * 100 + iy;
                    tempE9 += map_cellid[tmp];
                }
            }
            clusterE1E9 += tempE1 / tempE9;
        }
        clusterE1E9 /= nhits;
        return clusterE1E9;
    }, {"Hit_X", "Hit_Y", "layer", "Digi_Hit_Energy", "nhits"})
    .Define("shower_length", [] (vector<Double_t> layer_rms, Int_t shower_start)
    {
        Double_t shower_length = 0.0;
        Double_t startz = shower_start * thick;
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
    .Define("FD_2D", [] (RVec<Double_t> const& pos_x, RVec<Double_t> const& pos_y, RVec<Double_t> const& pos_z, Int_t nhits)
    {
        Double_t fd = 0;
        const Int_t num = 12;
        Int_t NResizeHit[num] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        Int_t scale[num] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30 };
        for (Int_t i = 0; i < num; i++)
        {
            NResizeHit[i] = NHScaleV2(pos_x, pos_y, pos_z, scale[i], scale[i], 1);
            fd += 0.1 * TMath::Log((Double_t) nhits / NResizeHit[i]) / TMath::Log((Double_t)scale[i]);
        }
        if (pos_x.size() == 0)
            fd = -1;
        return fd;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "nhits"})
    .Define("COG_X", [] (vector<Double_t> Hit_X, vector<Double_t> Digi_Hit_Energy, Double_t Edep, Int_t nhits)
    {
        Double_t cog_x = 0;
        for (Int_t i = 0; i < nhits; i++)
            cog_x += Hit_X.at(i) * Digi_Hit_Energy.at(i);
        cog_x /= Edep;
        return cog_x;
    }, {"Hit_X", "Digi_Hit_Energy", "Edep", "nhits"})
    .Define("COG_Y", [] (vector<Double_t> Hit_Y, vector<Double_t> Digi_Hit_Energy, Double_t Edep, Int_t nhits)
    {
        Double_t cog_y = 0;
        for (Int_t i = 0; i < nhits; i++)
            cog_y += Hit_Y.at(i) * Digi_Hit_Energy.at(i);
        cog_y /= Edep;
        return cog_y;
    }, {"Hit_Y", "Digi_Hit_Energy", "Edep", "nhits"})
    .Define("COG_Z", [] (vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy, Double_t Edep, Int_t nhits)
    {
        Double_t cog_z = 0;
        for (Int_t i = 0; i < nhits; i++)
            cog_z += Hit_Z.at(i) * Digi_Hit_Energy.at(i);
        cog_z /= Edep;
        return cog_z;
    }, {"Hit_Z", "Digi_Hit_Energy", "Edep", "nhits"})
    .Define("hclx", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy)
    {
        vector<Double_t> hclx;
        HTTool* httool = new HTTool(Hit_X, Hit_Y, Hit_Z, Digi_Hit_Energy);
        hclx = httool->GetHclX();
        delete httool;
        return hclx;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Digi_Hit_Energy"})
    .Define("hcly", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy)
    {
        vector<Double_t> hcly;
        HTTool* httool = new HTTool(Hit_X, Hit_Y, Hit_Z, Digi_Hit_Energy);
        hcly = httool->GetHclY();
        delete httool;
        return hcly;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Digi_Hit_Energy"})
    .Define("hclz", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy)
    {
        vector<Double_t> hclz;
        HTTool* httool = new HTTool(Hit_X, Hit_Y, Hit_Z, Digi_Hit_Energy);
        hclz = httool->GetHclZ();
        delete httool;
        return hclz;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Digi_Hit_Energy"})
    .Define("hcle", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy)
    {
        vector<Double_t> hcle;
        HTTool* httool = new HTTool(Hit_X, Hit_Y, Hit_Z, Digi_Hit_Energy);
        hcle = httool->GetHclE();
        delete httool;
        return hcle;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Digi_Hit_Energy"})
    .Define("ntrack", [] (vector<Double_t> Hit_X, vector<Double_t> Hit_Y, vector<Double_t> Hit_Z, vector<Double_t> Digi_Hit_Energy)
    {
        Int_t ntrack = 0;
        HTTool* httool = new HTTool(Hit_X, Hit_Y, Hit_Z, Digi_Hit_Energy);
        ntrack = httool->GetNtrack();
        delete httool;
        return ntrack;
    }, {"Hit_X", "Hit_Y", "Hit_Z", "Digi_Hit_Energy"})
    // Below are time analyses (relevant data not stored)
    /*
    .Define("time_mean_hit", [] (vector<Double_t> hit_time, Int_t nhits)
    {
        Double_t tot = 0;
        for (Double_t& i : hit_time)
            tot += i;
        return tot / nhits;
    }, {"Hit_Time", "nhits"})
    .Define("time_rms_hit", [] (vector<Double_t> hit_time, Int_t nhits)
    {
        Double_t tot2 = 0;
        for (Double_t& i : hit_time)
            tot2 += i * i;
        return TMath::Sqrt(tot2 / nhits);
    }, {"Hit_Time", "nhits"})
    .Define("time_mean_shower", [] (vector<Double_t> hit_time, vector<Int_t> layer, Int_t beginning, Int_t ending, Int_t nhits)
    {
        Double_t tot = 0;
        Int_t hits = 0;
        for (Int_t i = 0; i < n; i++)
        {
            if (layer.at(i) >= beginning && layer.at(i) < ending)
            {
                hits++;
                tot += hit_time.at(i), 2;
            }
        }
        if (hits == 0)
            return 0.0;
        else
            return tot / nhits;
    }, {"Hit_Time", "layer", "shower_start", "shower_end", "nhits"})
    .Define("time_rms_shower", [] (vector<Double_t> hit_time, vector<Int_t> layer, Int_t beginning, Int_t ending, Int_t nhits)
    {
        Double_t tot2 = 0;
        Int_t hits = 0;
        for (Int_t i = 0; i < n; i++)
        {
            if (layer.at(i) >= beginning && layer.at(i) < ending)
            {
                hits++;
                tot2 += TMath::Power(hit_time.at(i), 2);
            }
        }
        if (hits == 0)
            return 0.0;
        else
            return tot2 / nhits;
    }, {"Hit_Time", "layer", "shower_start", "shower_end", "nhits"})
    */
    //.Range(1)
    .Snapshot(tree, outname);
    delete dm;
    return 1;
}

Int_t PIDTool::TrainBDT()
{
    TMVA::Tools::Instance();

    unordered_map<TTree*, TString> trsig;
    unordered_map<TTree*, TString> trbkg;
    unordered_map<TTree*, TString> tesig;
    unordered_map<TTree*, TString> tebkg;

    for (auto i : train_sig)
    {
        TFile* f = TFile::Open(i.first.first, "READ");
        TTree* t = (TTree*) f->Get(i.first.second);
        trsig.insert(pair<TTree*, TString>(t, i.second));
    }

    for (auto j : train_bkg)
    {
        TFile* f = TFile::Open(j.first.first, "READ");
        TTree* t = (TTree*) f->Get(j.first.second);
        trbkg.insert(pair<TTree*, TString>(t, j.second));
    }

    for (auto k : test_sig)
    {
        TFile* f = TFile::Open(k.first.first, "READ");
        TTree* t = (TTree*) f->Get(k.first.second);
        tesig.insert(pair<TTree*, TString>(t, k.second));
    }

    for (auto l : test_bkg)
    {
        TFile* f = TFile::Open(l.first.first, "READ");
        TTree* t = (TTree*) f->Get(l.first.second);
        tebkg.insert(pair<TTree*, TString>(t, l.second));
    }

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName( "TMVAMulticlass.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory* factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
    for (auto i : var)
        dataloader->AddVariable(i.first, i.second);

    // You can add an arbitrary number of signal or background trees
    for (auto i : trsig)
        dataloader->AddSignalTree(i.first, 1.0, TMVA::Types::kTraining);
    for (auto j : trbkg)
        dataloader->AddBackgroundTree(j.first, 1.0, TMVA::Types::kTraining);
    for (auto k : tesig)
        dataloader->AddSignalTree(k.first, 1.0, TMVA::Types::kTesting);
    for (auto l : tebkg)
        dataloader->AddBackgroundTree(l.first, 1.0, TMVA::Types::kTesting);

    dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );

    factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
    // Now you can ask the factory to train, test, and evaluate the MVAs
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

    std::cout << "==> ROOT file written: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification finished!" << std::endl;

    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    //if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

    return 0;	
}

Int_t PIDTool::BDTNtuple(const string& fname, const string& tname)
{
    ROOT::EnableImplicitMT();
    string outname = fname;
    outname = outname.substr(outname.find_last_of('/') + 1);
    outname = "bdt_" + outname;
    // This loads the library
    TMVA::Tools::Instance();

    // Default MVA methods to be trained + tested
    std::map<std::string, Int_t> Use;

    // Cut optimisation
    Use["BDTG"] = 1;
    std::cout << "==> Start TMVAMulticlassApplication" << std::endl;
    TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );
    Float_t bdt_COG_X;
    Float_t bdt_COG_Y;
    Float_t bdt_COG_Z;
    Float_t bdt_Edep;
    Float_t bdt_Emean;
    Float_t bdt_FD_2D;
    Float_t bdt_clusterE1E9;
    Float_t bdt_hit_layer;
    Float_t bdt_nhits;
    Float_t	bdt_ntrack;
    Float_t bdt_shower_density;
    Float_t bdt_shower_end;
    Float_t bdt_shower_layer_ratio;
    Float_t bdt_shower_length;
    Float_t bdt_shower_radius;
    Float_t	bdt_shower_start;
    Float_t bdt_xwidth;
    Float_t bdt_ywidth;
    Float_t bdt_zwidth;

    reader->AddVariable("Edep",               &bdt_Edep);
    reader->AddVariable("COG_X",              &bdt_COG_X);
    reader->AddVariable("COG_Y",              &bdt_COG_Y);
    reader->AddVariable("COG_Z",              &bdt_COG_Z);
    reader->AddVariable("Emean",              &bdt_Emean);
    reader->AddVariable("FD_2D",              &bdt_FD_2D);
    reader->AddVariable("clusterE1E9",        &bdt_clusterE1E9);
    reader->AddVariable("hit_layer",          &bdt_hit_layer);
    reader->AddVariable("nhits",              &bdt_nhits);
    reader->AddVariable("ntrack",             &bdt_ntrack);
    reader->AddVariable("shower_density",     &bdt_shower_density);
    reader->AddVariable("shower_end",         &bdt_shower_end);
    reader->AddVariable("shower_layer_ratio", &bdt_shower_layer_ratio);
    reader->AddVariable("shower_length",      &bdt_shower_length);
    reader->AddVariable("shower_radius",      &bdt_shower_radius);
    reader->AddVariable("shower_start",       &bdt_shower_start);
    reader->AddVariable("xwidth",             &bdt_xwidth);
    reader->AddVariable("ywidth",             &bdt_ywidth);
    reader->AddVariable("zwidth",             &bdt_zwidth);

    reader->BookMVA("BDTG method", TString("dataset/weights/TMVAMulticlass_BDTG.weights.xml"));
    cout << "Booked" << endl;
    vector<string> rdf_input = { "COG_X", "COG_Y", "COG_Z", "Edep", "Emean", "FD_2D", "clusterE1E9", "hit_layer", "nhits", "ntrack", "shower_density", "shower_end", "shower_layer_ratio", "shower_length", "shower_radius", "shower_start", "xwidth", "ywidth", "zwidth"};

    ROOT::RDataFrame df(tname, fname);

    auto bdtout = df.Define("BDT_pi_plus", [&](Double_t cogx, Double_t cogy, Double_t cogz, Double_t e, Double_t em, Double_t fd, Double_t cl, Double_t hl, Int_t nh, Int_t n, Double_t d, Int_t se, Double_t lr, Double_t l, Double_t r, Int_t s, Double_t x, Double_t y, Double_t z)
    {
        bdt_COG_X              = cogx;
        bdt_COG_Y              = cogy;
        bdt_COG_Z              = cogz;
        bdt_Edep               = e;
        bdt_Emean              = em;
        bdt_FD_2D              = fd;
        bdt_clusterE1E9        = cl;
        bdt_hit_layer          = hl;
        bdt_nhits              = nh;
        bdt_ntrack             = n;
        bdt_shower_density     = d;
        bdt_shower_end         = se;
        bdt_shower_layer_ratio = lr;
        bdt_shower_length      = l;
        bdt_shower_radius      = r;
        bdt_shower_start       = s;
        bdt_xwidth             = x;
        bdt_ywidth             = y;
        bdt_zwidth             = z;
        return (reader->EvaluateMulticlass( "BDTG method" ))[0];
//        return (reader->EvaluateMulticlass( "BDTG method" ))[1];
    }, rdf_input)
    .Define("BDT_mu_plus", [&](Double_t cogx, Double_t cogy, Double_t cogz, Double_t e, Double_t em, Double_t fd, Double_t cl, Double_t hl, Int_t nh, Int_t n, Double_t d, Int_t se, Double_t lr, Double_t l, Double_t r, Int_t s, Double_t x, Double_t y, Double_t z)
    {
        bdt_COG_X              = cogx;
        bdt_COG_Y              = cogy;
        bdt_COG_Z              = cogz;
        bdt_Edep               = e;
        bdt_Emean              = em;
        bdt_FD_2D              = fd;
        bdt_clusterE1E9        = cl;
        bdt_hit_layer          = hl;
        bdt_nhits              = nh;
        bdt_ntrack             = n;
        bdt_shower_density     = d;
        bdt_shower_end         = se;
        bdt_shower_layer_ratio = lr;
        bdt_shower_length      = l;
        bdt_shower_radius      = r;
        bdt_shower_start       = s;
        bdt_xwidth             = x;
        bdt_ywidth             = y;
        bdt_zwidth             = z;
        return (reader->EvaluateMulticlass( "BDTG method" ))[1];
//        return (reader->EvaluateMulticlass( "BDTG method" ))[2];
    }, rdf_input)
    .Define("BDT_e_plus", [&](Double_t cogx, Double_t cogy, Double_t cogz, Double_t e, Double_t em, Double_t fd, Double_t cl, Double_t hl, Int_t nh, Int_t n, Double_t d, Int_t se, Double_t lr, Double_t l, Double_t r, Int_t s, Double_t x, Double_t y, Double_t z)
    {
        bdt_COG_X              = cogx;
        bdt_COG_Y              = cogy;
        bdt_COG_Z              = cogz;
        bdt_Edep               = e;
        bdt_Emean              = em;
        bdt_FD_2D              = fd;
        bdt_clusterE1E9        = cl;
        bdt_hit_layer          = hl;
        bdt_nhits              = nh;
        bdt_ntrack             = n;
        bdt_shower_density     = d;
        bdt_shower_end         = se;
        bdt_shower_layer_ratio = lr;
        bdt_shower_length      = l;
        bdt_shower_radius      = r;
        bdt_shower_start       = s;
        bdt_xwidth             = x;
        bdt_ywidth             = y;
        bdt_zwidth             = z;
        return (reader->EvaluateMulticlass( "BDTG method" ))[2];
//        return (reader->EvaluateMulticlass( "BDTG method" ))[3];
    }, rdf_input)
    /*
    .Define("bdt_proton", [&](Double_t cogx, Double_t cogy, Double_t cogz, Double_t e, Double_t em, Double_t fd, Double_t cl, Double_t hl, Int_t nh, Int_t n, Double_t d, Int_t se, Double_t lr, Double_t l, Double_t r, Int_t s, Double_t x, Double_t y, Double_t z)
    {
        bdt_COG_X              = cogx;
        bdt_COG_Y              = cogy;
        bdt_COG_Z              = cogz;
        bdt_Edep               = e;
        bdt_Emean              = em;
        bdt_FD_2D              = fd;
        bdt_clusterE1E9        = cl;
        bdt_hit_layer          = hl;
        bdt_nhits              = nh;
        bdt_ntrack             = n;
        bdt_shower_density     = d;
        bdt_shower_end         = se;
        bdt_shower_layer_ratio = lr;
        bdt_shower_length      = l;
        bdt_shower_radius      = r;
        bdt_shower_start       = s;
        bdt_xwidth             = x;
        bdt_ywidth             = y;
        bdt_zwidth             = z;
        return (reader->EvaluateMulticlass( "BDTG method" ))[0];
    }, rdf_input)
    */
    .Snapshot(tname, outname);
    return 1;
}
