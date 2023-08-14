#include "PIDTool.h"

using namespace std;

Int_t main(Int_t argc, char* argv[])
{
    string file = "", tree = "";
    string fntuple = "", tntuple = "";
    string help = "";
    Int_t train = 0, bdtntuple = 0;

    for (Int_t i = 1; i < argc; i++)
    {
        if (string(argv[i]) == string("-f"))
            file = string(argv[i + 1]);

        else if (string(argv[i]) == string("-t"))
            tree = string(argv[i + 1]);

        else if (string(argv[i]) == string("-train"))
            train = 1;

        else if (string(argv[i]) == string("-tag"))
        {
            bdtntuple = 1;
            fntuple = string(argv[i + 1]);
            tntuple = string(argv[i + 2]);
        }

        else if (string(argv[i]) == string("-h") || string(argv[i]) == string("-help"))
        {
            cout << "Reconstruct variables: PIDTool -f (file name) -t (tree name)" << endl;
            cout << "PID with BDT:          PIDTool -train" << endl;
            cout << "Classification:        PIDtool -tag (file name) (tree name)" << endl;
        }
    }

    PIDTool* pt = new PIDTool();

    if (file != "" && tree != "")
    {
        cout << "File: " << file << endl;
        cout << "Tree: " << tree << endl;

        pt->GenNtuple(file, tree);
    }

    if (train == 1)
    {
        cout << "Training and testing..." << endl;

        pt->AddVar("COG_X_overall",      'D');
        pt->AddVar("COG_Y_overall",      'D');
        pt->AddVar("COG_Z_overall",      'D');
        pt->AddVar("E1E9",               'D');
        pt->AddVar("E9E25",              'D');
        pt->AddVar("E9E49",              'D');
        pt->AddVar("Edep",               'D');
        pt->AddVar("Emean",              'D');
        pt->AddVar("FD_2D",              'D');
        pt->AddVar("FD_3D",              'D');
        pt->AddVar("hit_layer",          'D');
        pt->AddVar("hit_time_mean",      'D');
        pt->AddVar("hit_time_rms",       'D');
        pt->AddVar("nhits",              'I');
        pt->AddVar("ntrack",             'D');
        pt->AddVar("shower_density",     'D');
        pt->AddVar("shower_end",         'I');
        pt->AddVar("shower_layer",       'D');
        pt->AddVar("shower_layer_ratio", 'D');
        pt->AddVar("shower_length",      'D');
        pt->AddVar("shower_radius",      'D');
        pt->AddVar("shower_start",       'I');
        pt->AddVar("shower_time_mean",   'D');
        pt->AddVar("shower_time_rms",    'D');
        pt->AddVar("xwidth",             'D');
        pt->AddVar("ywidth",             'D');
        pt->AddVar("zwidth",             'D');

        pt->AddTrainSig("/lustre/collider/chenjiyuan/hbuana_data/build/run/pi-/train_pion.root", "Calib_Hit", "pion");
        pt->AddTestSig ("/lustre/collider/chenjiyuan/hbuana_data/build/run/pi-/test_pion.root",  "Calib_Hit", "pion");

        pt->AddTrainBkg("/lustre/collider/chenjiyuan/hbuana_data/build/run/mu-/train_muon.root", "Calib_Hit", "muon");
        pt->AddTrainBkg("/lustre/collider/chenjiyuan/hbuana_data/build/run/e-/train_e.root",     "Calib_Hit", "e");
        pt->AddTestBkg ("/lustre/collider/chenjiyuan/hbuana_data/build/run/mu-/test_muon.root",  "Calib_Hit", "muon");
        pt->AddTestBkg ("/lustre/collider/chenjiyuan/hbuana_data/build/run/e-/test_e.root",      "Calib_Hit", "e");

        pt->TrainBDT();
	}

    if (bdtntuple == 1)
    {
        cout << "Classifying..." << endl;
        pt->BDTNtuple(fntuple, tntuple);
    }

    delete pt;
    return 0;
}
