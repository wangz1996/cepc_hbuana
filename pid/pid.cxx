#include "PIDTool.h"

using namespace std;

Int_t main(Int_t argc, char* argv[])
{
    string file = "", tree = "";
    string fntuple = "", tntuple = "";
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
    }

    cout << "File: " << file << endl;
    cout << "Tree: " << tree << endl;

    PIDTool* pt = new PIDTool();

    if (file != "" && tree != "")
        pt->GenNtuple(file, tree);

    if (train == 1)
    {
//        pt->AddVar("E_dep",              'D');
        pt->AddVar("Hits_no",            'L');
        pt->AddVar("Shower_density",     'D');
        pt->AddVar("Shower_layer_ratio", 'D');
        pt->AddVar("Shower_length",      'D');
        pt->AddVar("Shower_radius",      'D');
        pt->AddVar("Shower_start",       'D');

        pt->AddTrainSignal("/lustre/collider/chenjiyuan/hbuana/build/e_signal/mc/e.root",               "Calib_Hit", "e");
        pt->AddTrainBkg   ("/lustre/collider/chenjiyuan/hbuana/build/e_signal/mc/background.root",      "Calib_Hit", "bkg");
        pt->AddTestSignal ("/lustre/collider/chenjiyuan/hbuana/build/e_signal/mc/test_e.root",          "Calib_Hit", "e");
        pt->AddTestBkg    ("/lustre/collider/chenjiyuan/hbuana/build/e_signal/mc/test_background.root", "Calib_Hit", "bkg");

        pt->TrainBDT();
	}

    if (bdtntuple == 1)
        pt->BDTNtuple(fntuple, tntuple);

    delete pt;
    return 0;
}
