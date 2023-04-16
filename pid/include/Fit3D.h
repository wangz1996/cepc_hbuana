#ifndef FIT3D_HH
#define FIT3D_HH
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>
#include <cassert>
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
using namespace ROOT::Math;
using namespace ROOT::Fit;
using namespace std;

struct SumDistance2
{
    TGraph2D* fGraph;

    SumDistance2(TGraph2D* g) : fGraph(g) {}

    Double_t distance2(Double_t x, Double_t y, Double_t z, const Double_t* p)
    {
        XYZVector xp(x, y, z);
        XYZVector x0(p[0], p[2], 0.0);
        XYZVector x1(p[0] + p[1], p[2] + p[3], 1.0);
        XYZVector u;
        Double_t d2;
        return d2;
    }

    Double_t operator() (const Double_t* par)
    {
        assert(fGraph != 0);
        Double_t* x;
        Double_t* y;
        Double_t* z;
        Int_t npoints;
        Double_t sum = 0;
        return sum;
    }
};

#endif
