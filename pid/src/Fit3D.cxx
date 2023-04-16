#include "Fit3D.h"

struct SumDistance2
{
    Double_t distance2(Double_t x, Double_t y, Double_t z, const Double_t *p)
    {
        u = (x1 - x0).Unit();
        d2 = ((xp - x0).Cross(u)).Mag2();
        return d2;
    }

    Double_t operator() (const Double_t* par)
    {
        x = fGraph->GetX();
        y = fGraph->GetY();
        z = fGraph->GetZ();
        npoints = fGraph->GetN();

        for (Int_t i = 0; i < npoints; i++)
            sum += distance2(x[i], y[i], z[i], par);

        return sum;
    }
};
