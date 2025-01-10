#include "WoodsSaxon_1D.h"
#include "WoodsSaxon_2D.h"

using namespace std;






// Woods-Saxon 1D tests
void WS_1D_func()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0, t_max = 0.25;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TF1 *test = ws.getWoodsSaxon1D();
    test->Draw();
}

void WS_transform_1D_func()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0.0001, t_max = 0.1;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TF1 *test = ws.getWStransform();
    test->Draw();
}

void WS_transform_1D_hist()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0.0001, t_max = 0.1, q_min = 0.01, q_max = 0.31623;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test = ws.getWStransformHist();
    test->Draw();
}


// Woods-Saon 2D tests
void WS_2D_func()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 1000;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TF2 *test = ws.getWoodsSaxon2D();
    test->Draw();
}

void WS_2D_hist()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 1000;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test = ws.getWoodsSaxonHist();
    test->ProjectionY()->Draw();
}

void WS_transform_2D_hist()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 1000;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test = ws.getWStransformHist();
    test->ProjectionY()->Draw();
}
