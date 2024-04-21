#ifndef TOPO_SPACEPARA_H_
#define TOPO_SPACEPARA_H_

#include <list>
#include <vector>
#include <Eigen/Core>
#include "filterOption.h"

using namespace std;
using namespace Eigen;

struct WeightPara {
    double weight;
    int position;
};

class StructureSpacePara {
private:

    VectorXi geometryPara;            //N dimension. N=number of dipoles. Each position stores the para index in VectorXi Para : 0->Para[0]...
    vector<vector<int>> Paratogeometry;
    VectorXd Para;                    // THESE ARE THE INPUTDIELS !!! P dimension. P=number of parameters. Same as Para_origin if Filter=False. Filtered and biased para if Filter=True.
    VectorXd Para_origin;             //Un-filtered, unbiased para. No use when Filter=False.
    VectorXd Para_filtered;           //Filtered but unbiased para. No use when Filter=False.
    bool Filter;                      //True for with Filter. False for without Filter. Defualt as False for most initilizations.
    FilterOption* Filterstats;        //Only used when Filter=True
    vector<vector<WeightPara>> FreeWeight;
    vector<int> ParaDividePos;
    bool Periodic;
    int Lx;
    int Ly;

    //From StructureAndSpace.h
    VectorXi geometry;
    int Nx, Ny, Nz, N;

public:
    StructureSpacePara(Vector3i bind_, VectorXi* geometry_, int Nx_, int Ny_, int Nz_, int N_, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_);
    void assignFreeWeightsForFilter( );
    void ChangeFilter( );
    VectorXi* get_geometryPara( );
    VectorXd* get_Para( );
    VectorXd* get_Para_origin( );
    VectorXd* get_Para_filtered( );
    Vector3i* get_bind( );
    bool get_Filter( );
    FilterOption* get_Filterstats( );
    vector<vector<WeightPara>>* get_FreeWeight( );
    vector<int>* get_ParaDividePos( );
    vector<vector<int>>* get_Paratogeometry( );

    //from StructureAndSpace.h
    VectorXi* get_geometry( );
    int get_geometry_size( );
    tuple<int, int, int, int> get_Ns( );

};

#endif


