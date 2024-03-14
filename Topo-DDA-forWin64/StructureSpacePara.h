#ifndef TOPO_STRUCTURESPACEPARA_H_
#define TOPO_STRUCTURESPACEPARA_H_

#include <vector>
#include <list>

#include "Eigen/Core"
#include "filterOption.h"

using namespace std;
using namespace Eigen;

struct WeightParas {
    double weight;
    int position;
};

class StructureSpacePara {
private:
    int Nx, Ny, Nz, N;
    VectorXi geometry;                //3N dimension
    VectorXi geometryPara;            //N dimension. N=number of dipoles. Each position stores the para index in VectorXi Para : 0->Para[0]...
    vector<vector<int>> Paratogeometry;
    VectorXd Para;                    // THESE ARE THE INPUTDIELS !!! P dimension. P=number of parameters. Same as Para_origin if Filter=False. Filtered and biased para if Filter=True.
    VectorXd Para_origin;             //Un-filtered, unbiased para. No use when Filter=False.
    VectorXd Para_filtered;           //Filtered but unbiased para. No use when Filter=False.
    VectorXi FreeparatoPara;          //Position of free parameters inside Para. dimension<=P. FreeparatoPara[i] is the index of a free parameter inside Para.
    //vector<list<int>> Paratogeometry;  //P dimension. Each position stores a list of corresponding dipole index for parameter for this specific position.
    bool Filter;                      //True for with Filter. False for without Filter. Defualt as False for most initilizations.
    FilterOption* Filterstats;        //Only used when Filter=True
    vector<vector<WeightParas>> FreeWeight;
    vector<int> ParaDividePos;
    bool Periodic;
    int Lx;
    int Ly;

public:

    StructureSpacePara(int Nx_, int Ny_, int Nz_, int N_, VectorXi* geometry_, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_);

    // from Structure and Space
    VectorXi* get_geometry( );
    int get_geometry_size( );
    tuple<int, int, int, int> get_Ns( );

    // from SpacePara
    void ChangeFilter( );
    //VectorXi get_geometry( );         getting rid of this one because already have pointer-based getter of this. might run into issues so check
    VectorXi* get_geometryPara( );
    VectorXd* get_Para( );
    VectorXd* get_Para_origin( );
    VectorXd* get_Para_filtered( );
    VectorXi* get_Free( );
    //vector<list<int>>* get_Paratogeometry();
    bool get_Filter( );
    FilterOption* get_Filterstats( );
    vector<vector<WeightParas>>* get_FreeWeight( );
    vector<int>* get_ParaDividePos( );
    vector<vector<int>>* get_Paratogeometry( );
};
#endif

