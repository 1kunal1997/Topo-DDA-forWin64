#ifndef TOPO_CORESTRUCTURE_H_
#define TOPO_CORESTRUCTURE_H_

//#include "Space.h"
//#include "SpacePara.h"
//#include "StructureSpacePara.h"

#include <list>
#include <vector>
#include "Eigen/Core"
#include "filterOption.h"

using namespace std;
using namespace Eigen;

struct WeightPara {
    double weight;
    int position;
};

class CoreStructure {
private:
    //---------------------------------Geometries, not related to wavelength-------------------------------
    //StructureSpacePara* structurespacepara;
    double d;
    VectorXi* R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.
    VectorXd diel_old;                //The 0~1 version of diel, 3*N
    VectorXd diel_old_max;

    // from StructureSpaceParahello
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

    VectorXi geometry;
    int Nx, Ny, Nz, N;

public:
    //CoreStructure(StructureSpacePara* structurespacepara_, double d_);
    CoreStructure(double d_, Vector3i bind_, VectorXi* geometry_, int Nx_, int Ny_, int Nz_, int N_, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_);
    void UpdateStr(VectorXd step, int current_it, int Max_it);
    //void UpdateStrCGD(VectorXd step, int current_it, int Max_it);
    void UpdateStrSingle(int idx, double value);
    void output_to_file();
    void output_to_file(string save_position, int iteration, string mode = "normal");

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    VectorXi* get_R();
    double get_d();
    //StructureSpacePara* get_structurespacepara( );
    VectorXd* get_diel_old();
    VectorXd* get_diel_old_max();
    double calculate_Penalty();

    // from StructureSpacePara
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