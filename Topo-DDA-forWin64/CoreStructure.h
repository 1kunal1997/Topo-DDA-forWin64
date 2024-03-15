#ifndef TOPO_CORESTRUCTURE_H_
#define TOPO_CORESTRUCTURE_H_

//#include "Space.h"
//#include "SpacePara.h"
#include "StructureSpacePara.h"

class CoreStructure {
private:
    //---------------------------------Geometries, not related to wavelength-------------------------------
    //Space* space;
    //SpacePara* spacepara;
    StructureSpacePara* structurespacepara;
    int N;                        //Number of dipoles
    int Nx;                       //scope of space. Nx*Ny*Nz!=N
    int Ny;
    int Nz;
    double d;
    VectorXi* R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.
    VectorXd diel_old;                //The 0~1 version of diel, 3*N
    VectorXd diel_old_max;

public:
    CoreStructure(StructureSpacePara* structurespacepara_, double d_);
    //CoreStructure(SpacePara* spacepara_, double d_);
    void UpdateStr(VectorXd step, int current_it, int Max_it);
    void UpdateStrCGD(VectorXd step, int current_it, int Max_it);
    //void UpdateStr(SpacePara* spacepara_);
    void UpdateStrSingle(int idx, double value);
    void output_to_file();
    void output_to_file(string save_position, int iteration, string mode = "normal");

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    VectorXi* get_R();
    double get_d();
    //SpacePara* get_spacepara();
    StructureSpacePara* get_structurespacepara( );
    VectorXd* get_diel_old();
    VectorXd* get_diel_old_max();
    double calculate_Penalty();

};

#endif