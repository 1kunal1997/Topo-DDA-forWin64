#ifndef TOPO_DDAMODEL_H_
#define TOPO_DDAMODEL_H_

#include "AProductCore.h"

class DDAModel {
private:
    //------------------------------------Get from AProductCore------------------------------ For the same AMatrix, these are all the same
    AProductCore* Core;
    //AProductCore OwnCore;

    //-----------------------------------Not from AProductCore------------------------------- For the same AMatrix, these can be diff for diff DDAModel
    bool RResultSwitch;               //0(false) for plot only E field on the structure points (fast), 1(true) for using RResult different from R to plot (slow but adjustable).
    VectorXi RResult;                //The position matrix for the EResult (where you want to plot the E field)
    double E0;
    Vector3d n_E0;
    Vector3d n_K;
    VectorXcd P;
    VectorXcd E;
    VectorXcd Einternal;              //E field on structure points
    VectorXcd EResult;                //E field on designated points
    VectorXcd al;                       // 1 over alpha instead of alpha.
    VectorXcd diel;                     //Real dielectric from diel_old. Needed to calculate the Lorentz factor.
    bool verbose;
    VectorXcd P_max;
    VectorXcd al_max;

    int N;
    int Nx;
    int Ny;
    int Nz;
    double d;
    VectorXi* R;
    VectorXd* diel_old;
    double K;
    double lam;
    VectorXcd* material;
    VectorXd* Para;

    //------------------different for different angles------------------
    int time;
    int ITERATION;
    double Error;


public:
    DDAModel(VectorXd* Para_, VectorXi* R_, VectorXd* diel_old_, int Nx_, int Ny_, int Nz_, int N_, Vector3d n_K_, double E0_, Vector3d n_E0_, double lam_, VectorXcd material_, double nback_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_, double d_, bool verbose_ = true);
    DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_, bool verbose_=true);
    DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_, VectorXi* RResult_, bool verbose_ = true);
    ~DDAModel( );
    void bicgstab(int MAX_ITERATION, double MAX_ERROR);
    void bicgstab(int MAX_ITERATION, double MAX_ERROR, int EVOITERATION);  //FOR DEBUG ONLY. OUTPUT SOME VALUE AT CERTAIN EVO ITERATION.
    void change_E(VectorXcd E_);
    void reset_E();             //reset E to E0                                
    void UpdateAlpha();                                //update alpha according to updated diel in AProductCore.
    void UpdateAlphaSingle(int idx);
    void solve_E();                                                        //update the result E field on each dipole or on a designated space
    void update_E_in_structure();                                          //update the result E field on each dipole 
    VectorXcd Aproductwithalb(VectorXcd& b);                    //add the al*b term on base of AproductCore
    void output_to_file();
    void output_to_file(string save_position, int iteration, int ModelLabel);              //especially used for EvoOptimization
    void output_to_file(string save_position, int iteration);             //For simplify output
    //void output_to_file(string save_position, double wavelength, int iteration);
    void InitializeP(VectorXcd& Initializer);
    VectorXcd* get_P();
    Vector3d get_nE0();
    Vector3d get_nK();
    double get_E0();
    VectorXcd* get_Einternal();
    AProductCore* get_Core();
    VectorXcd* get_al();
    VectorXcd* get_P_max();
    VectorXcd* get_al_max();
    int get_ITERATION();

    //-----------------From AProductCore-----------------------

    int get_N();
    int get_Nx( );
    int get_Ny( );
    int get_Nz( );
    double get_lam( );
    double get_d();
    VectorXi* get_R();
    VectorXd* get_diel_old();
    VectorXd* get_Para( );
};

#endif