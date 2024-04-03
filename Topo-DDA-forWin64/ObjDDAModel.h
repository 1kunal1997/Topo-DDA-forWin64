#ifndef TOPO_OBJ_H_
#define TOPO_OBJ_H_

#include "DDAModel.h"
#include <fstream>

class ObjDDAModel {
public:
    bool Have_Devx;
    //bool Have_Penalty;
    virtual void SingleResponse(int idx, bool deduction, bool hasPenalty = false) = 0;
    virtual void SingleResponseWithoutPenalty(int idx, bool deduction) = 0;
    virtual double GroupResponse() = 0;
    virtual void Reset() = 0;
    virtual double GetVal() = 0;
    //virtual double GetValWithoutPenalty() = 0;
    virtual double GetValWithPenalty(double coeff) = 0;
};


//Child classes for Obj function


class ObjPointEDDAModel : public ObjDDAModel {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    //EvoDDAModel* evomodel;
    Vector3cd E_sum;
    Vector3cd E_ext;
public:
    ObjPointEDDAModel(vector<double> parameters, DDAModel* model_);
    void SingleResponse(int idx, bool deduction, bool hasPenalty = false);
    void SingleResponseWithoutPenalty(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    double GetValWithPenalty(double coeff);
    void Reset();
};

class ObjIntegratedEDDAModel : public ObjDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    int powNum;        //Number on the exponential
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int zMin;
    int zMax;
    double ita;
    double beta;
    VectorXd* diel_old;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    VectorXd* Params;
    VectorXcd* al;
    double penalty;
    double E_int;

    string namedebugfile;
    
    
public:
    ObjIntegratedEDDAModel(vector<double> parameters, DDAModel* model_);
    void SingleResponse(int idx, bool deduction, bool hasPenalty = true);
    void SingleResponseWithoutPenalty(int idx, bool deduction);
    double GroupResponse();
    double GroupResponseWithPenalty(double coeff);
    double GetVal();
    double GetValWithPenalty(double coeff);
    void Reset();
};
#endif
