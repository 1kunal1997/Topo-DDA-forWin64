#ifndef TOPO_EVO_H_
#define TOPO_EVO_H_

#include "DDAModel.h"
#include "CoreStructure.h"

using namespace std;
using namespace Eigen;

class EvoDDAModel {
private:
    double output_time;
    CoreStructure* CStr;
    vector<DDAModel*> allModel;                    //List of DDA models sharing the same AProductCore : "Core"
    DDAModel* Model;
    int ModelNum = 1;                                 //number of DDA model
    string save_position;

    vector<VectorXcd> PforOrigin;
    vector<VectorXcd> PforAdjoint;
    vector<VectorXcd> PforOriginMax;
    vector<VectorXcd> PforAdjointMax;

    VectorXcd PolarizationforOrigin;
    VectorXcd PolarizationforAdjoint;
    VectorXcd PolarizationforOriginMax;
    VectorXcd PolarizationforAdjointMax;

    vector<double> objPara;
    string objName;

    VectorXd Originarray;                       //Record the Obj function for partial derivative (the value before change)   
    double originalObjValue;
    //bool HavePenalty;
    //double PenaltyFactor;

    double MaxObj;                                //Record the historical maximum obj func
    double PreviousObj;                            //The previous obj
    int CutoffHold;
    VectorXd MaxObjarray;                         //the individual objs for each model when the average obj is maximum(not necessaily the maximum individual objs)
    double maximumObjValue;
    double epsilon_fix;
    double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
    bool HavePathRecord;
    bool HaveOriginHeritage;
    bool HaveAdjointHeritage;
    int Stephold;

    VectorXd gradientsquare;                    //cumulative summation of gradients square. Used in Adagrad.
public:

    EvoDDAModel(string objName_, vector<double> objPara_, double epsilon_fix_, bool HavePathRecord_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, string save_position_, CoreStructure* CStr_, DDAModel* Model_);
   // VectorXd calculateGradient();

    //functions used to calculate partial derivatives                                                        
    tuple<VectorXd, VectorXcd> devx_and_Adevxp_stateless(double epsilon, double origin, VectorXd* para_, vector<vector<int>>* paratogeometry_);         //partial derivative of obj to parameter and A to x times p
    VectorXcd devp(double epsilon, DDAModel* CurrentModel, double origin);                       //partial derivative of obj to P. Size of P

    void EvoOptimizationQuick(double penaltyweight, string penaltytype, int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num = 0);

    double get_output_time();
    //double L1Norm();
    VectorXd gradients_filtered(VectorXd gradients, int current_it, int Max_it);



};

#endif