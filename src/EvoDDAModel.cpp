#include <chrono>
#include <iostream>
#include <fstream>

#include "EvoDDAModel.h"
#include "Tools.h"

using namespace std::chrono;

EvoDDAModel::EvoDDAModel(string objName_, vector<double> objPara_, double epsilon_fix_, bool HavePathRecord_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, string save_position_, DDAModel* Model_) {
    output_time = 0.0;
    objName = objName_;
    save_position = save_position_;
    objPara = objPara_;
    //HavePenalty = HavePenalty_;
    //PenaltyFactor = PenaltyFactor_;
    epsilon_fix = epsilon_fix_;
    epsilon_tmp = epsilon_fix;
    HavePathRecord = HavePathRecord_;
    HaveOriginHeritage = HaveOriginHeritage_;
    HaveAdjointHeritage = HaveAdjointHeritage_;
    MaxObj = 0.0;
    Stephold = 0;
    Model = Model_;

    double lam = (Model->get_Core( ))->get_lam( );
    cout << "lam inside of EvoDDAModel constructor is: " << lam << endl;

    MaxObjarray = VectorXd::Zero(ModelNum);
    Originarray = VectorXd::Zero(ModelNum);
    PreviousObj = 0.0;
    CutoffHold = 0;

    VectorXd* Para = Model->get_parameters();
    int n_para = (*Para).size();                     //Total number of parameters
    gradientsquare = VectorXd::Zero(n_para);

    int N = Model->get_N();
    VectorXcd Ptmp = VectorXcd::Zero(N * 3);
    PolarizationforOrigin = Ptmp;
    PolarizationforAdjoint = Ptmp;
    PolarizationforOriginMax = Ptmp;
    PolarizationforAdjointMax = Ptmp;

    //-----------------generate obj list----------------------------
    double origin = 0.0;
    //allObj.push_back(obj);

}

void EvoDDAModel::EvoOptimizationQuick(double penaltyweight, string penaltytype, int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num) {
    ofstream convergence;
    ofstream convergenceWithPenalty;
    ofstream Originiterations;
    ofstream Adjointiterations;
    int TotalOriginIt = 0;
    int TotalAdjointIt = 0;

    convergence.open(save_position + "convergence.txt");
    convergenceWithPenalty.open(save_position + "convergenceWithPenalty.txt");
    Originiterations.open(save_position + "Originiterations.txt");
    Adjointiterations.open(save_position + "Adjointiterations.txt");
    //Parameters for Adam Optimizer.
    //double beta1 = 0.9;
    //double beta2 = 0.99;
    //---------------new beta for THG test----------------
    double beta1 = 0.9;
    double beta2 = 1 - (1 - beta1) * (1 - beta1);
    VectorXd V;
    VectorXd S;


    high_resolution_clock::time_point TotalTime0 = high_resolution_clock::now();

    double epsilon_partial = 0.001;
    for (int iteration = 0; iteration <= MAX_ITERATION_EVO - 1; iteration++) {
        //solve DDA
        cout << "######################EVO ITERATION " << iteration << "#######################" << endl;
        //get object function value
        cout << "-----------------------------START ORIGINAL PROBLEM---------------------------" << endl;

        double obj;
        double objWithPenalty = 0.0;
        double coeff = penaltyweight;
        string coeff_type = penaltytype;
        // double coeff = 0.1;
        double penalty = 0.0;

        auto out_start = high_resolution_clock::now();
        Model->outputCStr_to_file(save_position + "CoreStructure\\", iteration + start_num, "simple");
        auto out_end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(out_end - out_start).count();
        output_time += duration;

        (*Model).InitializeP(PolarizationforOrigin);
        (*Model).bicgstab(MAX_ITERATION, MAX_ERROR);
        if (HaveOriginHeritage == true) {
            PolarizationforOrigin = *((*Model).get_P());
        }
        (*Model).update_E_in_structure();
        if (iteration == MAX_ITERATION_EVO - 1) {                                    //useless fix, not gonna to use RResultswithc = true feature in the future
            (*Model).solve_E();
        }

        (*Model).solve_E();
        out_start = high_resolution_clock::now();
        (*Model).output_to_file(save_position + "Model_output\\", iteration + start_num);
        out_end = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(out_end - out_start).count();
        output_time += duration;

        // penalty = 40*calculatePenalty(*Params);                 // mupltiplying by 40 becaue params is only 121 pixels (xy plane, 4 fold symmetry)

        cout << "about to calculate object function" << endl;
        obj = (*Model).calculateObjective();

        convergence << obj << " ";
        convergenceWithPenalty << objWithPenalty << " ";

        Originiterations << (*Model).get_ITERATION() << endl;
        TotalOriginIt += (*Model).get_ITERATION();

        cout << "Object function at iteration " << iteration << " is " << obj << endl;

        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        VectorXd* diel_old = Model->get_dielectric_old();
        VectorXd* diel_old_max = Model->get_diel_old_max();

        double epsilon = epsilon_fix;

        // CURRENTLY INITIALIZED TO FALSE AND NOT SURE WHY!
        if (HavePathRecord) {

            if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001 || epsilon <= 0.0001) {
                CutoffHold += 1;
            }
            else {
                if (CutoffHold > 0) {
                    CutoffHold -= 1;
                }
            }
            cout << "CutoffHold" << CutoffHold << endl;
            PreviousObj = obj;
            if (CutoffHold >= 5) {
                cout << "Five times with small change in obj, break the iterations" << endl;
                break;
            }

            if (obj < MaxObj) {
                epsilon_tmp = epsilon_tmp / 10;
                Stephold = 0;

                (*diel_old) = (*diel_old_max);
                VectorXcd* P = (*Model).get_P();
                VectorXcd* P_max = (*Model).get_P_max();
                VectorXcd* al = (*Model).get_al();
                VectorXcd* al_max = (*Model).get_al_max();
                (*P) = (*P_max);
                (*al) = (*al_max);
                obj = maximumObjValue;
                PolarizationforOrigin = PolarizationforOriginMax;
                PolarizationforAdjoint = PolarizationforAdjointMax;

                obj = MaxObj;
                cout << "New Obj smaller then Old One, back track to previous structure and search with new step size: " << epsilon_tmp << endl;
                /*
                if (obj != Obj->GetVal()) {
                    cout << "Reset failed, Obj is not equal to MaxObj" << endl;
                }
                */
            }
            else {

                if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001) {
                    CutoffHold += 1;
                }
                else {
                    if (CutoffHold > 0) {
                        CutoffHold -= 1;
                    }
                }

                (*diel_old_max) = (*diel_old);

                VectorXcd* P = (*Model).get_P();
                VectorXcd* P_max = (*Model).get_P_max();
                VectorXcd* al = (*Model).get_al();
                VectorXcd* al_max = (*Model).get_al_max();
                (*P_max) = (*P);
                (*al_max) = (*al);
                maximumObjValue = obj;
                PolarizationforOriginMax = PolarizationforOrigin;
                PolarizationforAdjointMax = PolarizationforAdjoint;

                MaxObj = obj;
                Stephold += 1;

                if (Stephold >= 2) {
                    epsilon_tmp = epsilon_tmp * 10;
                    cout << "Two times increase with previous step size, try with larger step size: " << epsilon_tmp << endl;
                    Stephold = 0;
                }
                else {
                    cout << "Not smaller obj nor three continuous increase. Current step size is: " << epsilon_tmp << endl;
                }
            }
            epsilon = epsilon_tmp;
        }

        originalObjValue = obj;

        convergence << "\n";
        convergenceWithPenalty << "\n";

        VectorXd* Para = Model->get_parameters();
        int n_para_all = (*Para).size();
        int N = Model->get_N();
        vector<vector<int>>* Paratogeometry = Model->get_Paratogeometry();

        cout << "n_para_all is: " << n_para_all << endl;

        VectorXd penaltygradients = VectorXd::Zero(n_para_all);
        VectorXd objgradients = Model->calculateGradients(epsilon_partial, originalObjValue, MAX_ITERATION, MAX_ERROR);
        VectorXd gradients = VectorXd::Zero(n_para_all);

        for (int i = 0; i <= n_para_all - 1; i++) {
            penaltygradients(i) = 1 - 2 * (*Para)(i);
        }

        const double coeff_min = 0.0;
        const double coeff_max = 0.5;


        if (coeff_type == "exp") {
            coeff = exp_update(iteration - 1, 299, coeff_min, coeff_max);
        }
        else if (coeff_type == "piecewise") {
            coeff = piecewise_update(iteration - 1, MAX_ITERATION_EVO - 1, coeff_min, coeff_max);
        }
        else if (coeff_type == "piecewise absolute") {
            coeff = piecewise_update_absolute(iteration - 1, MAX_ITERATION_EVO - 1, coeff_min, coeff_max);
        }
        else if (coeff_type == "linear") {
            coeff = linear_update(iteration - 1, MAX_ITERATION_EVO - 1, coeff_min, coeff_max);
        }
        else {
            cout << "ERROR: coeff_type not defined" << endl;
            throw 1;
            return;
        }


        gradients = objgradients - coeff * penaltygradients;

        if (Model->get_Filter() && (iteration >= 1)) {
            //For iteration=0, Para, Para_origin, Para_filtered are all the same. No need for updating gradients.
            //current_it is actually the it in current evo-1 as the str is updated in iteration-1.
            cout << "-----------ENTERED FILTER IF STATEMENT!!!----------------" << endl;
            gradients = Model->gradients_filtered(gradients, iteration - 1, MAX_ITERATION_EVO - 1);
        }

        double epsilon_final = epsilon;

        // HEEYOOO!!!!!
        if (method == "Adam") {
            cout << "Using Adam Optimizer." << endl;
            if (iteration == 0) {
                V = (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));

            }
            else {
                V = beta1 * V + (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = beta2 * S + (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            for (int i = 0; i <= n_para_all - 1; i++) {
                gradients(i) = V(i) / (sqrt(S(i)) + 0.00000001);
            }

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else {
                epsilon_final = epsilon;
            }
        }

        if (method == "Adamdecay") {
            cout << "Using Adam Optimizer With decay." << endl;
            if (iteration == 0) {
                V = (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            else {
                V = beta1 * V + (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = beta2 * S + (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            for (int i = 0; i <= n_para_all - 1; i++) {
                gradients(i) = V(i) / (sqrt(S(i)) + 0.00000001);
            }

            int decaybarrier = 150;

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else if (iteration >= decaybarrier) {
                double expconst = 10;
                epsilon_final = epsilon * exp(-(iteration - decaybarrier) / expconst);
                //epsilon_final = epsilon / exp(-(iteration - decaybarrier));
            }
            else {
                epsilon_final = epsilon;
            }
        }

        if (method == "Adagrad") {
            for (int i = 0; i <= n_para_all - 1; i++) {
                gradientsquare(i) += pow(gradients(i), 2) / 100000;
                gradients(i) = gradients(i) / sqrt(gradientsquare(i) + 1);
            }

        }
        cout << "abs(gradients.cwiseAbs().mean() before: " << abs(gradients.cwiseAbs().mean()) << endl;

        if (abs(gradients.cwiseAbs().mean()) < 0.1) {
            gradients /= abs(gradients.cwiseAbs().mean()) / 0.1;
        }

        cout << "abs(gradients.cwiseAbs().mean() after: " << abs(gradients.cwiseAbs().mean()) << endl;
        VectorXd step = epsilon_final * gradients;            //Find the maximum. If -1 find minimum

        cout << "epsilon = " << epsilon << endl;
        cout << "step = " << step.mean() << endl;

        if (Model->get_Filter()) {
            if ((*(Model->get_Filterstats())).filterchange(iteration)) {
                cout << "ABOUT TO CHANGE THE FILTER!!!" << endl;
                Model->assignFreeWeightsForFilter();
            }
        }
        cout << "right before spaceparams sentence" << endl;
        penalty = calculatePenalty(*Para);
        cout << "PENALTY AFTER GRADIENTS IS: " << 40 * penalty << endl;

        //Model->UpdateStr(step, iteration, MAX_ITERATION_EVO - 1);
        //(*Model).UpdateAlpha();                  //Dont forget this, otherwise bicgstab wont change

        Model->UpdateParameters(step, MAX_ITERATION_EVO - 1);

    }

    Originiterations << TotalOriginIt << endl;
    Adjointiterations << TotalAdjointIt << endl;

    convergence.close();
    convergenceWithPenalty.close();
    Originiterations.close();
    Adjointiterations.close();

}


