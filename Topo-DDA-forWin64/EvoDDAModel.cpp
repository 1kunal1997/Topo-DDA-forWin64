#include <chrono>
#include <iostream>
#include <fstream>

#include "EvoDDAModel.h"
#include "Tools.h"

using namespace std::chrono;

EvoDDAModel::EvoDDAModel(string objName_, vector<double> objPara_, double epsilon_fix_, bool HavePathRecord_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, string save_position_, CoreStructure* CStr_, DDAModel* Model_) {
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
    CStr = CStr_;
    Model = Model_;

    MaxObjarray = VectorXd::Zero(ModelNum);
    Originarray = VectorXd::Zero(ModelNum);
    PreviousObj = 0.0;
    CutoffHold = 0;

    StructureSpacePara* structurespacepara = (*CStr).get_structurespacepara();
    VectorXd* Para = (*structurespacepara).get_Para();
    int n_para = (*Para).size();                     //Total number of parameters
    gradientsquare = VectorXd::Zero(n_para);

    int N = (*CStr).get_N();
    VectorXcd Ptmp = VectorXcd::Zero(N * 3);
    PolarizationforOrigin = Ptmp;
    PolarizationforAdjoint = Ptmp;
    PolarizationforOriginMax = Ptmp;
    PolarizationforAdjointMax = Ptmp;

    //-----------------generate obj list----------------------------
    double origin = 0.0;
    objfunc = ObjFactory(objName, objPara, Model);
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
        double objWithPenalty;
        auto it_allModel = allModel.begin();
        auto it_allObj = allObj.begin();
        double coeff = penaltyweight;
        string coeff_type = penaltytype;
        // double coeff = 0.1;
        double penalty = 0.0;

        auto out_start = high_resolution_clock::now();
        (*CStr).output_to_file(save_position + "CoreStructure\\", iteration + start_num, "simple");
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
        objWithPenalty = (*objfunc).GetValWithPenalty(coeff);
        obj = (*objfunc).GetVal();

        convergence << obj << " ";
        convergenceWithPenalty << objWithPenalty << " ";

        Originiterations << (*Model).get_ITERATION() << endl;
        TotalOriginIt += (*Model).get_ITERATION();

        cout << "Object function at iteration " << iteration << " is " << obj << endl;

        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        VectorXd* diel_old = (*CStr).get_diel_old();
        VectorXd* diel_old_max = (*CStr).get_diel_old_max();

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



        /*
        if((*ObjectFunctionNames).size()>1){
            list<double> obj_minor =  this->MinorObj();
            list<double>::iterator it_obj_minor = obj_minor.begin();
            for(int i=0; i<=obj_minor.size()-1; i++){
                convergence << *it_obj_minor << " ";
                it_obj_minor++;
            }
        }
        */
        convergence << "\n";
        convergenceWithPenalty << "\n";

        StructureSpacePara* structurespacepara = (*CStr).get_structurespacepara();
        VectorXi* geometryPara = (*structurespacepara).get_geometryPara();
        VectorXd* Para = (*structurespacepara).get_Para();
        int n_para_all = (*Para).size();
        int N = (*CStr).get_N();
        vector<vector<int>>* Paratogeometry = (*structurespacepara).get_Paratogeometry();

        cout << "n_para_all is: " << n_para_all << endl;

        VectorXd penaltygradients = VectorXd::Zero(n_para_all);
        VectorXd objgradients = VectorXd::Zero(n_para_all);
        VectorXd gradients = VectorXd::Zero(n_para_all);

        cout << "about to start partial derivative part" << endl;

        //----------------------------------------get partial derivative of current model---------------------------
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        cout << "---------------------------START PARTIAL DERIVATIVE ----------------------" << endl;
        VectorXd devx;
        VectorXcd Adevxp;
        VectorXcd devp;
        tie(devx, Adevxp) = this->devx_and_Adevxp_stateless(epsilon_partial, objfunc, originalObjValue, Para, Paratogeometry);
        //tie(devx, Adevxp) = this->devx_and_Adevxp(epsilon_partial, Model, objfunc, originalObjValue);
        cout << "done with devx and adevxp, starting with devp" << endl;
        devp = this->devp(epsilon_partial, Model, objfunc, originalObjValue);
        cout << "done with devp, changing E using devp" << endl;
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(t2 - t1).count();
        cout << "------------------------PARTIAL DERIVATIVE finished in " << duration / 1000 << " s-------------------------" << endl;

        //------------------------------------Solving adjoint problem-----------------------------------------
        cout << "---------------------------START ADJOINT PROBLEM ----------------------" << endl;
        (*Model).change_E(devp);

        (*Model).InitializeP(PolarizationforAdjoint);

        (*Model).bicgstab(MAX_ITERATION, MAX_ERROR);

        if (HaveAdjointHeritage == true) {
            PolarizationforAdjoint = *((*Model).get_P());
        }

        VectorXcd* P = (*Model).get_P();
        VectorXcd lambdaT = (*P);
        (*Model).reset_E();                                  //reset E to initial value
        Adjointiterations << (*Model).get_ITERATION() << endl;
        TotalAdjointIt += (*Model).get_ITERATION();

        cout << "D O N E!" << endl;
        //times lambdaT and Adevxp together
        VectorXcd mult_result;
        mult_result = VectorXcd::Zero(n_para_all);               //multiplication result has the length of parameter

        for (int i = 0; i <= n_para_all - 1; i++) {
            int FreeParaPos = i;

            vector<int>::iterator it = (*Paratogeometry)[FreeParaPos].begin();
            for (int j = 0; j <= (*Paratogeometry)[FreeParaPos].size() - 1; j++) {
                int position = *it;
                mult_result(i) += lambdaT(3 * position) * Adevxp(3 * position);
                mult_result(i) += lambdaT(3 * position + 1) * Adevxp(3 * position + 1);
                mult_result(i) += lambdaT(3 * position + 2) * Adevxp(3 * position + 2);
                it++;
            }
        }
        cout << "D O N E!" << endl;
        VectorXd mult_result_real = VectorXd::Zero(n_para_all);
        for (int i = 0; i <= n_para_all - 1; i++) {
            complex<double> tmp = mult_result(i);
            mult_result_real(i) = tmp.real();
        }
        objgradients += devx - mult_result_real;              //What's the legitimacy in here to ignore the imag part?

        for (int i = 0; i <= n_para_all - 1; i++) {
            penaltygradients(i) = 1 - 2 * (*Para)(i);
        }

        const double coeff_min = 0.0;
        const double coeff_max = 0.0;

        if (iteration > 0) {

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
        }

        gradients = objgradients - coeff * penaltygradients;

        if ((*((*CStr).get_structurespacepara())).get_Filter() && (iteration >= 1)) {
            //For iteration=0, Para, Para_origin, Para_filtered are all the same. No need for updating gradients.
            //current_it is actually the it in current evo-1 as the str is updated in iteration-1.
            cout << "-----------ENTERED FILTER IF STATEMENT!!!----------------" << endl;
            gradients = gradients_filtered(gradients, iteration - 1, MAX_ITERATION_EVO - 1);
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

        if ((*structurespacepara).get_Filter()) {
            if ((*((*structurespacepara).get_Filterstats())).filterchange(iteration)) {
                (*structurespacepara).ChangeFilter();
            }
        }
        cout << "right before spaceparams sentence" << endl;
        StructureSpacePara* structurespaceparams = (*CStr).get_structurespacepara();
        VectorXd* Params = (*structurespaceparams).get_Para();

        penalty = calculatePenalty(*Params);

        cout << "PENALTY AFTER GRADIENTS IS: " << 40 * penalty << endl;

        /* if (penalty < 10 && (iteration + 1) % 10 == 0) {
             for (int i : step) {
                 step[i] = round(step[i]);
                 cout << step[i] << endl;
             }
         } */


        (*CStr).UpdateStr(step, iteration, MAX_ITERATION_EVO - 1);

        (*Model).UpdateAlpha();                  //Dont forget this, otherwise bicgstab wont change

    }


    Originiterations << TotalOriginIt << endl;
    Adjointiterations << TotalAdjointIt << endl;

    convergence.close();
    convergenceWithPenalty.close();
    Originiterations.close();
    Adjointiterations.close();

}

tuple<VectorXd, VectorXcd> EvoDDAModel::devx_and_Adevxp_stateless(double epsilon, ObjDDAModel* Obj, double origin, VectorXd* para_, vector<vector<int>>* paratogeometry_) {

    VectorXd* Para = para_;           // the values (0-1) of the free parameters ONLY (one quadrant, 2D)
    vector<vector<int>>* Paratogeometry = paratogeometry_;

    int N = (*Model).get_N();
    VectorXcd* al = (*Model).get_al();       // members from DDAModel
    VectorXcd* P = (*Model).get_P();

    int n_para_all = (*Para).size();
    VectorXcd Adevxp = VectorXcd::Zero(3 * N);
    VectorXd devx = VectorXd::Zero(n_para_all);

    cout << "n_para_all is: " << n_para_all << endl;

    for (int i = 0; i < n_para_all; i++) {
        int FreeParaPos = i;
        if (FreeParaPos != i) {
            cout << "----------------------------ERROR IN FREEPARAPOS!!--------------------------" << endl;
        }
        double diel_old_origin = (*Para)(FreeParaPos);
        double diel_old_tmp = diel_old_origin;
        int sign = 0;
        if (diel_old_origin >= epsilon) {
            sign = -1;
        }
        else {
            sign = 1;
        }
        diel_old_tmp += sign * epsilon;

        vector<int>::iterator it = (*Paratogeometry)[FreeParaPos].begin();
        for (int j = 0; j <= (*Paratogeometry)[FreeParaPos].size() - 1; j++) {
            //cout << (*it) << endl;
            int position = *it;
            complex<double> alphaorigin = (*al)(3 * position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, true);
            (*CStr).UpdateStrSingle(position, diel_old_tmp);
            (*Model).UpdateAlphaSingle(position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, false);
            complex<double> change = ((*al)(3 * position) - alphaorigin) / (sign * epsilon);
            Adevxp(3 * position) = change;
            Adevxp(3 * position + 1) = change;
            Adevxp(3 * position + 2) = change;

            it++;
        }

        devx(i) = (Obj->GroupResponse() - origin) / (sign * epsilon);  //If some obj has x dependency but you denote the havepenalty as false, it will still actually be calculated in an efficient way.
        it = (*Paratogeometry)[FreeParaPos].begin();
        for (int j = 0; j <= (*Paratogeometry)[FreeParaPos].size() - 1; j++) {
            int position = *it;
            if (Obj->Have_Devx) Obj->SingleResponse(position, true);
            (*CStr).UpdateStrSingle(position, diel_old_origin);
            (*Model).UpdateAlphaSingle(position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, false);
            it++;
        }

    }

    for (int i = 0; i <= 3 * N - 1; i++) {
        Adevxp(i) = Adevxp(i) * ((*P)(i));
    }

    return make_tuple(devx, Adevxp);
}

tuple<VectorXd, VectorXcd> EvoDDAModel::devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjDDAModel* Obj, double origin) {
    
    int N = (*CurrentModel).get_N();
    StructureSpacePara* spacepara = (*CurrentModel).get_structurespacepara();
    VectorXd* Para = (*spacepara).get_Para();           // the values (0-1) of the free parameters ONLY (one quadrant, 2D)
    vector<vector<int>>* Paratogeometry = (*spacepara).get_Paratogeometry();

    VectorXcd* al = (*CurrentModel).get_al();       // members from DDAModel
    VectorXcd* P = (*CurrentModel).get_P();

    int n_para_all = (*Para).size();
    VectorXcd Adevxp = VectorXcd::Zero(3 * N);
    VectorXd devx = VectorXd::Zero(n_para_all);

    cout << "n_para_all is: " << n_para_all << endl;

    for (int i = 0; i < n_para_all; i++) {
        int FreeParaPos = i;
        if (FreeParaPos != i) {
            cout << "----------------------------ERROR IN FREEPARAPOS!!--------------------------" << endl;
        }
        double diel_old_origin = (*Para)(FreeParaPos);
        double diel_old_tmp = diel_old_origin;
        int sign = 0;
        if (diel_old_origin >= epsilon) {
            sign = -1;
        }
        else {
            sign = 1;
        }
        diel_old_tmp += sign * epsilon;

        vector<int>::iterator it = (*Paratogeometry)[FreeParaPos].begin();
        for (int j = 0; j <= (*Paratogeometry)[FreeParaPos].size() - 1; j++) {
            //cout << (*it) << endl;
            int position = *it;
            complex<double> alphaorigin = (*al)(3 * position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, true);
            (*CStr).UpdateStrSingle(position, diel_old_tmp);
            (*CurrentModel).UpdateAlphaSingle(position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, false);
            complex<double> change = ((*al)(3 * position) - alphaorigin) / (sign * epsilon);
            Adevxp(3 * position) = change;
            Adevxp(3 * position + 1) = change;
            Adevxp(3 * position + 2) = change;

            it++;
        }

        devx(i) = (Obj->GroupResponse() - origin) / (sign * epsilon);  //If some obj has x dependency but you denote the havepenalty as false, it will still actually be calculated in an efficient way.
        it = (*Paratogeometry)[FreeParaPos].begin();
        for (int j = 0; j <= (*Paratogeometry)[FreeParaPos].size() - 1; j++) {
            int position = *it;
            if (Obj->Have_Devx) Obj->SingleResponse(position, true);
            (*CStr).UpdateStrSingle(position, diel_old_origin);
            (*CurrentModel).UpdateAlphaSingle(position);
            if (Obj->Have_Devx) Obj->SingleResponse(position, false);
            it++;
        }

    }
    
    for (int i = 0; i <= 3 * N - 1; i++) {
        Adevxp(i) = Adevxp(i) * ((*P)(i));
    }

    return make_tuple(devx, Adevxp);
}

VectorXcd EvoDDAModel::devp(double epsilon, DDAModel* CurrentModel, ObjDDAModel* Obj, double origin){
    //move origin=Obj0->GetVal() outside because it is the same for one partial derivative of the entire structure
    VectorXcd* P = (*CurrentModel).get_P();
    VectorXcd result=VectorXcd::Zero((*P).size());          // 3N dimension
    for(int i=0;i<= (*P).size() -1;i++){
        int position = i/3;
        
        Obj->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)+epsilon;
        
        Obj->SingleResponse(position, false);
        
        result(i)+=(Obj->GroupResponse()-origin)/epsilon;
        
        Obj->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)-epsilon;
        
        complex<double> epsilonimag=epsilon*1.0i;
        
        (*P)(i)= (*P)(i)+epsilonimag;
        
        Obj->SingleResponse(position, false);
        
        complex<double> tmpRes = (Obj->GroupResponse()-origin)/ epsilonimag;
        result(i)+=tmpRes;
        
        Obj->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)- epsilonimag;
        
        Obj->SingleResponse(position, false);
    }
    cout << "Devp_sum: " << result.sum() << endl;
    return result;
}


// for constrained gradient descent


// HEEYO !!

ObjDDAModel* EvoDDAModel::ObjFactory(string ObjectName, vector<double> ObjectParameters, DDAModel* ObjDDAModel){
    /*if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }*/
    if (objName == "PointE"){
        return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
    }
    //if (MajorObjectFunctionName == "PointEList") {
    //    return new ObjPointListEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "PointI") {
    //    return new ObjPointIDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    if (objName == "IntegratedE") {
        return new ObjIntegratedEDDAModel(ObjectParameters, ObjDDAModel);
    }
    //if (MajorObjectFunctionName == "MidAvgE") {
    //    return new ObjMidAvgEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "scattering0D") {
    //    return new Objscattering0D(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "scattering2D") {
    //    return new Objscattering2D(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "Abs") {
    //    return new ObjAbs(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "AbsPartial") {
    //    return new ObjAbsPartial(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "AbsPartialzslice") {
    //    return new ObjAbsPartialzslice(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "IntegratedEPartial") {
    //    return new ObjIntegrateEPartial(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}
    //if (MajorObjectFunctionName == "Absbyfar") {
    //    return new ObjAbsbyfar(ObjectParameters, ObjDDAModel, this, HavePenalty);
    //}

    // NOT FINALIZED. SHOULD RAISE AN EXCEPTION HERE.
    cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
    return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
}

//double EvoDDAModel::L1Norm(){
//    double Penalty = 0;
//    int N = (*CStr).get_N();
//    VectorXd* diel_old = (*CStr).get_diel_old();
//    for (int i=0;i<N;i++){
//        Penalty += 0.5-abs((*diel_old)(3*i)-0.5);
//    }
//    Penalty = Penalty * PenaltyFactor;
//    return Penalty;
//}

double EvoDDAModel::get_output_time() {
    return output_time / 1000;      //seconds
}

double PtoFderivative(const double input, const double beta, const double ita) {
    double result = 0.0;
    if (input <= ita && input >= 0.0) {
        return beta * exp(-beta * (1 - input / ita)) + exp(-beta);
    }
    else if (input > ita && input <= 1.0) {
        return beta * exp(-beta * (input - ita) / (1 - ita)) + exp(-beta);
    }
    else {
        cout << "ERROR: PtoFderivative(const double input, const double beta, const double ita)--input out of range" << endl;
        throw 1;
    }
}

VectorXd EvoDDAModel::gradients_filtered(VectorXd gradients, int current_it, int Max_it) {
    StructureSpacePara* sp = (*CStr).get_structurespacepara();
    FilterOption* fo = (*sp).get_Filterstats();
    const VectorXd* Para_filtered = (*sp).get_Para_filtered();
    (*fo).update_beta(current_it, Max_it);                     //current_it is actually the it in current evo-1 as the str is updated in iteration-1.
    const double gbeta = (*fo).get_beta();
    const double gita = (*fo).get_ita();
    

    int NFpara = gradients.size();
    VectorXd result = VectorXd::Zero(NFpara);
    const vector<vector<WeightPara>>* FreeWeight = (*sp).get_FreeWeight();
    if (NFpara != (*FreeWeight).size()) {
        cout << "ERROR: EvoDDAModel::gradients_filtered--NFpara != (*FreeWeight).size()" << endl;
        throw 1;
    }
    for (int i = 0; i <= NFpara - 1; i++) {
        //int position = (*Free)(i);
        //double pftmp = (*Para_filtered)(position);
        //double ptofd = PtoFderivative(pftmp, gbeta, gita);

        int num_weight = ((*FreeWeight)[i]).size();
        double sum_weight = 0.0;
        for (int j = 0; j <= num_weight - 1; j++) {
            double weight = (*FreeWeight)[i][j].weight;
            int weightpos = (*FreeWeight)[i][j].position;
            sum_weight += weight;
        }
        for (int j = 0; j <= num_weight - 1; j++) {
            double weight = (*FreeWeight)[i][j].weight;
            int weightpos = (*FreeWeight)[i][j].position;
            if (weightpos < NFpara) {                        //If there are non-paras also weighted but do not need to calculate gradients to.
                result(weightpos) += (weight / sum_weight);     
            }
            
        }
    }

    for (int i = 0; i <= NFpara - 1; i++) {
        int position = i;
        double pftmp = (*Para_filtered)(position);
        double ptofd = PtoFderivative(pftmp, gbeta, gita);
        result(position) = result(position) * ptofd * gradients(position);
    }

    return result;

}


