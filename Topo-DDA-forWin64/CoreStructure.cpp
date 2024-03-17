#include <iostream>
#include <fstream>

#include "CoreStructure.h"

CoreStructure::CoreStructure(StructureSpacePara* structurespacepara_, double d_) {

    structurespacepara = structurespacepara_; 
    d = d_;

    cout << "(d=" << d << ") " << endl;

    tie(Nx, Ny, Nz, N) = structurespacepara->get_Ns( );

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    R = structurespacepara->get_geometry( );
    VectorXi* geometryPara = structurespacepara->get_geometryPara( );
    VectorXd* Para = structurespacepara->get_Para( );
    //---------------------------------------------------initial diel------------------------------------
    diel_old = VectorXd::Zero(3 * N);
    diel_old_max = diel_old;
    for ( int i = 0; i <= N - 1; i++ ) {
        double dieltmp = ( *Para )( ( *geometryPara )( i ) );
        diel_old(3 * i) = dieltmp;
        diel_old(3 * i + 1) = dieltmp;
        diel_old(3 * i + 2) = dieltmp;
    }
} 

void CoreStructure::UpdateStr(VectorXd step, int current_it, int Max_it) {
    cout << "step in UpdateStr: " << step.mean() << endl;
    VectorXi* geometryPara = structurespacepara->get_geometryPara();
    VectorXd* Para = structurespacepara->get_Para();

    int Parasize = Para->size();
    if (Parasize != step.size()) {
        cout << "ERROR: In CoreStructure::UpdateStr(VectorXd step), step.size!=FreePara.size";
        throw 1;
    }
    

    if (structurespacepara->get_Filter()) {
        //When there is filter
        VectorXd* Para_origin = structurespacepara->get_Para_origin();
        VectorXd* Para_filtered = structurespacepara->get_Para_filtered();
        FilterOption* Filterstats = structurespacepara->get_Filterstats();
        (*Filterstats).update_beta(current_it, Max_it);                  //Update beta value according to current iteration
        
        for (int i = 0; i <= Parasize - 1; i++) {
            (*Para_origin)(i) += step(i);
            if ((*Para_origin)(i) >= 1) {
                (*Para_origin)(i) = 1;
            }
            if ((*Para_origin)(i) <= 0) {
                (*Para_origin)(i) = 0;
            }
        }
        
        cout << "Beta at iteration " << current_it << " is " << Filterstats->get_beta() << endl;

        vector<vector<WeightPara>>* FreeWeight = structurespacepara->get_FreeWeight();
        for (int i = 0; i <= Parasize - 1; i++) {
            int weightnum = ((*FreeWeight)[i]).size();
            double numerator = 0.0;
            double denominator = 0.0;
            for (int j = 0; j <= weightnum - 1; j++) {
                numerator += ((*FreeWeight)[i][j].weight) * (*Para_origin)((*FreeWeight)[i][j].position);
                denominator += ((*FreeWeight)[i][j].weight);

            }
            (*Para_filtered)(i) = numerator / denominator;

            double Para_physical = Filterstats->SmoothDensity((*Para_filtered)(i));     
            (*Para)(i) = Para_physical;

        }
    }
    else {//When there is no filter
        for (int i = 0; i <= Parasize - 1; i++) {
            (*Para)(i) += step(i);
            if ((*Para)(i) >= 1) {
                (*Para)(i) = 1;
            }
            if ((*Para)(i) <= 0) {
                (*Para)(i) = 0;
            }
        }
    }


    for (int i = 0; i <= N - 1; i++) {
        int position = (*geometryPara)(i);
        double value = (*Para)(position);
        diel_old(3 * i) = value;
        diel_old(3 * i + 1) = value;
        diel_old(3 * i + 2) = value;
    }
}


void CoreStructure::UpdateStrSingle(int idx, double value) {

    diel_old(3 * idx) = value;
    diel_old(3 * idx + 1) = value;
    diel_old(3 * idx + 2) = value;

}

void CoreStructure::output_to_file() {

    ofstream fout("CoreStructure.txt");
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout.close();
}

void CoreStructure::output_to_file(string save_position, int iteration, string mode) {

    if (mode == "normal") {
        string name;
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
        fout << R << endl;
        fout << diel_old << endl;
        fout << d << endl;
        fout.close();
    }
    else {
        string name;
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        fout << diel_old << endl;
        fout.close();
    }
}

int CoreStructure::get_N() {
    return N;
}
int CoreStructure::get_Nx() {
    return Nx;
}
int CoreStructure::get_Ny() {
    return Ny;
}
int CoreStructure::get_Nz() {
    return Nz;
}
VectorXi* CoreStructure::get_R() {
    return R;
}
double CoreStructure::get_d() {
    return d;
}

StructureSpacePara* CoreStructure::get_structurespacepara() {
    return structurespacepara;
}
VectorXd* CoreStructure::get_diel_old() {
    return &diel_old;
}
VectorXd* CoreStructure::get_diel_old_max() {
    return &diel_old_max;
}

double CoreStructure::calculate_Penalty() {
    double penalty = 0.0;

    // TODO: implement it!

    return penalty;
}

