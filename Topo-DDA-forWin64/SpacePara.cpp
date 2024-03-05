#include <iostream>
#include <map>
#include <set>

#include "SpacePara.h"
#include "Tools.h"

void FCurrentinsert(map<vector<int>, int>* FCurrent, vector<int> currentxy, int* currentpos, string insertmode, vector<double>* symaxis) {
    if (insertmode == "None") {
        (*FCurrent).insert(pair<vector<int>, int>(currentxy, *currentpos));
        (*currentpos) += 1;
    }
    else if(insertmode == "4fold") {
        vector<int> sym1{ int(round(2 * (*symaxis)[0] - currentxy[0])), currentxy[1] };
        vector<int> sym2{ int(round(2 * (*symaxis)[0] - currentxy[0])), int(round(2 * (*symaxis)[1] - currentxy[1])) };
        vector<int> sym3{ currentxy[0], int(round(2 * (*symaxis)[1] - currentxy[1])) };
        if ((*FCurrent).count(sym1)) {
            (*FCurrent).insert(pair<vector<int>, int>(currentxy, (*FCurrent)[sym1]));
        }
        else if ((*FCurrent).count(sym2)) {
            (*FCurrent).insert(pair<vector<int>, int>(currentxy, (*FCurrent)[sym2]));
        }
        else if ((*FCurrent).count(sym3)) {
            (*FCurrent).insert(pair<vector<int>, int>(currentxy, (*FCurrent)[sym3]));
        }
        else {
            (*FCurrent).insert(pair<vector<int>, int>(currentxy, *currentpos));
            (*currentpos) += 1;
        }
    }
    else {
        cout << "FCurrentinsert: This sym mode not supported yet" << endl;
        throw 1;
    }
}

double calweight(int xo, int yo, int x, int y, double r) {
    return r - sqrt(double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo));
}

bool circlerange(int xo, int yo, int x, int y, double r) {
    if (double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo) <= r * r - 0.01) {
        return true;
    }
    else {
        return false;
    }
}

int Get3divSize(VectorXi* geometry) {
    int N = (*geometry).size();
    if (N % 3 != 0) {
        cout << "int Get3divSize(VectorXi geometry):Geometry size must be time of 3" << endl;
        throw N;
    }
    return int(round(N / 3));
}

// return the geometry as a set of pixels, each pixel represented as a vector of 3 ints {[0,0,0], [1,0,0], ...}
set<vector<int>> Get3divSet(VectorXi* geometry) {
    set<vector<int>> result;
    int N = Get3divSize(geometry);
    for (int i = 0; i <= N - 1; i++) {
        vector<int> tmp{ (*geometry)(3 * i), (*geometry)(3 * i + 1), (*geometry)(3 * i + 2) };
       
        if (result.count(tmp)) {
            cout << "set<vector<int>> Get3divSet(VectorXi* geometry): Element is already present in the set" << endl;
            cout << "Which means grid number " << i << " has shown up at least twice which is not allowed" << endl;
            cout << "Corresponding grid is " << (*geometry)(3 * i) << " " << (*geometry)(3 * i + 1) << " " << (*geometry)(3 * i + 2) << endl;
            throw 1;
        }
        result.insert(tmp);
    }
    return result;
}

list<set<vector<int>>> Get3divSetList(vector<VectorXi*> geometry) {
    list<set<vector<int>>> result;
    for (int i = 0; i < geometry.size(); i++) {
        set<vector<int>> tmp = Get3divSet(geometry[i]);
        result.push_back(tmp);
    }
    return result;
}

bool CheckOverlap(VectorXi* geometry1, VectorXi* geometry2) {
    int n1 = Get3divSize(geometry1);
    int n2 = Get3divSize(geometry2);
    for (int i = 0; i <= n1 - 1; i++) {
        int x1 = (*geometry1)(3 * i);
        int y1 = (*geometry1)(3 * i + 1);
        int z1 = (*geometry1)(3 * i + 2);
        for (int j = 0; j <= n2 - 1; j++) {
            int x2 = (*geometry2)(3 * j);
            int y2 = (*geometry2)(3 * j + 1);
            int z2 = (*geometry2)(3 * j + 2);
            if (x1 == x2 && y1 == y2 && z1 == z2) {
                return false;
            }
        }
    }

    return true;
}

bool CheckOverlapList(vector<VectorXi*> geometry) {
    if (geometry.size() <= 1) {
        return true;
    }

    for (int i = 0; i <= int(geometry.size()) - 2; i++) {
        for (int j = i + 1; j <= int(geometry.size()) - 1;j++) {
            if (!CheckOverlap(geometry[i], geometry[j])) {
                cout << "geometry " << i << " and geometry " << j << " overlap" << endl;
                return false;
            }
        }
    }
    return true;
}

// return 3*N vector of geometry. if multiple structures, connects them into one geometry vector 
VectorXi ConnectGeometry(vector<VectorXi*> geometry) {
    if (!CheckOverlapList(geometry)) {
        cout << "CheckOverlapList Fail" << endl;
        throw 1;
    }

    vector<int> N_record;
    int N = 0;
    int Nprev;
    for (int i = 0; i < geometry.size(); i++) {
        cout << "i " << i << endl;
        Nprev = N;
        N += Get3divSize(geometry[i]);
        N_record.push_back(N - Nprev);
    }
    
    VectorXi result = VectorXi::Zero(3 * N);

    int pos = 0;
    for (int i = 0; i < geometry.size(); i++) {
        for (int j = 0; j <= N_record[i] - 1; j++) {
            result(3 * pos) = (*geometry[i])(3 * j);
            result(3 * pos + 1) = (*geometry[i])(3 * j + 1);
            result(3 * pos + 2) = (*geometry[i])(3 * j + 2);
            pos++;
        }
    }
    return result;

}

VectorXi SpacePara::cut(VectorXi* big, VectorXi* smalll) {

    int number_origin = round((*smalll).size() / 3);
    MatrixXi big_scope = find_scope_3_dim(big);
    //cout<<"big_scope "<<big_scope<<endl;
    list<int> positions_in;
    int number_out = 0;
    //cout<<"small_scope "<<find_scope_3_dim(small)<<endl;
    for (int i = 0; i <= number_origin - 1; i++) {
        if (((*smalll)(3 * i) < big_scope(0, 0)) || ((*smalll)(3 * i) > big_scope(0, 1)) ||
            ((*smalll)(3 * i + 1) < big_scope(1, 0)) || ((*smalll)(3 * i + 1) > big_scope(1, 1)) ||
            ((*smalll)(3 * i + 2) < big_scope(2, 0)) || ((*smalll)(3 * i + 2) > big_scope(2, 1))) {
            number_out += 1;
        }
        else {
            positions_in.push_back(i);
        }
    }
    int number_in = positions_in.size();
    VectorXi geometry = VectorXi::Zero(3 * number_in);
    for (int i = 0; i <= number_in - 1; i++) {
        int j = positions_in.front();
        positions_in.pop_front();
        geometry(3 * i) = (*smalll)(3 * j);
        geometry(3 * i + 1) = (*smalll)(3 * j + 1);
        geometry(3 * i + 2) = (*smalll)(3 * j + 2);
    }
    if (number_out == 0) {
        cout << "The geometry you built is entirely in the space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_real " << number_in << endl;
    }
    else {
        cout << "The geometry you built is at least partially outside of space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_out " << number_out << endl;
        cout << "number_real " << number_in << endl;
    }

    return geometry;
}


SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    vector<Structure>* ln = (*space).get_ln();

    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }


    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter;
    int lx, ly;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);

        cout << "xcenter=" << xcenter << " ycenter=" << ycenter << endl;
        cout << "lx=" << lx << " " << "ly=" << ly << endl;
        


        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }        
    }

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }
  
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    vector<Structure>* ln = (*space).get_ln();

    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter, zcenter;
    int lx, ly, lz;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        zcenter = round(((double)rand() / RAND_MAX) * (double(scope(2, 1) - scope(2, 0) + 1))) + scope(2, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);
        lz = round(((double)rand() / RAND_MAX) * (limitz2 - limitz1) + limitz1);


        cout << "xcenter=" << xcenter << " ycenter=" << ycenter << " zcenter=" << zcenter << endl;
        cout << "lx=" << lx << " " << "ly=" << ly << " " << "lz=" << lz << endl;



        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y, z;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    z = bind(2) * (2 * k + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly) && (abs(z - zcenter) <= lz)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }
    }

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, VectorXi* geometryPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    vector<Structure>* ln = (*space).get_ln();

    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;                //CAN IMPROVE
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter;
    int lx, ly;

    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);

        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, VectorXi* geometryPara_, double initial_num) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    vector<Structure>* ln = (*space).get_ln();

    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter, zcenter;
    int lx, ly, lz;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        zcenter = round(((double)rand() / RAND_MAX) * (double(scope(2, 1) - scope(2, 0) + 1))) + scope(2, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);
        lz = round(((double)rand() / RAND_MAX) * (limitz2 - limitz1) + limitz1);


        //cout << "xcenter=" << xcenter << " ycenter=" << ycenter << " zcenter=" << zcenter << endl;
        //cout << "lx=" << lx << " " << "ly=" << ly << " " << "lz=" << lz << endl;



        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y, z;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    z = bind(2) * (2 * k + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly) && (abs(z - zcenter) <= lz)) {
                        Para(pos) = initial_diel_func(initial_num);
                    }
                }
            }
        }
    }



}
/*
SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);
    
    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFparax, NFparay, NFparaz, NFpara;
    NFparax = ceil(double(FParascope(0, 1) - FParascope(0, 0) + 1) / bind(0));
    NFparay = ceil(double(FParascope(1, 1) - FParascope(1, 0) + 1) / bind(1));
    NFparaz = ceil(double(FParascope(2, 1) - FParascope(2, 0) + 1) / bind(2));
    NFpara = NFparax * NFparay * NFparaz;
    
    VectorXd Para1 = initial_diel_func(initial_diel, NFpara);
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }
    
    if (int(BParaGeometry_.size()) != int(BPara_.size())) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): BParaGeometry_.size() != BPara_.size()" << endl;
        throw 1;
    }
    list<VectorXi*>::iterator it1 = BParaGeometry_.begin();
    list<int> BParaDividePos;
    int Npara = NFpara;
    while (it1 != BParaGeometry_.end()) {
        BParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(*it1);
        it1++;
    }
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        Para(i) = Para1(i);
    }

    it1 = BParaGeometry_.begin();
    list<double>::iterator it2 = BPara_.begin();
    int Parapos = NFpara;
    while (it2 != BPara_.end()) {
        int tmpN = Get3divSize(*it1);
        for (int i = 0; i <= tmpN - 1; i++) {
            Para(Parapos) = (*it2);
            Parapos++;
        }
        it1++;
        it2++;
    }

    list<set<vector<int>>> BParaGeometrySetList=Get3divSetList(BParaGeometry_);
    
    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());
    list<int>::iterator it_tmp2 = BParaDividePos.begin();

    //cout << BParaGeometry_.size() - 1 << endl;
    //cout << int(BParaGeometry_.size()) - 1 << endl;                       These two are actually different when size=0

    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {
        BParaCurrentPos(i) = (*it_tmp2);
        it_tmp2++;
    }

    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {
            int parax = floor((double(x) - FParascope(0, 0)) / bind(0));
            int paray = floor((double(y) - FParascope(1, 0)) / bind(1));
            int paraz = floor((double(z) - FParascope(2, 0)) / bind(2));
            int pos = paraz + NFparaz * (paray + NFparay * parax);
            geometryPara(i) = pos;                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();
            
            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
            
        }

        
        
       
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, vector<string> initial_diel_list, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis) {
    Filter = Filter_;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    if (Filter && (bind(0) != 1.0) || bind(1) != 1.0) {
        cout << "SpacePara--If Filter=True, bind(0) and bind(1) must be 1.0" << endl;
        throw 1;
    }

    vector<Structure>* ln = (*space).get_ln();

    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFpara;
    //int NFparax, NFparay, NFparaz, NFpara;
    //NFparax = ceil(double(FParascope(0, 1) - FParascope(0, 0) + 1) / bind(0));
    //NFparay = ceil(double(FParascope(1, 1) - FParascope(1, 0) + 1) / bind(1));
    //NFparaz = ceil(double(FParascope(2, 1) - FParascope(2, 0) + 1) / bind(2));
    
    

    int dividesym;
    if (symmetry == "None") {
        dividesym = 1;
    }
    else if (symmetry == "4fold") {
        dividesym = 4;
    }
    else {
        cout << "SpacePara: not None nor 4 fold. Not supported" << endl;
        throw 1;
    }

    NFpara = int(round((FParaGeometry.size() / 3 / bind(2) / dividesym)));

    cout << "NFpara" << NFpara << endl;
    VectorXd Para1 = VectorXd::Zero(NFpara);
    list<VectorXd> Paratmplist;
    int initialcount = 0;
    if (FParaGeometry_.size() != initial_diel_list.size()) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): FParaGeometry_.size() != initial_diel_list.size()" << endl;
        throw 1;
    }

    int NFParacount = 0;
    for (list<VectorXi*>::iterator it = FParaGeometry_.begin(); it != FParaGeometry_.end(); it++) {
        int NFparatmp = int(round(((*it)->size() / 3 / bind(2) / dividesym)));
        cout << "NFparatmp" << NFparatmp << endl;
        ParaDividePos.push_back(NFParacount);
        NFParacount += NFparatmp;
        Paratmplist.push_back(initial_diel_func(initial_diel_list[initialcount], NFparatmp));
        initialcount += 1;
    }

    int Paracurrentpos = 0;
    for (list<VectorXd>::iterator it = Paratmplist.begin(); it != Paratmplist.end(); it++) {
        for (int i = 0; i <= (*it).size() - 1; i++) {
            Para1(Paracurrentpos) = (*it)(i);
            Paracurrentpos += 1;
        }
    }

    //VectorXd Para1 = initial_diel_func(initial_diel, NFpara);
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    if (int(BParaGeometry_.size()) != int(BPara_.size())) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): BParaGeometry_.size() != BPara_.size()" << endl;
        throw 1;
    }
    list<VectorXi*>::iterator it1 = BParaGeometry_.begin();
    //list<int> BParaDividePos;
    int Npara = NFpara;
    while (it1 != BParaGeometry_.end()) {
        ParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(*it1);
        it1++;
    }
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        Para(i) = Para1(i);
    }

    it1 = BParaGeometry_.begin();
    list<double>::iterator it2 = BPara_.begin();
    int Parapos = NFpara;
    while (it2 != BPara_.end()) {
        int tmpN = Get3divSize(*it1);
        for (int i = 0; i <= tmpN - 1; i++) {
            Para(Parapos) = (*it2);
            Parapos++;
        }
        it1++;
        it2++;
    }

    list<set<vector<int>>> BParaGeometrySetList = Get3divSetList(BParaGeometry_);

    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());
    //list<int>::iterator it_tmp2 = BParaDividePos.begin();

    //cout << BParaGeometry_.size() - 1 << endl;
    //cout << int(BParaGeometry_.size()) - 1 << endl;                       These two are actually different when size=0

    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {
        //BParaCurrentPos(i) = (*it_tmp2);
        //it_tmp2++;
        BParaCurrentPos(i) = ParaDividePos[i+(FParaGeometry_).size()];
        //cout << BParaDividePos[i] << endl;
    }
    //cout << "BParaDividePos.size()" << BParaDividePos.size() << endl;
    
    map<vector<int>, int> FCurrent;
    int currentpos = 0;
    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {
            //int parax = floor((double(x) - FParascope(0, 0)) / bind(0));
            //int paray = floor((double(y) - FParascope(1, 0)) / bind(1));
            //int paraz = floor((double(z) - FParascope(2, 0)) / bind(2));
            vector<int> currentxy{ x,y };
            if (!FCurrent.count(currentxy)) {
                //FCurrent.insert(pair<vector<int>, int>(currentxy, currentpos));//Change this for all different symmetries
                FCurrentinsert(&FCurrent, currentxy, &currentpos, symmetry, &symaxis);
                //cout << currentpos << endl;
            }
            geometryPara(i) = FCurrent[currentxy];                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();

            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
        }
    }
    if (Filter == true) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if (Filterstats == NULL) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }
        Paratogeometry = vector<vector<int>>(Npara);
        for (int i = 0; i <= N - 1; i++) {
            (Paratogeometry[geometryPara(i)]).push_back(i);
        }

        //Only works when freepara is 2D binding (Only do the filter in 2D)
        FreeWeight = vector<vector<WeightPara>>(NFpara);

        double rfilter = (*Filterstats).get_rfilter();
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);

            

            for (int j = 0; j <= Npara - 1; j++) {
                //bool inornot = false;
                //cout << j << endl;
                for (int k = 0; k <= Paratogeometry[j].size() - 1; k++) {
                    int posr = Paratogeometry[j][k];

                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }

                }

            }
        }



    }


}
*/




SpacePara::SpacePara(Space* space_, Vector3i bind_, vector<string> initial_diel_list, vector<double> BPara_, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis) {
    Filter = Filter_;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    if (Filter && (bind(0) != 1.0) || bind(1) != 1.0) {
        cout << "SpacePara--If Filter=True, bind(0) and bind(1) must be 1.0" << endl;
        throw 1;
    }

    vector<Structure>* ln = (*space).get_ln();
    
    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    vector<VectorXi*> FParaGeometry_, BParaGeometry_;

    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        if (true) {
            FParaGeometry_.push_back((*ln)[i].get_geometry());
        }
        else {
            BParaGeometry_.push_back((*ln)[i].get_geometry());
        }
    }

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    //MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFpara;

    int dividesym;
    if (symmetry == "None") {
        dividesym = 1;
    }
    else if (symmetry == "4fold") {
        dividesym = 4;
    }
    else {
        cout << "SpacePara: not None nor 4 fold. Not supported" << endl;
        throw 1;
    }

    NFpara = int(round((int(FParaGeometry.size()) / 3 / bind(2) / dividesym)));

    cout << "NFpara" << NFpara << endl;
    VectorXd Para1 = VectorXd::Zero(NFpara);
    vector<VectorXd> Paratmplist;
    
    if (FParaGeometry_.size() != initial_diel_list.size()) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): FParaGeometry_.size() != initial_diel_list.size()" << endl;
        throw 1;
    }

    int initialcount = 0;
    int NFParacount = 0;
    for (int i = 0; i < FParaGeometry_.size(); i++) {
        int NFparatmp = int(round((int((*FParaGeometry_[i]).size()) / 3 / bind(2) / dividesym)));
        cout << "NFparatmp" << NFparatmp << endl;
        ParaDividePos.push_back(NFParacount);
        NFParacount += NFparatmp;
        Paratmplist.push_back(initial_diel_func(initial_diel_list[initialcount], NFparatmp));
        initialcount += 1;
    }

    int Paracurrentpos = 0;
    for (int i = 0; i < Paratmplist.size(); i++) {
        for (int j = 0; j < Paratmplist[i].size(); j++) {
            Para1(Paracurrentpos) = (Paratmplist[i])(j);
            Paracurrentpos += 1;
        }
    }

    //VectorXd Para1 = initial_diel_func(initial_diel, NFpara);
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    if (int(BParaGeometry_.size()) != int(BPara_.size())) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): BParaGeometry_.size() != BPara_.size()" << endl;
        throw 1;
    }
    //list<int> BParaDividePos;
    int Npara = NFpara;
    for (int i = 0; i < BParaGeometry_.size(); i++) {
        ParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(BParaGeometry_[i]);
    }
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        Para(i) = Para1(i);
    }

    //it1 = BParaGeometry_.begin();
    //list<double>::iterator it2 = BPara_.begin();
    int Parapos = NFpara;
    for (int i = 0; i < BPara_.size(); i++) {
        int tmpN = Get3divSize(BParaGeometry_[i]);
        for (int j = 0; j <= tmpN - 1; j++) {
            Para(Parapos) = BPara_[i];
            Parapos++;
        }
    }

    list<set<vector<int>>> BParaGeometrySetList = Get3divSetList(BParaGeometry_);

    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());


    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {
        BParaCurrentPos(i) = ParaDividePos[i + int((FParaGeometry_).size())];
    }

    map<vector<int>, int> FCurrent;
    int currentpos = 0;
    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {
           
            vector<int> currentxy{ x,y };
            if (!FCurrent.count(currentxy)) {
                
                FCurrentinsert(&FCurrent, currentxy, &currentpos, symmetry, &symaxis);
            }
            geometryPara(i) = FCurrent[currentxy];                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();

            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
        }
    }

    Paratogeometry = vector<vector<int>>(Npara);
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
    if (Filter == true) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if (Filterstats == NULL) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }
        

        //Only works when freepara is 2D binding (Only do the filter in 2D)
        FreeWeight = vector<vector<WeightPara>>(NFpara);

        /*
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);

            double rfilter = (*Filterstats).get_rfilter();

            for (int j = 0; j <= NFpara - 1; j++) {
                //bool inornot = false;
                for (int k = 0; k <= Paratogeometry[FreeparatoPara(j)].size() - 1; k++) {
                    int posr = Paratogeometry[FreeparatoPara(j)][k];

                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        //1. Same xy plane 2. inside the circle in xy plane
                        //para>=2 wont be in NFpara
                        int posweight = FreeparatoPara(j);
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[i].push_back(WeightPara{ weight,posweight });
                        break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j.
                    }

                }

            }
        }
        */

        double rfilter = (*Filterstats).get_rfilter();
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);



            for (int j = 0; j <= Npara - 1; j++) {
                //bool inornot = false;
                //cout << j << endl;
                for (int k = 0; k < Paratogeometry[j].size() ; k++) {
                    int posr = Paratogeometry[j][k];

                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }

                }

            }
        }



    }


}

// vector of structures and BPara removed
SpacePara::SpacePara(Vector3i bind_, Space* space_, VectorXi* InputGeo, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_) {
    Filter = Filter_;
    space = space_;
    bind = bind_;
    VectorXi* total_space = InputGeo;
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = *InputGeo;
    Periodic = Periodic_;
    Lx = Lx_;
    Ly = Ly_;
    if (Filter && (bind(0) != 1.0) || bind(1) != 1.0) {
        cout << "SpacePara--If Filter=True, bind(0) and bind(1) must be 1.0" << endl;
        throw 1;
    }

    Structure* structure = (*space).get_structure();
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);
    
    // set of 1D vectors of size 3, representing each pixel. {[0,0,0], [1,0,0], ..., [Nx-1,Ny-1,Nz-1]}
    set<vector<int>> geometrySet = Get3divSet(&geometry);
    int NFpara;     // number of free parameters. for extruded, symmetric, take one quadrant of one xy-plane of geo. 121 in standard case

    int dividesym;
    if (symmetry == "None") {
        dividesym = 1;
    }
    else if (symmetry == "4fold") {
        dividesym = 4;
    }
    else {
        cout << "SpacePara: not None nor 4 fold. Not supported" << endl;
        throw 1;
    }

    NFpara = int(round((int(geometry.size()) / 3 / bind(2) / dividesym)));     // number of pixels in xy plane divided by symmetry 
    cout << "NFpara" << NFpara << endl;

    Para = VectorXd::Zero(NFpara);
    
    // seems useless because FreeparatoPara is just integers between 1 and NFPara (121), but it is used in devx as "get_Free()".
    // need to refactor that part before deleting this variable
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    // stores an index that can be used to find the free parameter index associated with that index.
    // uses symmetry and reflections in FCurrentInsert to keep the range [0,120]. for example, 
    // geometryPara(Nx-1) = 1, geometry(Nx) = 0, geometry(Nx*Ny+1)=0 because of symmetry/extrusions.
    geometryPara = VectorXi::Zero(N);
    map<vector<int>, int> FCurrent;         // only used to design geometryPara, so can be removed if geometryPara designed differently
    int currentpos = 0;
    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);

        vector<int> currentxy{ x,y };

        // for extruded geometries, same xy position can be used to refer to pixels that have different z-coord.
        if (!FCurrent.count(currentxy)) {

            FCurrentinsert(&FCurrent, currentxy, &currentpos, symmetry, &symaxis);
        }
        geometryPara(i) = FCurrent[currentxy];                                     
        //cout << "geometryPara(" << i << ") is: " << geometryPara(i) << endl;
    }

    // This is used to map a free parameter position (0 to NFPara, or 121) to a vector of
    // positions that this free parameter maps to, considering relfections and extrusions.
    // for example, the first vector would be (0, Nx, Nx*Ny-Nx, Nx*Ny, ...). if Nz is 10,
    // and you have symmetry, each position will map to 10*4, or 40 other pixels.

    Paratogeometry = vector<vector<int>>(NFpara);
    for (int i = 0; i <= N - 1; i++) {
        //cout << "dipole position : " << i << endl;

        (Paratogeometry[geometryPara(i)]).push_back(i);

        /*vector<int> currentvector = Paratogeometry[geometryPara(i)];
        for (int j = 0; j < currentvector.size(); j++) {
            cout << "Paratogeometry at index " << j << " is: " << currentvector[j] << endl;
        } */
        
    }
    
    // used to map a pixel vector to its inputdiel value (0-1)
    map<vector<int>, double> Inputmap;
    if ((*InputGeo).size() != (*Inputdiel).size()) {
        cout << "ERROR: SpacePara::SpacePara: Filter==(*InputGeo).size() != (*Inputdiel).size()" << endl;
        throw 1;
    }
    int Inputsize = int(round(int((*InputGeo).size()) / 3));
    for (int i = 0; i < Inputsize; i++) {
        Inputmap.insert(pair<vector<int>, double>(vector<int>{(*InputGeo)(3 * i), (*InputGeo)(3 * i + 1), (*InputGeo)(3 * i + 2)}, (*Inputdiel)(3 * i)));
    }

    // used to fill Para, which are the parameter values of the free indices (NFpara)
    // or one quadrant in a symmetric, extruded structure. using Paratogeometry here
    // to fetch each position, but if this is the only use of Paratogeometry, seems useless.
    for (int i = 0; i < NFpara; i++) {

        int pos = Paratogeometry[i][0];
        //cout << "Pos at position " << i << " is: " << pos << endl;
        vector<int> node{ geometry(3 * pos), geometry(3 * pos + 1), geometry(3 * pos + 2) };
        Para(i) = Inputmap[node];
    }
    cout << "para values:" << endl;

    for (int i = 0; i < Para.size(); i++) {
        cout << Para(i) << " ";
    }

    if (Filter == true) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if (Filterstats == NULL) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }
        //Only works when freepara is 2D binding (Only do the filter in 2D)
        FreeWeight = vector<vector<WeightPara>>(NFpara);

        double rfilter = (*Filterstats).get_rfilter();
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);



            for (int j = 0; j <= NFpara - 1; j++) {
                //bool inornot = false;
                //cout << j << endl;
                for (int k = 0; k < Paratogeometry[j].size(); k++) {
                    int posr = Paratogeometry[j][k];

                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if (Periodic == false) {
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                    }
                    else {
                        //Own cell
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //0, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr - Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr - Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr - Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr - Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, 0
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //0, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, 0
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr - Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr - Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                    }


                }

            }
        }
    }

}

// HEEYO!
SpacePara::SpacePara(Space* space_, Vector3i bind_, VectorXi* InputGeo, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_) {
    Filter = Filter_;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);
    Periodic = Periodic_;
    Lx = Lx_;
    Ly = Ly_;
    if (Filter && (bind(0) != 1.0) || bind(1) != 1.0) {
        cout << "SpacePara--If Filter=True, bind(0) and bind(1) must be 1.0" << endl;
        throw 1;
    }

    vector<Structure>* ln = (*space).get_ln();
    cout << "ln size is: " << (*ln).size() << endl;         // should be just 1
    cout << "ln geometry size is: " << (*ln)[0].get_geometry_size() << endl;        // number of pixels

    int n1 = 0;
    // assigns geometry to the structure's geometry in 'space'
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*ln)[i].get_geometry_size());        // 3*N
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*ln)[i].get_geometry()))(j);
        }
        n1 = n1 + n2;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    vector<VectorXi*> FParaGeometry_, BParaGeometry_;

    for (int i = 0; i < (*ln).size(); i++) {
        if (true) {
            FParaGeometry_.push_back((*ln)[i].get_geometry());
        }
        else {
            cout << "-------------------BPARAGEOMETRY ENTERED!!!!------------------------" << endl;
            BParaGeometry_.push_back((*ln)[i].get_geometry());
        }
    }

    cout << "FParaGeometry_ size: " << FParaGeometry_.size() << endl;       // size is 1 because ln.size is 1

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);       // this is just geometry of our structure. 3N

    cout << "FParaGeometry size: " << FParaGeometry.size() << endl;
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    //MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFpara;

    int dividesym;
    if (symmetry == "None") {
        dividesym = 1;
    }
    else if (symmetry == "4fold") {
        dividesym = 4;
    }
    else {
        cout << "SpacePara: not None nor 4 fold. Not supported" << endl;
        throw 1;
    }

    NFpara = int(round((int(FParaGeometry.size()) / 3 / bind(2) / dividesym)));     // number of pixels in xy plane divided by symmetry 
    cout << "NFpara" << NFpara << endl;

    VectorXd Para1 = VectorXd::Zero(NFpara);
    //vector<VectorXd> Paratmplist;
    int initialcount = 0;
    int NFParacount = 0;
    for (int i = 0; i < FParaGeometry_.size(); i++) {
        int NFparatmp = int(round((int((*FParaGeometry_[i]).size()) / 3 / bind(2) / dividesym)));
        cout << "NFparatmp" << NFparatmp << endl;
        ParaDividePos.push_back(NFParacount);               // vector of ints. defined in header file
        NFParacount += NFparatmp;
        //Paratmplist.push_back(initial_diel_func(initial_diel_list[initialcount], NFparatmp));
        //initialcount += 1;
    }

    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }
    cout << "FreeparaToParaDONE" << endl;
    int Npara = NFpara;

    // BPara is always empty as far as I can tell, so never enters this
    for (int i = 0; i < BParaGeometry_.size(); i++) {
        ParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(BParaGeometry_[i]);
    }
    cout << "Get3divSize" << endl;


    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        // cout << "Para1 at position " << i << "is: " << Para1[i] << endl;
        Para(i) = Para1(i);                     // but Para1 is empty...
    }

    list<set<vector<int>>> BParaGeometrySetList = Get3divSetList(BParaGeometry_);   // empty because BparaGeometry_ is empty
    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());
    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {                     // never enters this either
        cout << "-------------BARAGEOMTRY ENTERED!!!----------" << endl;
        BParaCurrentPos(i) = ParaDividePos[i + int((FParaGeometry_).size())];
    }

    map<vector<int>, int> FCurrent;
    int currentpos = 0;
    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {

            vector<int> currentxy{ x,y };
            if (!FCurrent.count(currentxy)) {

                FCurrentinsert(&FCurrent, currentxy, &currentpos, symmetry, &symaxis);
            }
            geometryPara(i) = FCurrent[currentxy];                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();

            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
        }
    }

    cout << "GEOMETRYPARA SIZE: " << geometryPara.size() << endl;
    cout << "NPara SIZE: " << Npara << endl;
    Paratogeometry = vector<vector<int>>(Npara);

    for (int i = 0; i <= N - 1; i++) {
        cout << "geometryPara(i) at " << i << " is: " << geometryPara(i) << endl;

        if (i == 6561) {
            cout << "for i:"<<geometryPara(i) << endl;
        }
        (Paratogeometry[geometryPara(i)]).push_back(i);
        
    }


    map<vector<int>, double> Inputmap;
    if ((*InputGeo).size() != (*Inputdiel).size()) {
        cout << "ERROR: SpacePara::SpacePara: Filter==(*InputGeo).size() != (*Inputdiel).size()" << endl;
        throw 1;
    }
    int Inputsize = int(round(int((*InputGeo).size()) / 3));
    for (int i = 0; i < Inputsize; i++) {
        Inputmap.insert(pair<vector<int>, double>(vector<int>{(*InputGeo)(3 * i), (*InputGeo)(3 * i + 1), (*InputGeo)(3 * i + 2)}, (*Inputdiel)(3 * i)));
    }
    for (int i = 0; i < Npara; i++) {
        /*cout << i << endl;*/
        int pos = Paratogeometry[i][0];
        vector<int> node{ geometry(3 * pos), geometry(3 * pos + 1), geometry(3 * pos + 2) };
        if (Inputmap.count(node)) {
            Para(i) = Inputmap[node];
        }
    }
    cout << "para values:" << endl;

     for (int i = 0; i < Para.size(); i++) {
        cout << Para(i) << " ";
     } 


    if (Filter == true) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if (Filterstats == NULL) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }
        //Only works when freepara is 2D binding (Only do the filter in 2D)
        FreeWeight = vector<vector<WeightPara>>(NFpara);

        double rfilter = (*Filterstats).get_rfilter();
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);



            for (int j = 0; j <= Npara - 1; j++) {
                //bool inornot = false;
                //cout << j << endl;
                for (int k = 0; k < Paratogeometry[j].size(); k++) {
                    int posr = Paratogeometry[j][k];

                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if (Periodic == false) {
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                    }
                    else {
                        //Own cell
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //0, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr-Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr-Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr - Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr - Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, 0
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //-1, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr - Lx, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //0, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, 1
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr + Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr + Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, 0
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                        //1, -1
                        if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr - Ly, rfilter))) {
                            //1. Same xy plane 2. inside the circle in xy plane 
                            //para>=2 wont be in NFpara
                            int posweight = j;
                            double weight = calweight(xo, yo, xr + Lx, yr - Ly, rfilter);
                            FreeWeight[i].push_back(WeightPara{ weight,posweight });
                            //break;
                            //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                        }
                    }
                    

                }

            }
        }
    }

    


}





void SpacePara::ChangeFilter() {
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    int Npara = Para.size();
    int NFpara = FreeparatoPara.size();

    double rfilter = (*Filterstats).get_rfilter();
    vector<vector<int>> Paratogeometry(Npara);
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
    vector<vector<WeightPara>> FreeWeight_tmp(NFpara);
    for (int i = 0; i <= NFpara - 1; i++) {
        int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
        int xo = geometry(3 * poso);
        int yo = geometry(3 * poso + 1);
        int zo = geometry(3 * poso + 2);

        

        for (int j = 0; j <= Npara - 1; j++) {
            //bool inornot = false;
            //cout << j << endl;
            for (int k = 0; k < Paratogeometry[j].size(); k++) {
                int posr = Paratogeometry[j][k];

                int xr = geometry(3 * posr);
                int yr = geometry(3 * posr + 1);
                int zr = geometry(3 * posr + 2);
                if (Periodic == false) {
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                }
                else {
                    
                    //Own cell
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //0, -1
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr - Ly, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr - Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, -1
                    if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr - Ly, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr - Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, 0
                    if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, 1
                    if ((zo == zr) && (circlerange(xo, yo, xr - Lx, yr + Ly, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr + Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //0, 1
                    if ((zo == zr) && (circlerange(xo, yo, xr, yr + Ly, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr + Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, 1
                    if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr + Ly, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr + Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, 0
                    if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr, rfilter))) {
                        /*if (i == 0) {
                            cout << xr << endl;
                            cout << yr << endl;
                            cout << zr << endl;
                        }*/
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, -1
                    if ((zo == zr) && (circlerange(xo, yo, xr + Lx, yr - Ly, rfilter))) {
                        //if (i == 0) {
                        //    cout << xr << endl;
                        //    cout << yr << endl;
                        //    cout << zr << endl;
                        //}
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr - Ly, rfilter);
                        FreeWeight_tmp[i].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                }

            }

        }
    }
    FreeWeight = FreeWeight_tmp;

}

void SpacePara::ChangeBind(Vector3i bind_) {
    bind = bind_;
    VectorXi geometryPara_before = geometryPara;
    VectorXd Para_before = Para;
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil((scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil((scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil((scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;

    vector<vector<int>> Paratogeometry(Npara, vector<int>(2, 0));
    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
        Paratogeometry[pos][0] += 1;
        Paratogeometry[pos][1] += Para_before(geometryPara_before(i));
    }

    for (int i = 0; i <= Npara - 1; i++) {
        Para(i) = Paratogeometry[i][1] / Paratogeometry[i][0];
    }
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }

    


}

Space* SpacePara::get_space() {
    return space;
}

VectorXi SpacePara::get_geometry() {
    return geometry;
}

VectorXi* SpacePara::get_geometryPara() {
    return &geometryPara;
}

VectorXd* SpacePara::get_Para() {
    return &Para;
}

VectorXd* SpacePara::get_Para_origin() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Para_origin()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_origin;
}

VectorXd* SpacePara::get_Para_filtered() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Para_filtered()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_filtered;
}

Vector3i* SpacePara::get_bind() {
    return &bind;
}

VectorXi* SpacePara::get_Free() {
    return &FreeparatoPara;
}

bool SpacePara::get_Filter() {
    return Filter;
}

FilterOption* SpacePara::get_Filterstats() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Filterstats()--Filter can not be false" << endl;
        throw 1;
    }
    return Filterstats;
}

vector<vector<WeightPara>>* SpacePara::get_FreeWeight() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_FreeWeight()---Filter can not be false" << endl;
        throw 1;
    }
    return &FreeWeight;
}

vector<int>* SpacePara::get_ParaDividePos() {
    return &ParaDividePos;
}

vector<vector<int>>* SpacePara::get_Paratogeometry() {
    return &Paratogeometry;
}
