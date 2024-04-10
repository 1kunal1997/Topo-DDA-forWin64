#define NUM_THREADS 6

#include <chrono>
#include <iostream>
#include <fstream>

#include "DDAModel.h"
#include "Tools.h"

using namespace std::chrono;

ObjDDAModel* DDAModel::ObjFactory(string ObjectName, vector<double> ObjectParameters) {
    /*if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }*/
    /*if ( objName == "PointE" ) {
        return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
    }*/

    if ( objName == "IntegratedE" ) {
        return new ObjIntegratedEDDAModel(ObjectParameters, N, &P, geometry, &al);
    }

    cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
    return new ObjIntegratedEDDAModel(ObjectParameters, N, &P, geometry, &al);
}

DDAModel::DDAModel(string objName_, vector<double> objPara_, VectorXd* Para_, VectorXi* geometry_, VectorXd* diel_old_, int Nx_, int Ny_, int Nz_, int N_, Vector3d n_K_, double E0_, Vector3d n_E0_, double lam_, VectorXcd material_, double nback_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_, double d_, bool verbose_) {
    Core = new AProductCore(Nx_, Ny_, Nz_, N_, d_, lam_, material_, nback_, MAXm_, MAXn_, Lm_ * d_, Ln_ * d_, AMatrixMethod_);
    time = 0;
    ITERATION = 0;
    Error = 0.0;
    E0 = E0_;
    n_K = n_K_;
    n_E0 = n_E0_;

    N = N_;
    Nx = Nx_;
    Ny = Ny_;
    Nz = Nz_;
    lam = lam_;
    cout << "lam in DDAModel is: " << lam << endl;
    double nback = nback_;
    K = Core->get_K( );
    d = d_;
    geometry = geometry_;
    diel_old = diel_old_;
    material = Core->get_material();
    Para = Para_;

    objDDAModel = ObjFactory(objName_, objPara_);
    RResultSwitch = false;
    RResult = *geometry;

    P = VectorXcd::Zero(N * 3);
    P_max = P;
    E = VectorXcd::Zero(N * 3);
    Einternal = VectorXcd::Zero(N * 3);
    EResult = VectorXcd::Zero(N * 3);
    for ( int i = 0; i < N; i++ ) {
        E(3 * i) = E0 * n_E0(0) * ( cos(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) + sin(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) * 1i );
        E(3 * i + 1) = E0 * n_E0(1) * ( cos(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) + sin(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) * 1i );
        E(3 * i + 2) = E0 * n_E0(2) * ( cos(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) + sin(K * d * ( n_K(0) * ( *geometry )( 3 * i ) + n_K(1) * ( *geometry )( 3 * i + 1 ) + n_K(2) * ( *geometry )( 3 * i + 2 ) )) * 1i );
    }
    al = VectorXcd::Zero(N * 3);
    diel = VectorXcd::Zero(N * 3);
    for ( int i = 0; i < N * 3; i++ ) {
        int labelfloor = int(floor(( *diel_old )( i )));
        int labelnext = labelfloor + 1;
        if ( labelfloor >= 1 ) {
            labelnext = labelfloor;
        }
        std::complex<double> diel_tmp = ( *material )( labelfloor ) + ( ( *diel_old )( i ) - double(labelfloor) ) * ( ( *material )( labelnext ) - ( *material )( labelfloor ) );
        diel(i) = diel_tmp;
        al(i) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    }

    //cout << "al" << al(0) << endl;

    al_max = al;
    verbose = verbose_;
}

DDAModel::~DDAModel( ) {
    delete Core;
    Core = nullptr;
}

double DDAModel::calculateObjective( ) {
    return objDDAModel->GetVal( );
}

bool DDAModel::get_HaveDevx( ) {
    return objDDAModel->Have_Devx;
}

void DDAModel::SingleResponse(int idx, bool deduction) {
    objDDAModel->SingleResponse(idx, deduction);
}

double DDAModel::GroupResponse( ) {
    return objDDAModel->GroupResponse( );
}

void DDAModel::bicgstab(int MAX_ITERATION,double MAX_ERROR){
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }
    
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    //fftw_init_threads();                                                       //////////Initialize the multi-thread
    //fftw_plan_with_nthreads(NUM_THREADS);
    //cout<<"Threads"<<NUM_THREADS<<endl;;
    VectorXcd p=VectorXcd::Zero(N*3); 
    VectorXcd t = VectorXcd::Zero(N*3);
    VectorXcd w = VectorXcd::Zero(N*3);
    VectorXcd r = VectorXcd::Zero(N*3);
    VectorXcd r0 = VectorXcd::Zero(N*3);
    VectorXcd rl = VectorXcd::Zero(N*3);
    VectorXcd y = VectorXcd::Zero(N*3);
    VectorXcd u = VectorXcd::Zero(N*3);
    VectorXcd z = VectorXcd::Zero(N*3);
    VectorXcd x = VectorXcd::Zero(N*3);
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> eta;
    std::complex<double> zeta;

    VectorXcd Ax0 = Aproductwithalb(P);
    
    r = E-Ax0;
    r0 = r;
    p = r;

    VectorXcd Ap0 = Aproductwithalb(p);
    alpha = r0.dot(r)/r0.dot(Ap0);
    t = r-alpha*Ap0;

    VectorXcd At0 = Aproductwithalb(t);
    zeta = At0.dot(t)/At0.dot(At0);
    u = zeta*Ap0;
    z = zeta*r-alpha*u;
    P = P+alpha*p+z;                                    //this will directly change P in this.
    rl = r;
    r = t-zeta*At0;

    VectorXcd Ap = VectorXcd::Zero(N*3);
    VectorXcd At = VectorXcd::Zero(N*3);
    for (int it=0;it<=MAX_ITERATION-1;it++) {
        if (verbose && (it+1)%1000==0) {
            //cout << "r.norm(): " << r.norm() << endl;
            //cout << "E.norm(): " << E.norm() << endl;
            cout << "                Iter " << it+1 << ", error=" << r.norm()/E.norm() << " MAX_ERROR="<<MAX_ERROR<<endl;
        }

        complex<double> r0dotrl = r0.dot(rl);
        if (r0dotrl.real() == 0 && r0dotrl.imag() == 0) {
            Error = r.norm() / E.norm();
            //cout << "r.norm(): " << r.norm() << endl;
            //cout << "E.norm(): " << E.norm() << endl;
            ITERATION = it + 1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end - t_start).count();
                time = duration_cast<milliseconds>(t_end - t_start).count();
                cout << "----------------------------r0dotrl==(0,0), nan is going to occur so stop now------------------------" << endl;
                cout << "--------------Calculation finished. Duration: " << duration / 1000.0 << "s.-------------------" << endl;

                cout << "              Error: " << Error << endl;
                cout << "              Iteration: " << ITERATION << endl;
                cout << endl;
            }
            return;
        }

        beta = (alpha/zeta)*r0.dot(r)/r0.dot(rl);
        p = r+beta*(p-u);
        Ap = Aproductwithalb(p);
        alpha = r0.dot(r)/r0.dot(Ap);
        t = r-alpha*Ap;
        At = Aproductwithalb(t);

        zeta = At.dot(t)/At.dot(At);
        u = zeta*Ap;
        z = zeta*r-alpha*u;
        P = P+alpha*p+z;
        rl = r;
        r = t-zeta*At;

        if (r.norm()/E.norm()<=MAX_ERROR) {
            Error = r.norm()/E.norm();
            
            ITERATION = it+1;
            if (verbose) {
                cout << "r.norm(): " << r.norm() << endl;
                cout << "E.norm(): " << E.norm() << endl;
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end-t_start).count();
                time = duration_cast<milliseconds>(t_end-t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration/1000.0 << "s.-------------------" << endl;

                cout << "              Error: "<<Error<<endl;
                cout << "              Iteration: "<<ITERATION<<endl;
                cout << endl;
            }
            return;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end-t_start).count();
    cout<<"                ERROR:does not converge in "<<MAX_ITERATION<<" iterations"<<endl;
    return;
}

void DDAModel::bicgstab(int MAX_ITERATION, double MAX_ERROR, int EVOITERATION) {
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }
    //int N = (*Core).get_N();

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    //fftw_init_threads();                                                       //////////Initialize the multi-thread
    //fftw_plan_with_nthreads(NUM_THREADS);
    //cout<<"Threads"<<NUM_THREADS<<endl;;
    VectorXcd p = VectorXcd::Zero(N * 3);
    VectorXcd t = VectorXcd::Zero(N * 3);
    VectorXcd w = VectorXcd::Zero(N * 3);
    VectorXcd r = VectorXcd::Zero(N * 3);
    VectorXcd r0 = VectorXcd::Zero(N * 3);
    VectorXcd rl = VectorXcd::Zero(N * 3);
    VectorXcd y = VectorXcd::Zero(N * 3);
    VectorXcd u = VectorXcd::Zero(N * 3);
    VectorXcd z = VectorXcd::Zero(N * 3);
    VectorXcd x = VectorXcd::Zero(N * 3);
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> eta;
    std::complex<double> zeta;

    //Always starts with P=0 to avoid strange behaviour
    //P = VectorXcd::Zero(N * 3);

    ofstream foutnew(".\\p330-lam542-beta8-TiO2-InE-circle-fordebug\\BUGINFO.txt");

    VectorXcd Ax0 = Aproductwithalb(P);
    r = E - Ax0;
    r0 = r;
    p = r;
    VectorXcd Ap0 = Aproductwithalb(p);
    alpha = r0.dot(r) / r0.dot(Ap0);
    t = r - alpha * Ap0;
    VectorXcd At0 = Aproductwithalb(t);
    zeta = At0.dot(t) / At0.dot(At0);
    u = zeta * Ap0;
    z = zeta * r - alpha * u;
    P = P + alpha * p + z;                                    //this will directly change P in this.
    rl = r;
    r = t - zeta * At0;

    VectorXcd Ap = VectorXcd::Zero(N * 3);
    VectorXcd At = VectorXcd::Zero(N * 3);
    for (int it = 0; it <= MAX_ITERATION - 1; it++) {
        if (verbose && (it + 1) % 1000 == 0) {
            cout << "r.norm(): " << r.norm() << endl;
            cout << "E.norm(): " << E.norm() << endl;
            cout << "                Iter " << it + 1 << ", error=" << r.norm() / E.norm() << " MAX_ERROR=" << MAX_ERROR << endl;
        }
        beta = (alpha / zeta) * r0.dot(r) / r0.dot(rl);

        p = r + beta * (p - u);
        Ap = Aproductwithalb(p);
        alpha = r0.dot(r) / r0.dot(Ap);
        t = r - alpha * Ap;
        At = Aproductwithalb(t);

        zeta = At.dot(t) / At.dot(At);
        u = zeta * Ap;
        z = zeta * r - alpha * u;
        P = P + alpha * p + z;
        rl = r;
        r = t - zeta * At;

        if (r.norm() / E.norm() <= MAX_ERROR) {
            Error = r.norm() / E.norm();
            
            ITERATION = it + 1;
            if (verbose) {
                cout << "r.norm(): " << r.norm() << endl;
                cout << "E.norm(): " << E.norm() << endl;
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end - t_start).count();
                time = duration_cast<milliseconds>(t_end - t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration / 1000.0 << "s.-------------------" << endl;
                //ofstream fout;
                //fout.open("DDATime.txt", fstream::app);
                //fout<<N<<" "<<duration/1000.0<<endl;
                //fout << ITERATION << " " << duration / 1000.0 / ITERATION << endl;
                //fout.close();

                cout << "              Error: " << Error << endl;
                cout << "              Iteration: " << ITERATION << endl;
                cout << endl;
            }
            return;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end - t_start).count();
    cout << "                ERROR:does not converge in " << MAX_ITERATION << " iterations" << endl;

    foutnew.close();
    return;
}

void DDAModel::change_E(VectorXcd E_){
    E=E_;
}

void DDAModel::reset_E(){

    for (int i=0;i<N;i++) {
        E(3 * i) = E0 * n_E0(0) * (cos(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) + sin(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) * 1i);
        E(3 * i + 1) = E0 * n_E0(1) * (cos(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) + sin(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) * 1i);
        E(3 * i + 2) = E0 * n_E0(2) * (cos(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) + sin(K * d * (n_K(0) * (*geometry)(3 * i) + n_K(1) * (*geometry)(3 * i + 1) + n_K(2) * (*geometry)(3 * i + 2))) * 1i);
    }
}

void DDAModel::UpdateAlpha() { 

    for (int i = 0; i <= (*diel_old).size() - 1; i++) {
        int labelfloor = int(floor((*diel_old)(i)));
        int labelnext = labelfloor + 1;
        if (labelfloor >= 1) {
            labelnext = labelfloor;
        }
        std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(i) - double(labelfloor)) * ((*material)(labelnext) - (*material)(labelfloor));
        diel(i) = diel_tmp;
        al(i) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    }
}

void DDAModel::UpdateAlphaSingle(int idx) {

    int labelfloor = int(floor((*diel_old)(3 * idx)));
    int labelnext = labelfloor + 1;
    if (labelfloor >= 1) {
        labelnext = labelfloor;
    }
    std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(3 * idx) - double(labelfloor)) * ((*material)(labelnext) - (*material)(labelfloor));
    diel(3 * idx) = diel_tmp;
    diel(3 * idx + 1) = diel_tmp;
    diel(3 * idx + 2) = diel_tmp;
    al(3 * idx) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    al(3 * idx + 1) = al(3 * idx);
    al(3 * idx + 2) = al(3 * idx);

}

void DDAModel::solve_E(){
    if(RResultSwitch == true){

        for (int j = 0;j<int(round(RResult.size()/3));j++){
            double x = d*RResult(3*j);
            double y = d*RResult(3*j+1);
            double z = d*RResult(3*j+2);
            Vector3cd sum=Vector3cd::Zero();
            Vector3cd E_ext=Vector3cd::Zero();
            E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            for (int i=0;i<N;i++){
                double rx = x - d * (*geometry)(3 * i);                  //R has no d in it, so needs to time d
                double ry = y - d * (*geometry)(3 * i + 1);
                double rz = z - d * (*geometry)(3 * i + 2);
                Matrix3cd A=(*Core).A_dic_generator(rx,ry,rz);
                sum(0)+=(A(0,0)*P(3*i)+A(0,1)*P(3*i+1)+A(0,2)*P(3*i+2));
                sum(1)+=(A(1,0)*P(3*i)+A(1,1)*P(3*i+1)+A(1,2)*P(3*i+2));
                sum(2)+=(A(2,0)*P(3*i)+A(2,1)*P(3*i+1)+A(2,2)*P(3*i+2));
            }
            EResult(3*j) = E_ext(0)-sum(0);
            EResult(3*j+1) = E_ext(1)-sum(1);
            EResult(3*j+2) = E_ext(2)-sum(2);
        }
    }

    else
        EResult = Einternal;
  
}

void DDAModel::update_E_in_structure(){
    
    for(int i=0;i<=3*N-1;i++){
        std::complex<double> lorentzfactor = 2.0 + diel(i);
        lorentzfactor = lorentzfactor / 3.0;
        Einternal(i) = al(i) * P(i)/lorentzfactor;
    }
}

VectorXcd DDAModel::Aproductwithalb(VectorXcd& b) {
    if (b.size() != al.size()) {
        cout << "In Aproductwithalb, bsize not equal to alsize" << endl;
    }
    VectorXcd result = VectorXcd::Zero(b.size());
    for (int i = 0; i <= al.size() - 1; i++) {
        result(i) = b(i) * al(i);
    }
    return (*Core).Aproduct(b, geometry) + result;
}

void DDAModel::output_to_file(string save_position, int iteration, int ModelLabel){
    
    string name;
    name = save_position + "Model_results" + to_string(ModelLabel) + "it" + to_string(iteration) + ".txt";
    //name = save_position + "Model_output_verify\\" + "Model_results" + "it" + to_string(iteration) + ".txt";
    ofstream fout(name);

    for (int i = 0; i <= EResult.size() - 1; i++) {

        if (EResult(i).imag() < 0) {
            fout << EResult(i).real() << EResult(i).imag() << "j" << endl;
        }
        else {
            fout << EResult(i).real() << "+" << EResult(i).imag() << "j" << endl;
        }
    }

    for (int i = 0; i <= P.size() - 1; i++) {

        if (P(i).imag() < 0) {
            fout << P(i).real() << P(i).imag() << "j" << endl;
        }
        else {
            fout << P(i).real() << "+" << P(i).imag() << "j" << endl;
        }
    }

    fout.close();
}

void DDAModel::output_to_file(string save_position, int iteration) {

    string name;
    name = save_position + "Model_results" + "it" + to_string(iteration) + ".txt";
    //name = save_position + "Model_output_verify\\" + "Model_results" + "it" + to_string(iteration) + ".txt";
    ofstream fout(name);
    for (int i = 0; i <= EResult.size() - 1; i++) {

        if (EResult(i).imag() < 0) {
            fout << EResult(i).real() << EResult(i).imag() << "j" << endl;
        }
        else {
            fout << EResult(i).real() << "+" << EResult(i).imag() << "j" << endl;
        }
    }
    
    fout.close();
}


void DDAModel::InitializeP(VectorXcd& Initializer) {
    P = Initializer;
}
VectorXcd* DDAModel::get_P() {
    return &P;
}
Vector3d DDAModel::get_nE0() {
    return n_E0;
}
Vector3d DDAModel::get_nK() {
    return n_K;
}
double DDAModel::get_E0() {
    return E0;
}
VectorXcd* DDAModel::get_Einternal() {
    return &Einternal;
}
AProductCore* DDAModel::get_Core() {
    return Core;
}

VectorXcd* DDAModel::get_al() {
    return &al;
}

VectorXcd* DDAModel::get_al_max() {
    return &al_max;
}

VectorXcd* DDAModel::get_P_max() {
    return &P_max;
}

int DDAModel::get_ITERATION() {
    return ITERATION;
}

//----------------------------------------heritage from AProductCore------------------------------------

int DDAModel::get_N() {
    return N;
}
int DDAModel::get_Nx( ) {
    return Nx;
}
int DDAModel::get_Ny( ) {
    return Ny;
}
int DDAModel::get_Nz( ) {
    return Nz;
}
double DDAModel::get_lam( ) {
    return lam;
}
double DDAModel::get_d( ) {
    return d;
}
VectorXi* DDAModel::get_geometry( ) {
    return geometry;
}
VectorXd* DDAModel::get_diel_old( ) {
    return diel_old;
}
VectorXd* DDAModel::get_Para( ) {
    return Para;
}