#ifndef TOPO_SPACE_H_
#define TOPO_SPACE_H_

#include "Eigen/Core"

#include "Structure.h"
using namespace std;
using namespace Eigen;

class Space {
private:
    VectorXi* total_space;
    int Nx, Ny, Nz;
    int N;                        //total size of all the geometry inside the list(dipole size which is actual size/3)
    vector<Structure>* ln;
    Structure* structure;
public:

    // new constructor to avoid creating total_space, and vector of structures (just using 1)
    Space(int Nx_, int Ny_, int Nz_, int N_, Structure* structure_);
    Space(VectorXi* total_space_, int Nx_, int Ny_, int Nz_, int N_, vector<Structure>* ln_);
    VectorXi* get_total_space();
    int get_ln_size();
    tuple<int, int, int, int> get_Ns();
    vector<Structure>* get_ln();
    Structure* get_structure();
    void show_something_about_Structures() const;
    friend Space operator+(const Space& s1, Structure& s2);


};

#endif

