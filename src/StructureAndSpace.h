#ifndef TOPO_STRUCTUREANDSPACE_H_
#define TOPO_STRUCTUREANDSPACE_H_

#include <vector>

#include <Eigen/Core>

using namespace std;
using namespace Eigen;
class StructureAndSpace {

private:
	VectorXi geometry;
	int Nx, Ny, Nz, N;

public:

	StructureAndSpace(VectorXi* geometry_, int Nx_, int Ny_, int Nz_, int N_);

	VectorXi* get_geometry();
	int get_geometry_size();
	tuple<int, int, int, int> get_Ns();

};
#endif

