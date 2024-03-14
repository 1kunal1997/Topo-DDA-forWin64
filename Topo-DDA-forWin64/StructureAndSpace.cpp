#include <iostream>
#include <set>

#include "StructureAndSpace.h"

StructureAndSpace::StructureAndSpace(VectorXi* geometry_, int Nx_, int Ny_, int Nz_, int N_) {
	geometry = *geometry_;
	Nx = Nx_;
	Ny = Ny_;
	Nz = Nz_;
	N = N_;
}
VectorXi* StructureAndSpace::get_geometry() {
	return &geometry;
}

int StructureAndSpace::get_geometry_size() {
	return round(geometry.size( ) / 3);
}
tuple<int, int, int, int> StructureAndSpace::get_Ns( ) {
	return make_tuple(Nx, Ny, Nz, N);
}