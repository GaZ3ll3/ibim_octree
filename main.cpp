#define DISP

#include <iostream>
#include "levelset.h"


int main() {


    Molecule mol;

//    mol.load("../data/2aid.pqr");
    mol.load("../data/1hje.pqr");
    scalar_t s = mol.centralize(200.0); mol.getCenter();
    scalar_t pr = 1.4 * s;

    scalar_t grid_lo = -300.0;
    scalar_t grid_hi = 300.0;

    index_t size = 100;
    scalar_t dx = (grid_hi - grid_lo) / (size);

    levelset ls(size, size, size, grid_lo, grid_lo, grid_lo, dx, s);

//    view& v= view::getInstance(2);

    RUN("OUTWARD", ls.expand(mol, pr));


    RUN("INWARD",ls.flow(ls.phi, 1.4, 20, 0.5));

    double* phi0 = (double*)malloc(ls.Nx * ls.Ny * ls.Nz * sizeof(double));

    for (int i = 0; i < ls.Nx * ls.Ny * ls.Nz; ++i) {
        phi0[i] = -ls.phi[i];
        ls.phi[i] = phi0[i];
    }

    RUN("REINIT", ls.reinit(ls.phi, phi0, 1.4, 20, 0.5));

//    v.loadLevelSet(ls);
//    v.run();


    free(phi0);



}