#define DISP

#include <iostream>
#include "view.h"


int main() {


    Molecule mol; mol.load("/home/lurker/ClionProjects/imbm/1hje.pqr");
    scalar_t s = mol.centralize(200.0); mol.getCenter();
    scalar_t pr = 1.4 * s;

    scalar_t grid_lo = -300.0;
    scalar_t grid_hi = 300.0;

    index_t size = 100;
    scalar_t dx = (grid_hi - grid_lo) / (size);

    levelset ls(size, size, size, grid_lo, grid_lo, grid_lo, dx, s);

    view& v= view::getInstance(2);

    std::chrono::steady_clock::time_point begin =std::chrono::steady_clock::now(); \

    for (index_t aId = 0; aId < mol.centers.size();++aId) {
        int r = (int)((mol.radii[aId] + pr)/dx) + ls.WENO;
        int icx = (int)((mol.centers[aId].data[0] - ls.sx) / dx);
        int icy = (int)((mol.centers[aId].data[1] - ls.sy) / dx);
        int icz = (int)((mol.centers[aId].data[2] - ls.sz) / dx);
        for (int a = icx - r; a <= icx + r; ++a) {
            for (int b = icy - r; b <= icy + r; ++b) {
                for (int c = icz - r; c<= icz + r; ++c) {
                    point grid_p = {ls.sx + a * dx, ls.sy + b * dx, ls.sz + c * dx};
                    scalar_t dist = (mol.radii[aId] + pr) - norm(grid_p - mol.centers[aId]) ;

                    dist = max(dist, ls.get(ls.phi, a, b,c));
                    ls.set(dist, ls.phi, a, b, c);
                }
            }
        }
    }

    std::chrono::steady_clock::time_point end =  end = std::chrono::steady_clock::now();\
std::cout << std::setw(15)<< "OUTWARD" << " "  << std::setprecision(5) << std::setw(8) <<\
 std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << " seconds"<<std::endl;


    RUN("INWARD",ls.flow(ls.phi, 1.4, 20, 0.5));

    double* phi0 = (double*)malloc(ls.Nx * ls.Ny * ls.Nz * sizeof(double));

    for (int i = 0; i < ls.Nx * ls.Ny * ls.Nz; ++i) {
        phi0[i] = -ls.phi[i];
        ls.phi[i] = phi0[i];
    }

    RUN("REINIT", ls.reinit(ls.phi, phi0, 1.4, 20, 0.5));

    v.loadLevelSet(ls);
    v.run();


    free(phi0);



}