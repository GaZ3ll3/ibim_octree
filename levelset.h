//
// Created by lurker on 2/18/17.
// for level set evolution.
//
// Grid cannot be too large.
//

#ifndef IBIM_OCTREE_LEVELSET_H
#define IBIM_OCTREE_LEVELSET_H


#include "point.h"
#include "molecule.h"

class levelset {
public:
    double* phi;
    double vel;
    double dx;
    double sx;
    double sy;
    double sz;
    index_t Nx;
    index_t Ny;
    index_t Nz;
    index_t WENO;
    index_t HALF_WENO;



    levelset(index_t X, index_t Y, index_t Z, double SX, double SY, double SZ, double DX, double VEL);
    ~levelset();
    double getNorm(point& Dun, point& Dup);
    double getUpwind(point& p, point& Dun, point& Dup);
    void getWindows(double* u, double* window, index_t i, index_t j, index_t k);
    void getGrad_WENO5(index_t dir, double* window, double& uxp, double& uxn);

    double get(double* u, index_t i, index_t j, index_t k);
    void set(double val, double* u, index_t i, index_t j, index_t k);

    void expand(Molecule& m, double pr);

    void flow(double* phi0, double final_t, index_t num_steps, double cfl);
    void reinit(double* u, double* phi0, double final_t, index_t num_steps, double cfl);

};


#endif //IBIM_OCTREE_LEVELSET_H
