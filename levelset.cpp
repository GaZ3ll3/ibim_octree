//
// Created by lurker on 2/18/17.
//

#include "levelset.h"

levelset::levelset(index_t X, index_t Y, index_t Z,
                   double SX, double SY, double SZ,
                   double DX, double VEL) {
    Nx = X; Ny = Y; Nz = Z; dx = DX; vel = VEL;
    sx = SX; sy = SY; sz = SZ;
    phi = (double*)malloc(Nx * Ny * Nz * sizeof(double));
    WENO = 7; HALF_WENO = WENO / 2;

    std::fill(phi, phi + Nx * Ny * Nz, -512.0);
}

levelset::~levelset() {
    free(phi);
}


void levelset::expand(Molecule &mol, double pr) {
    for (index_t aId = 0; aId < mol.centers.size();++aId) {
        int r = (int)((mol.radii[aId] + pr)/dx) + WENO;
        int icx = (int)((mol.centers[aId].data[0] - sx) / dx);
        int icy = (int)((mol.centers[aId].data[1] - sy) / dx);
        int icz = (int)((mol.centers[aId].data[2] - sz) / dx);
        for (int a = icx - r; a <= icx + r; ++a) {
            for (int b = icy - r; b <= icy + r; ++b) {
                for (int c = icz - r; c<= icz + r; ++c) {
                    point grid_p = {sx + a * dx, sy + b * dx, sz + c * dx};
                    scalar_t dist = (mol.radii[aId] + pr) - norm(grid_p - mol.centers[aId]) ;

                    dist = max(dist, get(phi, a, b,c));
                    set(dist, phi, a, b, c);
                }
            }
        }
    }
}

double levelset::getNorm(point &Dun, point &Dup) {
    for (int i = 0; i < DIM; ++i) {
        Dun.data[i] = max((Dun.data[i] > 0. ? Dun.data[i] : 0.), (Dup.data[i] < 0. ? -Dup.data[i] : 0.));
    }
    return norm(Dun);
}

double levelset::getUpwind(point &p, point &Dun, point &Dup) {
    double val = 0.;
    for (int i = 0; i < DIM; ++i) {
        if (p.data[i] >= 0) {
            val += p.data[i] * Dun.data[i];
        }
        else {
            val += p.data[i] * Dup.data[i];
        }
    }
    return val;
}

double levelset::get(double* u, index_t i, index_t j, index_t k) {
    i = min(max(i, 0), Nx - 1);
    j = min(max(j, 0), Ny - 1);
    k = min(max(k, 0), Nz - 1);
    return u[i * Ny * Nz + j * Nz + k];
}

void levelset::set(double val , double *u, index_t i, index_t j, index_t k) {
    u[i * Ny * Nz + j * Nz + k] = val;
}


void levelset::getWindows(double* u, double* window, index_t i, index_t j, index_t k) {
    index_t _loc = 0;
    for (index_t loc = 0; loc < WENO; ++loc) {
        window[_loc++] = get(u, i - HALF_WENO + loc, j, k);
    }

    for (index_t loc = 0; loc < WENO; ++loc) {
        window[_loc++] = get(u, i, j - HALF_WENO + loc, k);
    }

    for (index_t loc = 0; loc < WENO; ++loc) {
        window[_loc++] = get(u, i, j, k - HALF_WENO + loc);
    }
}

void levelset::getGrad_WENO5(index_t dir, double* window, double& uxp, double& uxn) {
    index_t start = dir * WENO;
    double diff[WENO - 1];
    double df, is0, is1, is2, a0, a1, a2, w0, w1, w2;

    for (index_t i = 0; i < WENO - 1; ++i) {
        diff[i] = window[start + i + 1] - window[start + i];
    }

    df = (-diff[1] + 7.0 * (diff[2] + diff[3]) - diff[4]) / 12.0;

    for (index_t i = 0; i < WENO - 2; ++i) {
        diff[i] = diff[i + 1] - diff[i];
    }

    is0 = 13.0 * SQR(diff[4] - diff[3]) + 3.0 * SQR(diff[4] - 3.0 * diff[3]);
    is1 = 13.0 * SQR(diff[3] - diff[2]) + 3.0 * SQR(diff[3] + diff[2]);
    is2 = 13.0 * SQR(diff[2] - diff[1]) + 3.0 * SQR(3.0 * diff[2] - diff[1]);

    a0 = 1.0 / (SQR(EPS + is0));
    a1 = 6.0 / (SQR(EPS + is1));
    a2 = 3.0 / (SQR(EPS + is2));

    w1 = 1.0 / (a0 + a1 + a2);
    w0 = a0 * w1;
    w2 = a2 * w1;

    uxp = df +
            w0 * ( diff[4] - 2.0 * diff[3] + diff[2]) / 3.0 +
            (w2 - 0.5) * (diff[3] - 2.0 * diff[2] + diff[1]) / 6.0;

    is0 = 13.0 * SQR(diff[0] - diff[1]) + 3.0 * (diff[0] - 3.0 * diff[1]);
    is1 = 13.0 * SQR(diff[1] - diff[2]) + 3.0 * (diff[1] + diff[2]);
    is2 = 13.0 * SQR(diff[2] - diff[3]) + 3.0 * (3.0 * diff[2] - diff[3]);

    a0 = 1.0 / (SQR(EPS + is0));
    a1 = 6.0 / (SQR(EPS + is1));
    a2 = 3.0 / (SQR(EPS + is2));

    w1 = 1.0 / (a0 + a1 + a2);
    w0 = a0 * w1;
    w2 = a2 * w1;

    uxn = df -
            w0 * (diff[0] - 2 * diff[1] + diff[2]) / 3.0 -
            (w2 - 0.5) * (diff[1] - 2.0 * diff[2] + diff[3]) / 6.0;

    uxp /= dx;
    uxn /= dx;
}

void levelset::flow(double *u0, double final_t, index_t num_steps, double cfl_thres) {
    assert(num_steps > 0);
    double dt = final_t / (double)num_steps;

    bool valid = false;
    while (dt > cfl_thres * dx / vel) {
        num_steps *= 2;
        dt *= 0.5;
        valid = true;
    }

    if (valid) {std::cout << "adapt dt and steps to dt : " << dt << " and steps : "  << num_steps << std::endl;}

    /*
     * RK3, careful with the memory management.
     */
    double* u1 = (double*)malloc(Nx * Ny * Nz * sizeof(double));
    double* u2 = (double*)malloc(Nx * Ny * Nz * sizeof(double));

    index_t step = 0;
    point Dup, Dun;

    const int core = omp_get_max_threads();

    omp_set_num_threads(core);
    double** window = (double**)malloc(core * sizeof(double*));

    for (int i = 0; i < core; ++i) {
        window[i] = (double*)malloc(DIM * WENO * sizeof(double));
    }

    for (step = 0; step < num_steps; ++step) {
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;

                    index_t tid = omp_get_thread_num();
                    getWindows(u0, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }
                    u1[I] = u0[I] - dt * vel * getNorm(Dun, Dup);
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;

                    index_t tid = omp_get_thread_num();
                    getWindows(u1, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }
                    u2[I] = (3 * u0[I] + u1[I] - dt * vel * getNorm(Dun, Dup)) / 4.0;
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;

                    index_t tid = omp_get_thread_num();
                    getWindows(u2, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }
                    u0[I] = (u0[I] + 2 * (u2[I] - dt * vel * getNorm(Dun, Dup))) / 3.0;
                }
            }
        }
#pragma omp barrier
    }

    free(u1);free(u2);


    for (int i = 0; i < core; ++i) {
        free(window[i]);
    }
    free(window);
}


void levelset::reinit(double*u, double *u0, double final_t, index_t num_steps, double cfl_thres) {
    double const eps = dx * dx;

    assert(num_steps > 0);
    double dt = final_t / (double)num_steps;

    bool valid = false;
    while (dt > cfl_thres * dx / vel) {
        num_steps *= 2;
        dt *= 0.5;
        valid = true;
    }

    if (valid) {std::cout << "adapt dt and steps to dt : " << dt << " and steps : "  << num_steps << std::endl;}

    point Dup, Dun;

    /*
     * RK3, careful with the memory management.
     */
    double* u1 = (double*)malloc(Nx * Ny * Nz * sizeof(double));
    double* u2 = (double*)malloc(Nx * Ny * Nz * sizeof(double));

    index_t step = 0;

    const int core = omp_get_max_threads();

    omp_set_num_threads(core);

    double** window = (double**)malloc(core * sizeof(double*));

    for (int i = 0; i < core; ++i) {
        window[i] = (double*)malloc(DIM * WENO * sizeof(double));
    }

    for (step = 0; step < num_steps; ++step) {
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3) num_threads(4)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    /*
                     * get windows for all directions.
                     */
                    index_t tid = omp_get_thread_num();
                    getWindows(u0, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }

                    double sign = u0[I] / (sqrt(SQR(u0[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));
                    u1[I] = u[I] - dt * sign * (normDu - 1.0);
                }
            }
        }
        // second stage
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3)  num_threads(4)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    /*
                     * get windows for all directions.
                     */
                    index_t tid = omp_get_thread_num();
                    getWindows(u1, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }

                    double sign = u0[I] / (sqrt(SQR(u0[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));
                    u2[I] = (3 * u[I] + u1[I] - dt * sign * (normDu - 1.0)) / 4.0;
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static) collapse(3)  num_threads(4)
        // third stage
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    /*
                     * get windows for all directions.
                     */
                    index_t tid = omp_get_thread_num();
                    getWindows(u2, window[tid], i, j, k);


                    for (index_t dir = 0; dir < DIM; ++dir) {
                        getGrad_WENO5(dir, window[tid], Dup.data[dir], Dun.data[dir]);
                    }

                    double sign = u0[I] / (sqrt(SQR(u0[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));
                    u[I] = (u[I] + 2 * (u2[I] - dt * sign * (normDu - 1.0))) / 3.0;
                }
            }
        }
#pragma omp barrier
    }

    free(u1);free(u2);

    for (int i = 0; i < core; ++i) {
        free(window[i]);
    }
    free(window);
}