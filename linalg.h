//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_LINALG_H
#define IBIM_OCTREE_LINALG_H

#include "Vector.h"
#include "Matrix.h"

void dscal(scalar_t alpha, Vector &v);
void dscal(scalar_t alpha, Matrix &v);
void daxpy(scalar_t a, Vector &x, Vector &y);
void daxpy(scalar_t a, Vector &x, scalar_t *y);
void daxpy(scalar_t a, Matrix &x, Matrix &y);
void dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C);
void dgemm_t(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C);
void t_dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C);
void dger(scalar_t alpha, Vector &X, Vector &Y, Matrix &A);
void dgemv(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y);
void dgemv_t(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y);
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, Vector &y);
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, scalar_t *y);
scalar_t nrm2(Vector &x);
scalar_t ddot(Vector &x, Vector &y);


#endif //IBIM_OCTREE_LINALG_H