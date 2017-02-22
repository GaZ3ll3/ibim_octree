//
// Created by lurker on 2/14/17.
//

#ifndef IBIM_OCTREE_GRID_H
#define IBIM_OCTREE_GRID_H

#include "utils.h"
#include "molecule.h"

#define MAX_SIZE 512
#define HALF_SIZE 256

/*
 * a uniform grid to capture the detailed (but coarse)
 * structure information of the molecule.
 *
 * all grid are integer based to save memory.
 *
 * grid size is pre-defined.
 *
 */
typedef struct nb
{
    short x;
    short y;
    short z;
    struct nb* next;
} nbNode;


typedef struct p
{
    short x;
    short y;
    short z;
} PPoint;

typedef struct gp
{
    PPoint point;
    int phi;
    int from;
    float dist;
} GridPoint;

typedef struct g
{
    short N;
    GridPoint*** matrix;
    short stepSize;
} *Grid;

Grid createGrid(short);
void signDistanceGridMol(Grid, Molecule&, double);
void shrink(Grid, double);
int fastMarching(Grid, char);
int convexity(Grid,Molecule&);
void probing(Grid, Molecule&);




#endif //IBIM_OCTREE_GRID_H
