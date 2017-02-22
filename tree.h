//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_OCTREE_H
#define IBIM_OCTREE_OCTREE_H

#include "node.h"
#include "molecule.h"
#include "grid.h"

class tree {
public:
    index_t root;
    index_t level;
    index_t minLevel;
    scalar_t resolution;
    Molecule* molecule;
    Grid g;



    /*
     * initialization of tree. No nodes assigned.
     *
     *  _maxId : current node numbering.
     *  root : root of the tree. if exists, it is 0, otherwise it is -1.
     *  level: depth of the current tree, if tree is null, it is -1.
     *  resolution: it is assigned to be any number now. should come from user input.
     *  molecule: pointer of the molecule instance.
     *
     */
    tree() {
        _maxId = -1;
        root = -1;
        level = -1;
        resolution = std::numeric_limits<scalar_t >::infinity();
        molecule = nullptr;
    }

    ~tree() {}

    /*
     *  public methods
     *
     *  1. "build", construction of tree.
     *  2. "coarse_evaluate", calculate iso-value.
     *
     *
     *  protected methods
     *
     *  1. "_split", split current node.
     *  2. "_mark", mark the nodes to be _split.
     *  3. "_intersect", decide whether box is intersected with atom.
     *  4. "_neighbor", locate all neighbors (of same level) of a node.
     *  5. "_find", find node id containing the point, iteratively.
     *
     */

    void build(Molecule& m, Grid grid, scalar_t dmin, scalar_t eps, scalar_t iso_val, scalar_t pr);
    scalar_t coarse_evaluate(index_t id);


    vector<node> _dict;
    index_t _maxId;
    point _center;
    scalar_t _radius;

    unordered_set<index_t> _newNodes;
    unordered_set<index_t> _markedNodes;

    bool _mark(scalar_t iso_val);
    void _split(index_t id);
    bool _intersect(index_t id, point &c, scalar_t r);
    void _neighbor(index_t id);
    index_t _find(point& p, index_t _parent);
    void _refine();



};



#endif //IBIM_OCTREE_OCTREE_H
