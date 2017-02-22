//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_NODE_H
#define IBIM_OCTREE_NODE_H

#include "point.h"

class node {
public:
    index_t parent;
    index_t child[8];
    vector<index_t> neighbors;
    vector<index_t> atoms;
    index_t nLevel;
    index_t nodeIndex;


    point center;
    scalar_t  radius;
    scalar_t  value;

    bool onBoundary;
    bool isLeaf;

    node(index_t level, index_t index) {
        parent = -1;
        for (index_t i = 0; i < 8; ++i) {
            child[i] = -1;
        }
        nLevel = level;
        nodeIndex = index;
        isLeaf = true;
        onBoundary = false;
    }

    virtual ~node() {}
};

#endif //IBIM_OCTREE_NODE_H
