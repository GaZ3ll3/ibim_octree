//
// Created by lurker on 2/11/17.
//

#include "tree.h"

#define SQR(X) ((X)*(X))

/*
 * build tree.
 *
 * root level node. the value of F at center should not be necessary, since
 * there is no neighbor.
 *
 *
 * usage:
 *
 *   tree t;
 *   RUN("build tree", t.build(mol, 0.05, 0.7, 0.));
 *   view& v= view::getInstance(true);
 *   v.loadTree(t);
 *   v.run();
 *
 *
 */
void tree::build(Molecule& m, Grid grid, scalar_t dmin, scalar_t eps, scalar_t iso_val, scalar_t pr) {
    molecule = &m; molecule->getCenter();
    resolution = dmin;
    g = grid;

    _center = molecule->center;
    _radius = molecule->radius;

    minLevel = index_t(log(_radius / eps) / log(2));

    level = 0;
    root = 0;

    _dict.push_back(node(0, 0));
    _maxId = root;

    _dict[root].center = _center;
    _dict[root].radius = HALF_SIZE;
    _dict[root].atoms.resize(molecule->centers.size());

    // assign value to root
    _dict[root].value = coarse_evaluate(root);

    for (index_t i = 0; i < molecule->centers.size(); ++i) {
        _dict[root].atoms[i] = i;
    }

    _newNodes.insert(root);
    _markedNodes.clear();

    while(_mark(iso_val)) {

        _newNodes.clear();



        for (auto marked_itr = _markedNodes.begin(); marked_itr != _markedNodes.end(); ++marked_itr) {
            /*
             * *marked_itr is the iterator of marked nodes.
             *
             * after splitting, marked nodes must be cleared.
             *
             */
            _split(*marked_itr);

        }

        /*
         * get neighbors information for new nodes.
         */
        for (auto new_iter = _newNodes.begin(); new_iter != _newNodes.end(); ++new_iter) {
            _neighbor(*new_iter);
        }

        _markedNodes.clear();

    }
}

/*
 * mark all nodes to be _split, each time we mark, we go one level down of the tree.
 *
 * neighbors are only valid for the same level.
 */
bool tree::_mark(scalar_t iso_val) {
    bool marked = false;

    for (auto itr = _newNodes.begin(); itr != _newNodes.end(); ++itr) {
        if (_dict[*itr].nLevel < minLevel) {
            _markedNodes.insert(*itr);
            marked = true;
        }
        /*
         * consider if current node and its neighbors are separated by contour.
         */
        for (auto neighbor_itr =  _dict[*itr].neighbors.begin();
             neighbor_itr != _dict[*itr].neighbors.end(); ++neighbor_itr) {
            if ((_dict[*itr].value - iso_val) * (_dict[*neighbor_itr].value - iso_val) <= 0) {
                _markedNodes.insert(*itr);
                _markedNodes.insert(*neighbor_itr);
                _dict[*itr].onBoundary = true;
                _dict[*neighbor_itr].onBoundary = true;
            }
            marked = true;
        }

    }

    return marked;
}

/*
 * split node into 8 subnodes.
 */
void tree::_split(index_t id) {
    if (_dict[id].radius < 2 * resolution){
        return;
    }
    else {
        for (index_t i = 0; i < 8 ; ++i) {
            _maxId += 1;
            _dict[id].child[i] = _maxId;
            /*
             * new node is automatically leaf.
             */
            _dict.push_back(node(_dict[id].nLevel + 1, i));

            if (level < _dict[_maxId].nLevel) level = _dict[_maxId].nLevel;

            _dict[_maxId].parent = id;
            _dict[_maxId].center.data[0] = _dict[id].center.data[0] + ((i & 1) - 0.5) * _dict[id].radius;
            _dict[_maxId].center.data[1] = _dict[id].center.data[1] + (((i >> 1) & 1) - 0.5) * _dict[id].radius;
            _dict[_maxId].center.data[2] = _dict[id].center.data[2] + ((i >> 2) - 0.5) * _dict[id].radius;
            _dict[_maxId].radius = _dict[id].radius * 0.5;

            for (index_t atomId : _dict[id].atoms) {
                if (_intersect(_maxId, molecule->centers[atomId], molecule->radii[atomId])) {
                    _dict[_maxId].atoms.push_back(atomId);
                }
            }

            /*
             * assign value to current node
             */
            _dict[_maxId].value = coarse_evaluate(_maxId);

            _newNodes.insert(_maxId);
        }

        _dict[id].atoms.clear();
        _dict[id].isLeaf = false;

    }
}

/*
 * calculate shortest distance between box and sphere.
 */
bool tree::_intersect(index_t id, point &c, scalar_t r) {
    scalar_t r2 = r * r;
    scalar_t dmin = 0.;
    point bmax = _dict[id].center + _dict[id].radius;
    point bmin = _dict[id].center - _dict[id].radius;

    for (index_t i = 0; i < DIM; ++i) {
        if (bmax.data[i] < c.data[i]) dmin += (bmax.data[i] - c.data[i]) * (bmax.data[i] - c.data[i]);
        else if (bmin.data[i] > c.data[i]) dmin += (bmin.data[i] - c.data[i]) * (bmin.data[i] - c.data[i]);
    }
    return dmin <= r2;
}

/*
 * locate neighbors of node. in O(26logN) time.
 */
void tree::_neighbor(index_t id) {
    index_t _level = _dict[id].nLevel;
    scalar_t r = _dict[id].radius;

    for (index_t i = -1; i < 2; ++i) {
        for (index_t j = -1; j < 2; ++j) {
            for (index_t k = -1; k < 2; ++k) {
                point cur_p = {
                        _dict[id].center.data[0] + 2 * i * r,
                        _dict[id].center.data[1] + 2 * j * r,
                        _dict[id].center.data[2] + 2 * k * r
                };

                index_t ret = _find(cur_p, 0);

                if (ret != -1 && ret != id) {
                    if (_dict[ret].nLevel == _level) {
                        _dict[id].neighbors.push_back(ret);
                    }
                }
            }
        }
    }
}

/*
 * find the location of node containing point p.
 *
 * it is iterative scheme, with time cost logN.
 *
 * constant time scheme also exists by looking for parent's neighbors' children.
 */
index_t tree::_find(point &p, index_t _parent) {
    node& n = _dict[_parent];
    if (n.center == p) return _parent;
    else {
        if (n.isLeaf) {
            /*
             * if it is inside this node
             */
            if ((n.center + n.radius) >= p && (n.center - n.radius) <= p) {
                return _parent;
            }
            else {
                /*
                 * otherwise it is outside of the box
                 */
                return -1;
            }
        }
        else {
            index_t x_bit = n.center.data[0] > p.data[0] ? 0 : 1;
            index_t y_bit = n.center.data[1] > p.data[1] ? 0 : 1;
            index_t z_bit = n.center.data[2] > p.data[2] ? 0 : 1;
            index_t id = 4 * z_bit + 2 *y_bit + x_bit;

            return _find(p, n.child[id]);
        }
    }
}

/*
 * calculate the value of a point.
 */
scalar_t tree::coarse_evaluate(index_t id) {
    auto s = g->stepSize;

    auto icx = (int) ((_dict[id].center.data[0] + HALF_SIZE) / s);
    auto icy = (int) ((_dict[id].center.data[1] + HALF_SIZE) / s);
    auto icz = (int) ((_dict[id].center.data[2] + HALF_SIZE) / s);

    return g->matrix[icx][icy][icz].phi;
}


