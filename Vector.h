//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_VECTOR_H
#define IBIM_OCTREE_VECTOR_H

#include "utils.h"

class Vector {
public:
    /*
     * data members
     */
    index_t _row;
    bool_t _ownership;
    scalar_t *_data;

    Vector(int row = 0) : _row(row), _ownership(true) {
        if (_row > 0) {
            _data = new scalar_t[_row];
            assert(_data != nullptr);
            memset(_data, 0, _row * sizeof(scalar_t));
        }
        else { _data = nullptr; }
    }

    Vector(int row, bool_t ownership, scalar_t *data) : _row(row), _ownership(ownership) {
        if (_ownership) {
            if (_row > 0) {
                _data = new scalar_t[_row];
                assert(_data != nullptr);
                memcpy(_data, data, _row * sizeof(scalar_t));
            }
            else {
                _data = nullptr;
            }
        } else {
            /*
             * has to assure data holds exact "row" elements.
             */
            _data = data;
        }
    }


    Vector(const Vector &v) : _row(v._row), _ownership(v._ownership) {
        if (_ownership) {
            if (_row > 0) {
                _data = new scalar_t[_row];
                assert(_data != nullptr);
                memcpy(_data, v._data, _row * sizeof(scalar_t));
            }
            else {
                _data = nullptr;
            }
        } else {
            _data = v._data;
        }
    }

    ~Vector();

    Vector &operator=(const Vector &v);

    void resize(int row);

    const scalar_t operator()(int i) const;

    scalar_t &operator()(int i);

    scalar_t *data();

    int row();
};

inline std::ostream &operator<<(std::ostream &os, Vector &v);
inline void setValue(Vector &v, scalar_t val);
inline void clear(Vector &v);

#endif //IBIM_OCTREE_VECTOR_H
