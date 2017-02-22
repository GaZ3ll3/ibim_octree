//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_MATRIX_H
#define IBIM_OCTREE_MATRIX_H

#include "utils.h"
#include "Vector.h"

class Matrix {
public:
    int _row;
    int _col;
    bool _ownership;
    scalar_t *_data;


    Matrix(int row = 0, int col = 0) : _row(row), _col(col), _ownership(true) {
        if (_row > 0 && _col > 0) {
            _data = new scalar_t[_row * _col];
            assert(_data != nullptr);
            memset(_data, 0, _row * _col * sizeof(scalar_t));
        } else { _data = nullptr; }
    }

    Matrix(int row, int col, bool_t ownership, scalar_t *data) : _row(row), _col(col), _ownership(ownership) {
        if (_ownership) {
            if (_row > 0 && _col > 0) {
                _data = new scalar_t[_row * _col];
                assert(_data != nullptr);
                memcpy(_data, data, _row * _col * sizeof(scalar_t));
            } else {
                _data = nullptr;
            }
        } else {
            /*
             * has to assure data holds exact "row" elements.
             */
            _data = data;
        }
    }


    Matrix(const Matrix &v) : _row(v._row), _col(v._col), _ownership(v._ownership) {
        if (_ownership) {
            if (_row > 0 && _col > 0) {
                _data = new scalar_t[_row * _col];
                assert(_data != nullptr);
                memcpy(_data, v._data, _row * _col * sizeof(scalar_t));
            } else {
                _data = nullptr;
            }
        } else {
            _data = v._data;
        }
    }

    ~Matrix();

    Matrix &operator=(const Matrix &v);

    void resize(int row, int col);

    const scalar_t operator()(int i, int j) const;

    scalar_t &operator()(int i, int j);

    scalar_t *data();

    scalar_t *column(int j);

    int row();

    int col();

    void setColumn(int j, Vector &v);
};

inline std::ostream &operator<<(std::ostream &os, Matrix &v);
inline void setValue(Matrix &v, scalar_t val);
inline void clear(Matrix &v);
inline void setBlock(Matrix &A, Matrix &B, int r, int c, int rs, int cs);

#endif //IBIM_OCTREE_MATRIX_H
