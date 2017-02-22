//
// Created by lurker on 2/10/17.
//

#include "Matrix.h"


Matrix::~Matrix() {
    if (_ownership) {
        if (_row > 0 && _col > 0) {
            delete[] _data;
            _data = nullptr;
        }
    }
}

Matrix &Matrix::operator=(const Matrix &v) {
    if (_ownership) {
        if (_row > 0 && _col > 0) {
            delete[] _data;
            _data = nullptr;
        }
    }

    _row = v._row;
    _col = v._col;
    _ownership = v._ownership;

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
    return *this;
}


void Matrix::resize(int row, int col) {
    assert(_ownership);
    if (row != _row || col != _col) {
        if (_row > 0 && _col > 0) {
            delete[] _data;
            _data = nullptr;
        }
        _row = row;
        _col = col;
        if (_row > 0 && _col > 0) {
            _data = new scalar_t[_row * _col];
            assert(_data != nullptr);
            memset(_data, 0, _row * _col * sizeof(scalar_t));
        } else { _data = nullptr; }
    }
}

const scalar_t Matrix::operator()(int i, int j) const {
    assert(i >= 0 && i < _row);
    assert(j >= 0 && j <= _col);
    return _data[i + j * _row];
}

scalar_t &Matrix::operator()(int i, int j) {
    assert(i >= 0 && i < _row);
    assert(j >= 0 && j < _col);
    return _data[i + j * _row];
}

scalar_t *Matrix::data() { return _data; }

scalar_t *Matrix::column(int j) {
    assert(j >= 0 && j < _col);
    return &(_data[j * _row]);
}

int Matrix::row() { return _row; }

int Matrix::col() { return _col; }

void Matrix::setColumn(int j, Vector &v) {
    assert(j >= 0 && j < _col);
    assert(v.row() == _row);
    memcpy(&(_data[j * _row]), v.data(), sizeof(scalar_t) * _row);
}

inline std::ostream &operator<<(std::ostream &os, Matrix &v) {
    os << "rows: " << v._row << " cols: " << v._col << std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < v._row; ++i) {
        for (int j = 0; j < v._col; ++j) {
            os << " " << v(i, j);
        }
        os << std::endl;
    }
    return os;
}

inline void setValue(Matrix &v, scalar_t val) {
    std::fill_n(v.data(), v._row * v._col, val);
}

inline void clear(Matrix &v) {
    memset(v._data, 0, v._row * v._col * sizeof(scalar_t));
}

inline void setBlock(Matrix &A, Matrix &B, int r, int c, int rs, int cs) {
    for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < cs; ++j) {
            A(i, j) = B(r + i, c + j);
        }
    }
}
