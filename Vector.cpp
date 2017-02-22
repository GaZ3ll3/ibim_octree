//
// Created by lurker on 2/10/17.
//

#include "Vector.h"

Vector::~Vector() {
    if (_ownership) {
        if (_row > 0) {
            delete[] _data;
            _data = nullptr;
        }
    }
}

Vector &Vector::operator=(const Vector &v) {
    if (_ownership) {
        if (_row > 0) {
            delete[] _data;
            _data = nullptr;
        }
    }

    _row = v._row;
    _ownership = v._ownership;

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
    return *this;
}

void Vector::resize(int row) {
    assert(_ownership);
    if (row != _row) {
        if (_row > 0) {
            delete[] _data;
            _data = nullptr;
        }
        _row = row;
        if (_row > 0) {
            _data = new scalar_t[_row];
            assert(_data != nullptr);
            memset(_data, 0, _row * sizeof(scalar_t));
        }
        else { _data = nullptr; }
    }
}

const scalar_t Vector::operator()(int i) const {
    assert(i >= 0 && i < _row);
    return _data[i];
}

scalar_t &Vector::operator()(int i) {
    assert(i >= 0 && i < _row);
    return _data[i];
}

scalar_t *Vector::data() {
    return _data;
}

int Vector::row() {
    return _row;
}

inline std::ostream &operator<<(std::ostream &os, Vector &v) {
    os << "rows: " << v._row << std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < v._row; ++i) {
        os << " " << v(i);
    }
    os << std::endl;
    return os;
}

inline void setValue(Vector &v, scalar_t val) {
    std::fill_n(v.data(), v._row, val);
}

inline void clear(Vector &v) {
    memset(v._data, 0, v._row * sizeof(scalar_t));
}