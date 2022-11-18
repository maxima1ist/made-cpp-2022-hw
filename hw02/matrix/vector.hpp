#pragma once

#include <cstddef>
#include <array>
#include <stdexcept>
#include <cmath>

template <typename MatrixValueType, size_t m, size_t n>
class Matrix;

template <class T, size_t n>
class Vector
{
    template <typename MatrixValueType, size_t matrix_m, size_t m>
    friend class Matrix;

public:
    Vector();
    Vector(const std::array<T, n> &values);
    Vector(const Vector<T, n> &other);

    Vector<T, n> &operator=(const Vector<T, n> &other);

    size_t size() const;

    bool operator==(const Vector<T, n> &other) const;
    bool operator!=(const Vector<T, n> &other) const;

    const T &operator[](size_t i) const;
    T &operator[](size_t i);

    Vector<T, n> operator+(const Vector<T, n> &other) const;
    Vector<T, n> &operator+=(const Vector<T, n> &other);
    Vector<T, n> &operator+=(const T &value);

    Vector<T, n> operator-(const Vector<T, n> &other) const;
    Vector<T, n> &operator-=(const Vector<T, n> &other);
    Vector<T, n> &operator-=(const T &value);

    Vector<T, n> &operator*=(const T &value);
    Vector<T, n> operator*(const Vector<T, n> &other) const;
    Vector<T, n> &operator*=(const Vector<T, n> &other);

    T dot(const Vector<T, n> &other) const;

    template <size_t m>
    Vector<T, m> dot(const Matrix<T, n, m> &matrix) const;

    template <ssize_t l = 0, ssize_t r = n, ssize_t step = 1, ssize_t slice_size = (r - l) / step>
    Vector<T, slice_size> slice();

private:
    std::array<T, n> arr_;

    static constexpr double kRealNumberPrecision = 1e-6;
};

template <class T, size_t n>
Vector<T, n>::Vector() {}

template <class T, size_t n>
Vector<T, n>::Vector(const std::array<T, n> &values)
{
    std::copy(values.begin(), values.end(), arr_.begin());
}

template <class T, size_t n>
Vector<T, n>::Vector(const Vector<T, n> &other)
{
    std::copy(other.arr_.begin(), other.arr_.end(), arr_.begin());
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator=(const Vector<T, n> &other)
{
    if (this == &other)
    {
        return *this;
    }

    std::copy(other.arr_.begin(), other.arr_.end(), arr_.begin());
    return *this;
}

template <class T, size_t n>
size_t Vector<T, n>::size() const
{
    return n;
}

template <class T, size_t n>
bool Vector<T, n>::operator==(const Vector<T, n> &other) const
{
    for (size_t i = 0; i < n; ++i)
    {
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
        {
            if (std::abs(arr_[i] - other.arr_[i]) >= kRealNumberPrecision)
            {
                return false;
            }
        }
        else
        {
            if (arr_[i] != other.arr_[i])
            {
                return false;
            }
        }
    }
    return true;
}

template <class T, size_t n>
bool Vector<T, n>::operator!=(const Vector<T, n> &other) const
{
    return !(*this == other);
}

template <class T, size_t n>
const T &Vector<T, n>::operator[](size_t i) const
{
    if (i >= n)
    {
        throw std::runtime_error("out of the range");
    }

    return arr_[i];
}

template <class T, size_t n>
T &Vector<T, n>::operator[](size_t i)
{
    if (i >= n)
    {
        throw std::runtime_error("out of the range");
    }

    return arr_[i];
}

template <class T, size_t n>
Vector<T, n> operator+(const Vector<T, n> &vector, const T &value)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = vector[i] + value;
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> operator+(const T &value, const Vector<T, n> &vector)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = value + vector[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> Vector<T, n>::operator+(const Vector<T, n> &other) const
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res.arr_[i] = arr_[i] + other.arr_[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator+=(const Vector<T, n> &other)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] += other.arr_[i];
    }
    return *this;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator+=(const T &value)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] += value;
    }
    return *this;
}

template <class T, size_t n>
Vector<T, n> operator-(const Vector<T, n> &vector, const T &value)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = vector[i] - value;
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> operator-(const T &value, const Vector<T, n> &vector)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = value - vector[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> Vector<T, n>::operator-(const Vector<T, n> &other) const
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res.arr_[i] = arr_[i] - other.arr_[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator-=(const Vector<T, n> &other)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] -= other.arr_[i];
    }
    return *this;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator-=(const T &value)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] -= value;
    }
    return *this;
}

template <class T, size_t n>
Vector<T, n> operator*(const Vector<T, n> &vector, const T &value)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = vector[i] * value;
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> operator*(const T &value, const Vector<T, n> &vector)
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res[i] = value * vector[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator*=(const T &value)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] *= value;
    }
    return *this;
}

template <class T, size_t n>
Vector<T, n> Vector<T, n>::operator*(const Vector<T, n> &other) const
{
    Vector<T, n> res;
    for (size_t i = 0; i < n; ++i)
    {
        res.arr_[i] = arr_[i] * other.arr_[i];
    }
    return res;
}

template <class T, size_t n>
Vector<T, n> &Vector<T, n>::operator*=(const Vector<T, n> &other)
{
    for (size_t i = 0; i < n; ++i)
    {
        arr_[i] *= other.arr_[i];
    }
    return *this;
}

template <class T, size_t n>
T Vector<T, n>::dot(const Vector<T, n> &other) const
{
    T res = 0;
    for (size_t i = 0; i < n; ++i)
    {
        res += arr_[i] * other.arr_[i];
    }
    return res;
}

template <class T, size_t n>
template <size_t m>
Vector<T, m> Vector<T, n>::dot(const Matrix<T, n, m> &matrix) const
{
    Vector<T, m> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = 0;
        for (size_t j = 0; j < n; ++j)
        {
            res[i] += arr_[j] * matrix[j][i];
        }
    }
    return res;
}

template <class T, size_t n>
template <ssize_t l, ssize_t r, ssize_t step, ssize_t slice_size>
Vector<T, slice_size> Vector<T, n>::slice()
{
    static_assert(l <= r);
    static_assert(r <= n);
    static_assert(0 <= l);

    if (l == r)
    {
        return {};
    }

    std::array<T, slice_size> tmp;
    ssize_t i = 0;
    for (ssize_t j = l; j < r; j += step, ++i)
    {
        tmp[i] = arr_[j];
    }
    return {tmp};
}
