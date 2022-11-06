#pragma once

#include <cstddef>
#include <array>
#include <stdexcept>
#include <cmath>

#define REAL_NUMBER_PRECISION 1e-6

template <typename MatrixValueType, size_t m, size_t n>
class Matrix;

template <class ValueType, size_t n = 1>
class Vector
{
    template <typename MatrixValueType, size_t matrix_m, size_t matrix_n>
    friend class Matrix;

public:
    Vector() {}

    Vector(std::array<ValueType, n> values)
    {
        std::copy(values.begin(), values.end(), arr_.begin());
    }

    Vector(const Vector<ValueType, n> &other)
    {
        std::copy(other.arr_.begin(), other.arr_.end(), arr_.begin());
    }

    Vector<ValueType, n> &operator=(const Vector<ValueType, n> &other)
    {
        if (this == &other)
        {
            return *this;
        }

        std::copy(other.arr_.begin(), other.arr_.end(), arr_.begin());
        return *this;
    }

    size_t size() const
    {
        return n;
    }

    bool operator==(const Vector<ValueType, n> &other) const
    {
        for (size_t i = 0; i < n; ++i)
        {
            if constexpr (std::is_same_v<ValueType, float> || std::is_same_v<ValueType, double>)
            {
                if (std::abs(arr_[i] - other.arr_[i]) >= REAL_NUMBER_PRECISION)
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

    bool operator!=(const Vector<ValueType, n> &other) const
    {
        return !(*this == other);
    }

    const ValueType &operator[](size_t i) const
    {
        if (i > n)
        {
            throw std::runtime_error("out of the range");
        }

        return arr_[i];
    }

    ValueType &operator[](size_t i)
    {
        if (i > n)
        {
            throw std::runtime_error("out of the range");
        }

        return arr_[i];
    }

    Vector<ValueType, n> operator+(const ValueType &value) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] + value;
        }
        return res;
    }

    friend Vector<ValueType, n> operator+(const ValueType &value, const Vector<ValueType, n> &vector)
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value + vector.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator+=(const ValueType &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] += value;
        }
        return *this;
    }

    Vector<ValueType, n> operator-(const ValueType &value) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] - value;
        }
        return res;
    }

    friend Vector<ValueType, n> operator-(const ValueType &value, const Vector<ValueType, n> &vector)
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value - vector.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator-=(const ValueType &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] -= value;
        }
        return *this;
    }

    Vector<ValueType, n> operator*(const ValueType &value) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] * value;
        }
        return res;
    }

    friend Vector<ValueType, n> operator*(const ValueType &value, const Vector<ValueType, n> &vector)
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value * vector.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator*=(const ValueType &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] *= value;
        }
        return *this;
    }

    Vector<ValueType, n> operator+(const Vector<ValueType, n> &other) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] + other.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator+=(const Vector<ValueType, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] += other.arr_[i];
        }
        return *this;
    }

    Vector<ValueType, n> operator-(const Vector<ValueType, n> &other) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] - other.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator-=(const Vector<ValueType, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] -= other.arr_[i];
        }
        return *this;
    }

    Vector<ValueType, n> operator*(const Vector<ValueType, n> &other) const
    {
        Vector<ValueType, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] * other.arr_[i];
        }
        return res;
    }

    Vector<ValueType, n> &operator*=(const Vector<ValueType, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] *= other.arr_[i];
        }
        return *this;
    }

    ValueType dot(const Vector<ValueType, n> &other) const
    {
        ValueType res = 0;
        for (size_t i = 0; i < n; ++i)
        {
            res += arr_[i] * other.arr_[i];
        }
        return res;
    }

    template <size_t matrix_n>
    Vector<ValueType, matrix_n> dot(const Matrix<ValueType, n, matrix_n> &matrix) const
    {
        Vector<ValueType, matrix_n> res;
        for (size_t i = 0; i < matrix_n; ++i)
        {
            res[i] = 0;
            for (size_t j = 0; j < n; ++j)
            {
                res[i] += arr_[j] * matrix[j][i];
            }
        }
        return res;
    }

    template <ssize_t l = 0, ssize_t r = n, ssize_t step = 1, ssize_t slice_size = (r - l) / step>
    Vector<ValueType, slice_size> slice()
    {
        static_assert(l <= r);
        static_assert(r <= n);
        static_assert(0 <= l);

        if (l == r)
        {
            return {};
        }

        std::array<ValueType, slice_size> tmp;
        ssize_t i = 0;
        for (ssize_t j = l; j < r; j += step, ++i)
        {
            tmp[i] = arr_[j];
        }
        return {tmp};
    }

private:
    std::array<ValueType, n> arr_;
};