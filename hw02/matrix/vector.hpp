#pragma once

#include <cstddef>
#include <array>
#include <stdexcept>
#include <cmath>

template <typename MatrixValueType, size_t m, size_t n>
class Matrix;

template <class Т, size_t n = 1>
class Vector
{
    template <typename MatrixValueType, size_t matrix_m, size_t m>
    friend class Matrix;

public:
    Vector() {}

    Vector(const std::array<Т, n> &values)
    {
        std::copy(values.begin(), values.end(), arr_.begin());
    }

    Vector(const Vector<Т, n> &other)
    {
        std::copy(other.arr_.begin(), other.arr_.end(), arr_.begin());
    }

    Vector<Т, n> &operator=(const Vector<Т, n> &other)
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

    bool operator==(const Vector<Т, n> &other) const
    {
        for (size_t i = 0; i < n; ++i)
        {
            if constexpr (std::is_same_v<Т, float> || std::is_same_v<Т, double>)
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

    bool operator!=(const Vector<Т, n> &other) const
    {
        return !(*this == other);
    }

    const Т &operator[](size_t i) const
    {
        if (i >= n)
        {
            throw std::runtime_error("out of the range");
        }

        return arr_[i];
    }

    Т &operator[](size_t i)
    {
        if (i >= n)
        {
            throw std::runtime_error("out of the range");
        }

        return arr_[i];
    }

    Vector<Т, n> operator+(const Т &value) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] + value;
        }
        return res;
    }

    friend Vector<Т, n> operator+(const Т &value, const Vector<Т, n> &vector)
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value + vector.arr_[i];
        }
        return res;
    }

    Vector<Т, n> operator+(const Vector<Т, n> &other) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] + other.arr_[i];
        }
        return res;
    }

    Vector<Т, n> &operator+=(const Vector<Т, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] += other.arr_[i];
        }
        return *this;
    }

    Vector<Т, n> &operator+=(const Т &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] += value;
        }
        return *this;
    }

    Vector<Т, n> operator-(const Т &value) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] - value;
        }
        return res;
    }

    friend Vector<Т, n> operator-(const Т &value, const Vector<Т, n> &vector)
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value - vector.arr_[i];
        }
        return res;
    }

    Vector<Т, n> operator-(const Vector<Т, n> &other) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] - other.arr_[i];
        }
        return res;
    }

    Vector<Т, n> &operator-=(const Vector<Т, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] -= other.arr_[i];
        }
        return *this;
    }

    Vector<Т, n> &operator-=(const Т &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] -= value;
        }
        return *this;
    }

    Vector<Т, n> operator*(const Т &value) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] * value;
        }
        return res;
    }

    friend Vector<Т, n> operator*(const Т &value, const Vector<Т, n> &vector)
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = value * vector.arr_[i];
        }
        return res;
    }

    Vector<Т, n> &operator*=(const Т &value)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] *= value;
        }
        return *this;
    }

    Vector<Т, n> operator*(const Vector<Т, n> &other) const
    {
        Vector<Т, n> res;
        for (size_t i = 0; i < n; ++i)
        {
            res.arr_[i] = arr_[i] * other.arr_[i];
        }
        return res;
    }

    Vector<Т, n> &operator*=(const Vector<Т, n> &other)
    {
        for (size_t i = 0; i < n; ++i)
        {
            arr_[i] *= other.arr_[i];
        }
        return *this;
    }

    Т dot(const Vector<Т, n> &other) const
    {
        Т res = 0;
        for (size_t i = 0; i < n; ++i)
        {
            res += arr_[i] * other.arr_[i];
        }
        return res;
    }

    template <size_t m>
    Vector<Т, m> dot(const Matrix<Т, n, m> &matrix) const
    {
        Vector<Т, m> res;
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

    template <ssize_t l = 0, ssize_t r = n, ssize_t step = 1, ssize_t slice_size = (r - l) / step>
    Vector<Т, slice_size> slice()
    {
        static_assert(l <= r);
        static_assert(r <= n);
        static_assert(0 <= l);

        if (l == r)
        {
            return {};
        }

        std::array<Т, slice_size> tmp;
        ssize_t i = 0;
        for (ssize_t j = l; j < r; j += step, ++i)
        {
            tmp[i] = arr_[j];
        }
        return {tmp};
    }

private:
    std::array<Т, n> arr_;

    static constexpr double kRealNumberPrecision = 1e-6;
};