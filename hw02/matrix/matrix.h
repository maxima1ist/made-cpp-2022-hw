#pragma once

#include "vector.h"

template <class ValueType, size_t m = 1, size_t n = 1>
class Matrix
{
public:
    Matrix() {}

    Matrix(std::array<ValueType, m * n> values)
    {
        for (size_t i = 0; i < m; ++i)
        {
            std::copy(values.begin() + i * n, values.begin() + n + i * n, matrix_[i].arr_.begin());
        }
    }

    Matrix(std::array<Vector<ValueType, n>, m> vectors)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] = vectors[i];
        }
    }

    Matrix(const Matrix<ValueType, m, n> &other)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] = other[i];
        }
    }

    size_t size() const
    {
        return m;
    }

    bool operator==(const Matrix<ValueType, m, n> &other) const
    {
        for (size_t i = 0; i < m; ++i)
        {
            if (matrix_[i] != other.matrix_[i])
            {
                return false;
            }
        }
        return true;
    }

    const Vector<ValueType, n> &operator[](size_t i) const
    {
        if (i > m)
        {
            throw std::runtime_error("out of the range");
        }

        return matrix_[i];
    }

    Vector<ValueType, n> &operator[](size_t i)
    {
        if (i > m)
        {
            throw std::runtime_error("out of the range");
        }

        return matrix_[i];
    }

    template <size_t diag_size = std::min(m, n)>
    Vector<ValueType, diag_size> get_diag() const
    {
        static_assert(diag_size == std::min(m, n));

        std::array<ValueType, diag_size> diag;
        for (size_t i = 0; i < diag_size; ++i)
        {
            diag[i] = matrix_[i][i];
        }
        return {std::move(diag)};
    }

    Vector<ValueType, n> get_row(size_t i) const
    {
        if (i > m)
        {
            throw std::runtime_error("out of the range");
        }

        return {matrix_[i].arr_};
    }

    Vector<ValueType, m> get_column(size_t i) const
    {
        if (i > n)
        {
            throw std::runtime_error("out of the range");
        }

        std::array<ValueType, m> column;
        for (size_t j = 0; j < m; ++j)
        {
            column[j] = matrix_[j][i];
        }
        return {std::move(column)};
    }

    Matrix<ValueType, m, n> operator+(const ValueType &value) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] + value;
        }
        return res;
    }

    friend Matrix<ValueType, m, n> operator+(const ValueType &value, const Matrix<ValueType, m, n> &matrix)
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix.matrix_[i] + value;
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator+=(const ValueType &value)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] += value;
        }
        return *this;
    }

    Matrix<ValueType, m, n> operator-(const ValueType &value) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] - value;
        }
        return res;
    }

    friend Matrix<ValueType, m, n> operator-(const ValueType &value, const Matrix<ValueType, m, n> &matrix)
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = value - matrix.matrix_[i];
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator-=(const ValueType &value)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] -= value;
        }
        return *this;
    }

    Matrix<ValueType, m, n> operator*(const ValueType &value) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] * value;
        }
        return res;
    }

    friend Matrix<ValueType, m, n> operator*(const ValueType &value, const Matrix<ValueType, m, n> &matrix)
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix.matrix_[i] * value;
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator*=(const ValueType &value)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] *= value;
        }
        return *this;
    }

    Matrix<ValueType, m, n> add_by_column(const Vector<ValueType, m> &vector)
    {
        Matrix<ValueType, m, n> res = *this;
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                res[i][j] += vector[i];
            }
        }
        return res;
    }

    Matrix<ValueType, m, n> add_by_row(const Vector<ValueType, n> &vector)
    {
        Matrix<ValueType, m, n> res = *this;
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                res[i][j] += vector[j];
            }
        }
        return res;
    }

    Matrix<ValueType, m, n> sub_by_column(const Vector<ValueType, m> &vector)
    {
        Matrix<ValueType, m, n> res = *this;
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                res[i][j] -= vector[i];
            }
        }
        return res;
    }

    Matrix<ValueType, m, n> sub_by_row(const Vector<ValueType, n> &vector)
    {
        Matrix<ValueType, m, n> res = *this;
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                res[i][j] -= vector[j];
            }
        }
        return res;
    }

    Matrix<ValueType, m, n> operator+(const Matrix<ValueType, m, n> &other) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] + other.matrix_[i];
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator+=(const Matrix<ValueType, m, n> &other)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] += other.matrix_[i];
        }
        return *this;
    }

    Matrix<ValueType, m, n> operator-(const Matrix<ValueType, m, n> &other) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] - other.matrix_[i];
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator-=(const Matrix<ValueType, m, n> &other)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] -= other.matrix_[i];
        }
        return *this;
    }

    Matrix<ValueType, m, n> operator*(const Matrix<ValueType, m, n> &other) const
    {
        Matrix<ValueType, m, n> res;
        for (size_t i = 0; i < m; ++i)
        {
            res.matrix_[i] = matrix_[i] * other.matrix_[i];
        }
        return res;
    }

    Matrix<ValueType, m, n> &operator*=(const Matrix<ValueType, m, n> &other)
    {
        for (size_t i = 0; i < m; ++i)
        {
            matrix_[i] *= other.matrix_[i];
        }
        return *this;
    }

    Vector<ValueType, m> dot(const Vector<ValueType, n> &vector) const
    {
        Vector<ValueType, m> res;
        for (size_t i = 0; i < m; ++i)
        {
            res[i] = 0;
            for (size_t j = 0; j < n; ++j)
            {
                res[i] += matrix_[i][j] * vector[j];
            }
        }
        return res;
    }

    template <size_t k>
    Matrix<ValueType, m, k> dot(const Matrix<ValueType, n, k> &other) const
    {
        Matrix<ValueType, m, k> res;
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < k; ++j)
            {
                res[i][j] = 0;
                for (size_t l = 0; l < n; ++l)
                {
                    res[i][j] += matrix_[i][l] * other[l][j];
                }
            }
        }
        return res;
    }

    Matrix<ValueType, n, m> get_transpose() const
    {
        Matrix<ValueType, n, m> res;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                res[i][j] = matrix_[j][i];
            }
        }
        return res;
    }

    template <size_t min_size = std::min(m, n)>
    ValueType get_determinant() const
    {
        return find_determinant<min_size>(*this, min_size);
    }

    template <size_t min_size = std::min(m, n)>
    Matrix<ValueType, min_size, min_size> get_inversed() const
    {
        return find_inversed<min_size>(*this);
    }

private:
    template <size_t size>
    static ValueType find_determinant(const Matrix<ValueType, size, size> &matrix, size_t dimension)
    {
        ValueType res = 0;

        if (dimension == 1)
        {
            return matrix[0][0];
        }

        if (dimension == 2)
        {
            return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
        }

        Matrix<ValueType, size, size> tmp;
        for (size_t x = 0; x < dimension; ++x)
        {
            size_t subi = 0;
            for (size_t i = 1; i < dimension; ++i)
            {
                size_t subj = 0;
                for (size_t j = 0; j < dimension; ++j)
                {
                    if (j == x)
                    {
                        continue;
                    }

                    tmp[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            res = res + (pow(-1, x) * matrix[0][x] * find_determinant(tmp, dimension - 1));
        }

        return res;
    }

    template <size_t size>
    static void get_cofactor(const Matrix<ValueType, size, size> &matrix, Matrix<ValueType, size, size> &tmp,
                             size_t p, size_t q, size_t tmp_size)
    {
        size_t i = 0, j = 0;
        for (size_t row = 0; row < tmp_size; ++row)
        {
            for (size_t col = 0; col < tmp_size; ++col)
            {
                if (row != p && col != q)
                {
                    tmp[i][j] = matrix[row][col];
                    if (++j == tmp_size - 1)
                    {
                        j = 0;
                        ++i;
                    }
                }
            }
        }
    }

    template <size_t size>
    static void adjoint(const Matrix<ValueType, size, size> &matrix, Matrix<ValueType, size, size> &adj)
    {
        if (size == 1)
        {
            adj[0][0] = 1;
            return;
        }

        ssize_t sign = 1;
        Matrix<ValueType, size, size> tmp;
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < size; ++j)
            {
                get_cofactor(matrix, tmp, i, j, size);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                adj[j][i] = (sign) * (find_determinant<size>(tmp, size - 1));
            }
        }
    }

    template <size_t size>
    static Matrix<ValueType, size, size> find_inversed(const Matrix<ValueType, size, size> &matrix)
    {
        ValueType det = find_determinant<size>(matrix, size);
        if (det == 0)
        {
            throw std::runtime_error("nonsingular matrix");
        }

        Matrix<ValueType, size, size> adj;
        adjoint<size>(matrix, adj);

        Matrix<ValueType, size, size> res;
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < size; j++)
            {
                res[i][j] = adj[i][j] / det;
            }
        }

        return res;
    }

    Vector<Vector<ValueType, n>, m> matrix_;
};