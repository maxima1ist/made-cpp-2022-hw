#pragma once

#include "vector.hpp"

template <class T, size_t m, size_t n>
class Matrix
{
public:
    Matrix();
    Matrix(const std::array<T, m * n> &values);
    Matrix(const std::array<Vector<T, n>, m> &vectors);
    Matrix(const Matrix<T, m, n> &other);

    std::pair<size_t, size_t> size() const;

    bool operator==(const Matrix<T, m, n> &other) const;

    const Vector<T, n> &operator[](size_t i) const;
    Vector<T, n> &operator[](size_t i);

    Vector<T, std::min(m, n)> get_diag() const;
    Vector<T, n> get_row(size_t i) const;
    Vector<T, m> get_column(size_t i) const;

    Matrix<T, m, n> operator+(const Matrix<T, m, n> &other) const;
    Matrix<T, m, n> &operator+=(const Matrix<T, m, n> &other);
    Matrix<T, m, n> &operator+=(const T &value);

    Matrix<T, m, n> operator-(const Matrix<T, m, n> &other) const;
    Matrix<T, m, n> &operator-=(const Matrix<T, m, n> &other);
    Matrix<T, m, n> &operator-=(const T &value);

    Matrix<T, m, n> &operator*=(const T &value);
    Matrix<T, m, n> operator*(const Matrix<T, m, n> &other) const;
    Matrix<T, m, n> &operator*=(const Matrix<T, m, n> &other);

    Matrix<T, m, n> add_by_column(const Vector<T, m> &vector);
    Matrix<T, m, n> add_by_row(const Vector<T, n> &vector);
    Matrix<T, m, n> sub_by_column(const Vector<T, m> &vector);
    Matrix<T, m, n> sub_by_row(const Vector<T, n> &vector);

    template <size_t k>
    Matrix<T, m, k> dot(const Matrix<T, n, k> &other) const;
    Vector<T, m> dot(const Vector<T, n> &vector) const;

    Matrix<T, n, m> get_transpose() const;

    T get_determinant() const;

    template <size_t min_size = std::min(m, n)>
    Matrix<T, min_size, min_size> get_inversed() const;

private:
    template <size_t size>
    static T find_determinant(const Matrix<T, size, size> &matrix, size_t dimension)
    {
        T res = 0;

        if (dimension == 1)
        {
            return matrix[0][0];
        }

        if (dimension == 2)
        {
            return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
        }

        Matrix<T, size, size> tmp;
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
    static void get_cofactor(const Matrix<T, size, size> &matrix, Matrix<T, size, size> &tmp,
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
    static void adjoint(const Matrix<T, size, size> &matrix, Matrix<T, size, size> &adj)
    {
        if (size == 1)
        {
            adj[0][0] = 1;
            return;
        }

        ssize_t sign = 1;
        Matrix<T, size, size> tmp;
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
    static Matrix<T, size, size> find_inversed(const Matrix<T, size, size> &matrix)
    {
        T det = find_determinant<size>(matrix, size);
        if (det == 0)
        {
            throw std::runtime_error("nonsingular matrix");
        }

        Matrix<T, size, size> adj;
        adjoint<size>(matrix, adj);

        Matrix<T, size, size> res;
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < size; j++)
            {
                res[i][j] = adj[i][j] / det;
            }
        }

        return res;
    }

    Vector<Vector<T, n>, m> matrix_;
};

template <class T, size_t m, size_t n>
Matrix<T, m, n>::Matrix() {}

template <class T, size_t m, size_t n>
Matrix<T, m, n>::Matrix(const std::array<T, m * n> &values)
{
    for (size_t i = 0; i < m; ++i)
    {
        std::copy(values.begin() + i * n, values.begin() + n + i * n, matrix_[i].arr_.begin());
    }
}

template <class T, size_t m, size_t n>
Matrix<T, m, n>::Matrix(const std::array<Vector<T, n>, m> &vectors)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] = vectors[i];
    }
}

template <class T, size_t m, size_t n>
Matrix<T, m, n>::Matrix(const Matrix<T, m, n> &other)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] = other[i];
    }
}

template <class T, size_t m, size_t n>
std::pair<size_t, size_t> Matrix<T, m, n>::size() const
{
    return std::make_pair(m, n);
}

template <class T, size_t m, size_t n>
bool Matrix<T, m, n>::operator==(const Matrix<T, m, n> &other) const
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

template <class T, size_t m, size_t n>
const Vector<T, n> &Matrix<T, m, n>::operator[](size_t i) const
{
    if (i >= m)
    {
        throw std::runtime_error("out of the range");
    }

    return matrix_[i];
}

template <class T, size_t m, size_t n>
Vector<T, n> &Matrix<T, m, n>::operator[](size_t i)
{
    if (i >= m)
    {
        throw std::runtime_error("out of the range");
    }

    return matrix_[i];
}

template <class T, size_t m, size_t n>
Vector<T, std::min(m, n)> Matrix<T, m, n>::get_diag() const
{
    static constexpr size_t diag_size = std::min(m, n);
    Vector<T, diag_size> diag;
    for (size_t i = 0; i < diag_size; ++i)
    {
        diag[i] = matrix_[i][i];
    }
    return diag;
}

template <class T, size_t m, size_t n>
Vector<T, n> Matrix<T, m, n>::get_row(size_t i) const
{
    if (i >= m)
    {
        throw std::runtime_error("out of the range");
    }

    return matrix_[i];
}

template <class T, size_t m, size_t n>
Vector<T, m> Matrix<T, m, n>::get_column(size_t i) const
{
    if (i >= n)
    {
        throw std::runtime_error("out of the range");
    }

    Vector<T, m> column;
    for (size_t j = 0; j < m; ++j)
    {
        column[j] = matrix_[j][i];
    }
    return column;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator+(const Matrix<T, m, n> &matrix, const T &value)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = matrix[i] + value;
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator+(const T &value, const Matrix<T, m, n> &matrix)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = matrix[i] + value;
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::operator+(const Matrix<T, m, n> &other) const
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res.matrix_[i] = matrix_[i] + other.matrix_[i];
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator+=(const Matrix<T, m, n> &other)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] += other.matrix_[i];
    }
    return *this;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator+=(const T &value)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] += value;
    }
    return *this;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator-(const Matrix<T, m, n> &matrix, const T &value)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = matrix[i] - value;
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator-(const T &value, const Matrix<T, m, n> &matrix)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = value - matrix[i];
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::operator-(const Matrix<T, m, n> &other) const
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res.matrix_[i] = matrix_[i] - other.matrix_[i];
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator-=(const Matrix<T, m, n> &other)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] -= other.matrix_[i];
    }
    return *this;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator-=(const T &value)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] -= value;
    }
    return *this;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator*(const Matrix<T, m, n> &matrix, const T &value)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = matrix[i] * value;
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> operator*(const T &value, const Matrix<T, m, n> &matrix)
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res[i] = matrix[i] * value;
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator*=(const T &value)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] *= value;
    }
    return *this;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::add_by_column(const Vector<T, m> &vector)
{
    Matrix<T, m, n> res = *this;
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            res[i][j] += vector[i];
        }
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::add_by_row(const Vector<T, n> &vector)
{
    Matrix<T, m, n> res = *this;
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            res[i][j] += vector[j];
        }
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::sub_by_column(const Vector<T, m> &vector)
{
    Matrix<T, m, n> res = *this;
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            res[i][j] -= vector[i];
        }
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::sub_by_row(const Vector<T, n> &vector)
{
    Matrix<T, m, n> res = *this;
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            res[i][j] -= vector[j];
        }
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> Matrix<T, m, n>::operator*(const Matrix<T, m, n> &other) const
{
    Matrix<T, m, n> res;
    for (size_t i = 0; i < m; ++i)
    {
        res.matrix_[i] = matrix_[i] * other.matrix_[i];
    }
    return res;
}

template <class T, size_t m, size_t n>
Matrix<T, m, n> &Matrix<T, m, n>::operator*=(const Matrix<T, m, n> &other)
{
    for (size_t i = 0; i < m; ++i)
    {
        matrix_[i] *= other.matrix_[i];
    }
    return *this;
}

template <class T, size_t m, size_t n>
Vector<T, m> Matrix<T, m, n>::dot(const Vector<T, n> &vector) const
{
    Vector<T, m> res;
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

template <class T, size_t m, size_t n>
template <size_t k>
Matrix<T, m, k> Matrix<T, m, n>::dot(const Matrix<T, n, k> &other) const
{
    Matrix<T, m, k> res;
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

template <class T, size_t m, size_t n>
Matrix<T, n, m> Matrix<T, m, n>::get_transpose() const
{
    Matrix<T, n, m> res;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            res[i][j] = matrix_[j][i];
        }
    }
    return res;
}

template <class T, size_t m, size_t n>
T Matrix<T, m, n>::get_determinant() const
{
    static_assert(m == n, "determinant computes only for square matrix");
    return find_determinant<m>(*this, m);
}

template <class T, size_t m, size_t n>
template <size_t min_size>
Matrix<T, min_size, min_size> Matrix<T, m, n>::get_inversed() const
{
    return find_inversed<min_size>(*this);
}
