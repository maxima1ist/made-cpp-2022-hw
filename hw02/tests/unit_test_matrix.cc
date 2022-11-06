#include "matrix.h"

#include <gtest/gtest.h>

TEST(VectorTest, Size)
{
    Vector<double, 5> vector;
    ASSERT_EQ(vector.size(), 5);
}

TEST(MatrixTest, TestSize)
{
    Matrix<double, 5, 10> matrix;
    ASSERT_EQ(matrix.size(), 5);
    for (size_t i; i < matrix.size(); ++i)
    {
        ASSERT_EQ(matrix[i].size(), 10);
    }
}

TEST(VectorTest, ByIndex)
{
    Vector<double, 3> vector({0, 1, 2});
    for (size_t i = 0; i < vector.size(); ++i)
    {
        ASSERT_EQ(vector[i], i);

        auto prev_value = vector[i];
        vector[i] = i + 1;
        ASSERT_EQ(vector[i], i + 1);
        ASSERT_EQ(prev_value, i);
    }
}

TEST(MatrixTest, ByIndex)
{
    Matrix<double, 2, 2> matrix({Vector<double, 2>({0, 1}),
                                 Vector<double, 2>({0, 1})});
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {
            ASSERT_EQ(matrix[i][j], j);

            auto prev_value = matrix[i][j];
            matrix[i][j] = j + 1;
            ASSERT_EQ(matrix[i][j], j + 1);
            ASSERT_EQ(prev_value, j);
        }
    }
}

TEST(MatrixTest, Diag)
{
    Matrix<double, 1, 1> digit({1});
    ASSERT_TRUE((bool)(Vector<double, 1>({1}) == digit.get_diag()));

    Matrix<double, 2, 2> square({1, 2, 3, 4});
    ASSERT_TRUE((bool)(Vector<double, 2>({1, 4}) == square.get_diag()));

    Matrix<double, 3, 2> vertical({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 2>({1, 4}) == vertical.get_diag()));

    Matrix<double, 2, 3> horizontal({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 2>({1, 5}) == horizontal.get_diag()));
}

TEST(MatrixTest, Row)
{
    Matrix<double, 1, 1> digit({1});
    ASSERT_TRUE((bool)(Vector<double, 1>({1}) == digit.get_row(0)));

    Matrix<double, 2, 2> square({1, 2, 3, 4});
    ASSERT_TRUE((bool)(Vector<double, 2>({3, 4}) == square.get_row(1)));

    Matrix<double, 3, 2> vertical({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 2>({3, 4}) == vertical.get_row(1)));

    Matrix<double, 2, 3> horizontal({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 3>({4, 5, 6}) == horizontal.get_row(1)));
}

TEST(MatrixTest, Column)
{
    Matrix<double, 1, 1> digit({1});
    ASSERT_TRUE((bool)(Vector<double, 1>({1}) == digit.get_column(0)));

    Matrix<double, 2, 2> square({1, 2, 3, 4});
    ASSERT_TRUE((bool)(Vector<double, 2>({2, 4}) == square.get_column(1)));

    Matrix<double, 3, 2> vertical({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 3>({2, 4, 6}) == vertical.get_column(1)));

    Matrix<double, 2, 3> horizontal({1, 2, 3, 4, 5, 6});
    ASSERT_TRUE((bool)(Vector<double, 2>({2, 5}) == horizontal.get_column(1)));
}

TEST(VectorTest, AddSubMultWithNumber)
{
    ASSERT_TRUE((bool)(Vector<double, 3>({11, 12, 13}) == Vector<double, 3>({1, 2, 3}) + 10));
    ASSERT_TRUE((bool)(Vector<double, 3>({11, 12, 13}) == (Vector<double, 3>({1, 2, 3}) += 10)));
    ASSERT_TRUE((bool)(Vector<double, 3>({6, 7, 8}) == 5 + Vector<double, 3>({1, 2, 3})));

    ASSERT_TRUE((bool)(Vector<double, 3>({0, 1, 2}) == Vector<double, 3>({1, 2, 3}) - 1));
    ASSERT_TRUE((bool)(Vector<double, 3>({0, 1, 2}) == (Vector<double, 3>({1, 2, 3}) -= 1)));
    ASSERT_TRUE((bool)(Vector<double, 3>({1, 0, -1}) == 2 - Vector<double, 3>({1, 2, 3})));

    ASSERT_TRUE((bool)(Vector<double, 3>({3, 6, 9}) == Vector<double, 3>({1, 2, 3}) * 3));
    ASSERT_TRUE((bool)(Vector<double, 3>({3, 6, 9}) == (Vector<double, 3>({1, 2, 3}) *= 3)));
    ASSERT_TRUE((bool)(Vector<double, 3>({2, 4, 6}) == 2 * Vector<double, 3>({1, 2, 3})));
}

TEST(VectorTest, AddSubMultWithVector)
{
    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}) + Vector<double, 3>({3, 2, 1}) == Vector<double, 3>({4, 4, 4})));
    ASSERT_TRUE((bool)((Vector<double, 3>({1, 2, 3}) += Vector<double, 3>({3, 2, 1})) == Vector<double, 3>({4, 4, 4})));

    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}) - Vector<double, 3>({3, 2, 1}) == Vector<double, 3>({-2, 0, 2})));
    ASSERT_TRUE((bool)((Vector<double, 3>({1, 2, 3}) -= Vector<double, 3>({3, 2, 1})) == Vector<double, 3>({-2, 0, 2})));

    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}) * Vector<double, 3>({3, 2, 1}) == Vector<double, 3>({3, 4, 3})));
    ASSERT_TRUE((bool)((Vector<double, 3>({1, 2, 3}) *= Vector<double, 3>({3, 2, 1})) == Vector<double, 3>({3, 4, 3})));
}

TEST(MatrixTest, AddSubMultWithNumber)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({3, 4, 5, 6}) == Matrix<double, 2, 2>({2, 3, 4, 5}) + 1));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({3, 4, 5, 6}) == (Matrix<double, 2, 2>({2, 3, 4, 5}) += 1)));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({3, 4, 5, 6}) == 1 + Matrix<double, 2, 2>({2, 3, 4, 5})));

    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}) == Matrix<double, 2, 2>({2, 3, 4, 5}) - 1));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}) == (Matrix<double, 2, 2>({2, 3, 4, 5}) -= 1)));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({-1, -2, -3, -4}) == 1 - Matrix<double, 2, 2>({2, 3, 4, 5})));

    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({4, 6, 8, 10}) == Matrix<double, 2, 2>({2, 3, 4, 5}) * 2));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({4, 6, 8, 10}) == (Matrix<double, 2, 2>({2, 3, 4, 5}) *= 2)));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({4, 6, 8, 10}) == 2 * Matrix<double, 2, 2>({2, 3, 4, 5})));
}

TEST(MatrixTest, AddSubVector)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).add_by_column(Vector<double, 2>({1, 2})) == Matrix<double, 2, 3>({2, 3, 4, 6, 7, 8})));
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).add_by_row(Vector<double, 3>({1, 2, 3})) == Matrix<double, 2, 3>({2, 4, 6, 5, 7, 9})));

    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).sub_by_column(Vector<double, 2>({1, 2})) == Matrix<double, 2, 3>({0, 1, 2, 2, 3, 4})));
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).sub_by_row(Vector<double, 3>({1, 2, 3})) == Matrix<double, 2, 3>({0, 0, 0, 3, 3, 3})));
}

TEST(MatrixTest, AddSubMultWithMatrix)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}) + Matrix<double, 2, 2>({4, 3, 2, 1}) == Matrix<double, 2, 2>({5, 5, 5, 5})));
    ASSERT_TRUE((bool)((Matrix<double, 2, 2>({1, 2, 3, 4}) += Matrix<double, 2, 2>({4, 3, 2, 1})) == Matrix<double, 2, 2>({5, 5, 5, 5})));

    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}) - Matrix<double, 2, 2>({4, 3, 2, 1}) == Matrix<double, 2, 2>({-3, -1, 1, 3})));
    ASSERT_TRUE((bool)((Matrix<double, 2, 2>({1, 2, 3, 4}) -= Matrix<double, 2, 2>({4, 3, 2, 1})) == Matrix<double, 2, 2>({-3, -1, 1, 3})));

    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}) * Matrix<double, 2, 2>({4, 3, 2, 1}) == Matrix<double, 2, 2>({4, 6, 6, 4})));
    ASSERT_TRUE((bool)((Matrix<double, 2, 2>({1, 2, 3, 4}) *= Matrix<double, 2, 2>({4, 3, 2, 1})) == Matrix<double, 2, 2>({4, 6, 6, 4})));
}

TEST(VectorTest, DotVector)
{
    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}).dot(Vector<double, 3>({3, 2, 1})) == 3 + 4 + 3));
}

TEST(MatrixTest, DotVector)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).dot(Vector<double, 3>({1, 2, 3})) == Vector<double, 2>({1 + 4 + 9, 4 + 10 + 18})));
    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}).dot(Matrix<double, 3, 2>({1, 2, 3, 4, 5, 6})) == Vector<double, 2>({1 + 6 + 15, 2 + 8 + 18})));
}

TEST(MatrixTest, DotMatrix)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).dot(Matrix<double, 3, 2>({1, 2, 3, 4, 5, 6})) == Matrix<double, 2, 2>({1 + 6 + 15, 2 + 8 + 18, 4 + 15 + 30, 8 + 20 + 36})));
}

TEST(MatrixTest, Transpose)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}).get_transpose() == Matrix<double, 2, 2>({1, 3, 2, 4})));
    ASSERT_TRUE((bool)(Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6}).get_transpose() == Matrix<double, 3, 2>({1, 4, 2, 5, 3, 6})));
    ASSERT_TRUE((bool)(Matrix<double, 3, 2>({1, 4, 2, 5, 3, 6}).get_transpose() == Matrix<double, 2, 3>({1, 2, 3, 4, 5, 6})));
}

TEST(MatrixTest, Determinant)
{
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}).get_determinant() == -2));
    ASSERT_TRUE((bool)(Matrix<double, 3, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9}).get_determinant() == 0));
    ASSERT_TRUE((bool)(Matrix<double, 3, 3>({7, 1, 3, 2, 4, 1, 1, 5, 1}).get_determinant() == 10));
}

TEST(MatrixTest, Inversed)
{
    ASSERT_TRUE((bool)(Matrix<double, 1, 1>({100}).get_inversed() == Matrix<double, 1>({0.01})));
    ASSERT_TRUE((bool)(Matrix<double, 2, 2>({1, 2, 3, 4}).get_inversed() == Matrix<double, 2, 2>({-2, 1, 1.5, -0.5})));
    ASSERT_TRUE((bool)(Matrix<double, 3, 3>({1, 1, 4, 1, 2, 4, 1, 2, 2}).get_inversed() == Matrix<double, 3, 3>({2, -3, 2, -1, 1, 0, 0, 0.5, -0.5})));
}

TEST(VectorTest, Slice)
{
    ASSERT_TRUE((bool)(Vector<double, 4>({0, 1, 2, 3}).slice() == Vector<double, 4>({0, 1, 2, 3})));
    ASSERT_TRUE((bool)(Vector<double, 4>({0, 1, 2, 3}).slice<1, 1>() == Vector<double, 0>()));
    ASSERT_TRUE((bool)(Vector<double, 4>({0, 1, 2, 3}).slice<1, 3>() == Vector<double, 2>({1, 2})));
    ASSERT_TRUE((bool)(Vector<double, 4>({0, 1, 2, 3}).slice<0, 4, 2>() == Vector<double, 2>({0, 2})));
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}