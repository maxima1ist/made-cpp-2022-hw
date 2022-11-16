#include "matrix.hpp"

#include <gtest/gtest.h>

TEST(MatrixTest, TestSize)
{
    Matrix<double, 5, 10> matrix;
    ASSERT_EQ(matrix.size().first, 5);
    ASSERT_EQ(matrix.size().second, 10);
    for (size_t i; i < matrix.size().first; ++i)
    {
        ASSERT_EQ(matrix[i].size(), 10);
    }
}

TEST(MatrixTest, ByIndex)
{
    Matrix<double, 2, 2> matrix({Vector<double, 2>({0, 1}),
                                 Vector<double, 2>({0, 1})});
    for (size_t i = 0; i < matrix.size().first; ++i)
    {
        for (size_t j = 0; j < matrix.size().second; ++j)
        {
            ASSERT_EQ(matrix[i][j], j);

            auto prev_value = matrix[i][j];
            matrix[i][j] = j + 1;
            ASSERT_EQ(matrix[i][j], j + 1);
            ASSERT_EQ(prev_value, j);
        }
    }
}

TEST(MatrixTest, CopyConstructor)
{
    Matrix<double, 2, 2> matrix({Vector<double, 2>({0, 1}),
                                 Vector<double, 2>({0, 1})});
    Matrix<double, 2, 2> other(matrix);
    ASSERT_TRUE((bool)(matrix == other));
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

    Matrix<double, 1, 1>({100}).get_diag<1>();
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}