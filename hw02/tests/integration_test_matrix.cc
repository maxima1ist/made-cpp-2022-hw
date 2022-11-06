#include "matrix.h"

#include <gtest/gtest.h>

TEST(TestMatrixVector, FullTest)
{
    Matrix<double, 2, 2> one({1, 2, 3, 4});
    Matrix<double, 2, 2> two({Vector<double, 2>({1, 2}),
                              Vector<double, 2>({3, 4})});
    ASSERT_NEAR(one[0][0], two[0][0], REAL_NUMBER_PRECISION);
    ASSERT_NEAR(one[0][1], two[0][1], REAL_NUMBER_PRECISION);
    ASSERT_NEAR(one[1][0], two[1][0], REAL_NUMBER_PRECISION);
    ASSERT_NEAR(one[1][1], two[1][1], REAL_NUMBER_PRECISION);

    ASSERT_TRUE((bool)(one.get_diag().slice<0>() == two.get_diag()));
    ASSERT_TRUE((bool)(one.get_row(0).slice<0, 1>() == two.get_row(0).slice<0, 1>()));
    ASSERT_TRUE((bool)(one.get_column(1).slice<0, 2, 2>() == two.get_column(1).slice<0, 2, 2>()));

    ASSERT_TRUE((bool)(3 * one == one + two + one));
    ASSERT_TRUE((bool)(3 * (one + two) == 6 * two));

    Matrix<double, 2, 2> zero({0, 0, 0, 0});
    ASSERT_TRUE((bool)(zero == one - two));
    ASSERT_TRUE((bool)(zero == one + two - 2 * one));

    Matrix<double, 2, 2> mult({1 + 10, 4 + 10, 9 + 10, 16 + 10});
    ASSERT_TRUE((bool)(mult == one * two + 10));

    Matrix<double, 2, 2> dot_product({2 * 7, 2 * 10, 2 * 15, 2 * 22});
    ASSERT_TRUE((bool)(dot_product == 2 * one.dot(two)));

    Matrix<double, 2, 2> transposed({1, 3, 2, 4});
    ASSERT_TRUE((bool)(one.get_transpose() == transposed));
    ASSERT_NEAR(one.get_transpose().get_determinant(), -2, REAL_NUMBER_PRECISION);

    Matrix<double, 2, 2> inversed({-2, 1, 1.5, -0.5});
    ASSERT_TRUE((bool)(one.get_inversed() == inversed));
    ASSERT_NEAR(one.get_inversed().get_determinant(), 1 / (double)-2, REAL_NUMBER_PRECISION);
}