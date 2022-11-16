#include "vector.hpp"

#include <gtest/gtest.h>

TEST(VectorTest, Size)
{
    Vector<double, 5> vector;
    ASSERT_EQ(vector.size(), 5);
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

TEST(VectorTest, CopyConstructor)
{
    Vector<double, 3> vector({0, 1, 2});
    Vector<double, 3> other(vector);
    ASSERT_TRUE((bool)(vector == other));
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

TEST(VectorTest, DotVector)
{
    ASSERT_TRUE((bool)(Vector<double, 3>({1, 2, 3}).dot(Vector<double, 3>({3, 2, 1})) == 3 + 4 + 3));
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