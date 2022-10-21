extern "C"
{
#include "tag.h"
}

#include <string>
#include <iostream>

#include <gtest/gtest.h>

TEST(TestTag, TestEmpty)
{
    testing::internal::CaptureStdout();

    parse_tag("");

    const std::string real_str = testing::internal::GetCapturedStdout();
    ASSERT_EQ(real_str, "Cann't parse empty string.\n");
}

TEST(TestTag, TestInvalidInput)
{
    testing::internal::CaptureStdout();

    parse_tag("abc");

    const std::string real_str = testing::internal::GetCapturedStdout();
    ASSERT_EQ(real_str, "Invalid tag structure. It must start with '<' and end with '>'.\n");
}

TEST(TestTag, TestOpenTag)
{
    Tag *tag = parse_tag("<br>");

    ASSERT_STREQ(tag->name, "br");
    ASSERT_EQ(tag->tage_type, OPEN);

    destroy_tag(tag);
}

TEST(TestTag, TestCloseTag)
{
    Tag *tag = parse_tag("</br>");

    ASSERT_STREQ(tag->name, "br");
    ASSERT_EQ(tag->tage_type, CLOSE);

    destroy_tag(tag);
}

TEST(TestTag, TestTagWithInvalidAttributeValue)
{
    testing::internal::CaptureStdout();

    Tag *tag = parse_tag("<p atr1=  >");

    const std::string real_str = testing::internal::GetCapturedStdout();
    ASSERT_EQ(real_str, "Invalid attribute value for 'atr1='. It will be skipped...\n");

    ASSERT_STREQ(tag->name, "p");
    ASSERT_EQ(tag->tage_type, OPEN);
    ASSERT_FALSE(tag->attributes);

    destroy_tag(tag);
}

TEST(TestTag, TestTagWithOneAttribute)
{
    Tag *tag = parse_tag("<p atr1=\"value1\">");

    ASSERT_STREQ(tag->name, "p");
    ASSERT_EQ(tag->tage_type, OPEN);

    std::vector<std::pair<std::string, std::string>> atr_to_values{{"atr1", "\"value1\""}};
    size_t i = 0;
    while (tag->attributes[i])
    {
        ASSERT_LE(i, atr_to_values.size());
        ASSERT_EQ(tag->attributes[i]->name, atr_to_values[i].first);
        ASSERT_EQ(tag->attributes[i]->value, atr_to_values[i].second);
        ++i;
    }
    ASSERT_EQ(i, atr_to_values.size());

    destroy_tag(tag);
}

TEST(TestTag, TestTagWithMultipleAttribute)
{
    Tag *tag = parse_tag("<xxx atr1=\"value1\"    atr2=\"value2\"   some_big_attribure="
                         "\"https://coderefinery.github.io/cmake-workshop/testing/\">");

    ASSERT_STREQ(tag->name, "xxx");
    ASSERT_EQ(tag->tage_type, OPEN);

    std::vector<std::pair<std::string, std::string>>
        atr_to_values{{"atr1", "\"value1\""},
                      {"atr2", "\"value2\""},
                      {"some_big_attribure",
                       "\"https://coderefinery.github.io/cmake-workshop/testing/\""}};
    size_t i = 0;
    while (tag->attributes[i])
    {
        ASSERT_LE(i, atr_to_values.size());
        ASSERT_EQ(tag->attributes[i]->name, atr_to_values[i].first);
        ASSERT_EQ(tag->attributes[i]->value, atr_to_values[i].second);
        ++i;
    }
    ASSERT_EQ(i, atr_to_values.size());

    destroy_tag(tag);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}