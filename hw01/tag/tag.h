#ifndef TAG_H
#define TAG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ATTRIBURE_COUNT 100
#define MAX_NAME_LEN 20
#define MAX_VALUE_LEN 100

typedef enum TagType
{
    OPEN,
    CLOSE
} TagType;

typedef struct Attribute
{
    char name[MAX_NAME_LEN];
    char value[MAX_VALUE_LEN];
} Attribute;

typedef struct Tag
{
    char name[MAX_NAME_LEN];
    TagType tage_type;
    Attribute **attributes;
} Tag;

Tag *parse_tag(const char *str);

void print_tag(Tag *tag);

void destroy_tag(Tag *tag);

#endif // TAG_H