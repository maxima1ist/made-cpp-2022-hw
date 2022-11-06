#include "tag.h"

Tag *parse_tag(const char *str)
{
    if (!strlen(str))
    {
        printf("Cann't parse empty string.\n");
        return NULL;
    }

    if (str[0] != '<' || str[strlen(str) - 1] != '>')
    {
        printf("Invalid tag structure. It must start with '<' and end with '>'.\n");
        return NULL;
    }

    char tag_body[strlen(str) - 1];
    memcpy(tag_body, str + 1, strlen(str) - 2);
    tag_body[strlen(str) - 2] = '\0';

    Tag *tag = malloc(sizeof(Tag));

    char *body_token = strtok(tag_body, " ");
    memset(&tag->name[0], 0, sizeof(tag->name));
    if (body_token)
    {
        if (body_token[0] == '/')
        {
            tag->tage_type = CLOSE;
            memcpy(tag->name, body_token + 1, strlen(body_token) - 1);
        }
        else
        {
            tag->tage_type = OPEN;
            memcpy(tag->name, body_token, strlen(body_token));
        }
        body_token = strtok(NULL, " ");
    }

    char **attribute_strs = malloc(MAX_ATTRIBURE_COUNT * sizeof(void *));
    size_t attribute_num = 0;
    while (body_token)
    {
        if (!strlen(body_token))
        {
            continue;
        }

        attribute_strs[attribute_num] = malloc(strlen(body_token));
        strcpy(attribute_strs[attribute_num], body_token);

        body_token = strtok(NULL, " ");
        ++attribute_num;
    }

    if (attribute_num)
    {
        tag->attributes = malloc((attribute_num + 1) * sizeof(Attribute *));
        tag->attributes[attribute_num] = NULL;
    }

    size_t valid_attribute_num = 0;
    for (size_t i = 0; i < attribute_num; ++i)
    {
        char attribute_strs_copy[strlen(attribute_strs[i])];
        strcpy(attribute_strs_copy, attribute_strs[i]);

        char *attribute_name = strtok(attribute_strs_copy, "=");
        if (!attribute_name || strlen(attribute_name) == 0)
        {
            printf("Invalid attribute name for '%s'. It will be skipped...\n", attribute_strs[i]);
            continue;
        }

        char *attribut_value = strtok(NULL, "=");
        if (!attribut_value || strlen(attribut_value) == 0)
        {
            printf("Invalid attribute value for '%s'. It will be skipped...\n", attribute_strs[i]);
            continue;
        }

        tag->attributes[valid_attribute_num] = malloc(sizeof(Attribute));

        memset(&tag->attributes[valid_attribute_num]->name[0], 0, sizeof(tag->attributes[valid_attribute_num]->name));
        memcpy(tag->attributes[valid_attribute_num]->name, attribute_name, strlen(attribute_name));

        memset(&tag->attributes[valid_attribute_num]->value[0], 0, sizeof(tag->attributes[valid_attribute_num]->value));
        memcpy(tag->attributes[valid_attribute_num]->value, attribut_value, strlen(attribut_value));

        ++valid_attribute_num;
    }

    if (valid_attribute_num)
    {
        tag->attributes = realloc(tag->attributes, valid_attribute_num * sizeof(tag->attributes[0]));
    }
    else
    {
        tag->attributes = NULL;
    }

    for (size_t i = 0; i < attribute_num; ++i)
    {
        free(attribute_strs[i]);
    }
    free(attribute_strs);

    return tag;
}

void print_tag(Tag *tag)
{
    if (!tag)
    {
        return;
    }

    printf("Tag name: '%s', tag type: '%s'.",
           tag->name,
           tag->tage_type == OPEN ? "open" : "close");

    if (tag->attributes)
    {
        printf(" List of tag attributes and values: ");
        for (size_t i = 0; tag->attributes[i]; ++i)
        {
            printf("%s : %s", tag->attributes[i]->name, tag->attributes[i]->value);
            if (tag->attributes[i + 1])
            {
                printf(", ");
            }
        }
        printf(".\n");
    }
}

void destroy_tag(Tag *tag)
{
    if (!tag)
    {
        return;
    }

    if (tag->attributes)
    {
        for (size_t i = 0; tag->attributes[i]; ++i)
        {
            free(tag->attributes[i]);
        }
        free(tag->attributes);
    }

    free(tag);
}