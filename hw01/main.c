#include "tag.h"

typedef enum ReturnCode
{
    OK = 0,
    INVALID_INPUT = 1
} ReturnCode;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s 'html tag with attributes'\n", argv[0]);
        return INVALID_INPUT;
    }

    Tag *tag = parse_tag(argv[1]);

    print_tag(tag);

    destroy_tag(tag);

    return OK;
}