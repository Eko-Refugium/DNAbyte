#include "hedges_wrapper.h"
#include "DNAcode.h"

#include <vector>
#include <string>
#include <cstring>


int hedges_encode(
    const unsigned char* data,
    int length,
    char*** sequences,
    int* count
)
{
    /*
       Here we call:
       
       messtodna()
       
       from HEDGES
    */

    return 0;
}


void hedges_free_sequences(
    char** sequences,
    int count
)
{
    for(int i=0;i<count;i++)
    {
        delete[] sequences[i];
    }

    delete[] sequences;
}


int hedges_decode(
    char** sequences,
    int count,
    unsigned char** output,
    int* length
)
{
    /*
       Here we call:
       
       dnatomess()
       
       and extractplaintext()
    */

    return 0;
}