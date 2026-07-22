#pragma once

#ifdef _WIN32
#define HEDGES_API __declspec(dllexport)
#else
#define HEDGES_API
#endif

extern "C" {

HEDGES_API int hedges_encode(
    const unsigned char* data,
    int length,
    char*** sequences,
    int* count
);

HEDGES_API void hedges_free_sequences(
    char** sequences,
    int count
);

HEDGES_API int hedges_decode(
    char** sequences,
    int count,
    unsigned char** output,
    int* length
);

}