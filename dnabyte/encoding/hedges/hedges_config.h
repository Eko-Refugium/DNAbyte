#pragma once

#include <string>

struct HedgesConfig
{
    int coderate = 3;

    int total_strand_length = 300;

    int strand_id_bytes = 2;
    int strand_runout_bytes = 2;

    int gc_window = 12;
    int gc_max = 8;
    int max_homopolymer = 4;

    int heap_limit = 1000000;

    std::string left_primer =
        "TCGAAGTCAGCGTGTATTGTATG";

    std::string right_primer =
        "TAGTGAGTGCGATTAAGCGTGTT";
};