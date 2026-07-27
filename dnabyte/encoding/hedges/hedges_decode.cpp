#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "HEDGES-main/nr3b.h"
#include "HEDGES-main/RSecc.h"
#include "HEDGES-main/DNAcode.h"
#include "hedges_wrapper.h"
#include "dna_export.h"


static unsigned char dna_value(char c)
{
    switch(c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:
            throw std::runtime_error("Invalid DNA base");
    }
}


std::vector<uint8_t> hedges_decode(
    const std::vector<std::string>& strands
)
{
    MatUchar dna(
        strands.size(),
        strands[0].size()
    );

    for(size_t i = 0; i < strands.size(); i++)
    {
        for(size_t j = 0; j < strands[i].size(); j++)
        {
            dna[i][j] = dna_value(strands[i][j]);
        }
    }

    // Decode DNA -> packets
    Dnatomess_out decoded = dnatomess(dna);

    // Error correction
    Correctmesspacket_out corrected =
        correctmesspacket(
            decoded.mpacket,
            decoded.epacket
        );

    // Extract payload
    VecUchar plaintext =
        extractplaintext(corrected.packet);


    // Remove HEDGES length prefix
    uint16_t len =
        ((uint16_t)plaintext[0] << 8) |
        plaintext[1];


    std::vector<uint8_t> result;

    for(int i = 0; i < len; i++)
    {
        result.push_back(plaintext[i + 2]);
    }

    return result;
}