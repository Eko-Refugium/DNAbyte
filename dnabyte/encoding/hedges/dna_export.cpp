#include "dna_export.h"


std::vector<std::string> dna_to_strings(
    const MatUchar& dna
)
{
    const char bases[] = {'A','C','G','T'};

    std::vector<std::string> result;


    for(Int i = 0; i < dna.nrows(); i++)
    {
        std::string strand;

        for(Int j = 0; j < dna.ncols(); j++)
        {
            strand += bases[dna[i][j]];
        }

        result.push_back(strand);
    }

    return result;
}
