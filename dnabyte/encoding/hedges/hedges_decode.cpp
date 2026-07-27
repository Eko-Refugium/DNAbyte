#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "HEDGES-main/nr3b.h"
#include "HEDGES-main/RSecc.h"
#include "HEDGES-main/DNAcode.h"
#include "hedges_wrapper.h"
#include "dna_export.h"

extern Int bytesperstrand;
extern Int strandsperpacketmessage;
extern Int messbytesperstrand;
extern Int messbytesperpacket;
extern Int strandIDbytes;
extern Int strandrunoutbytes;
extern Int totstrandlen;
extern Int leftlen;
extern Int rightlen;
extern Int strandlen;
// globals not normally user settable because assumed by Reed - Solomon outer code :
extern Int strandsperpacket;
extern Int strandsperpacketcheck;

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
std::cout << "Before init" << std::endl;

std::cout << "totstrandlen: " << totstrandlen << std::endl;
std::cout << "strandIDbytes: " << strandIDbytes << std::endl;
std::cout << "strandrunoutbytes: " << strandrunoutbytes << std::endl;

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

    std::cout << "Calling dnatomess..." << std::endl;

    Dnatomess_out decoded = dnatomess(dna);

    std::cout << "dnatomess finished" << std::endl;

    std::cout << "mpacket rows: "
            << decoded.mpacket.nrows()
            << " cols: "
            << decoded.mpacket.ncols()
            << std::endl;

    std::cout << "epacket rows: "
            << decoded.epacket.nrows()
            << " cols: "
            << decoded.epacket.ncols()
            << std::endl;


    // Error correction
    std::cout << "Calling correctmesspacket..." << std::endl;

    Correctmesspacket_out corrected =
        correctmesspacket(
            decoded.mpacket,
            decoded.epacket
        );

    std::cout << "correctmesspacket finished" << std::endl;
    std::cout << "Corrected packet first 5 rows:\n";

    for(int r=0;r<5;r++)
    {
        for(int c=0;c<10;c++)
        {
            std::cout << (int)corrected.packet[r][c] << " ";
        }
        std::cout << "\n";
    }

    // Extract payload
    VecUchar plaintext =
        extractplaintext(corrected.packet);

    std::cout << "corrected.packet rows: "
          << corrected.packet.nrows()
          << " cols: "
          << corrected.packet.ncols()
          << std::endl;

    std::cout << "Plaintext first 20:\n";

    for(int i=0;i<20;i++)
    {
        std::cout << (int)plaintext[i] << " ";
    }

    std::cout << "\n";

    std::cout << std::endl;

    std::cout << "plaintext size: "
          << plaintext.size()
          << std::endl;

    std::cout << "First bytes: ";

    for(int i = 0; i < 10 && i < plaintext.size(); i++)
    {
        std::cout << (int)plaintext[i] << " ";
    }

    std::vector<uint8_t> result;

    for(size_t i = 0; i < plaintext.size(); i++)
    {
        result.push_back(plaintext[i]);
    }

    return result;
}