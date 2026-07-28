#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "HEDGES-main/nr3b.h"
#include "HEDGES-main/RSecc.h"
#include "HEDGES-main/DNAcode.h"
#include "hedges_wrapper.h"
#include "dna_export.h"
#include "hedges_config.h"

extern Int strandsperpacket;


static unsigned char dna_value(char c)
{
    switch(c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }

    throw std::runtime_error("Invalid DNA base");
}


static int decode_base(char c)
{
    switch(c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }

    return 0;
}


// first 4 DNA bases = packet id byte
static uint32_t get_packet_id(const std::string& s)
{
    uint32_t id = 0;

    for(int i=0;i<4;i++)
    {
        id <<= 2;
        id |= decode_base(s[i]);
    }

    return id;
}


// next 4 DNA bases = strand number
static uint32_t get_strand_id(const std::string& s)
{
    uint32_t id = 0;

    for(int i=4;i<8;i++)
    {
        id <<= 2;
        id |= decode_base(s[i]);
    }

    return id;
}



std::vector<uint8_t> hedges_decode(
    const std::vector<std::string>& input,
    const HedgesConfig& cfg
)
{

    std::vector<std::string> strands = input;


    // sort by packet, then strand number
    std::sort(
        strands.begin(),
        strands.end(),
        [](const std::string& a,
           const std::string& b)
        {

            uint32_t pa = get_packet_id(a);
            uint32_t pb = get_packet_id(b);

            if(pa != pb)
                return pa < pb;


            return get_strand_id(a)
                 < get_strand_id(b);
        }
    );



    std::vector<uint8_t> result;



    size_t pos = 0;


    while(pos < strands.size())
    {

        std::vector<std::string> packet;


        uint32_t packet_id =
            get_packet_id(strands[pos]);


        while(pos < strands.size() &&
            packet.size() < strandsperpacket)
        {

            if(get_packet_id(strands[pos]) != packet_id)
                break;


            packet.push_back(
                strands[pos]
            );

            pos++;
        }


    std::cout
        << "Decoding packet "
        << packet_id
        << " strands="
        << packet.size()
        << std::endl;


        std::cout
            << "Decoding packet "
            << packet_id
            << " strands="
            << packet.size()
            << std::endl;



        if(packet.size() != strandsperpacket)
        {
            std::cerr
                << "WARNING: incomplete packet "
                << packet_id
                << std::endl;
        }



        MatUchar dna(
            packet.size(),
            packet[0].size()
        );


        for(size_t r=0;r<packet.size();r++)
        {
            for(size_t c=0;c<packet[r].size();c++)
            {
                dna[r][c] =
                    dna_value(packet[r][c]);
            }
        }



        Dnatomess_out decoded =
            dnatomess(dna);

        for (int r = 0; r < 10; r++)
        {
            std::cout
                << "row " << r
                << " packet=" << (int)decoded.mpacket[r][0]
                << " strand=" << (int)decoded.mpacket[r][1]
                << '\n';
        }

        Correctmesspacket_out corrected =
            correctmesspacket(
                decoded.mpacket,
                decoded.epacket
            );



        VecUchar plain =
            extractplaintext(
                corrected.packet
            );



        for(int i=0;i<plain.size();i++)
        {
            result.push_back(
                plain[i]
            );
        }

    }



    std::cout
        << "Before trim: "
        << result.size()
        << std::endl;



    if(result.size() >= 4)
    {

        uint32_t original_size =
            (uint32_t(result[0]) << 24) |
            (uint32_t(result[1]) << 16) |
            (uint32_t(result[2]) << 8) |
             uint32_t(result[3]);


        std::cout
            << "Original size="
            << original_size
            << std::endl;



        result.erase(
            result.begin(),
            result.begin()+4
        );



        if(original_size <= result.size())
        {
            result.resize(
                original_size
            );
        }
    }



    std::cout
        << "Final decoded size="
        << result.size()
        << std::endl;


    return result;
}
// std::vector<uint8_t> hedges_decode(
//     const std::vector<std::string>& strands
// )
// {
// std::cout << "Before init" << std::endl;

// std::cout << "totstrandlen: " << totstrandlen << std::endl;
// std::cout << "strandIDbytes: " << strandIDbytes << std::endl;
// std::cout << "strandrunoutbytes: " << strandrunoutbytes << std::endl;

//     MatUchar dna(
//         strands.size(),
//         strands[0].size()
//     );

//     for(size_t i = 0; i < strands.size(); i++)
//     {
//         for(size_t j = 0; j < strands[i].size(); j++)
//         {
//             dna[i][j] = dna_value(strands[i][j]);
//         }
//     }

//     std::cout << "Calling dnatomess..." << std::endl;

//     Dnatomess_out decoded = dnatomess(dna);

//     std::cout << "dnatomess finished" << std::endl;

//     std::cout << "mpacket rows: "
//             << decoded.mpacket.nrows()
//             << " cols: "
//             << decoded.mpacket.ncols()
//             << std::endl;

//     std::cout << "epacket rows: "
//             << decoded.epacket.nrows()
//             << " cols: "
//             << decoded.epacket.ncols()
//             << std::endl;


//     // Error correction
//     std::cout << "Calling correctmesspacket..." << std::endl;

//     Correctmesspacket_out corrected =
//         correctmesspacket(
//             decoded.mpacket,
//             decoded.epacket
//         );

//     std::cout << "Decoded mpacket first 5 rows:\n";

//     for (int r = 0; r < 5; r++)
//     {
//         for (int c = 0; c < 10; c++)
//         {
//             std::cout << (int)decoded.mpacket[r][c] << " ";
//         }
//         std::cout << "\n";
//     }

//     // Extract payload
//     VecUchar plaintext =
//         extractplaintext(corrected.packet);

//     std::cout << "corrected.packet rows: "
//           << corrected.packet.nrows()
//           << " cols: "
//           << corrected.packet.ncols()
//           << std::endl;

//     std::cout << "Plaintext first 20:\n";

//     for(int i=0;i<20;i++)
//     {
//         std::cout << (int)plaintext[i] << " ";
//     }

//     std::cout << "\n";

//     std::cout << std::endl;

//     std::cout << "plaintext size: "
//           << plaintext.size()
//           << std::endl;

//     std::cout << "First bytes: ";

//     for(int i = 0; i < 10 && i < plaintext.size(); i++)
//     {
//         std::cout << (int)plaintext[i] << " ";
//     }

//     std::vector<uint8_t> result;

//     for(size_t i = 0; i < plaintext.size(); i++)
//     {
//         result.push_back(plaintext[i]);
//     }

//     return result;
// }