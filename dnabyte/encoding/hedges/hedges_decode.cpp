#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "HEDGES-main/nr3b.h"
#include "HEDGES-main/RSecc.h"
#include "HEDGES-main/DNAcode.h"

#include "hedges_decode.h"
#include "hedges_config.h"
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

extern Int strandsperpacket;
extern Int strandsperpacketcheck;


static unsigned char dna_value(char c)
{
    switch (c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }

    throw std::runtime_error("Invalid DNA base");
}


static void initialize_hedges(const HedgesConfig& cfg)
{
    std::cerr << ">>> INIT 1: entered" << std::endl;

    Doub coderates_[] =
    {
        0.,
        0.75,
        0.6,
        0.5,
        1. / 3.,
        0.25,
        1. / 6.
    };

    VecDoub coderates(7, coderates_);

    std::cerr << ">>> INIT 2: coderates created" << std::endl;

    std::vector<char> leftprimer(
        cfg.left_primer.begin(),
        cfg.left_primer.end()
    );
    leftprimer.push_back('\0');

    std::vector<char> rightprimer(
        cfg.right_primer.begin(),
        cfg.right_primer.end()
    );
    rightprimer.push_back('\0');

    std::cerr << ">>> INIT 3: primers created" << std::endl;

    totstrandlen = cfg.total_strand_length;
    strandIDbytes = cfg.strand_id_bytes;
    strandrunoutbytes = cfg.strand_runout_bytes;

    std::cerr << ">>> INIT 4: globals set" << std::endl;

    Getparams_out params = getparams();

    std::cerr << ">>> INIT 5: getparams returned" << std::endl;

    setparams(
        8 * strandIDbytes,
        params.MAXSEQ,
        params.NSTAK,
        cfg.heap_limit
    );

    std::cerr << ">>> INIT 6: setparams returned" << std::endl;

    setcoderate(
        cfg.coderate,
        leftprimer.data(),
        rightprimer.data()
    );

    std::cerr << ">>> INIT 7: setcoderate returned" << std::endl;

    Int min_GC = cfg.gc_window - cfg.gc_max;

    setdnaconstraints(
        cfg.gc_window,
        cfg.gc_max,
        min_GC,
        cfg.max_homopolymer
    );

    std::cerr << ">>> INIT 8: setdnaconstraints returned" << std::endl;

    leftlen = Int(cfg.left_primer.size());
    rightlen = Int(cfg.right_primer.size());

    strandlen = totstrandlen - leftlen - rightlen;

    strandsperpacketmessage =
        strandsperpacket - strandsperpacketcheck;

    bytesperstrand =
        Int(
            strandlen
            * coderates[cfg.coderate]
            / 4.
        );

    messbytesperstrand =
        bytesperstrand
        - strandIDbytes
        - strandrunoutbytes;

    messbytesperpacket =
        strandsperpacketmessage
        * messbytesperstrand;

    std::cerr << ">>> INIT 9: initialization complete" << std::endl;
}


std::vector<uint8_t> hedges_decode(
    const std::vector<std::string>& input,
    const HedgesConfig& cfg
)
{
    std::cerr << ">>> DECODE 1: entered hedges_decode" << std::endl;

    std::cerr << ">>> DECODE 2: before initialize_hedges" << std::endl;
    initialize_hedges(cfg);
    std::cerr << ">>> DECODE 3: after initialize_hedges" << std::endl;

    std::vector<uint8_t> result;

    size_t offset = 0;

    std::cerr << ">>> DECODE 4: input.size = "
              << input.size() << std::endl;

    while (offset < input.size())
    {
        std::cerr << ">>> DECODE 5: packet offset = "
                  << offset << std::endl;

        size_t remaining = input.size() - offset;

        size_t count =
            std::min(
                static_cast<size_t>(strandsperpacket),
                remaining
            );

        std::cerr << ">>> DECODE 6: count = "
                  << count
                  << ", strandsperpacket = "
                  << strandsperpacket
                  << std::endl;


        if (count != strandsperpacket)
        {
            std::cerr
                << "Incomplete HEDGES packet: "
                << count
                << " / "
                << strandsperpacket
                << " strands"
                << std::endl;

            break;
        }


        //
        // Build ONE DNA packet
        //
        MatUchar dna(
            strandsperpacket,
            totstrandlen,
            Uchar(0)
        );


        for (size_t r = 0; r < count; ++r)
        {
            const std::string& strand =
                input[offset + r];

            size_t len =
                std::min(
                    strand.size(),
                    static_cast<size_t>(totstrandlen)
                );

            for (size_t c = 0; c < len; ++c)
            {
                dna[r][c] =
                    dna_value(strand[c]);
            }
        }


        //
        // HEDGES decoding.
        //
        // This calls the ORIGINAL DNAcode.cpp implementation.
        //
        Dnatomess_out decoded =
            dnatomess(dna);


        std::cout
            << "HEDGES bad decodes: "
            << decoded.baddecodes
            << std::endl;


        std::cout
            << "HEDGES erasures: "
            << decoded.erasures
            << std::endl;


        std::cout << ">>> BEFORE correctmesspacket" << std::endl;

        //
        // Reed-Solomon correction
        //
        Correctmesspacket_out corrected =
            correctmesspacket(
                decoded.mpacket,
                decoded.epacket
            );

        std::cout << ">>> AFTER correctmesspacket" << std::endl;
        //
        // Extract the actual message bytes
        //
        VecUchar plain =
            extractplaintext(
                corrected.packet
            );

        std::cout << ">>> AFTER extractplaintext" << std::endl;
        //
        // Append this packet to the FINAL result.
        //
        for (Int i = 0; i < plain.size(); ++i)
        {
            result.push_back(
                plain[i]
            );
        }

        std::cout << ">>> HEDGES PACKET LOOP FINISHED" << std::endl;


        //
        // Move to next packet.
        //
        offset += count;
    }


    //
    // Remove 4-byte original-size header.
    //
    if (result.size() < 4)
    {
        throw std::runtime_error(
            "Decoded data is too short"
        );
    }


    uint32_t original_size =
        (uint32_t(result[0]) << 24) |
        (uint32_t(result[1]) << 16) |
        (uint32_t(result[2]) << 8) |
        uint32_t(result[3]);


    result.erase(
        result.begin(),
        result.begin() + 4
    );


    //
    // Remove packet padding.
    //
    if (original_size > result.size())
    {
        throw std::runtime_error(
            "Decoded data is shorter than original size"
        );
    }


    result.resize(
        original_size
    );

    std::cout
        << ">>> HEDGES DECODE RETURNING "
        << result.size()
        << " BYTES"
        << std::endl;

    return result;
}