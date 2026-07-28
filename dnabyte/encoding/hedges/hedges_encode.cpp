#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "HEDGES-main/nr3b.h"
#include "HEDGES-main/RSecc.h"
#include "HEDGES-main/DNAcode.h"
#include "hedges_wrapper.h"
#include "dna_export.h"
#include "hedges_config.h"

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


std::vector<std::string> hedges_encode(
    const std::vector<uint8_t>& data,
    const HedgesConfig& cfg
)
{
    // if (argc != 3)
    // {
    //     std::cerr << "Usage:\n";
    //     std::cerr << "hedges_encode input.bin output.json\n";
    //     return 1;
    // }

    // std::string input_file = argv[1];
    // std::string output_file = argv[2];

    // std::ifstream file(input_file, std::ios::binary);

    // if (!file)
    // {
    //     std::cerr << "Cannot open " << input_file << std::endl;
    //     return 1;
    // }

    // std::vector<uint8_t> data(
    //     (std::istreambuf_iterator<char>(file)),
    //     std::istreambuf_iterator<char>()
    // );

    // std::cout << "Read " << data.size() << " bytes\n";
    // std::cout << "Output will be " << output_file << "\n";

            
    Doub coderates_[] = { 0., 0.75, 0.6, 0.5, 1. / 3., 0.25, 1. / 6. }; // table of coderates 1..6
    VecDoub coderates(7, coderates_);

    // see the PNAS paper to understand the intended use of primers
    // primers
    std::vector<char> leftprimer_s(
        cfg.left_primer.begin(),
        cfg.left_primer.end()
    );

    leftprimer_s.push_back('\0');


    std::vector<char> rightprimer_s(
        cfg.right_primer.begin(),
        cfg.right_primer.end()
    );

    rightprimer_s.push_back('\0');


    // user parameters
    Int coderatecode =
        cfg.coderate;

    Int hlimit =
        cfg.heap_limit;


    // strand parameters
    totstrandlen =
        cfg.total_strand_length;

    strandIDbytes =
        cfg.strand_id_bytes;

    strandrunoutbytes =
        cfg.strand_runout_bytes;

    // sub,del,ins rates to simulate in this test (as multiple of our experimentally observed values):
    Doub ratefac = 1.5;
    Doub srate = ratefac * 0.0238;
    Doub drate = ratefac * 0.0082;
    Doub irate = ratefac * 0.0039;

    // set parameters for DNA constrants (normally not changed, except for no constraint)
    Int max_hpoly_run =
        cfg.max_homopolymer;

    Int GC_window =
        cfg.gc_window;

    Int max_GC =
        cfg.gc_max;
    Int min_GC = GC_window - max_GC;

    // get and reset parameters per above
    Getparams_out params = getparams();
    Int NSALT = params.NSALT;
    Int MAXSEQ = params.MAXSEQ;
    Int NSTAK = params.NSTAK;
    Int HLIMIT = params.HLIMIT;
    setparams(8 * strandIDbytes, MAXSEQ, NSTAK, hlimit); // change NSALT and HLIMIT
    setcoderate(
    coderatecode,
    leftprimer_s.data(),
    rightprimer_s.data()
    );
    setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run); // set DNA constraints (see paper)

    // set values of global variables derived from above
    leftlen = Int(strlen(leftprimer_s.data()));
    rightlen = Int(strlen(rightprimer_s.data()));
    strandlen = totstrandlen - leftlen - rightlen;
    strandsperpacketmessage = strandsperpacket - strandsperpacketcheck;
    bytesperstrand = Int(strandlen * coderates[coderatecode] / 4.);
    messbytesperstrand = bytesperstrand - strandIDbytes - strandrunoutbytes; // payload bytes per strand
    messbytesperpacket = strandsperpacketmessage * messbytesperstrand; // payload bytes per packet of 255 strands

    std::cout << "HEDGES initialized.\n";
    std::cout << "Message bytes per packet: "
            << messbytesperpacket << std::endl;


    std::vector<std::string> all_strands;

    size_t offset = 0;
    uint32_t packet_id = 0;

    std::cout << "data.size() = "
          << data.size()
          << std::endl;

    while(offset < data.size())
    {
        size_t remaining = data.size() - offset;

        size_t capacity = messbytesperpacket;

        if(packet_id == 0)
        {
            capacity -= 4;
        }

        size_t chunk_size = std::min(
            remaining,
            capacity
        );
        std::cout << "remaining = " << remaining << '\n';
        std::cout << "chunk_size = " << chunk_size << '\n';

        std::vector<uint8_t> chunk(
            data.begin() + offset,
            data.begin() + offset + chunk_size
        );

        // Add length only if decoder needs to know final size
        // (usually yes)
        std::vector<uint8_t> payload;


        // Add original size ONLY on first packet
        if(packet_id == 0)
        {
            uint32_t size = data.size();

            payload.push_back((size >> 24) & 0xff);
            payload.push_back((size >> 16) & 0xff);
            payload.push_back((size >> 8) & 0xff);
            payload.push_back(size & 0xff);
        }


        // Add actual data
        payload.insert(
            payload.end(),
            chunk.begin(),
            chunk.end()
        );


        // pad packet to exact message size
        while(payload.size() < messbytesperpacket)
        {
            payload.push_back(0);
}

        VecUchar input(payload.size());

        for(size_t i=0;i<payload.size();i++)
            input[i]=payload[i];

        std::cout << "input.size() = "
          << input.size()
          << std::endl;


        MatUchar packet = build_packet(
            input,
            packet_id
        );

        std::cout << "packet rows: "
                << packet.nrows()
                << " cols: "
                << packet.ncols()
                << std::endl;

        std::cout << "packet first row: ";

        for(int i=0;i<10;i++)
        {
            std::cout << (int)packet[0][i] << " ";
        }

        std::cout << std::endl;

        MatUchar protected_packet =
            protectmesspacket(packet);

        std::cout << "Protected packet first 5 rows:\n";

        for (int r = 0; r < 5; r++)
        {
            for (int c = 0; c < 10; c++)
            {
                std::cout << (int)protected_packet[r][c] << " ";
            }
            std::cout << "\n";
        }


        MatUchar dna =
            messtodna(protected_packet);




        auto strands = dna_to_strings(dna);


        all_strands.insert(
            all_strands.end(),
            strands.begin(),
            strands.end()
        );


        offset += chunk_size;
        packet_id++;
    }

    return all_strands;
}