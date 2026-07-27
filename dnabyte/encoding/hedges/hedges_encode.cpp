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


std::vector<std::string> hedges_encode(
    const std::vector<uint8_t>& data
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
    char leftprimer_s[] = "TCGAAGTCAGCGTGTATTGTATG"; // _s means "as a string"
    char rightprimer_s[] = "TAGTGAGTGCGATTAAGCGTGTT"; // for direct right appending (no reverse complement)

    // user-settable parameters for this test
    Int coderatecode = 3; // test this coderate in coderates table above
    Int npackets = 20; // number of packets (of 255 strands each) to generate and test
    Int hlimit = 1000000; // maximum size of decode heap, see paper

    // these lines are setting global variables in DNAcode.cpp
    totstrandlen = 300; // total length of DNA strand
    strandIDbytes = 2; // ID bytes each strand for packet and sequence number (is global)
    strandrunoutbytes = 2; // confirming bytes end of each strand (see paper)

    // sub,del,ins rates to simulate in this test (as multiple of our experimentally observed values):
    Doub ratefac = 1.5;
    Doub srate = ratefac * 0.0238;
    Doub drate = ratefac * 0.0082;
    Doub irate = ratefac * 0.0039;

    // set parameters for DNA constrants (normally not changed, except for no constraint)
    Int max_hpoly_run = 4; // max homopolymer length allowed (0 for no constraint)
    Int GC_window = 12; // window for GC count (0 for no constraint)
    Int max_GC = 8; // max GC allowed in window (0 for no constraint)
    Int min_GC = GC_window - max_GC;

    // get and reset parameters per above
    Getparams_out params = getparams();
    Int NSALT = params.NSALT;
    Int MAXSEQ = params.MAXSEQ;
    Int NSTAK = params.NSTAK;
    Int HLIMIT = params.HLIMIT;
    setparams(8 * strandIDbytes, MAXSEQ, NSTAK, hlimit); // change NSALT and HLIMIT
    setcoderate(coderatecode, leftprimer_s, rightprimer_s); // set code rate with left and right primers
    setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run); // set DNA constraints (see paper)

    // set values of global variables derived from above
    leftlen = Int(strlen(leftprimer_s));
    rightlen = Int(strlen(rightprimer_s));
    strandlen = totstrandlen - leftlen - rightlen;
    strandsperpacketmessage = strandsperpacket - strandsperpacketcheck;
    bytesperstrand = Int(strandlen * coderates[coderatecode] / 4.);
    messbytesperstrand = bytesperstrand - strandIDbytes - strandrunoutbytes; // payload bytes per strand
    messbytesperpacket = strandsperpacket * messbytesperstrand; // payload bytes per packet of 255 strands

    std::cout << "HEDGES initialized.\n";
    std::cout << "Message bytes per packet: "
            << messbytesperpacket << std::endl;


    std::vector<std::string> all_strands;

    size_t offset = 0;
    uint32_t packet_id = 0;

    while(offset < data.size())
    {
        size_t remaining = data.size() - offset;

        size_t chunk_size = std::min(
            remaining,
            static_cast<size_t>(messbytesperpacket)
        );

        std::vector<uint8_t> chunk(
            data.begin() + offset,
            data.begin() + offset + chunk_size
        );

        // Add length only if decoder needs to know final size
        // (usually yes)
        std::vector<uint8_t> payload;

        

        payload.insert(
            payload.end(),
            chunk.begin(),
            chunk.end()
        );


        VecUchar input(payload.size());

        for(size_t i=0;i<payload.size();i++)
            input[i]=payload[i];


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