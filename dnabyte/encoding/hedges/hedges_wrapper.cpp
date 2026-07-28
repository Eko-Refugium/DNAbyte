#include "hedges_wrapper.h"

extern Int bytesperstrand;
extern Int messbytesperpacket;
extern Int strandsperpacket;
extern Int messbytesperstrand;
extern Int strandIDbytes;
extern Int strandrunoutbytes;
extern Int strandsperpacketmessage;

MatUchar build_packet(
    const VecUchar& input,
    Int packet_number
)
{
MatUchar packet(strandsperpacket, bytesperstrand, Uchar(0));

for (Int strand = 0; strand < strandsperpacketmessage; strand++)
{
    packet[strand][0] = static_cast<Uchar>(packet_number);
    packet[strand][1] = static_cast<Uchar>(strand);

    for (Int b = 0; b < messbytesperstrand; b++)
    {
        Int index = strand * messbytesperstrand + b;

        if (index < input.size())
            packet[strand][strandIDbytes + b] = input[index];
    }

    packet[strand][bytesperstrand - 2] = 0;
    packet[strand][bytesperstrand - 1] = 0;
}

    return packet;
}