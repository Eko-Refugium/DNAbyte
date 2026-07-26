#include "hedges_wrapper.h"


extern Int messbytesperpacket;
extern Int strandsperpacket;
extern Int messbytesperstrand;
extern Int strandIDbytes;


// Create a HEDGES message packet from arbitrary bytes
MatUchar build_packet(
    const VecUchar& input,
    Int packet_number
)
{
    MatUchar packet(strandsperpacket, messbytesperstrand);

    Int offset = packet_number * messbytesperpacket;


    for (Int strand = 0; strand < strandsperpacket; strand++)
    {
        for (Int byte = 0; byte < messbytesperstrand; byte++)
        {
            Int index =
                offset +
                strand * messbytesperstrand +
                byte;


            if (index < input.size())
            {
                packet[strand][byte] = input[index];
            }
            else
            {
                // padding
                packet[strand][byte] = 0;
            }
        }
    }


    return packet;
}