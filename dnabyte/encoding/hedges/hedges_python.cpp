#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hedges_wrapper.h"
#include "hedges_encode.h"
#include "hedges_decode.h"

namespace py = pybind11;


PYBIND11_MODULE(hedges_python, m)
{
    m.def(
        "encode_hedges",
        [](py::bytes input)
        {
            std::string buffer = input;

            std::vector<uint8_t> data(
                buffer.begin(),
                buffer.end()
            );

            return hedges_encode(data);
        }
    );
    
    m.def(
        "decode_hedges",
        [](std::vector<std::string> strands)
        {
            return hedges_decode(strands);
        }
    );
}


// cl /EHsc /LD /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\include" /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\Lib\site-packages\pybind11\include" hedges_python.cpp HEDGES-main\DNAcode.cpp HEDGES-main\RSecc.cpp hedges_encode.cpp hedges_wrapper.cpp dna_export.cpp /link /LIBPATH:"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\libs" python311.lib /OUT:hedges_python.pyd

// cl /EHsc /LD /I"C:\Users\kaya\AppData\Local\Programs\Python\Python313\include" /I"C:\Users\kaya\AppData\Local\Programs\Python\Python313\Lib\site-packages\pybind11\include" hedges_python.cpp HEDGES-main\DNAcode.cpp HEDGES-main\RSecc.cpp hedges_encode.cpp hedges_wrapper.cpp dna_export.cpp /link /LIBPATH:"C:\Users\kaya\AppData\Local\Programs\Python\Python313\libs" python313.lib /OUT:hedges_python.pyd