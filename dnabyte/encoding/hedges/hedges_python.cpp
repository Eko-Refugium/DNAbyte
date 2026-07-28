#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hedges_wrapper.h"
#include "hedges_encode.h"
#include "hedges_decode.h"
#include "hedges_config.h"

namespace py = pybind11;


PYBIND11_MODULE(hedges_python, m)
{
    py::class_<HedgesConfig>(m, "HedgesConfig")
        .def(py::init<>())

        .def_readwrite("coderate",
            &HedgesConfig::coderate)

        .def_readwrite("total_strand_length",
            &HedgesConfig::total_strand_length)

        .def_readwrite("strand_id_bytes",
            &HedgesConfig::strand_id_bytes)

        .def_readwrite("strand_runout_bytes",
            &HedgesConfig::strand_runout_bytes)

        .def_readwrite("gc_window",
            &HedgesConfig::gc_window)

        .def_readwrite("gc_max",
            &HedgesConfig::gc_max)

        .def_readwrite("max_homopolymer",
            &HedgesConfig::max_homopolymer)

        .def_readwrite("heap_limit",
            &HedgesConfig::heap_limit)

        .def_readwrite("left_primer",
            &HedgesConfig::left_primer)

        .def_readwrite("right_primer",
            &HedgesConfig::right_primer);
    m.def(
        "encode_hedges",
        [](py::bytes input, const HedgesConfig& cfg)
        {
            std::string buffer = input;

            std::vector<uint8_t> data(
                buffer.begin(),
                buffer.end()
            );

            return hedges_encode(data, cfg);
        }
    );
    
    m.def(
        "decode_hedges",
        [](std::vector<std::string> strands, const HedgesConfig& cfg)
        {
            return hedges_decode(strands, cfg);
        }
    );
}


// cl /EHsc /LD /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\include" /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\Lib\site-packages\pybind11\include" hedges_python.cpp HEDGES-main\DNAcode.cpp HEDGES-main\RSecc.cpp hedges_encode.cpp hedges_wrapper.cpp dna_export.cpp hedges_decode.cpp /link /LIBPATH:"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\libs" python311.lib /OUT:hedges_python.pyd

// cl /EHsc /LD /I"C:\Users\kaya\AppData\Local\Programs\Python\Python313\include" /I"C:\Users\kaya\AppData\Local\Programs\Python\Python313\Lib\site-packages\pybind11\include" hedges_python.cpp HEDGES-main\DNAcode.cpp HEDGES-main\RSecc.cpp hedges_encode.cpp hedges_wrapper.cpp dna_export.cpp /link /LIBPATH:"C:\Users\kaya\AppData\Local\Programs\Python\Python313\libs" python313.lib /OUT:hedges_python.pyd