#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hedges_wrapper.h"
#include "hedges_encode.h"


namespace py = pybind11;


PYBIND11_MODULE(hedges_cpp, m)
{
    m.def(
        "encode_hedges",
        &hedges_encode,
        "Encode binary data into DNA"
    );
}


cl /EHsc /LD /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\include" /I"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\Lib\site-packages\pybind11\include" hedges_python.cpp HEDGES-main\DNAcode.cpp HEDGES-main\RSecc.cpp hedges_encode.cpp hedges_wrapper.cpp dna_export.cpp /link /LIBPATH:"C:\Users\Kaya\AppData\Local\Programs\Python\Python311\libs" python311.lib /OUT:hedges_python.pyd