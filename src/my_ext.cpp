#include <nanobind/nanobind.h>

#include "ext.h"

/*
// Minimal example
NB_MODULE(my_ext, m) {
    m.def("add", &add);
}
*/

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(my_ext, m) {
    m.doc() = "A simple example python extension";

    m.def("add", &add, "a"_a, "b"_a = 1,
          "This function adds two numbers and increments if only one is provided.");
}
