#include <nanobind/nanobind.h>

#include <nanobind/ndarray.h>

#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include "ext.h"

/*
// Minimal example
NB_MODULE(my_ext, m) {
    m.def("add", &add);
}
*/

namespace nb = nanobind;
using namespace nb::literals;
using RGBImage = nb::ndarray<uint8_t, nb::shape<-1, -1, 3>, nb::device::cpu>;

// Nonstandard arithmetic types
// This doesn't seem to work...
// namespace nanobind::detail {
//     template <> struct dtype_traits<std::float16_t> {
//         static constexpr dlpack::dtype value {
//             (uint8_t) dlpack::dtype_code::Float, // type code
//             16, // size in bits
//             1   // lanes (simd), usually set to 1
//         };
//         static constexpr auto name = const_name("float16");
//     };
// }


// Fast array views
// You create a view, and work with that instead of doing arg(i,j)
// NOTE - you need to know array's scalar type and dimensions to use arg(i,j) or v(i,j)
// void fill(nb::ndarray<float, nb::ndim<2>, nb::c_contig, nb::device::cpu> arg) {
//     auto v = arg.view(); // <-- new!

//     for (size_t i = 0; i < v.shape(0); ++i) // Important; use 'v' instead of 'arg' everywhere in loop
//         for (size_t j = 0; j < v.shape(1); ++j)
//             v(i, j) = /* ... */;
// }


// Specializing views at runtime (Fast array views)
// This accepts contiguous CPU arrays, then does runtime checking for specific functionality.
// Useful for when we will know the type of the void * at run-time.
void fill(nb::ndarray<nb::c_contig, nb::device::cpu> arg) {
    if (arg.dtype() == nb::dtype<float>() && arg.ndim() == 2) {
        auto v = arg.view<float, nb::ndim<2>>(); // <-- new!

        for (size_t i = 0; i < v.shape(0); ++i) {
            for (size_t j = 0; j < v.shape(1); ++j) {
                v(i, j) = 42.0;
            }
        }
     } else { /* ... */ }
}


struct Matrix4f { float m[4][4] { }; };

using Array = nb::ndarray<float, nb::numpy, nb::shape<4, 4>, nb::f_contig>;

NB_MODULE(my_ext, m) {
    m.doc() = "A simple example python extension";

    m.def("add", &add, "a"_a, "b"_a = 1,
          "This function adds two numbers and increments if only one is provided.");

    /*************************************************************************/
    // The nb::ndarray<..> class examples

    // Array input arguments
    m.def("inspect", [](const nb::ndarray<>& a) {
        printf("Array data pointer : %p\n", a.data());
        printf("Array dimension : %zu\n", a.ndim());
        for (size_t i = 0; i < a.ndim(); ++i) {
            printf("Array dimension [%zu] : %zu\n", i, a.shape(i));
            printf("Array stride    [%zu] : %zd\n", i, a.stride(i));
        }
        printf("Device ID = %u (cpu=%i, cuda=%i)\n", a.device_id(),
            int(a.device_type() == nb::device::cpu::value),
            int(a.device_type() == nb::device::cuda::value)
        );
        printf("Array dtype: int16=%i, uint32=%i, float32=%i\n",
            a.dtype() == nb::dtype<int16_t>(),
            a.dtype() == nb::dtype<uint32_t>(),
            a.dtype() == nb::dtype<float>()
        );
    });


    // Array constraints (trying to point out RGBImage typedef)
    m.def("process", [](RGBImage data) {
        // Double brightness of the MxNx3 RGB image
        for (size_t y = 0; y < data.shape(0); ++y)
            for (size_t x = 0; x < data.shape(1); ++x)
                for (size_t ch = 0; ch < 3; ++ch)
                    data(y, x, ch) = (uint8_t) std::min(255, data(y, x, ch) * 2);
    });


    // Notably, this doesn't do anything.
    m.def("modify", [](nb::ndarray<uint8_t, nb::shape<-1, -1>>& a) {
        a(0,0) = (uint8_t)1776;
    });


    // Returning arrays from C++ to Python Example
    // - Note the policy
    // - compile-time "Array"/ndarray spec, which isn't useful for us
    nb::class_<Matrix4f>(m, "Matrix4f")
        .def(nb::init<>())
        .def("view",
             [](Matrix4f &data){ return Array(data.m); },
             nb::rv_policy::reference_internal);


    // Data Ownership (Dynamic array configurations)
    m.def("create_2d",
          [](size_t rows, size_t cols) {
              // Allocate a memory region and initialize it
              float *data = new float[rows * cols];
              for (size_t i = 0; i < rows * cols; ++i)
              {
                  data[i] = (float) i;                
              }

              // Delete 'data' when the 'owner' capsule expires
              nb::capsule owner(data, [](void *p) noexcept {
                 delete[] (float *) p;
              });

              return nb::ndarray<nb::numpy, float, nb::ndim<2>>(
                  /* data = */ data,
                  /* shape = */ { rows, cols },
                  /* owner = */ owner
              );
    });


    // Capsule with data structure with multiple regions (Dynamic array configurations)
    m.def("return_multiple", []() {
        struct Temp {
            std::vector<float> vec_1;
            std::vector<float> vec_2;
        };

        Temp *temp = new Temp();
        temp->vec_1 = std::vector<float>{1.0, 2.0, 3.0};
        temp->vec_2 = std::vector<float>{4.0, 5.0, 6.0};

        nb::capsule deleter(temp, [](void *p) noexcept {
            delete (Temp *) p;
        });

        size_t size_1 = temp->vec_1.size();
        size_t size_2 = temp->vec_2.size();

        return std::make_pair(
            nb::ndarray<nb::numpy, float>(temp->vec_1.data(), { size_1 }, deleter),
            nb::ndarray<nb::numpy, float>(temp->vec_2.data(), { size_2 }, deleter)
        );
    });


    // What NOT to do - Returning temporaries
    // (What happens - you get garbage instead of expected data - yikes)
    using Vector3f_Bad = nb::ndarray<float, nb::numpy, nb::shape<3>>;
    m.def("return_vec3_bad", []{
        float data[] { 1, 2, 3 };
        // !!! BAD don't do this !!!
        return Vector3f_Bad(data);
    },
    nb::rv_policy::automatic // default
    // nb::rv_policy::reference_internal // should fail, doesn't - what is "self" here?
    );


    // What you CAN do - Returning temporaries
    // (What happens - calling "help(my_ext.return_vec3_better)" to get
    // more information about the function has a non-descriptive "object" return
    // type instead of numpy array.)
    using Vector3f_Better = nb::ndarray<float, nb::numpy, nb::shape<3>>;
    m.def("return_vec3_better", []{
        float data[] { 1, 2, 3 };
        // OK.
        return nb::cast(Vector3f_Better(data));
    });


    // What is BEST - Returning temporaries
    // (Fixes the "better" problem)
    using Vector3f_Best = nb::ndarray<float, nb::numpy, nb::shape<3>>;
    // Tried some stuff:
    // shape<-1> - not allowed
    // shape<4> - 4th element is garbage
    m.def("return_vec3", []{
        float data[] { 1, 2, 3 };
        // Perfect.
        return Vector3f_Best(data).cast();
    });


    // Array Libraries
    // Section is for if you want to create your own custom Array.
}
