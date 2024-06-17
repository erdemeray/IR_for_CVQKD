#include "../cpp/include/h_files/reconciliation.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void init_reconciliation(py::module &m){

    m.def("reconcile", &reconciliation::reconcile, "Performs reconciliation.",
          py::arg("alice_states"),
          py::arg("bob_states"),
          py::arg("rate"),
          py::arg("noise_variance"),
          py::arg("NoI") = 500,
          py::arg("MDR_dim") = 8,
          py::arg("layered_decoding") = 1,
          py::arg("fast_decoding") = 1,
          py::arg("print_flag") = 0,
          py::arg("H_file_name") = "R_0p2_R_0p01_RA",
          py::arg("lifting_factor") = 5000);

}
