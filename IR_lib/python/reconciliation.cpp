/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: reconciliation.cpp
 * Description:
 *    This file contains the implementation of the Python bindings for the "reconciliation" namespace.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../cpp/include/h_files/reconciliation.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void init_reconciliation(py::module &m){

    py::class_<reconciliation::decoder>(m, "Decoder")
        .def(py::init<size_t, bool, bool, const std::string&, int, bool>(),
             py::arg("NoI") = 500,
             py::arg("layered_flag") = true,
             py::arg("fast_flag") = true,
             py::arg("H_file_name") = PCM_NAME,
             py::arg("lifting_factor") = 5000,
             py::arg("print_flag") = false)
        .def("get_N", &reconciliation::decoder::get_N)
        .def("get_M", &reconciliation::decoder::get_M)
        .def("set_rate", &reconciliation::decoder::set_rate)
        .def("get_rate", &reconciliation::decoder::get_rate)
        .def("get_syndrome", &reconciliation::decoder::get_syndrome<std::vector<uint8_t>>)
        .def("SPA_decoder", &reconciliation::decoder::SPA_decoder<std::vector<double>, std::vector<uint8_t>>)
        .def("get_decoding_error_count", &reconciliation::decoder::get_decoding_error_count<std::vector<uint8_t>, std::vector<uint8_t>>);

    py::class_<reconciliation::MDR>(m, "MDR")
        .def(py::init<int, int>(),
             py::arg("dim") = 8,
             py::arg("print_flag") = false)
        .def("normalize", &reconciliation::MDR::normalize)
                .def("multiplication_Alice", 
             (std::vector<double> (reconciliation::MDR::*)(const std::vector<double>&, const std::vector<double>&) const)
             &reconciliation::MDR::multiplication_Alice)
        .def("multiplication_Bob", 
             (std::vector<double> (reconciliation::MDR::*)(const std::vector<double>&, const std::vector<double>&) const)
             &reconciliation::MDR::multiplication_Bob);

    py::class_<reconciliation::QRNG>(m, "QRNG")
        .def(py::init<size_t, size_t, int>(),
             py::arg("seq_length"),
             py::arg("NoF"),
             py::arg("print_flag") = false)
        .def("generate_random_sequence", &reconciliation::QRNG::generate_random_sequence);

    py::class_<reconciliation::CRC>(m, "CRC")
        .def(py::init<uint32_t>(),
             py::arg("polynomial") = 0xEDB88320)
        .def("get_CRC_checksum", &reconciliation::CRC::get_CRC_checksum<std::vector<uint8_t>>);

    m.def("binary_to_bipolar", &reconciliation::binary_to_bipolar,
          py::arg("QRNG_output"), py::arg("N"), py::arg("num_of_decoding_frames"), py::arg("MDR_dim") = 8);

    m.def("get_LLR", &reconciliation::get_LLR,
          py::arg("synthetic_channel_output"), py::arg("alice_normalization_vector"),
          py::arg("bob_normalization_vector"), py::arg("noise_variance"),
          py::arg("MDR_dim"), py::arg("num_of_decoding_frames"),
          py::arg("N"), py::arg("print_flag") = false);

    m.def("get_CRC_values", &reconciliation::get_CRC_values<uint8_t>,
         py::arg("decoded_frame"), py::arg("crc"), py::arg("num_of_decoding_frames"));

    m.def("check_CRC", &reconciliation::check_CRC,
          py::arg("CRC_decoded"), py::arg("CRC_QRNG"), py::arg("syndrome_flags"),
          py::arg("stats"), py::arg("num_of_decoding_frames"), py::arg("print_flag") = false);

    m.def("calculate_syndrome", &reconciliation::calculate_syndrome,
          py::arg("LDPC_decoder"), py::arg("QRNG_output"),
          py::arg("num_of_decoding_frames"), py::arg("print_flag") = false); 
    
    m.def("reconcile", &reconciliation::reconcile, "Performs reconciliation.",
          py::arg("alice_states"),
          py::arg("bob_states"),
          py::arg("rate"),
          py::arg("noise_variance"),
          py::arg("NoI") = 500,
          py::arg("MDR_dim") = 8,
          py::arg("layered_decoding") = true,
          py::arg("fast_decoding") = true,
          py::arg("print_flag") = false,
          py::arg("H_file_name") = PCM_NAME,
          py::arg("lifting_factor") = 5000);

}
