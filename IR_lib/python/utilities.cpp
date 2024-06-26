#include "../cpp/include/h_files/utilities.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

void init_statistics(py::module &m) {

    py::class_<utilities::statistics>(m, "statistics")
        .def(py::init<size_t, size_t,size_t,size_t>(), py::arg("blocklength") = 20000, py::arg("num_of_frames") = 1000, py::arg("MDR_dim") = 8, py::arg("crc_length") = 32)
        .def("set_rate", &utilities::statistics::set_rate)
        .def("set_num_of_quantum_states", &utilities::statistics::set_num_of_quantum_states)
        .def("set_number_of_unused_states", &utilities::statistics::set_number_of_unused_states)
        .def("set_bit_error_count", &utilities::statistics::set_bit_error_count)
        .def("set_frame_error_count", &utilities::statistics::set_frame_error_count)
        .def("set_num_of_iterations", &utilities::statistics::set_num_of_iterations)
        .def("mark_codeword_as_wrong", &utilities::statistics::mark_codeword_as_wrong)
        .def("set_count_of_undetectable_error_after_CRC", &utilities::statistics::set_count_of_undetectable_error_after_CRC)
        .def("set_count_of_wrong_codewords_detected_by_CRC", &utilities::statistics::set_count_of_wrong_codewords_detected_by_CRC)
        .def("start_timing", &utilities::statistics::start_timing)
        .def("end_timing", &utilities::statistics::end_timing)
        .def("get_count_of_undetectable_error_after_CRC", &utilities::statistics::get_count_of_undetectable_error_after_CRC)
        .def("get_count_of_wrong_codewords_detected_by_CRC", &utilities::statistics::get_count_of_wrong_codewords_detected_by_CRC)
        .def("get_elapsed_time", &utilities::statistics::get_elapsed_time)
        .def("get_bit_error_count", &utilities::statistics::get_bit_error_count)
        .def("get_frame_error_count", &utilities::statistics::get_frame_error_count)
        .def("get_num_of_iterations", &utilities::statistics::get_num_of_iterations)
        .def("get_bit_error_rate", &utilities::statistics::get_bit_error_rate)
        .def("get_frame_error_rate", &utilities::statistics::get_frame_error_rate)
        .def("get_average_num_of_iterations", &utilities::statistics::get_average_num_of_iterations)
        .def("get_total_frame_error_count", &utilities::statistics::get_total_frame_error_count)
        .def("get_num_of_decoding_frames", &utilities::statistics::get_num_of_decoding_frames)
        .def("get_count_of_wrong_codewords", &utilities::statistics::get_count_of_wrong_codewords)
        .def("get_total_number_of_states", &utilities::statistics::get_total_number_of_states)
        .def("get_number_of_unused_states", &utilities::statistics::get_number_of_unused_states)
        .def("get_rate", &utilities::statistics::get_rate)
        .def("get_blocklength", &utilities::statistics::get_blocklength)
        .def("get_average_statistics", &utilities::statistics::get_average_statistics)
        .def("get_crc_length", &utilities::statistics::get_crc_length)
        .def("get_MDR_dim", &utilities::statistics::get_MDR_dim)
        .def("get_noise_variance", &utilities::statistics::get_noise_variance)
        .def("set_noise_variance", &utilities::statistics::set_noise_variance)
        .def("get_layered_flag", &utilities::statistics::get_layered_flag);

    py::class_<utilities::print>(m, "print")
        .def(py::init<int>(), py::arg("flag") = 1)
        .def("print_statistics_report", &utilities::print::print_statistics_report);

    m.def("read_quantum_states", &utilities::read_quantum_states, "Reads quantum states from a file.");
    m.def("capacity", &utilities::capacity, "Calculates channel capacity.", py::arg("snr"));
    m.def("snr_from_capacity", &utilities::snr_from_capacity, "Calculates SNR from capacity.", py::arg("capacity"));
}