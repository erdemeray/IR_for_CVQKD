/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: utilities.cpp
 * Description:
 *    This file contains the implementation of the classes/methods described in utilities.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../include/h_files/utilities.hpp"
#include "utilities/statistics.cpp"
#include "utilities/print.cpp"

namespace utilities
{

    print::print(bool flag) : print_flag(flag)
    {
    }

    template<typename T>
    void print::print_line(const std::string& label, const T& value, const std::string& suffix) const {
        std::cout << std::left << std::setw(35) << label << ": " << value << suffix << '\n';
    }

    void print::print_statistics_report(statistics stats) const
    {
        if (!print_flag)
            return;

                    std::cout << "Input Parameters:\n"
                  << "----------------------------------------------------\n";
        print_line("Total number of states read", stats.get_total_number_of_states());
        print_line("Number of unused states", stats.get_number_of_unused_states());
        print_line("Total number of used states", stats.get_total_number_of_states() - stats.get_number_of_unused_states());
        print_line("SNR", 10 * std::log10(1 / stats.get_noise_variance()), " dB");
        std::cout << '\n';

        std::cout << "Multi-dimensional Reconciliation Parameters:\n"
                  << "----------------------------------------------------\n";
        print_line("MDR dimension", stats.get_MDR_dim());
        std::cout << '\n';

        std::cout << "Decoding Statistics:\n"
                  << "----------------------------------------------------\n";
        print_line("Code rate", stats.get_rate());
        print_line("Layered decoding", stats.get_layered_flag());
        print_line("Blocklength", stats.get_blocklength());
        print_line("Total number of frames", stats.get_num_of_decoding_frames());
        print_line("Bit error rate", stats.get_bit_error_rate());
        print_line("Frame error rate", stats.get_frame_error_rate());
        print_line("Avg. number of iterations", stats.get_average_num_of_iterations());
        print_line("Frame error count/total NoF", std::to_string(stats.get_total_frame_error_count()) + " / " + std::to_string(stats.get_num_of_decoding_frames()));
        print_line("Total wrong codeword count", stats.get_count_of_wrong_codewords());
        print_line("Decoding - Elapsed time", stats.get_elapsed_time(), " seconds");

        int num_of_threads = 1;
        #pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();
        }

        print_line("Decoding duration", stats.get_elapsed_time() / static_cast<double>(stats.get_num_of_decoding_frames()) * num_of_threads, " seconds/frame per thread");
        std::cout << '\n';

        std::cout << "CRC Statistics:\n"
                  << "----------------------------------------------------\n";
        print_line("CRC length", stats.get_crc_length());
        print_line("Wrong codewords detected by CRC", stats.get_count_of_wrong_codewords_detected_by_CRC());
        print_line("Undetected frame errors after CRC", stats.get_count_of_undetectable_error_after_CRC());
        std::cout << '\n';

        std::cout << "Output Statistics:\n"
                  << "----------------------------------------------------\n";
        print_line("Total number output bits", (stats.get_num_of_decoding_frames() - stats.get_total_frame_error_count()) * stats.get_blocklength());
        
    }

    double capacity(double snr)
    {
        return 0.5 * std::log2(1 + std::pow(10.0, snr / 10));
    }

    double snr_from_capacity(double capacity)
    {
        return 10 * std::log10(std::pow(2.0, 2 * capacity) - 1);
    }

    std::vector<double> read_quantum_states(const std::string &filename)
    {
        std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open the file: " + filename);
    }

    std::vector<double> states;
    double state;
    while (file >> state) {
        states.push_back(state);
    }

    return states;
    }

} // namespace utilities