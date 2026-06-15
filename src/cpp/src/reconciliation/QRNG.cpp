/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: QRNG.cpp
 * Description:
 *    This file contains the implementation of the class "QRNG" described in QRNG.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../../include/h_files/reconciliation.hpp"

namespace reconciliation
{
    
    QRNG::QRNG(size_t seq_length, size_t NoF, bool print_flag)
        : sequence_length(seq_length), num_of_frames(NoF), uniform_dist(0, 1), print_flag(print_flag)
    {
    }

    std::vector<std::vector<int8_t>> QRNG::generate_random_sequence()
    {

        std::vector<std::vector<int8_t>> QRNG_output(num_of_frames, std::vector<int8_t>(sequence_length));

#pragma omp parallel
        {
            thread_local auto rng = generate_randomly_seeded_engine();

#pragma omp parallel for
            for (size_t i = 0; i < num_of_frames; i++)
            {
                for (size_t j = 0; j < sequence_length; j++)
                    QRNG_output[i][j] = uniform_dist(rng);
            }
        }

        if (print_flag)
        {
            std::cout << "QRNG output is generated." << std::endl;
        }

        return QRNG_output;
    }

    std::mt19937 QRNG::generate_randomly_seeded_engine()
    {
        std::vector<unsigned int> random_data(std::mt19937::state_size);
        std::random_device source;
        std::generate(random_data.begin(), random_data.end(), [&source]() { return source(); });
        std::seed_seq seeds(random_data.begin(), random_data.end());
        return std::mt19937(seeds);
    }
}
