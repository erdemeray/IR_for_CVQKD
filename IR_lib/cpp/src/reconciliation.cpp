/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: reconciliation.cpp
 * Description:
 *    This file contains the implementation of the classes/methods described in reconciliation.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../include/h_files/reconciliation.hpp"

#include "reconciliation/decoder.cpp"
#include "reconciliation/MDR.cpp"
#include "reconciliation/QRNG.cpp"
#include "reconciliation/CRC.cpp"

namespace reconciliation
{
    std::vector<double> binary_to_bipolar(const std::vector<std::vector<int8_t>> &QRNG_output, const size_t N, const size_t num_of_decoding_frames, const int MDR_dim)
    {

        std::vector<double> QRNG_output_bipolar(num_of_decoding_frames * N, 0);

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                QRNG_output_bipolar[i * N + j] = (-2 * QRNG_output[i][j] + 1) / std::sqrt(MDR_dim);
            }
        }

        return QRNG_output_bipolar;
    }

    std::vector<std::vector<double>> get_LLR(std::vector<double> synthetic_channel_output, const std::vector<double> &alice_normalization_vector, const std::vector<double> &bob_normalization_vector, double noise_variance, int MDR_dim, size_t num_of_decoding_frames, int N, bool print_flag)
    {

#pragma omp parallel for
        for (size_t i = 0; i < synthetic_channel_output.size(); i++)
        {
            double Lc = 2 * alice_normalization_vector[i] * bob_normalization_vector[i] / noise_variance / std::sqrt(MDR_dim);

            synthetic_channel_output[i] = Lc * synthetic_channel_output[i];
        }

        // Reshape the LLRs to the LDPC frames

        std::vector<std::vector<double>> LLR(num_of_decoding_frames, std::vector<double>(N, 0));

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            std::copy(synthetic_channel_output.begin() + i * N, synthetic_channel_output.begin() + (i + 1) * N,  LLR[i].begin());
        }

        if (print_flag)
        {
            std::cout << "LLR values are calculated." << std::endl;
        }

        return LLR;
    }

    template <typename T>
    std::vector<std::vector<uint8_t>> get_CRC_values(const std::vector<std::vector<T>> &decoded_frame, CRC &crc, size_t num_of_decoding_frames)
    {

        std::vector<std::vector<uint8_t>> CRC_decoded(num_of_decoding_frames, std::vector<uint8_t>(32, 0));

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            auto temp = crc.get_CRC_checksum(decoded_frame[i]);
            std::copy(temp.begin(), temp.end(), CRC_decoded[i].begin());
        }

        return CRC_decoded;
    }

    std::tuple<double, int> check_CRC(const std::vector<std::vector<uint8_t>> &CRC_decoded, const std::vector<std::vector<uint8_t>> &CRC_QRNG, const std::vector<int> &syndrome_flags, utilities::statistics &stats, size_t num_of_decoding_frames, bool print_flag)
    {
        std::vector<double> CRC_mismatch(num_of_decoding_frames, 0);
        int undetected_error = 0;
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {

            if (syndrome_flags[i] == 1)
            {
                if (!std::equal(CRC_decoded[i].begin(), CRC_decoded[i].end(), CRC_QRNG[i].begin()))
                {
                    CRC_mismatch[i] = 1;
                }
                else if (stats.get_frame_error_count(i) == 1)
                {
                    ++undetected_error;
                }
            }
        }

        double No_CRC_mismatch = std::accumulate(CRC_mismatch.begin(), CRC_mismatch.end(), 0);

        if (print_flag)
        {
            std::cout << "Number of decoded codeword CRC mismatches: " << No_CRC_mismatch << std::endl;
            std::cout << "Number of undetected erroneous codewords: " << undetected_error << std::endl;
        }

        return std::make_tuple(No_CRC_mismatch, undetected_error);
    }

    std::vector<std::vector<uint8_t>> calculate_syndrome(const decoder LDPC_decoder, const std::vector<std::vector<int8_t>>& QRNG_output, const size_t num_of_decoding_frames, bool print_flag)
    {

        std::vector<std::vector<uint8_t>> syndrome(num_of_decoding_frames, std::vector<uint8_t>(LDPC_decoder.get_M(), 0));
#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            std::vector<uint8_t> temp_syndrome = LDPC_decoder.get_syndrome(QRNG_output[i]);
            std::copy(temp_syndrome.begin(), temp_syndrome.end(), syndrome[i].begin());
        }

        if (print_flag)
            std::cout << "Syndromes are calculated." << std::endl;

        return syndrome;
    }

    utilities::statistics reconcile(const std::vector<double> &alice_states_input, const std::vector<double> &bob_states_input, const double rate, const double noise_variance, const size_t NoI, const int MDR_dim, bool layered_flag, bool fast_flag, bool print_flag, const std::string H_file_name, int lifting_factor)
    {
        // make a copy of alice and bob states
        std::vector<double> alice_states = alice_states_input;
        std::vector<double> bob_states = bob_states_input;

        // Initialize the reconciliation parameters and objects
        size_t num_of_total_states = alice_states.size();

        reconciliation::decoder LDPC_decoder(NoI, layered_flag, fast_flag, H_file_name, lifting_factor);

        LDPC_decoder.set_rate(rate);

        double actual_rate = LDPC_decoder.get_rate();

        size_t N = LDPC_decoder.get_N();

        size_t num_of_decoding_frames = static_cast<size_t>(std::floor((double)(num_of_total_states) / N));

        // Make sure we don't use more resources than we need

        size_t num_of_threads;
#pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();    
        }

        num_of_threads = std::min(num_of_threads, num_of_decoding_frames);
        
        omp_set_num_threads(num_of_threads);

        utilities::statistics stats(N, num_of_decoding_frames, MDR_dim, layered_flag);

        stats.set_noise_variance(noise_variance);

        if (print_flag)
        {
            std::cout << "Number of quantum states: " << num_of_total_states << std::endl;
            std::cout << "Number of decoding frames: " << num_of_decoding_frames << std::endl;
            std::cout << "Blocklength: " << N << std::endl;
            std::cout << "Number of unused quantum states:" << num_of_total_states - num_of_decoding_frames * N << std::endl;
            std::cout << std::endl << std::endl;
            std::cout << "Bob :" << std::endl;
        }

        stats.set_rate(actual_rate);
        stats.set_num_of_quantum_states(num_of_total_states);
        stats.set_number_of_unused_states(num_of_total_states - num_of_decoding_frames * N);

        QRNG qrng(N, num_of_decoding_frames);
        MDR mdr(MDR_dim);
        CRC crc;

        // Generate the random raw key material using QRNG and convert it to bipolar form
        auto QRNG_output = qrng.generate_random_sequence();
        auto QRNG_output_bipolar = binary_to_bipolar(QRNG_output, N, num_of_decoding_frames, MDR_dim);

        // Normalize the quantum states and perform the MDR multiplication at Bobs's side
        auto bob_normalization_vector = mdr.normalize(bob_states);
        auto channel_message = mdr.multiplication_Bob(bob_states, QRNG_output_bipolar);

        // Calculate the syndromes of the QRNG output
        auto syndrome = calculate_syndrome(LDPC_decoder, QRNG_output, num_of_decoding_frames);

        // Perform the MDR multiplication at Alice's side
        if (print_flag)
        {
            std::cout << std::endl;
            std::cout << "Alice :" << std::endl;
        }

        auto alice_normalization_vector = mdr.normalize(alice_states);
        auto synthetic_channel_output = mdr.multiplication_Alice(channel_message, alice_states);

        // Calculate the LLR values
        auto LLR = get_LLR(synthetic_channel_output, alice_normalization_vector, bob_normalization_vector, noise_variance, MDR_dim, num_of_decoding_frames, N);

        if (print_flag)
        {
            std::cout << "Decoding:" << std::endl;
        }

        // Perform the LDPC decoding
        std::vector<std::vector<uint8_t>> decoded_frame(num_of_decoding_frames, std::vector<uint8_t>(N, 0));
        std::vector<int> syndrome_flags(num_of_decoding_frames, 0);

        stats.start_timing();

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            auto result = LDPC_decoder.SPA_decoder(LLR[i], syndrome[i]);
            auto temp = std::get<0>(result);
            auto iterations = std::get<1>(result);
            auto syndrome_flag = std::get<2>(result);

            std::copy(temp.begin(), temp.end(), decoded_frame[i].begin());

            stats.set_bit_error_count(i, LDPC_decoder.get_decoding_error_count(decoded_frame[i], QRNG_output[i]));
            if (stats.get_bit_error_count(i) > 0)
                stats.set_frame_error_count(i);

            stats.set_num_of_iterations(i, iterations);
            syndrome_flags[i] = syndrome_flag;

            if ((syndrome_flags[i] == 1) && (stats.get_frame_error_count(i) == 1))
            {
                stats.mark_codeword_as_wrong(i);
            }
        }

        stats.end_timing();

        // Calculate the CRC of the decoded frames:
        auto CRC_decoded = get_CRC_values(decoded_frame, crc, num_of_decoding_frames);

        // Calculate the CRC of the QRNG output
        auto CRC_QRNG = get_CRC_values(QRNG_output, crc, num_of_decoding_frames);

        // Calculate the CRC mismatchs and compare it to the frame errors
        auto crc_result = check_CRC(CRC_decoded, CRC_QRNG, syndrome_flags, stats, num_of_decoding_frames);
        auto wrong_codewords_detected = std::get<0>(crc_result);
        auto undetectable_errors = std::get<1>(crc_result);

        stats.set_count_of_undetectable_error_after_CRC(undetectable_errors);
        stats.set_count_of_wrong_codewords_detected_by_CRC(wrong_codewords_detected);

        if (print_flag)
        {
            utilities::print display;
            display.print_statistics_report(stats);
        }

        return stats;
    }

    template std::vector<double> MDR::multiplication_Alice<std::vector<double>, std::vector<double>>(
        const std::vector<double> &, const std::vector<double> &) const;

    template std::vector<double> MDR::multiplication_Bob<std::vector<double>, std::vector<double>>(
        const std::vector<double> &, const std::vector<double> &) const;

    template std::vector<uint8_t> decoder::get_syndrome<std::vector<uint8_t>>( const std::vector<uint8_t> &) const;

    template std::tuple<std::vector<uint8_t>, int, bool> decoder::SPA_decoder<std::vector<double>, std::vector<uint8_t>>(const std::vector<double> &, const std::vector<uint8_t> &) const;

    template int decoder::get_decoding_error_count<std::vector<uint8_t>, std::vector<uint8_t>>(const std::vector<uint8_t> &, const std::vector<uint8_t> &) const;

    template std::vector<double> MDR::compute_higher_dimension_conjugate<std::vector<double>>(const std::vector<double> &) const;

    template std::vector<int> CRC::get_CRC_checksum<std::vector<uint8_t>>(const std::vector<uint8_t> &);

    template std::vector<std::vector<uint8_t>> get_CRC_values<uint8_t>(const std::vector<std::vector<uint8_t>> &, CRC &, size_t);

}
