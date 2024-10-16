/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: reconciliation.hpp
 * Description: 
 *    This file contains the classes and functions to perform information reconciliation.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <sstream>
#include <omp.h>
#include <cstdint>

#include "utilities.hpp"
#pragma once

#include "reconciliation/decoder.hpp"
#include "reconciliation/MDR.hpp"
#include "reconciliation/QRNG.hpp"
#include "reconciliation/CRC.hpp"

namespace reconciliation
{

  /**
     * @brief Convert a binary vector to a bipolar vector.
     *
     * @param QRNG_output Binary vector
     * @param N Blocklength of the codeword
     * @param num_of_decoding_frames Number of frames
     * @param MDR_dim Dimension of the MDR, default = 8
     * @return std::vector<double> Bipolar vector
     */
    std::vector<double> binary_to_bipolar(const std::vector<std::vector<int8_t>> &QRNG_output, const size_t N, const size_t num_of_decoding_frames, const int MDR_dim = 8);

    /**
     * @brief Calculate the log-likelihood ratio (LLR) values from the synthetic channel output.
     *
     * @param synthetic_channel_output Synthetic channel output
     * @param alice_normalization_vector Normalization vector at Alice's side
     * @param bob_normalization_vector Normalization vector at Bob's side
     * @param noise_variance Noise variance of the channel
     * @param MDR_dim Dimension of the MDR
     * @param num_of_decoding_frames Number of frames
     * @param N Blocklength of the codeword
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
     * @return std::vector<std::vector<double>> LLR values
     */
    std::vector<std::vector<double>> get_LLR(std::vector<double> synthetic_channel_output, const std::vector<double> &alice_normalization_vector, const std::vector<double> &bob_normalization_vector, double noise_variance, int MDR_dim, size_t num_of_decoding_frames, int N, bool print_flag = false);

    /**
     * @brief Calculate the CRC values of the decoded frames.
     *
     * @param decoded_frame Decoded frames
     * @param crc CRC object
     * @param num_of_decoding_frames Number of frames
     * @return std::vector<std::vector<uint8_t>> CRC values of the decoded frames
     */
    template <typename T>
    std::vector<std::vector<uint8_t>> get_CRC_values(const std::vector<std::vector<T>> &decoded_frame, CRC &crc, size_t num_of_decoding_frames);

    /**
     * @brief Check the CRC values of the decoded frames against the CRC values of the QRNG output.
     *
     * @param CRC_decoded CRC values of the decoded frames
     * @param CRC_QRNG CRC values of the QRNG output
     * @param syndrome_flags Flags indicating whether the decoder converged to a codeword or not
     * @param stats Statistics object
     * @param num_of_decoding_frames Number of frames
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
     * @return std::tuple<double, int> Number of CRC mismatches and the number of undetected errors
     */
    std::tuple<double, int> check_CRC(const std::vector<std::vector<uint8_t>> &CRC_decoded, const std::vector<std::vector<uint8_t>> &CRC_QRNG, const std::vector<int> &syndrome_flags, utilities::statistics &stats, size_t num_of_decoding_frames, bool print_flag = false);

    /**
     * @brief Calculate the syndromes of the QRNG output.
     *
     * @param LDPC_decoder LDPC decoder object
     * @param QRNG_output QRNG output
     * @param num_of_decoding_frames Number of frames
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
     * @return auto Syndromes of the QRNG output
     */
    std::vector<std::vector<uint8_t>> calculate_syndrome(const decoder LDPC_decoder, const std::vector<std::vector<int8_t>> & QRNG_output, const size_t num_of_decoding_frames, bool print_flag = false);

    /**
     * @brief Sample uniform numbers to be used as raw key material, perform multi-dimensional reconciliation (MDR) and error correction, and check the CRC values of the decoded frames.
     *
     * @param alice_states_input Alice's quantum states
     * @param bob_states_input Bob's quantum states
     * @param rate Code rate
     * @param noise_variance Noise variance of the channel
     * @param NoI Maximum number of iterations for the SPA decoder, default = 500
     * @param MDR_dim Dimension of the MDR, default = 8
     * @param layered_flag Flag to indicate whether to use layered decoding or not, default = true
     * @param fast_flag Flag to indicate whether to use the fast decoder or not, default = true
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
     * @param H_file_name Name of the file containing the LDPC code parity-check matrix, default is defined in CMakelists.txt file
     * @param lifting_factor Lifting factor of the LDPC code, default = 5000
     * @return utilities::statistics Statistics object containing the simulation results
     */
    utilities::statistics reconcile(const std::vector<double> &alice_states_input, const std::vector<double> &bob_states_input, const double rate, const double noise_variance, const size_t NoI = 500, const int MDR_dim = 8, bool layered_flag = true, bool fast_flag = true, bool print_flag = false, const std::string H_file_name = PCM_NAME, int lifting_factor = 5000);

     /**
     * @brief Samples uniform numbers to be used as raw key material, perform the operations required at Bob's side for multi-dimensional reconciliation (MDR) and output the messages to be sent through the classical channel and the raw keys.
     *
     * @param bob_states_input Bob's quantum states
     * @param beta reconciliation efficiency
     * @param SNR  Signal-to-noise ratio of the channel
     * @param MDR_dim Dimension of the MDR, default = 8
     * @param print_flag Flag to indicate whether to print information, default = false
     * @return the classical channel message, syndromes of the decoded frames, the normalization factors of Bob's measurements and the raw keys
     */
    std::tuple<std::vector<double>,std::vector<std::vector<uint8_t>>,std::vector<double>,std::vector<std::vector<int8_t>>> reconcile_Bob( const std::vector<double> &bob_states_input, const double beta, const double SNR,  const int MDR_dim = 8, bool print_flag = false);

         /**
     * @brief Perform the operations required at Alice's side for multi-dimensional reconciliation (MDR) using the classical channel message and the transmitted states. 
     *
     * @param alice_states_input Alice's quantum states
     * @param channel_message Classical channel message from Bob
     * @param syndrome Syndromes of the QRNG output
     * @param bob_normalization_vector Normalization vector for Bob's measurements
     * @param SNR Signal-to-noise ratio of the channel
     * @param MDR_dim Dimension of the MDR, default = 8
     * @param NoI Maximum number of iterations for the SPA decoder, default = 500
     * @param layered_flag Flag to indicate whether to use layered decoding or not, default = true
     * @param fast_flag Flag to indicate whether to use the fast decoder or not, default = true
     * @param print_flag Flag to indicate whether to print information, default = false
     * @param H_file_name Name of the file containing the LDPC code parity-check matrix, default is defined in CMakelists.txt file
     * @param lifting_factor Lifting factor of the LDPC code, default = 5000
     * @return the CRCs of the decoded frames, flags for the frames that could not be decoded, and the decoded frames
     */
    std::tuple<std::vector<std::vector<uint8_t>>, std::vector<bool>, std::vector<std::vector<uint8_t>>> reconcile_Alice( const std::vector<double> &alice_states_input, std::vector<double> channel_message, std::vector<std::vector<uint8_t>> syndrome,   std::vector<double> bob_normalization_vector, const double SNR,  const int MDR_dim = 8, size_t NoI = 500, bool layered_flag = true, bool fast_flag = true, bool print_flag = false, const std::string H_file_name = PCM_NAME, int lifting_factor = 5000);

    /** 
     * @brief Compare the CRC values of the decoded frames and the raw keys.
     * 
     * @param CRC_decoded CRC values of the decoded frames
     * @param CRC_QRNG CRC values of the raw keys
     * @param not_a_CW Flags indicating whether the decoder converged to a codeword or not
     * 
     * @return Flags indicating whether the frames are to be discarded or not
     */
    std::vector<bool> get_discard_flags(const std::vector<std::vector<uint8_t>> &CRC_decoded, const std::vector<std::vector<uint8_t>> &CRC_QRNG, const std::vector<bool> &not_a_CW);

    /**
     * @brief Check the CRC values of the decoded frames of Alice against the CRC values of the raw keys. Then discard the frames that do not match the CRC values of the raw keys or could not be decoded. 
     * 
     * @param CRC_decoded CRC values of the decoded frames
     * @param CRC_QRNG CRC values of the raw keys
     * @param not_a_CW Flags indicating whether the decoder converged to a codeword or not
     * 
     * @return std::tuple<std::vector<bool>, std::vector<std::vector<int8_t>>> Tuple containing the flags indicating whether the frames are to be discarded or not and the reconciled raw keys
     */
    std::tuple<std::vector<bool>, std::vector<std::vector<int8_t>>> Bob_CRC_check(std::vector<std::vector<int8_t>> QRNG_output, std::vector<std::vector<uint8_t>> CRC_decoded, std::vector<bool> frame_not_a_CW);

    extern template std::vector<double> MDR::multiplication_Alice<std::vector<double>, std::vector<double>>(
        const std::vector<double>&, const std::vector<double>&) const;
    
    extern template std::vector<double> MDR::multiplication_Bob<std::vector<double>, std::vector<double>>(
        const std::vector<double>&, const std::vector<double>&) const;

    extern template std::vector<uint8_t> decoder::get_syndrome<std::vector<uint8_t>>(const std::vector<uint8_t> &) const;

    extern template std::tuple<std::vector<uint8_t>, int, bool> decoder::SPA_decoder<std::vector<double>,std::vector<uint8_t>>(const std::vector<double>&, const std::vector<uint8_t>&) const;    

    extern template int decoder::get_decoding_error_count<std::vector<uint8_t>,std::vector<uint8_t>>(const std::vector<uint8_t> &, const std::vector<uint8_t> &) const;

    extern template std::vector<double> MDR::compute_higher_dimension_conjugate<std::vector<double>>(const std::vector<double> &) const;

    extern template std::vector<int> CRC::get_CRC_checksum<std::vector<uint8_t>>(const std::vector<uint8_t> &);    

    extern template std::vector<std::vector<uint8_t>> get_CRC_values<uint8_t>(const std::vector<std::vector<uint8_t>> &, CRC &, size_t);    

}
