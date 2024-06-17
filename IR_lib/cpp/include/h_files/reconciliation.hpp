/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: reconciliation.hpp
 * Description: 
 *    This file contains the classes and functions to perform information reconciliation.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
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

namespace reconciliation
{
    /**
     * @brief This class is used to decode LDPC codes using the sum-product algorithm (SPA).
     *
     */
    class decoder
    {
    private:

        ///> Vector of the degrees of each check node
        std::vector<size_t> CN_degrees;    
        
        ///> Vector of the variable node indices connected to each check node
        std::vector<size_t> Connected_VNs; 
        
        ///> Maximum degree of the check nodes in the LDPC code
        size_t max_CN_degree;              
        
        ///> Number of variable nodes in the LDPC code
        size_t number_of_VNs;              
        
        ///> Lookup table for the hyperbolic tangent function
        std::vector<double> tanh_table;    
        
        ///> Lookup table for the inverse hyperbolic tangent function
        std::vector<double> atanh_table;   
        
        ///> Blocklength of the codeword
        size_t N;                          
        
        ///> Number of check nodes in the LDPC code
        size_t M;                     

        ///> Maximum number of iterations for the SPA decoder
        size_t max_iter;                   
        
        ///> Code rate. It must be between 0.2 and 0.01.
        double rate;                       
        
        ///> Lifting factor of the LDPC code
        int lifting_factor;                
        
        ///> Inverse of the step size used for the tanh lookup table
        double inverse_stepsize_tanh;      
        
        ///> Inverse of the step size used for the atanh lookup table
        double inverse_stepsize_atanh;     
        
        ///> Flag to indicate whether to use the LUT for the decoding or not
        int fast_decoder;                  
        
        ///> Flag to indicate whether to print the decoding results or not
        int print_flag;                    
        
        ///> Flag to indicate whether to use the layered decoding or not
        int layered_flag;                  

    public:
        /**
         * @brief Construct a new decoder object.
         *
         * @param NoI Maximum number of iterations for the SPA decoder, default = 500
         * @param layered_flag Flag to indicate whether to use layered decoding or not, default = 1
         * @param fast_flag Flag to indicate whether to use the fast decoder or not, default = 0
         * @param H_file_name Name of the file containing the LDPC code parity-check matrix, default = "R_0p2_R_0p01_RA"
         * @param lifting_factor Lifting factor of the LDPC code, default = 5000
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
         */
        decoder(size_t NoI = 500, int layered_flag = 1, int fast_flag = 0, const std::string H_file_name = "R_0p2_R_0p01_RA", int lifting_factor = 5000, int print_flag = 0);

        /**
         * @brief Get the Blocklength of the codeword.
         *
         * @return size_t Blocklength of the codeword
         */
        size_t get_N();

        /**
         * @brief Get the number of check nodes in the LDPC code.
         *
         * @return size_t Number of check nodes
         */
        size_t get_M();

        /**
         * @brief Set the code rate.
         *
         * @param rate Code rate
         */
        void set_rate(double rate);

        /**
         * @brief Get the code rate.
         *
         * @return rate Code rate
         */
        double get_rate();

        /**
         * @brief Calculate the syndrome of a given word.
         *
         * @param word Codeword
         * @return std::vector<uint8_t> Syndrome of the word
         */
        template <typename T>
        std::vector<uint8_t> get_syndrome(T word) const;

        /**
         * @brief Decode the received LLRs using the SPA.
         *
         * @param llr_in Received LLRs
         * @param syndrome Syndrome of the received codeword
         * @return std::tuple<std::vector<uint8_t>, int, int> Tuple containing the decoded codeword, the number of iterations taken to decode the codeword, and a flag indicating whether the decoder converged to a codeword or not
         */
        template <typename T1, typename T2>
        std::tuple<std::vector<uint8_t>, int, int> SPA_decoder(const T1 &llr_in, const T2 &syndrome) const;

        /**
         * @brief Calculate the number of bit errors between the decoded word and the true codeword.
         *
         * @param decoded_message Decoded word
         * @param true_message True codeword
         * @return int Number of bit errors
         */
        template <typename T1, typename T2>
        int get_decoding_error_count(const T1 &decoded_message, const T2 &true_message) const;

    private:
        /**
         * @brief Lookup table implementation of the hyperbolic tangent function.
         *
         * @param a Input to the hyperbolic tangent function
         * @return double Output of the hyperbolic tangent function
         */
        double tanh_lookup(double a) const;

        /**
         * @brief Lookup table implementation of the inverse hyperbolic tangent function.
         *
         * @param a Input to the inverse hyperbolic tangent function
         * @return double Output of the inverse hyperbolic tangent function
         */
        double atanh_lookup(double a) const;

        /**
         * @brief Read the LDPC code parity-check matrix from a file.
         *
         * @param H_file_name Name of the file containing the LDPC code parity-check matrix
         */
        void read_H(const std::string H_file_name = "R_0p2_R_0p01_RA");
    };

    /**
     * @brief This class is used to perform multi-dimensional reconciliation (MDR).
     *
     */
    class MDR
    {
    private:

        ///> Dimension of the MDR
        size_t dimension; 

        ///> Flag to indicate whether to print the statistics report or not    
        int print_flag;   

    public:
        /**
         * @brief Construct a new MDR object.
         *
         * @param dim Dimension of the MDR, default = 8
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
         */
        MDR(int dim = 8, int print_flag = 0);

        /**
         * @brief Normalize the input vector.
         *
         * @param x Input vector
         * @return std::vector<double> Normalization vector
         */
        std::vector<double> normalize(std::vector<double> &x) const;

        /**
         * @brief Perform the MDR multiplication at Alice's side.
         *
         * @param channel_message Channel message
         * @param q_states Quantum states
         * @return std::vector<double> MDR output at Alice's side
         */
        template <typename T1, typename T2>
        std::vector<double> multiplication_Alice(const T1 &channel_message, const T2 &q_states) const;

        /**
         * @brief Perform the MDR multiplication at Bob's side.
         *
         * @param received_quantum_states Received quantum states
         * @param QRNG_output_bipolar QRNG output in bipolar form
         * @return std::vector<double> MDR output at Bob's side
         */
        template <typename T1, typename T2>
        std::vector<double> multiplication_Bob(const T1 &received_quantum_states, const T2 &QRNG_output_bipolar) const;

    private:
        /**
         * @brief Compute the higher dimension conjugate of a vector.
         *
         * @param x Input vector
         * @return std::vector<double> Higher dimension conjugate of the input vector
         */
        template <typename T>
        std::vector<double> compute_higher_dimension_conjugate(const T &x) const;

        /**
         * @brief Compute the Cayley-Dickson construction of two vectors.
         *
         * @param x First input vector
         * @param y Second input vector
         * @return std::vector<double> Cayley-Dickson construction of the two input vectors
         */
        std::vector<double> compute_cayley_dickson_construction(const std::vector<double> &x, const std::vector<double> &y) const;
    };

    /**
     * @brief This class is used to simulate quantum random number generator (QRNG) output.
     *
     */
    class QRNG
    {
    public:

        ///> Length of the QRNG sequence
        size_t sequence_length;
    
        ///> Number of frames
        size_t num_of_frames;   
    
        ///> Flag to indicate whether to print the statistics report or not
        int print_flag;        

        /**
         * @brief Construct a new QRNG object.
         *
         * @param seq_length Length of the QRNG sequence
         * @param NoF Number of frames
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
         */
        QRNG(size_t seq_length, size_t NoF, int print_flag = 0);

        /**
         * @brief Generate a random sequence of bits using the QRNG.
         *
         * @return std::vector<std::vector<int8_t>> QRNG output
         */
        std::vector<std::vector<int8_t>> generate_random_sequence();

    private:
        /**
         * @brief Generate a randomly seeded Mersenne Twister engine.
         *
         * @return std::mt19937 Randomly seeded Mersenne Twister engine
         */
        std::mt19937 generate_randomly_seeded_engine();

        std::uniform_int_distribution<int8_t> uniform_dist; ///> Uniform distribution used to generate random bits
    };

    /**
     * @brief This class is used to calculate the cyclic redundancy check (CRC) of a message.
     *
     */
    class CRC
    {
    private:

        ///> Generator polynomial of the CRC
        uint32_t poly; 

    public:
        /**
         * @brief Construct a new CRC object.
         *
         * @param polynomial Generator polynomial of the CRC, default = 0xEDB88320
         */
        CRC(uint32_t polynomial = 0xEDB88320);

        /**
         * @brief Calculate the CRC checksum of a message.
         *
         * @param data Message
         * @return std::vector<int> CRC checksum of the message
         */
        template <typename T>
        std::vector<int> get_CRC_checksum(const T &data);
    };

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
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
     * @return std::vector<std::vector<double>> LLR values
     */
    std::vector<std::vector<double>> get_LLR(std::vector<double> synthetic_channel_output, const std::vector<double> &alice_normalization_vector, const std::vector<double> &bob_normalization_vector, double noise_variance, int MDR_dim, size_t num_of_decoding_frames, int N, int print_flag = 0);

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
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
     * @return std::tuple<double, int> Number of CRC mismatches and the number of undetected errors
     */
    std::tuple<double, int> check_CRC(const std::vector<std::vector<uint8_t>> &CRC_decoded, const std::vector<std::vector<uint8_t>> &CRC_QRNG, const std::vector<int> &syndrome_flags, utilities::statistics &stats, size_t num_of_decoding_frames, int print_flag = 0);

    /**
     * @brief Calculate the syndromes of the QRNG output.
     *
     * @param LDPC_decoder LDPC decoder object
     * @param QRNG_output QRNG output
     * @param num_of_decoding_frames Number of frames
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
     * @return auto Syndromes of the QRNG output
     */
    auto calculate_syndrome(decoder LDPC_decoder, std::vector<std::vector<int8_t>> QRNG_output, size_t num_of_decoding_frames, int print_flag = 0);

    /**
     * @brief Sample uniform numbers to be used as raw key material, perform multi-dimensional reconciliation (MDR) and error correction, and check the CRC values of the decoded frames.
     *
     * @param alice_states_input Alice's quantum states
     * @param bob_states_input Bob's quantum states
     * @param rate Code rate
     * @param noise_variance Noise variance of the channel
     * @param NoI Maximum number of iterations for the SPA decoder, default = 500
     * @param MDR_dim Dimension of the MDR, default = 8
     * @param layered_flag Flag to indicate whether to use layered decoding or not, default = 1
     * @param fast_flag Flag to indicate whether to use the fast decoder or not, default = 1
     * @param print_flag Flag to indicate whether to print the statistics report or not, default = 0
     * @param H_file_name Name of the file containing the LDPC code parity-check matrix, default = "R_0p2_R_0p01_RA"
     * @param lifting_factor Lifting factor of the LDPC code, default = 5000
     * @return utilities::statistics Statistics object containing the simulation results
     */
    utilities::statistics reconcile(const std::vector<double> &alice_states_input, const std::vector<double> &bob_states_input, const double rate, const double noise_variance, const size_t NoI = 500, const int MDR_dim = 8, int layered_flag = 1, int fast_flag = 1, int print_flag = 0, const std::string H_file_name = "R_0p2_R_0p01_RA", int lifting_factor = 5000);

}
