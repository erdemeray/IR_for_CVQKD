/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: QRNG.hpp
 * Description: 
 *    This file contains the class "QRNG" to simulate quantum random number generation. It uses std::mt19937 engine to generate random bits.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#pragma once

namespace reconciliation
{

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
        bool print_flag;        

        /**
         * @brief Construct a new QRNG object.
         *
         * @param seq_length Length of the QRNG sequence
         * @param NoF Number of frames
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
         */
        QRNG(size_t seq_length, size_t NoF, bool print_flag = false);

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



}
