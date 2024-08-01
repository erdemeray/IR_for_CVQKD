/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: MDR.hpp
 * Description: 
 *    This file contains the class "MDR" to multidimensional reconciliation as described in DOI: 10.1103/PhysRevA.77.042325.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#pragma once

namespace reconciliation
{

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
        bool print_flag;   

    public:
        /**
         * @brief Construct a new MDR object.
         *
         * @param dim Dimension of the MDR, default = 8
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
         */
        MDR(int dim = 8, bool print_flag = false);

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
   

}
