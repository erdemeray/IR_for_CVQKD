/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: decoder.hpp
 * Description: 
 *    This file contains the class "decoder" to perform BP decoding. It uses SPA algorithm to decode the received message.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

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
        bool fast_decoder;                  
        
        ///> Flag to indicate whether to print the decoding results or not
        bool print_flag;                    
        
        ///> Flag to indicate whether to use the layered decoding or not
        bool layered_flag;                  

    public:
        /**
         * @brief Construct a new decoder object.
         *
         * @param NoI Maximum number of iterations for the SPA decoder, default = 500
         * @param layered_flag Flag to indicate whether to use layered decoding or not, default = true
         * @param fast_flag Flag to indicate whether to use the fast decoder or not, default = true
         * @param H_file_name Name of the file containing the LDPC code parity-check matrix, default is defined in CMakelists.txt file
         * @param lifting_factor Lifting factor of the LDPC code, default = 5000
         * @param print_flag Flag to indicate whether to print the statistics report or not, default = false
         */
        decoder(size_t NoI = 500, bool layered_flag = true, bool fast_flag = true, const std::string& H_file_name = PCM_NAME, int lifting_factor = 5000, bool print_flag = false);

        /**
         * @brief Get the Blocklength of the codeword.
         *
         * @return size_t Blocklength of the codeword
         */
        size_t get_N() const;

        /**
         * @brief Get the number of check nodes in the LDPC code.
         *
         * @return size_t Number of check nodes
         */
        size_t get_M() const;

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
        double get_rate() const;

        /**
         * @brief Calculate the syndrome of a given word.
         *
         * @param word Codeword
         * @return std::vector<uint8_t> Syndrome of the word
         */
        template <typename T>
        std::vector<uint8_t> get_syndrome(const T& word) const;

        /**
         * @brief Decode the received LLRs using the SPA.
         *
         * @param llr_in Received LLRs
         * @param syndrome Syndrome of the received codeword
         * @return std::tuple<std::vector<uint8_t>, int, bool> Tuple containing the decoded codeword, the number of iterations taken to decode the codeword, and a flag indicating whether the decoder converged to a codeword or not
         */
        template <typename T1, typename T2>
        std::tuple<std::vector<uint8_t>, int, bool> SPA_decoder(const T1 &llr_in, const T2 &syndrome) const;

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
        void read_H(const std::string& H_file_name);
    };   
}
