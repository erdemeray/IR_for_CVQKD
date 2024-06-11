/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: utilities.cpp
 * Description: 
 *    This file contains the implementation of the classes/methods described in utilities.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 ********************************************************************/

#include "../include/h_files/utilities.hpp"

namespace utilities
{
    statistics::statistics(size_t blocklength, size_t num_of_frames, size_t MDR_dim, int layered_flag, size_t crc_length)
    {
        bit_error_counts.resize(num_of_frames);
        frame_error_counts.resize(num_of_frames);
        num_of_iterations.resize(num_of_frames);
        converged_wrong_code.resize(num_of_frames);
        N = blocklength;
        this->num_of_frames = num_of_frames;
        elapsed_time = 0.0;
        No_undetected_error_after_CRC = 0;
        Wrong_codewords_detected_by_CRC = 0;
        this->MDR_dim = MDR_dim;
        this->crc_length = crc_length;
        this->layered_flag = layered_flag; 

    }

    void statistics::set_noise_variance(double noise_variance)
    {
        this->noise_variance = noise_variance;
    }

    void statistics::set_rate(double rate)
    {
        if (rate < 0 || rate > 1)
            throw std::invalid_argument("statistics::set_rate(): Rate must be between 0 and 1");
        
        this->rate = rate;
    }

    void statistics::set_num_of_quantum_states(size_t total_number_of_states)
    {
        this->total_number_of_states = total_number_of_states;
    }

    void statistics::set_number_of_unused_states(size_t number_of_unused_states)
    {
        this->number_of_unused_states = number_of_unused_states;
    }

    void statistics::set_bit_error_count(size_t frame_index, double bit_error_count)
    {
        
         if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::set_bit_error_count(): Invalid frame index.");
    
        bit_error_counts[frame_index] = bit_error_count;
    }

    void statistics::set_frame_error_count(size_t frame_index)
    {
        if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::set_frame_error_count(): Invalid frame index.");

        frame_error_counts[frame_index] = 1;
    }

    void statistics::set_num_of_iterations(size_t frame_index, double num_of_iteration)
    {
        if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::set_num_of_iterations(): Invalid frame index.");

        num_of_iterations[frame_index] = num_of_iteration;
    }

    void statistics::mark_codeword_as_wrong(size_t frame_index)
    {
                 if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::mark_codeword_as_wrong(): Invalid frame index.");

        
        this->converged_wrong_code[frame_index] = 1;
    }

    void statistics::set_count_of_undetectable_error_after_CRC(size_t No_undetected_errors)
    {
        No_undetected_error_after_CRC += No_undetected_errors;
    }

    void statistics::set_count_of_wrong_codewords_detected_by_CRC(size_t No_wrong_codewords)
    {
        Wrong_codewords_detected_by_CRC += No_wrong_codewords;
    }

    void statistics::start_timing()
    {
        start_time = std::chrono::system_clock::now();
    }

    void statistics::end_timing()
    {
        end_time = std::chrono::system_clock::now();
        elapsed_time = (double)(std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());
    }

    int statistics::get_count_of_undetectable_error_after_CRC()
    {
        return No_undetected_error_after_CRC;
    }

    int statistics::get_count_of_wrong_codewords_detected_by_CRC()
    {
        return Wrong_codewords_detected_by_CRC;
    }

    double statistics::get_elapsed_time()
    {
        return (elapsed_time );
    }

    double statistics::get_bit_error_count(size_t frame_index)
    {
                 if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::get_bit_error_count(): Invalid frame index.");
        
        return bit_error_counts[frame_index];
    }

    int statistics::get_frame_error_count(size_t frame_index)
    {
                 if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::get_frame_error_count(): Invalid frame index.");
        
        return frame_error_counts[frame_index];
    }

    double statistics::get_num_of_iterations(size_t frame_index)
    {
        
        if (frame_index < 0 || frame_index > num_of_frames-1) 
            throw std::invalid_argument("statistics::get_num_of_iterations(): Invalid frame index.");
        
        return num_of_iterations[frame_index];
    }

    double statistics::get_bit_error_rate()
    {
        return (std::accumulate(bit_error_counts.begin(), bit_error_counts.end(), 0.0)) / num_of_frames / N;
    }

    double statistics::get_frame_error_rate()
    {
        return (std::accumulate(frame_error_counts.begin(), frame_error_counts.end(), 0.0)) / num_of_frames;
    }

    double statistics::get_average_num_of_iterations()
    {
        return (std::accumulate(num_of_iterations.begin(), num_of_iterations.end(), 0.0)) / num_of_frames;
    }

    int statistics::get_total_frame_error_count()
    {
        return std::accumulate(frame_error_counts.begin(), frame_error_counts.end(), 0.0);
    }

    int statistics::get_num_of_decoding_frames()
    {
        return num_of_frames;
    }

    int statistics::get_count_of_wrong_codewords()
    {
        return (std::accumulate(converged_wrong_code.begin(), converged_wrong_code.end(), 0));
    }

    int statistics::get_total_number_of_states()
    {
        return total_number_of_states;
    }

    int statistics::get_number_of_unused_states()
    {
        return number_of_unused_states;
    }

    double statistics::get_rate()
    {
        return rate;
    }

    int statistics::get_blocklength()
    {
        return N;
    }

    int statistics::get_crc_length()
    {
        return crc_length;
    }

    int statistics::get_MDR_dim()
    {
        return MDR_dim;
    }

    double statistics::get_noise_variance()
    {
        return noise_variance;
    }

    int statistics::get_layered_flag()
    {
        return layered_flag;
    }

    std::tuple<double, double, int, int, int, double, double, double> statistics::get_average_statistics()
    {

        int num_of_threads = omp_get_num_threads();
#pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();
        }

        double elapsed_time_per_thread_per_frame = get_elapsed_time() / (double)(num_of_frames)*num_of_threads;

        return std::make_tuple(get_frame_error_rate(), get_bit_error_rate(), get_total_frame_error_count(), get_count_of_wrong_codewords(), get_num_of_decoding_frames(), get_average_num_of_iterations(), get_elapsed_time(), elapsed_time_per_thread_per_frame);
    }

    print::print(int flag)
    {
        this->print_flag = flag;
    }

    void print::print_statistics_report(statistics stats)
    {

        if (print_flag == 1)
        {
            std::cout << "Input Parameters: " << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Total number of states read" <<":"  << stats.get_total_number_of_states() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Number of unused states" <<":"  << stats.get_number_of_unused_states() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Total number of used states" <<":"  << stats.get_total_number_of_states() - stats.get_number_of_unused_states() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "SNR" <<":"  << 10*std::log10(1/stats.get_noise_variance()) << " dB" << std::endl;
            std::cout << std::endl;

            std::cout << "Multi-dimensional Reconciliation Parameters: " << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout << std::left <<std::setw(35) <<  "MDR dimension" <<":"  << stats.get_MDR_dim() << std::endl;
            std::cout << std::endl;

            std::cout << "Decoding Statistics: " << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout << std::left <<std::setw(35) << "Code rate"  << ":" << stats.get_rate() << std::endl;
            std::cout << std::left <<std::setw(35) << "Layered decoding" <<":"  << stats.get_layered_flag() << std::endl;
            std::cout << std::left <<std::setw(35) << "Blocklength" <<":"  << stats.get_blocklength() << std::endl;
            std::cout << std::left <<std::setw(35) << "Total number of frames" <<":"  << stats.get_num_of_decoding_frames() << std::endl;
            std::cout << std::left <<std::setw(35) << "Bit error rate" <<":"  << stats.get_bit_error_rate() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Frame error rate" <<":"  << stats.get_frame_error_rate() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Avg. number of iterations" <<":"  << stats.get_average_num_of_iterations() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Frame error count/total NoF" <<":"  << stats.get_total_frame_error_count() << " / " << stats.get_num_of_decoding_frames() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Total wrong codeword count" <<":" << stats.get_count_of_wrong_codewords() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Decoding - Elapsed time" <<":"  << stats.get_elapsed_time() << " seconds" << std::endl;

                        int num_of_threads = omp_get_num_threads();
#pragma omp parallel
            {
                num_of_threads = omp_get_num_threads();
            }

            std::cout << std::left <<std::setw(35) << "Decoding duration" <<":"  << stats.get_elapsed_time() / (double)(stats.get_num_of_decoding_frames()) * num_of_threads << " seconds/frame per thread" << std::endl;
            std::cout << std::endl;    

            std::cout << "CRC Statistics: " << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout  << std::left <<std::setw(35) <<   "CRC length" <<":" << stats.get_crc_length() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Wrong codewords detected by CRC" <<":"  << stats.get_count_of_wrong_codewords_detected_by_CRC() << std::endl;
            std::cout << std::left <<std::setw(35) <<  "Undetected frame errors after CRC" <<":" << stats.get_count_of_undetectable_error_after_CRC() << std::endl;
            std::cout << std::endl;

            std::cout << "Output Statistics: " << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout << std::left <<std::setw(35) << "Total number output bits" <<":"  << (stats.get_num_of_decoding_frames() - stats.get_total_frame_error_count()) * stats.get_blocklength() << std::endl;


        }
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
        std::vector<double> states;
        double state;

        if (file.is_open())
        {
            while (file >> state)
            {
                states.push_back(state);
            }
            file.close();
        }
        else
        {
            throw std::runtime_error("Unable to open the file: " + filename);
            // return 0;
        }

        return states;
    }
} // namespace utilities