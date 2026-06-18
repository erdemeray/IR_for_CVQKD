/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: statistics.cpp
 * Description:
 *    This file contains the implementation of the class "statistics" described in statistics.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First release version 
 *    16/08/2024 - v1.1 - Added the get_average_time() method 
 ********************************************************************/

#include "../../include/h_files/utilities.hpp"

namespace utilities
{
    statistics::statistics(std::size_t blocklength, std::size_t num_of_frames, std::size_t MDR_dim, bool layered_flag, size_t crc_length) : N(blocklength), num_of_frames(num_of_frames), MDR_dim(MDR_dim), layered_flag(layered_flag), crc_length(crc_length), total_num_of_frames(num_of_frames), elapsed_time(0.0), No_undetected_error_after_CRC(0), Wrong_codewords_detected_by_CRC(0)
    {
        bit_error_counts.resize(num_of_frames);
        frame_error_counts.resize(num_of_frames);
        num_of_iterations.resize(num_of_frames);
        converged_wrong_code.resize(num_of_frames);
    }

    statistics::statistics() : N(0), num_of_frames(0), MDR_dim(0), layered_flag(0), crc_length(0), rate(0), noise_variance(0), total_number_of_states(0), number_of_unused_states(0), total_num_of_frames(0), elapsed_time(0), No_undetected_error_after_CRC(0), Wrong_codewords_detected_by_CRC(0)
    {

    }

    statistics statistics::operator+(const statistics &rhs) const
    {

        if ((num_of_frames != 0) && (num_of_frames != rhs.num_of_frames))
        {
            throw std::invalid_argument("Cannot add statistics objects with different number of frames");
        }
        if ((rate != 0) && (rate != rhs.rate))
        {
            throw std::invalid_argument("Cannot add statistics objects with different rates");
        }
        if ((N != 0) && (N != rhs.N))
        {
            throw std::invalid_argument("Cannot add statistics objects with different blocklengths");
        }
        if ((MDR_dim != 0) && (MDR_dim != rhs.MDR_dim))
        {
            throw std::invalid_argument("Cannot add statistics objects with different MDR dimensions");
        }

        if ((crc_length != 0) && (crc_length != rhs.crc_length))
        {
            throw std::invalid_argument("Cannot add statistics objects with different CRC lengths");
        }

        if ((noise_variance != 0) && (noise_variance != rhs.noise_variance))
        {
            throw std::invalid_argument("Cannot add statistics objects with different noise variances");
        }

        statistics result(std::max(N, rhs.N), std::max(num_of_frames, rhs.num_of_frames), std::max(MDR_dim, rhs.MDR_dim), std::max(layered_flag, rhs.layered_flag), std::max(crc_length, rhs.crc_length));

        // Use the larger of the two for each property
        result.rate = std::max(rate, rhs.rate);
        result.noise_variance = std::max(noise_variance, rhs.noise_variance);

        for (size_t i = 0; i < result.num_of_frames; ++i)
        {
            if (num_of_frames > 0)
            {
                result.bit_error_counts[i] = bit_error_counts[i] + rhs.bit_error_counts[i];
                result.frame_error_counts[i] = frame_error_counts[i] + rhs.frame_error_counts[i];
                result.num_of_iterations[i] = num_of_iterations[i] + rhs.num_of_iterations[i];
                result.converged_wrong_code[i] = converged_wrong_code[i] + rhs.converged_wrong_code[i];
            }
            else
            {
                result.bit_error_counts[i] = rhs.bit_error_counts[i];
                result.frame_error_counts[i] = rhs.frame_error_counts[i];
                result.num_of_iterations[i] = rhs.num_of_iterations[i];
                result.converged_wrong_code[i] = rhs.converged_wrong_code[i];
            }
        }

        result.total_num_of_frames = total_num_of_frames + rhs.total_num_of_frames;
        result.number_of_unused_states = number_of_unused_states + rhs.number_of_unused_states;
        result.total_number_of_states = total_number_of_states + rhs.total_number_of_states;

        result.elapsed_time = (elapsed_time + rhs.elapsed_time );
        result.No_undetected_error_after_CRC = No_undetected_error_after_CRC + rhs.No_undetected_error_after_CRC;
        result.Wrong_codewords_detected_by_CRC = Wrong_codewords_detected_by_CRC + rhs.Wrong_codewords_detected_by_CRC;

        return result;
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

    void statistics::set_num_of_quantum_states(std::size_t total_number_of_states)
    {
        this->total_number_of_states = total_number_of_states;
    }

    void statistics::set_number_of_unused_states(std::size_t number_of_unused_states)
    {
        this->number_of_unused_states = number_of_unused_states;
    }

    void statistics::set_bit_error_count(std::size_t frame_index, double bit_error_count)
    {

        if (frame_index >= num_of_frames)
            throw std::out_of_range("statistics::set_bit_error_count(): Invalid frame index.");

        bit_error_counts[frame_index] = bit_error_count;
    }

    void statistics::set_frame_error_count(std::size_t frame_index)
    {
        if (frame_index >= num_of_frames)
            throw std::out_of_range("statistics::set_frame_error_count(): Invalid frame index.");

        frame_error_counts[frame_index] = 1;
    }

    void statistics::set_num_of_iterations(std::size_t frame_index, double num_of_iteration)
    {
        if (frame_index >= num_of_frames)
            throw std::out_of_range("Statistics::set_num_of_iterations(): Invalid frame index.");

        num_of_iterations[frame_index] = num_of_iteration;
    }

    void statistics::mark_codeword_as_wrong(std::size_t frame_index)
    {
        if (frame_index >= num_of_frames)
            throw std::out_of_range("Statistics::mark_codeword_as_wrong(): Invalid frame index.");

        this->converged_wrong_code[frame_index] = 1;
    }

    void statistics::set_count_of_undetectable_error_after_CRC(std::size_t No_undetected_errors)
    {
        No_undetected_error_after_CRC += No_undetected_errors;
    }

    void statistics::set_count_of_wrong_codewords_detected_by_CRC(std::size_t No_wrong_codewords)
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
        elapsed_time = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());
    }

    int statistics::get_count_of_undetectable_error_after_CRC() const
    {
        return No_undetected_error_after_CRC;
    }

    int statistics::get_count_of_wrong_codewords_detected_by_CRC() const
    {
        return Wrong_codewords_detected_by_CRC;
    }

    double statistics::get_elapsed_time() const
    {
        return (elapsed_time);
    }

    double statistics::get_average_time() const
    {

        int num_of_threads = 1;
        #pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();
        }

        return (elapsed_time / static_cast<double>(total_num_of_frames) * num_of_threads);
    }

    double statistics::get_bit_error_count(std::size_t frame_index) const
    {
        if (frame_index >= num_of_frames)
            throw std::out_of_range("statistics::get_bit_error_count(): Invalid frame index.");

        return bit_error_counts[frame_index];
    }

    int statistics::get_frame_error_count(std::size_t frame_index) const
    {
        if (frame_index >= num_of_frames)
            throw std::out_of_range("statistics::get_frame_error_count(): Invalid frame index.");

        return frame_error_counts[frame_index];
    }

    double statistics::get_num_of_iterations(std::size_t frame_index) const
    {

        if (frame_index >= num_of_frames)
            throw std::out_of_range("statistics::get_num_of_iterations(): Invalid frame index.");

        return num_of_iterations[frame_index];
    }

    double statistics::get_bit_error_rate() const
    {
        return (std::accumulate(bit_error_counts.begin(), bit_error_counts.end(), 0.0)) / (total_num_of_frames * N);
    }

    double statistics::get_frame_error_rate() const
    {
        return (std::accumulate(frame_error_counts.begin(), frame_error_counts.end(), 0.0)) / total_num_of_frames;
    }

    double statistics::get_average_num_of_iterations() const
    {
        return (std::accumulate(num_of_iterations.begin(), num_of_iterations.end(), 0.0)) / total_num_of_frames;
    }

    int statistics::get_total_frame_error_count() const
    {
        return std::accumulate(frame_error_counts.begin(), frame_error_counts.end(), 0.0);
    }

    int statistics::get_num_of_decoding_frames() const
    {
        return total_num_of_frames;
    }

    int statistics::get_count_of_wrong_codewords() const
    {
        return std::accumulate(converged_wrong_code.begin(), converged_wrong_code.end(), 0);
    }

    int statistics::get_total_number_of_states() const
    {
        return total_number_of_states;
    }

    int statistics::get_number_of_unused_states() const
    {
        return number_of_unused_states;
    }

    double statistics::get_rate() const
    {
        return rate;
    }

    int statistics::get_blocklength() const
    {
        return N;
    }

    int statistics::get_crc_length() const
    {
        return crc_length;
    }

    int statistics::get_MDR_dim() const
    {
        return MDR_dim;
    }

    double statistics::get_noise_variance() const
    {
        return noise_variance;
    }

    int statistics::get_layered_flag() const
    {
        return layered_flag;
    }

    std::tuple<double, double, int, int, int, double, double, double> statistics::get_average_statistics() const
    {

        int num_of_threads = 1;
#pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();
        }

        double elapsed_time_per_thread_per_frame = get_elapsed_time() / (static_cast<double>(total_num_of_frames))*num_of_threads;

        return std::make_tuple(get_frame_error_rate(), get_bit_error_rate(), get_total_frame_error_count(), get_count_of_wrong_codewords(), get_num_of_decoding_frames(), get_average_num_of_iterations(), get_elapsed_time(), elapsed_time_per_thread_per_frame);
    }

} // namespace utilities