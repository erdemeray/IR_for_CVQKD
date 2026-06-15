/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: statistics.hpp
 * Description: This file defines the class "statistics" used for monitoring the statistics of information reconciliation.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First stable release
 *    16/08/2024 - v1.1 - Added get_average_time() method 
 ********************************************************************/

#pragma once

namespace utilities
{
  /**
   * @brief This class is used to store and manage the statistics of the simulation.
   *
   */
  class statistics
  {
  private:
    
    ///> Vector to store the number of bit errors for each frame
    std::vector<double> bit_error_counts;                  
    
    ///> Vector to store whether the frame is decoded correctly or not for each frame
    std::vector<int> frame_error_counts;                  
    
    ///> Vector to store the number of iterations it takes to decode each frame
    std::vector<double> num_of_iterations;                 
    
    ///> Vector to store whether the decoder converged to a wrong codeword for each frame
    std::vector<int> converged_wrong_code;                 
    
    ///> Variable to store the start time of the simulation
    std::chrono::system_clock::time_point start_time;
    
    ///> Variable to store the end time of the simulation
    std::chrono::system_clock::time_point end_time;   
    
    ///> Variable to store the elapsed time of the simulation
    double elapsed_time;                                   
    
    ///> Blocklength of the codeword
    std::size_t N;                                              
    
    ///> Number of frames to be simulated
    std::size_t num_of_frames;  

    ///> Total Number of simulated frames
    std::size_t total_num_of_frames;                                                                  
    
    ///> CRC length
    std::size_t crc_length;                                     
    
    ///> Dimension of the multi-dimensional reconciliation
    std::size_t MDR_dim;                                        
    
    ///> Noise variance of the channel
    double noise_variance;            
    
    ///> Number of undetected frame errors after CRC
    std::size_t No_undetected_error_after_CRC;   
    
    ///> Number of wrong codewords detected by CRC
    std::size_t Wrong_codewords_detected_by_CRC;
    
    ///> Total number of quantum states
    std::size_t total_number_of_states; 
    
    ///> Total number of quantum states
    std::size_t number_of_unused_states; 
    
    ///> Rate of the code
    double rate;                    
    
    ///> Flag to indicate whether the decoder uses layered decoding or not
    bool layered_flag;               
  
  public:
    /**
     * @brief Construct a new statistics object and initialize the variables.
     *
     * @param blocklength Codeword blocklength
     * @param num_of_frames Number of frames to be simulated
     * @param MDR_dim Dimension of the MDR
     * @param layered_flag Flag to indicate whether the decoder uses layered decoding or not
     * @param crc_length CRC length
     */
    statistics(std::size_t blocklength, std::size_t num_of_frames, std::size_t MDR_dim, bool layered_flag, std::size_t crc_length = 32);

    /**
     * @brief Construct a new statistics object.
     *
     */
    statistics();

/**
 * @brief Overload the addition operator to add two statistics objects.
 * 
 * @param rhs Right-hand side statistics object
 * @return statistics Sum of the two statistics objects
 */
    statistics operator+(const statistics& rhs) const;

    /**
     * @brief Set the rate of the code.
     *
     * @param rate Rate of the code
     */
    void set_rate(double rate);

    /**
     * @brief Set the total number of quantum states.
     *
     * @param total_number_of_states Total number of quantum states
     */
    void set_num_of_quantum_states(std::size_t total_number_of_states);

    /**
     * @brief Set the number of unused quantum states.
     *
     * @param number_of_unused_states Number of unused quantum states
     */
    void set_number_of_unused_states(std::size_t number_of_unused_states);

    /**
     * @brief Set the bit error count for a given frame index.
     *
     * @param frame_index Frame index
     * @param bit_error_count Bit error count
     */
    void set_bit_error_count(std::size_t frame_index, double bit_error_count);

    /**
     * @brief Set the frame error count for a given frame index.
     *
     * @param frame_index Frame index
     */
    void set_frame_error_count(std::size_t frame_index);

    /**
     * @brief Set the number of iterations for a given frame index.
     *
     * @param frame_index Frame index
     * @param num_of_iteration Number of iterations
     */
    void set_num_of_iterations(std::size_t frame_index, double num_of_iteration);

    /**
     * @brief Mark a given frame as wrongly decoded codeword.
     *
     * @param frame_index Frame index
     */
    void mark_codeword_as_wrong(std::size_t frame_index);

    /**
     * @brief Set the number of undetected frame errors after CRC check.
     *
     * @param No_undetected_errors Number of undetected errors
     */
    void set_count_of_undetectable_error_after_CRC(std::size_t No_undetected_errors);

    /**
     * @brief Set the number of wrong codewords detected by CRC.
     *
     * @param No_wrong_codewords Number of wrong codewords
     */
    void set_count_of_wrong_codewords_detected_by_CRC(std::size_t No_wrong_codewords);

        /**
     * @brief Set the noise variance of the channel.
     *
     * @param noise_variance Noise variance of the channel
     */
    void set_noise_variance(double noise_variance);

    /**
     * @brief Start the timing of the simulation.
     *
     */
    void start_timing();

    /**
     * @brief End the timing of the simulation.
     *
     */
    void end_timing();

    /**
     * @brief Get the number of undetected frame errors after CRC check.
     *
     * @return int Number of undetected frame errors
     */
    int get_count_of_undetectable_error_after_CRC()  const;

    /**
     * @brief Get the number of wrong codewords detected by CRC.
     *
     * @return int Number of wrong codewords
     */
    int get_count_of_wrong_codewords_detected_by_CRC()  const;

    /**
     * @brief Get the elapsed time of the simulation.
     *
     * @return double Elapsed time of the simulation
     */
    double get_elapsed_time()  const;

        /**
     * @brief Get the average decoding duration per frame per thread of the simulation.
     *
     * @return double Average decoding duration per frame per thread in seconds
     */
    double get_average_time()  const;

    /**
     * @brief Get the bit error count for a given frame index.
     *
     * @param frame_index Frame index
     * @return double Bit error count
     */
    double get_bit_error_count(size_t frame_index)  const;

    /**
     * @brief Get the frame error count for a given frame index.
     *
     * @param frame_index Frame index
     * @return int Frame error count
     */
    int get_frame_error_count(size_t frame_index)  const;

    /**
     * @brief Get the number of decoding iterations for a given frame index.
     *
     * @param frame_index Frame index
     * @return double Number of decoding iterations
     */
    double get_num_of_iterations(size_t frame_index)   const;

    /**
     * @brief Get the bit error rate.
     *
     * @return double Bit error rate
     */
    double get_bit_error_rate()  const;

    /**
     * @brief Get the frame error rate.
     *
     * @return double Frame error rate
     */
    double get_frame_error_rate()  const;

    /**
     * @brief Get the average number of decoding iterations.
     *
     * @return double Average number of decoding iterations
     */
    double get_average_num_of_iterations()  const;

    /**
     * @brief Get the total number of frame errors.
     *
     * @return int Total number of frame errors
     */
    int get_total_frame_error_count()  const;

    /**
     * @brief Get the number of decoding frames.
     *
     * @return int Number of decoding frames
     */
    int get_num_of_decoding_frames() const;

    /**
     * @brief Get the number of wrong codewords.
     *
     * @return int Number of wrong codewords
     */
    int get_count_of_wrong_codewords()  const;

    /**
     * @brief Get the total number of quantum states.
     *
     * @return int Total number of quantum states
     */
    int get_total_number_of_states()  const;

    /**
     * @brief Get the number of unused quantum states.
     *
     * @return int Number of unused quantum states
     */
    int get_number_of_unused_states()  const;

    /**
     * @brief Get the rate of the code.
     *
     * @return double Rate of the code
     */
    double get_rate()  const;

    /**
     * @brief Get the blocklength of the codeword.
     *
     * @return int Blocklength of the codeword
     */
    int get_blocklength()  const;

    /**
     * @brief Get the CRC length.
     *
     * @return int CRC length
     */
    int get_crc_length()  const;

    /**
     * @brief Get the dimension of the multi-dimensional reconciliation.
     *
     * @return int Dimension of the multi-dimensional reconciliation
     */
    int get_MDR_dim()  const;

    /**
     * @brief Get the noise variance of the channel.
     *
     * @return double Noise variance of the channel
     */
    double get_noise_variance()  const;

    /**
     * @brief Get the layered decoding flag.
     *
     * @return int Layered decoding flag
     */
    int get_layered_flag()  const;

    /**
     * @brief Get the average statistics of the simulation.
     *
     *
     * The tuple includes:
     * - Frame error rate
     * - Bit error rate
     * - Total number of frame errors
     * - Total number of wrong codewords
     * - Total number of frames
     * - Average number of iterations
     * - Total simulation time
     * - Simulation time per frame, per thread
     *
     * @return std::tuple<double, double, int, int, int, double, double, double> Tuple of average statistics
     */
    std::tuple<double, double, int, int, int, double, double, double> get_average_statistics() const;
  };

}
